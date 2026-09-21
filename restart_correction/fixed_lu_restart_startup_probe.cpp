#include "config.h"
#include "guess.h"
#include "guessserializer.h"
#include "inputmodule.h"
#include "outputmodule.h"
#include "parameters.h"
#include "shell.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <vector>

xtring file_log = "/tmp/fixed_lu_restart_startup_probe.log";

int main(int argc, char** argv) {
    if (argc != 7) {
        std::fprintf(stderr, "usage: %s INS STATE_DIR YEAR SHARD SHARDS OUTPUT\n", argv[0]);
        return 2;
    }
    const int year = std::atoi(argv[3]);
    const int shard = std::atoi(argv[4]);
    const int shards = std::atoi(argv[5]);
    if (shard < 0 || shards <= 0 || shard >= shards) return 2;
    try {
        set_shell(new CommandLineShell(file_log));
        std::auto_ptr<InputModule> input(
            InputModuleRegistry::get_instance().create_input_module("ece"));
        GuessOutput::OutputModuleContainer outputs;
        GuessOutput::OutputModuleRegistry::get_instance().create_all_modules(outputs);
        read_instruction_file(argv[1]);
        ecearth.fixedLUafter = 1850;
        input->init();
        date.set_first_calendar_year(year);
        date.year = 0;
        date.day = 0;

        GuessDeserializer deserializer(argv[2]);
        const std::vector<std::pair<double, double> > coords =
            deserializer.get_source_coordinates();
        size_t checked = 0;
        double maximum_st_difference = 0.0;
        double maximum_lc_difference = 0.0;
        for (size_t index = static_cast<size_t>(shard); index < coords.size();
                index += static_cast<size_t>(shards)) {
            Gridcell cell;
            cell.set_coordinates(coords[index].first, coords[index].second);
            deserializer.deserialize_gridcell(cell);
            cell.isspinup = false;
            if (!input->getgridcell(cell))
                throw std::runtime_error("could not load grid-cell forcing");

            std::vector<double> physical_st(cell.st.nobj, 0.0);
            double physical_lc[NLANDCOVERTYPES] = {0.0};
            for (unsigned int s = 0; s < cell.nbr_stands(); ++s) {
                const Stand& stand = cell[s];
                physical_st[stand.stid] += stand.get_gridcell_fraction();
                physical_lc[stand.landcover] += stand.get_gridcell_fraction();
            }
            input->getlandcover(cell);
            input->validate_landcover_state(cell, "full fixed-LU startup audit");
            for (unsigned int st = 0; st < cell.st.nobj; ++st)
                maximum_st_difference = std::max(maximum_st_difference,
                    std::fabs(cell.st[st].frac - physical_st[st]));
            for (int lc = 0; lc < NLANDCOVERTYPES; ++lc)
                maximum_lc_difference = std::max(maximum_lc_difference,
                    std::fabs(cell.landcover.frac[lc] - physical_lc[lc]));
            ++checked;
        }
        FILE* output = std::fopen(argv[6], "w");
        if (!output) throw std::runtime_error("cannot create audit output");
        std::fprintf(output,
            "shard\tshards\tcells\tmaximum_st_difference\tmaximum_lc_difference\n"
            "%d\t%d\t%zu\t%.17g\t%.17g\n",
            shard, shards, checked, maximum_st_difference, maximum_lc_difference);
        std::fclose(output);
    }
    catch (const std::exception& error) {
        std::fprintf(stderr, "STARTUP AUDIT ERROR shard=%d: %s\n", shard, error.what());
        return 1;
    }
    return 0;
}
