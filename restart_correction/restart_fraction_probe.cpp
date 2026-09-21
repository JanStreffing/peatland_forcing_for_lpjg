#include "config.h"
#include "guess.h"
#include "guessserializer.h"
#include "inputmodule.h"
#include "outputmodule.h"
#include "parameters.h"
#include "shell.h"

#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

xtring file_log = "/tmp/restart_fraction_probe.log";

static int crop_index(const char* name) {
	const char* names[10] = {"CC3ann", "CC3per", "CC3nfx", "CC4ann", "CC4per",
		"CC3anni", "CC3peri", "CC3nfxi", "CC4anni", "CC4peri"};
	for (int i = 0; i < 10; ++i)
		if (std::strcmp(name, names[i]) == 0) return i;
	return -1;
}

int main(int argc, char** argv) {
	if (argc != 6) {
		std::fprintf(stderr,
			"usage: %s INS STATE_DIR SHARD SHARDS OUTPUT\n", argv[0]);
		return 2;
	}
	const int shard = std::atoi(argv[3]);
	const int shards = std::atoi(argv[4]);
	if (shard < 0 || shards <= 0 || shard >= shards)
		return 2;

	try {
		set_shell(new CommandLineShell(file_log));
		std::auto_ptr<InputModule> input_module(
			InputModuleRegistry::get_instance().create_input_module("ece"));
		GuessOutput::OutputModuleContainer output_modules;
		GuessOutput::OutputModuleRegistry::get_instance().create_all_modules(output_modules);
		read_instruction_file(argv[1]);

		GuessDeserializer deserializer(argv[2]);
		const std::vector<std::pair<double, double> > coords =
			deserializer.get_source_coordinates();
		FILE* output = std::fopen(argv[5], "w");
		if (!output)
			throw std::runtime_error("cannot create fraction-probe output");
		std::fprintf(output,
			"index\tlon\tlat\tnstands\tphysical_total\tphysical_urban\t"
			"physical_crop\tphysical_pasture\tphysical_forest\tphysical_natural\t"
			"physical_peat\tphysical_barren\tmetadata_total\tmetadata_urban\t"
			"metadata_crop\tmetadata_pasture\tmetadata_forest\tmetadata_natural\t"
			"metadata_peat\tmetadata_barren\t"
			"physical_CC3ann\tphysical_CC3per\tphysical_CC3nfx\tphysical_CC4ann\tphysical_CC4per\t"
			"physical_CC3anni\tphysical_CC3peri\tphysical_CC3nfxi\tphysical_CC4anni\tphysical_CC4peri\t"
			"metadata_CC3ann\tmetadata_CC3per\tmetadata_CC3nfx\tmetadata_CC4ann\tmetadata_CC4per\t"
			"metadata_CC3anni\tmetadata_CC3peri\tmetadata_CC3nfxi\tmetadata_CC4anni\tmetadata_CC4peri\n");

		for (size_t index = static_cast<size_t>(shard); index < coords.size();
				index += static_cast<size_t>(shards)) {
			Gridcell gridcell;
			gridcell.set_coordinates(coords[index].first, coords[index].second);
			deserializer.deserialize_gridcell(gridcell);
			double physical[NLANDCOVERTYPES] = {0.0};
			double physical_crop[10] = {0.0};
			double metadata_crop[10] = {0.0};
			long double physical_total = 0.0L;
			for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s) {
				Stand& stand = gridcell[s];
				const double area = stand.get_gridcell_fraction();
				if (!std::isfinite(area) || area < 0.0 ||
						stand.landcover < 0 || stand.landcover >= NLANDCOVERTYPES)
					throw std::runtime_error("invalid physical stand in restart");
				physical[stand.landcover] += area;
				if (stand.landcover == CROPLAND && stand.stid >= 0 &&
						stand.stid < static_cast<int>(stlist.nobj)) {
					const int crop = crop_index(stlist[stand.stid].name);
					if (crop >= 0) physical_crop[crop] += area;
				}
				physical_total += static_cast<long double>(area);
			}
			long double metadata_total = 0.0L;
			for (unsigned int st = 0; st < gridcell.st.nobj; ++st) {
				if (stlist[st].landcover == CROPLAND) {
					const int crop = crop_index(stlist[st].name);
					if (crop >= 0) metadata_crop[crop] += gridcell.st[st].frac;
				}
			}
			for (int lc = 0; lc < NLANDCOVERTYPES; ++lc) {
				if (!std::isfinite(gridcell.landcover.frac[lc]))
					throw std::runtime_error("non-finite land-cover metadata");
				metadata_total += static_cast<long double>(gridcell.landcover.frac[lc]);
			}
			std::fprintf(output,
				"%zu\t%.17g\t%.17g\t%u\t%.17Lg\t"
				"%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17Lg\t"
				"%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t"
				"%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t"
				"%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
				index, coords[index].first, coords[index].second,
				gridcell.nbr_stands(), physical_total,
				physical[URBAN], physical[CROPLAND], physical[PASTURE],
				physical[FOREST], physical[NATURAL], physical[PEATLAND], physical[BARREN],
				metadata_total,
				gridcell.landcover.frac[URBAN], gridcell.landcover.frac[CROPLAND],
				gridcell.landcover.frac[PASTURE], gridcell.landcover.frac[FOREST],
				gridcell.landcover.frac[NATURAL], gridcell.landcover.frac[PEATLAND],
				gridcell.landcover.frac[BARREN],
				physical_crop[0], physical_crop[1], physical_crop[2], physical_crop[3], physical_crop[4],
				physical_crop[5], physical_crop[6], physical_crop[7], physical_crop[8], physical_crop[9],
				metadata_crop[0], metadata_crop[1], metadata_crop[2], metadata_crop[3], metadata_crop[4],
				metadata_crop[5], metadata_crop[6], metadata_crop[7], metadata_crop[8], metadata_crop[9]);
		}
		std::fclose(output);
	}
	catch (const std::exception& error) {
		std::fprintf(stderr, "FRACTION PROBE ERROR shard=%d: %s\n",
			shard, error.what());
		return 1;
	}
	return 0;
}
