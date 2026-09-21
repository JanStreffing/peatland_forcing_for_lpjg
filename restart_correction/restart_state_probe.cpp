#include "config.h"
#include "externalinput.h"
#include "guess.h"
#include "guessserializer.h"
#include "inputmodule.h"
#include "outputmodule.h"
#include "parameters.h"
#include "shell.h"
#include "soilinput.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

xtring file_log = "/tmp/restart_state_probe.log";

namespace {

typedef std::pair<double, double> Coordinate;

class SoilCodes {
public:
	explicit SoilCodes(const char* path) {
		std::ifstream input(path);
		std::string header;
		if (!input || !std::getline(input, header) || header != "Lon Lat SoilCode")
			throw std::runtime_error("cannot read soil-code map");
		double lon = 0.0;
		double lat = 0.0;
		int code = 0;
		while (input >> lon >> lat >> code) {
			if (code < 1 || code > 7 || !codes.insert(
					std::make_pair(Coordinate(lon, lat), code)).second)
				throw std::runtime_error("invalid or duplicate soil-code map row");
		}
		if (codes.empty())
			throw std::runtime_error("soil-code map contains no rows");
	}

	int get(double lon, double lat) const {
		std::map<Coordinate, int>::const_iterator found = codes.find(Coordinate(lon, lat));
		if (found == codes.end())
			throw std::runtime_error("restart coordinate is absent from soil-code map");
		return found->second;
	}

private:
	std::map<Coordinate, int> codes;
};

double physical_water(Stand& stand) {
	double result = 0.0;
	for (unsigned int p = 0; p < stand.npatch(); ++p) {
		Soil& soil = stand[p].soil;
		double patch_water = soil.snowpack;
		for (int layer = 0; layer < NSOILLAYER; ++layer) {
			const int physical_layer = IDX_STD + layer;
			const double pwp = stand.is_highlatitude_peatland_stand() ?
				peat_wp : stand.get_gridcell().soiltype.water_below_wp;
			patch_water += (soil.Frac_water[physical_layer] +
				soil.Frac_ice_yesterday[physical_layer] + pwp) *
				soil.Dz[physical_layer];
		}
		result += patch_water / stand.npatch();
	}
	return result;
}

} // namespace

int main(int argc, char** argv) {
	if (argc != 7) {
		std::fprintf(stderr,
			"usage: %s INS STATE_DIR SOIL_CODES SHARD SHARDS OUTPUT\n", argv[0]);
		return 2;
	}
	const int shard = std::atoi(argv[4]);
	const int shards = std::atoi(argv[5]);
	if (shard < 0 || shards <= 0 || shard >= shards)
		return 2;

	FILE* output = std::fopen(argv[6], "w");
	if (!output) {
		std::perror(argv[5]);
		return 2;
	}

	set_shell(new CommandLineShell(file_log));
	std::auto_ptr<InputModule> input_module(
		InputModuleRegistry::get_instance().create_input_module("ece"));
	GuessOutput::OutputModuleContainer output_modules;
	GuessOutput::OutputModuleRegistry::get_instance().create_all_modules(output_modules);
	read_instruction_file(argv[1]);

	GuessDeserializer deserializer(argv[2]);
	SoilCodes soil_codes(argv[3]);
	const std::vector<std::pair<double, double> > coords =
		deserializer.get_source_coordinates();
	std::fprintf(output,
		"index\tlon\tlat\tnstands\tphysical_total\tphysical_peat\t"
		"c_total\tn_total\twater_total\tc_peat\tn_peat\twater_peat\n");

	for (size_t index = static_cast<size_t>(shard); index < coords.size();
			index += static_cast<size_t>(shards)) {
		Gridcell gridcell;
		gridcell.set_coordinates(coords[index].first, coords[index].second);
		deserializer.deserialize_gridcell(gridcell);
		if (!ifs_soil_parameters(gridcell,
				soil_codes.get(coords[index].first, coords[index].second)))
			throw std::runtime_error("IFS soil-code initialization failed");
		validate_fixed_peat_restart_soil_state(gridcell, "independent read-back");

		double physical_total = 0.0;
		double physical_peat = 0.0;
		double c_peat = 0.0;
		double n_peat = 0.0;
		double water_peat = 0.0;
		for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s) {
			Stand& stand = gridcell[s];
			const double area = stand.get_gridcell_fraction();
			physical_total += area;
			if (stand.landcover == PEATLAND) {
				physical_peat += area;
				c_peat += area * stand.ccont();
				n_peat += area * stand.ncont();
				water_peat += area * physical_water(stand);
			}
		}
		double water_total = 0.0;
		for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s)
			water_total += gridcell[s].get_gridcell_fraction() *
				physical_water(gridcell[s]);

		std::fprintf(output,
			"%zu\t%.17g\t%.17g\t%u\t%.17g\t%.17g\t%.17g\t%.17g\t"
			"%.17g\t%.17g\t%.17g\t%.17g\n",
			index, coords[index].first, coords[index].second,
			gridcell.nbr_stands(), physical_total, physical_peat,
			gridcell.ccont(), gridcell.ncont(), water_total,
			c_peat, n_peat, water_peat);
	}

	std::fclose(output);
	return 0;
}
