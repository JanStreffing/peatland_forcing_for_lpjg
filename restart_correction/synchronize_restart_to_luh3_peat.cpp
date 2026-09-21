#include "archive.h"
#include "config.h"
#include "externalinput.h"
#include "guess.h"
#include "guessserializer.h"
#include "inputmodule.h"
#include "management.h"
#include "outputmodule.h"
#include "parameters.h"
#include "shell.h"
#include "soilinput.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

xtring file_log = "/tmp/migrate_fixed_peat_restart.log";

namespace {

const double AREA_TOLERANCE = 1.0e-12;
const double CONSERVATIVE_CLOSURE_TOLERANCE = 1.0e-10;
// The TCO319 restart coordinates are stored as float-derived values, whereas
// the regridded LUH3 coordinates and text peat map use different decimal
// precision. The audited maximum mismatch is 5e-4 degree.
const double MAP_TOLERANCE = 1.0e-3;
const double BUCKET_SCALE = 1.0e4;

struct Pools {
	double c;
	double n;
	double water;
	Pools() : c(0.0), n(0.0), water(0.0) {}
};

struct MigrationAccounting {
	double recovered_1878;
	double recovered_1854;
	double recovered_synthetic;
	double removed_area;
	Pools recovered;
	Pools removed;
	MigrationAccounting()
		: recovered_1878(0.0), recovered_1854(0.0),
		  recovered_synthetic(0.0), removed_area(0.0) {}
};

struct MapRow {
	double lon;
	double lat;
	double peat;
};

struct TargetRow {
	double lon;
	double lat;
	bool luh3_valid;
	double natural;
	double cropland;
	double pasture;
	double urban;
	double peat;
	double crop_st[10];
};

typedef std::pair<long long, long long> BucketKey;

double periodic_longitude(double lon) {
	double result = std::fmod(lon + 180.0, 360.0);
	if (result < 0.0)
		result += 360.0;
	return result - 180.0;
}

double longitude_distance(double a, double b) {
	double difference = std::fabs(periodic_longitude(a) - periodic_longitude(b));
	return std::min(difference, 360.0 - difference);
}

BucketKey bucket_key(double lon, double lat) {
	return BucketKey(
		static_cast<long long>(std::llround(periodic_longitude(lon) * BUCKET_SCALE)),
		static_cast<long long>(std::llround(lat * BUCKET_SCALE)));
}

class LandcoverTargetMap {
public:
	explicit LandcoverTargetMap(const char* path) {
		std::ifstream input(path);
		if (!input)
			throw std::runtime_error(std::string("cannot open LUH3/peat target: ") + path);
		std::string header;
		if (!std::getline(input, header) || header !=
				"Lon Lat LUH3Valid RawNatural RawCropland RawPasture RawUrban CorrectedPeat "
				"Crop_CC3ann Crop_CC3per Crop_CC3nfx Crop_CC4ann Crop_CC4per "
				"Crop_CC3anni Crop_CC3peri Crop_CC3nfxi Crop_CC4anni Crop_CC4peri")
			throw std::runtime_error("invalid LUH3/peat/crop target header");
		double lon = 0.0;
		double lat = 0.0;
		int valid = 0;
		double natural = 0.0;
		double cropland = 0.0;
		double pasture = 0.0;
		double urban = 0.0;
		double peat = 0.0;
		while (input >> lon >> lat >> valid >> natural >> cropland >> pasture >> urban >> peat) {
			double crop_st[10] = {0.0};
			for (int crop = 0; crop < 10; ++crop)
				if (!(input >> crop_st[crop]))
					throw std::runtime_error("incomplete LUH3 crop target row");
			if (!std::isfinite(lon) || !std::isfinite(lat) || !std::isfinite(peat) ||
					!std::isfinite(natural) || !std::isfinite(cropland) ||
					!std::isfinite(pasture) || !std::isfinite(urban) ||
					lat < -90.0 || lat > 90.0 || (valid != 0 && valid != 1) ||
					natural < 0.0 || cropland < 0.0 || pasture < 0.0 || urban < 0.0 ||
					peat < 0.0 || peat > 1.0)
				throw std::runtime_error("invalid LUH3/peat target row");
			TargetRow row = {periodic_longitude(lon), lat, valid != 0,
				natural, cropland, pasture, urban, peat, {0.0}};
			for (int crop = 0; crop < 10; ++crop) {
				if (!std::isfinite(crop_st[crop]) || crop_st[crop] < 0.0)
					throw std::runtime_error("invalid LUH3 crop target");
				row.crop_st[crop] = crop_st[crop];
			}
			const size_t index = rows.size();
			rows.push_back(row);
			buckets[bucket_key(row.lon, row.lat)].push_back(index);
		}
		if (rows.empty())
			throw std::runtime_error("LUH3/peat target contains no data rows");
	}

	const TargetRow& lookup(double lon, double lat) const {
		const BucketKey center = bucket_key(lon, lat);
		int found = -1;
		const long long bucket_radius = static_cast<long long>(
			std::ceil(MAP_TOLERANCE * BUCKET_SCALE)) + 1;
		for (long long di = -bucket_radius; di <= bucket_radius; ++di) {
			for (long long dj = -bucket_radius; dj <= bucket_radius; ++dj) {
				const BucketKey candidate_key(center.first + di, center.second + dj);
				std::map<BucketKey, std::vector<size_t> >::const_iterator bucket =
					buckets.find(candidate_key);
				if (bucket == buckets.end())
					continue;
				for (size_t j = 0; j < bucket->second.size(); ++j) {
					const size_t index = bucket->second[j];
					const TargetRow& row = rows[index];
					if (longitude_distance(row.lon, lon) <= MAP_TOLERANCE &&
							std::fabs(row.lat - lat) <= MAP_TOLERANCE) {
						if (found >= 0 && found != static_cast<int>(index))
							throw std::runtime_error("ambiguous target coordinate lookup");
						found = static_cast<int>(index);
					}
				}
			}
		}
		if (found < 0)
			throw std::runtime_error("no target row within coordinate tolerance");
		return rows[static_cast<size_t>(found)];
	}

private:
	std::vector<TargetRow> rows;
	std::map<BucketKey, std::vector<size_t> > buckets;
};

class SoilCodeMap {
public:
	explicit SoilCodeMap(const char* path) {
		std::ifstream input(path);
		if (!input)
			throw std::runtime_error(std::string("cannot open soil-code map: ") + path);
		std::string header;
		if (!std::getline(input, header) || header != "Lon Lat SoilCode")
			throw std::runtime_error("invalid soil-code map header");
		double lon = 0.0;
		double lat = 0.0;
		int code = 0;
		while (input >> lon >> lat >> code) {
			if (!std::isfinite(lon) || !std::isfinite(lat) || lat < -90.0 ||
					lat > 90.0 || code < 1 || code > 7)
				throw std::runtime_error("invalid soil-code map row");
			MapRow row = {periodic_longitude(lon), lat, static_cast<double>(code)};
			const size_t index = rows.size();
			rows.push_back(row);
			buckets[bucket_key(row.lon, row.lat)].push_back(index);
		}
		if (rows.empty())
			throw std::runtime_error("soil-code map contains no rows");
	}

	int lookup(double lon, double lat) const {
		const BucketKey center = bucket_key(lon, lat);
		int found = -1;
		for (long long di = -1; di <= 1; ++di) {
			for (long long dj = -1; dj <= 1; ++dj) {
				const BucketKey candidate_key(center.first + di, center.second + dj);
				std::map<BucketKey, std::vector<size_t> >::const_iterator bucket =
					buckets.find(candidate_key);
				if (bucket == buckets.end()) continue;
				for (size_t j = 0; j < bucket->second.size(); ++j) {
					const size_t index = bucket->second[j];
					const MapRow& row = rows[index];
					if (longitude_distance(row.lon, lon) <= MAP_TOLERANCE &&
							std::fabs(row.lat - lat) <= MAP_TOLERANCE) {
						if (found >= 0 && found != static_cast<int>(index))
							throw std::runtime_error("ambiguous soil-code coordinate lookup");
						found = static_cast<int>(index);
					}
				}
			}
		}
		if (found < 0)
			throw std::runtime_error("no soil-code row within coordinate tolerance");
		return static_cast<int>(rows[static_cast<size_t>(found)].peat);
	}

private:
	std::vector<MapRow> rows;
	std::map<BucketKey, std::vector<size_t> > buckets;
};

struct LayerGeometry {
	double dz;
	double porosity;
	double pwp;
	double available_capacity;
};

struct LayerWater {
	double liquid;
	double ice;
	double below_pwp_liquid;
};

typedef std::vector<std::vector<LayerWater> > StandWater;

LayerGeometry layer_geometry(Stand& stand, Patch& patch, int layer) {
	LayerGeometry result;
	Soil& soil = patch.soil;
	Soiltype& soiltype = stand.get_gridcell().soiltype;
	if (stand.is_highlatitude_peatland_stand()) {
		const bool acrotelm = layer < NACROTELM;
		result.dz = acrotelm ? Dz_acro : Dz_cato;
		result.porosity = acrotelm ? soil.acro_por : soil.cato_por;
		result.pwp = peat_wp;
		result.available_capacity = result.porosity - result.pwp;
	} else {
		result.dz = Dz_soil;
		result.porosity = soiltype.porosity;
		result.pwp = soiltype.water_below_wp;
		result.available_capacity = soiltype.awc[layer] / result.dz;
	}
	result.available_capacity = std::min(result.available_capacity,
		result.porosity - result.pwp);
	if (!std::isfinite(result.dz) || result.dz <= 0.0 ||
			!std::isfinite(result.porosity) || !std::isfinite(result.pwp) ||
			!std::isfinite(result.available_capacity) || result.pwp < 0.0 ||
			result.available_capacity < 0.0 ||
			result.pwp + result.available_capacity > result.porosity + 1.0e-12)
		throw std::runtime_error("invalid destination soil geometry");
	return result;
}

void materialize_geometry(Stand& stand) {
	for (unsigned int p = 0; p < stand.npatch(); ++p) {
		Patch& patch = stand[p];
		Soil& soil = patch.soil;
		soil.IDX = IDX_STD;
		soil.ngroundl = NSOILLAYER;
		for (int layer = 0; layer < NSOILLAYER; ++layer) {
			const int physical_layer = IDX_STD + layer;
			const LayerGeometry geometry = layer_geometry(stand, patch, layer);
			if (!std::isfinite(soil.Dz[physical_layer]) ||
					std::fabs(soil.Dz[physical_layer] - geometry.dz) > 1.0e-8)
				throw std::runtime_error("restart layer thickness differs from model geometry");
			soil.por[physical_layer] = geometry.porosity;
			soil.Fpwp_ref[physical_layer] = geometry.pwp;
			soil.Frac_peat[physical_layer] = stand.is_highlatitude_peatland_stand() ?
				1.0 - (geometry.porosity + Fgas) : 0.0;
			if (stand.is_highlatitude_peatland_stand()) {
				soil.Frac_org[physical_layer] = 0.0;
				soil.Frac_min[physical_layer] = 0.0;
				stand.get_gridcell().soiltype.awc_peat[layer] =
					geometry.available_capacity * geometry.dz;
			} else {
				soil.Frac_org[physical_layer] = stand.get_gridcell().soiltype.organic_frac;
				soil.Frac_min[physical_layer] = stand.get_gridcell().soiltype.mineral_frac;
			}
		}
	}
}

StandWater capture_water(Stand& stand) {
	// Fpwp_ref/porosity are reconstructed at startup rather than serialized.
	// Materialize the source stand's own geometry before interpreting its saved
	// liquid/ice fields; retyping to the destination happens only afterwards.
	materialize_geometry(stand);
	StandWater result(stand.npatch(), std::vector<LayerWater>(NSOILLAYER));
	for (unsigned int p = 0; p < stand.npatch(); ++p) {
		Soil& soil = stand[p].soil;
		for (int layer = 0; layer < NSOILLAYER; ++layer) {
			const int physical_layer = IDX_STD + layer;
			LayerWater& state = result[p][layer];
			state.below_pwp_liquid = soil.Frac_water_belowpwp[physical_layer];
			state.liquid = soil.Frac_water[physical_layer] + state.below_pwp_liquid;
			// Gridcell/Patch water accounting uses the serialized current ice
			// fraction. The migration writes this normalized value to both current
			// and yesterday arrays so restart initialization cannot diverge.
			state.ice = soil.Frac_ice[physical_layer] +
				soil.Fpwp_ref[physical_layer] -
				state.below_pwp_liquid;
			if (!std::isfinite(state.liquid) || !std::isfinite(state.ice) ||
					state.liquid < -1.0e-10 || state.ice < -1.0e-10) {
				std::fprintf(stderr,
					"SOURCE WATER DETAIL lon=%.17g lat=%.17g lc=%d stid=%d patch=%u "
					"layer=%d Frac_water=%.17g belowPWP=%.17g Frac_ice=%.17g "
					"Frac_ice_yesterday=%.17g source_PWP=%.17g liquid=%.17g ice=%.17g\n",
					stand.get_gridcell().get_lon(), stand.get_gridcell().get_lat(),
					static_cast<int>(stand.landcover), stand.stid, p, layer,
					soil.Frac_water[physical_layer], state.below_pwp_liquid,
					soil.Frac_ice[physical_layer], soil.Frac_ice_yesterday[physical_layer],
					soil.Fpwp_ref[physical_layer], state.liquid, state.ice);
				throw std::runtime_error("restart contains invalid liquid/ice state");
			}
			state.liquid = std::max(0.0, state.liquid);
			state.ice = std::max(0.0, state.ice);
		}
	}
	return result;
}

void apply_water_to_current_geometry(Stand& stand, const StandWater& source) {
	if (source.size() != stand.npatch())
		throw std::runtime_error("patch count changed while retyping a stand");
	materialize_geometry(stand);
	for (unsigned int p = 0; p < stand.npatch(); ++p) {
		Patch& patch = stand[p];
		Soil& soil = patch.soil;
		double surface_transfer = 0.0;
		double evaporation_liquid = 0.0;
		double evaporation_capacity = 0.0;
		for (int layer = 0; layer < NSOILLAYER; ++layer) {
			const int physical_layer = IDX_STD + layer;
			const LayerGeometry geometry = layer_geometry(stand, patch, layer);
			const LayerWater& old = source[p][layer];
			const double old_total = old.liquid + old.ice;
			const double minimum = geometry.pwp;
			const double maximum = geometry.pwp + geometry.available_capacity;
			const double soil_total = std::max(minimum, std::min(old_total, maximum));
			// Positive transfer is pore water displaced to the surface; negative
			// transfer is existing surface water needed to fill the destination
			// soil's irreducible below-wilting-point reservoir.
			surface_transfer += (old_total - soil_total) * geometry.dz;
			const double phase_scale = old_total > 0.0 ? soil_total / old_total : 0.0;
			const double liquid = old_total > 0.0 ?
				old.liquid * phase_scale : soil_total;
			const double ice = old_total > 0.0 ? old.ice * phase_scale : 0.0;
			const double below_min = std::max(0.0, geometry.pwp - ice);
			const double below_max = std::min(geometry.pwp, liquid);
			if (below_min > below_max + 1.0e-10)
				throw std::runtime_error("water phase cannot be represented in destination soil");
			const double below = std::max(below_min,
				std::min(below_max, old.below_pwp_liquid));
			const double new_liquid = std::max(0.0, liquid - below);
			const double new_ice = std::max(0.0, ice - (geometry.pwp - below));

			soil.Frac_water_belowpwp[physical_layer] = below;
			soil.Frac_water[physical_layer] = new_liquid;
			soil.Frac_ice[physical_layer] = new_ice;
			soil.Frac_ice_yesterday[physical_layer] = new_ice;
			soil.Frac_air[physical_layer] = geometry.porosity - soil_total;
			soil.alwhc_init[layer] = geometry.available_capacity;
			soil.alwhc[layer] = std::max(0.0, geometry.available_capacity - new_ice);
			soil.aw_max[layer] = geometry.available_capacity * geometry.dz;
			soil.whc[layer] = std::max(0.0, soil.aw_max[layer] - new_ice * geometry.dz);
			const double relative_liquid = geometry.available_capacity > 0.0 ?
				new_liquid / geometry.available_capacity : 0.0;
			soil.set_layer_soil_water(layer, relative_liquid);
			if (layer < 2) {
				evaporation_liquid += new_liquid * geometry.dz;
				evaporation_capacity += geometry.available_capacity * geometry.dz;
			}
			if (soil.Frac_air[physical_layer] < -1.0e-10 ||
					relative_liquid > 1.0 + 1.0e-10)
				throw std::runtime_error("normalized destination soil is not physical");
		}
		soil.set_layer_soil_water_evap(evaporation_capacity > 0.0 ?
			evaporation_liquid / evaporation_capacity : 0.0);
		// Exchange water with the already serialized surface snow/water reservoir.
		// This conserves the patch total exactly for both capacity excess and a
		// destination PWP deficit; never manufacture water to satisfy geometry.
		// A temporarily negative value records water required by this patch. The
		// grid-cell rebalancer below draws it from positive serialized surface
		// reservoirs in other stands using exact area weights.
		soil.snowpack += surface_transfer;
		soil.snow_active = soil.snowpack > 1.0;
	}
}

void rebalance_gridcell_surface_water(Gridcell& gridcell) {
	long double positive = 0.0L;
	long double deficit = 0.0L;
	long double movable_soil = 0.0L;
	for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s) {
		Stand& stand = gridcell[s];
		const long double patch_weight =
			static_cast<long double>(stand.get_gridcell_fraction()) / stand.npatch();
		for (unsigned int p = 0; p < stand.npatch(); ++p) {
			Soil& soil = stand[p].soil;
			const double surface = soil.snowpack;
			if (!std::isfinite(surface))
				throw std::runtime_error("non-finite surface water after soil conversion");
			if (surface >= 0.0)
				positive += patch_weight * surface;
			else
				deficit -= patch_weight * surface;
			for (int layer = 0; layer < NSOILLAYER; ++layer) {
				const int physical_layer = IDX_STD + layer;
				const LayerGeometry geometry = layer_geometry(stand, stand[p], layer);
				movable_soil += patch_weight * geometry.dz *
					(soil.Frac_water[physical_layer] + soil.Frac_ice[physical_layer]);
			}
		}
	}
	const long double soil_draw = std::max(0.0L, deficit - positive);
	if (soil_draw > movable_soil + 1.0e-9L) {
		std::fprintf(stderr,
			"GRIDCELL WATER RESERVOIR DETAIL lon=%.17g lat=%.17g "
			"positive=%.21Lg deficit=%.21Lg movable_soil=%.21Lg "
			"shortfall=%.21Lg\n",
			gridcell.get_lon(), gridcell.get_lat(), positive, deficit,
			movable_soil, soil_draw - movable_soil);
		throw std::runtime_error(
			"destination soil PWP requires more water than the grid cell contains");
	}
	const long double surface_scale = positive > 0.0L ?
		std::max(0.0L, (positive - std::min(positive, deficit)) / positive) : 0.0L;
	const long double soil_scale = movable_soil > 0.0L ?
		std::max(0.0L, (movable_soil - soil_draw) / movable_soil) : 1.0L;
	for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s) {
		Stand& stand = gridcell[s];
		for (unsigned int p = 0; p < stand.npatch(); ++p) {
			Soil& soil = stand[p].soil;
			soil.snowpack = soil.snowpack > 0.0 ?
				static_cast<double>(soil.snowpack * surface_scale) : 0.0;
			soil.snow_active = soil.snowpack > 1.0;

			double evaporation_liquid = 0.0;
			double evaporation_capacity = 0.0;
			for (int layer = 0; layer < NSOILLAYER; ++layer) {
				const int physical_layer = IDX_STD + layer;
				const LayerGeometry geometry = layer_geometry(stand, stand[p], layer);
				soil.Frac_water[physical_layer] *= static_cast<double>(soil_scale);
				soil.Frac_ice[physical_layer] *= static_cast<double>(soil_scale);
				soil.Frac_ice_yesterday[physical_layer] = soil.Frac_ice[physical_layer];
				const double above_pwp = soil.Frac_water[physical_layer] +
					soil.Frac_ice[physical_layer];
				soil.Frac_air[physical_layer] = geometry.porosity -
					geometry.pwp - above_pwp;
				soil.alwhc_init[layer] = geometry.available_capacity;
				soil.alwhc[layer] = std::max(0.0,
					geometry.available_capacity - soil.Frac_ice[physical_layer]);
				soil.aw_max[layer] = geometry.available_capacity * geometry.dz;
				soil.whc[layer] = std::max(0.0, soil.aw_max[layer] -
					soil.Frac_ice[physical_layer] * geometry.dz);
				const double relative_liquid = geometry.available_capacity > 0.0 ?
					soil.Frac_water[physical_layer] / geometry.available_capacity : 0.0;
				soil.set_layer_soil_water(layer, relative_liquid);
				if (layer < 2) {
					evaporation_liquid += soil.Frac_water[physical_layer] * geometry.dz;
					evaporation_capacity += geometry.available_capacity * geometry.dz;
				}
			}
			soil.set_layer_soil_water_evap(evaporation_capacity > 0.0 ?
				evaporation_liquid / evaporation_capacity : 0.0);
		}
	}
}

void normalize_gridcell_hydrology(Gridcell& gridcell) {
	for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s) {
		Stand& stand = gridcell[s];
		const StandWater water = capture_water(stand);
		apply_water_to_current_geometry(stand, water);
	}
	rebalance_gridcell_surface_water(gridcell);
}

Pools gridcell_pools(Gridcell& gridcell) {
	Pools result;
	result.c = gridcell.ccont();
	result.n = gridcell.ncont();
	result.water = gridcell.water_content();
	return result;
}

void add_weighted_pools(Pools& total, Stand& stand, double area) {
	total.c += stand.ccont() * area;
	total.n += stand.ncont() * area;
	total.water += stand.water_content() * area;
}

std::vector<double> physical_st(Gridcell& gridcell) {
	std::vector<double> result(gridcell.st.nobj, 0.0);
	for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s) {
		Stand& stand = gridcell[s];
		if (stand.stid < 0 || stand.stid >= static_cast<int>(result.size()))
			throw std::runtime_error("stand has invalid stand-type id");
		const double area = stand.get_gridcell_fraction();
		if (!std::isfinite(area) || area < 0.0 || stand.npatch() == 0 ||
				stand.landcover != stlist[stand.stid].landcover)
			throw std::runtime_error("stand has invalid physical structure");
		result[stand.stid] += area;
	}
	return result;
}

double vector_sum(const std::vector<double>& values) {
	long double result = 0.0L;
	for (size_t i = 0; i < values.size(); ++i)
		result += static_cast<long double>(values[i]);
	return static_cast<double>(result);
}

int largest_nonpeat_st(const std::vector<double>& fractions) {
	int selected = -1;
	for (size_t st = 0; st < fractions.size(); ++st) {
		if (stlist[st].landcover == PEATLAND)
			continue;
		if (selected < 0 || fractions[st] > fractions[static_cast<size_t>(selected)])
			selected = static_cast<int>(st);
	}
	return selected;
}

void close_fraction_vector(std::vector<double>& fractions) {
	const int closure = largest_nonpeat_st(fractions);
	const double residual = 1.0 - vector_sum(fractions);
	if (closure < 0) {
		if (std::fabs(residual) > AREA_TOLERANCE)
			throw std::runtime_error("no non-peat stand type available for closure");
		return;
	}
	if (fractions[static_cast<size_t>(closure)] + residual < -AREA_TOLERANCE)
		throw std::runtime_error("fraction closure would create a negative stand type");
	fractions[static_cast<size_t>(closure)] =
		std::max(0.0, fractions[static_cast<size_t>(closure)] + residual);
}

std::vector<double> restart_target(Gridcell& gridcell) {
	std::vector<double> target(gridcell.st.nobj, 0.0);
	for (unsigned int st = 0; st < gridcell.st.nobj; ++st) {
		target[st] = gridcell.st[st].frac;
		if (!std::isfinite(target[st]) || target[st] < 0.0)
			throw std::runtime_error("invalid serialized stand-type fraction");
	}
	close_fraction_vector(target);
	return target;
}

void production_landcover_target(const TargetRow& row,
		double target[NLANDCOVERTYPES]) {
	Landcover lc;
	for (int i = 0; i < NLANDCOVERTYPES; ++i)
		lc.frac[i] = 0.0;
	lc.frac[NATURAL] = row.natural;
	lc.frac[CROPLAND] = row.cropland;
	lc.frac[PASTURE] = row.pasture;
	lc.frac[URBAN] = row.urban;
	close_luh3_base_fractions(lc);
	reserve_fixed_peatland(lc, row.peat);
	filter_subresolution_luh3_landcovers(lc);
	for (int i = 0; i < NLANDCOVERTYPES; ++i)
		target[i] = lc.frac[i];

	long double total = 0.0L;
	for (int i = 0; i < NLANDCOVERTYPES; ++i) {
		if (!std::isfinite(target[i]) || target[i] < 0.0)
			throw std::runtime_error("production target contains invalid land fraction");
		total += static_cast<long double>(target[i]);
	}
	if (std::fabs(static_cast<double>(total - 1.0L)) > AREA_TOLERANCE ||
			std::fabs(target[PEATLAND] - row.peat) > AREA_TOLERANCE)
		throw std::runtime_error("production LUH3/peat target is not closed");
}

int crop_target_index(const char* name) {
	const char* names[10] = {"CC3ann", "CC3per", "CC3nfx", "CC4ann", "CC4per",
		"CC3anni", "CC3peri", "CC3nfxi", "CC4anni", "CC4peri"};
	for (int i = 0; i < 10; ++i)
		if (std::strcmp(name, names[i]) == 0)
			return i;
	return -1;
}

std::vector<double> prescribed_landcover_target(
		const std::vector<double>& old_target, const TargetRow& row,
		double target_lc[NLANDCOVERTYPES]) {
	production_landcover_target(row, target_lc);

	double old_lc[NLANDCOVERTYPES] = {0.0};
	for (size_t st = 0; st < old_target.size(); ++st)
		old_lc[stlist[st].landcover] += old_target[st];

	std::vector<double> result(old_target.size(), 0.0);
	for (int lc = 0; lc < NLANDCOVERTYPES; ++lc) {
		if (target_lc[lc] <= 0.0 || lc == CROPLAND)
			continue;
		if (old_lc[lc] > 0.0) {
			for (size_t st = 0; st < result.size(); ++st)
				if (stlist[st].landcover == lc)
					result[st] = target_lc[lc] * old_target[st] / old_lc[lc];
		}
		else {
			int receiver = -1;
			for (size_t st = 0; st < result.size(); ++st)
				if (stlist[st].landcover == lc) { receiver = static_cast<int>(st); break; }
			if (receiver < 0)
				throw std::runtime_error("positive land-cover target has no configured stand type");
			result[static_cast<size_t>(receiver)] = target_lc[lc];
		}
	}

	if (target_lc[CROPLAND] > 0.0) {
		double raw_sum = 0.0;
		for (int crop = 0; crop < 10; ++crop) raw_sum += row.crop_st[crop];
		if (!(raw_sum > 0.0))
			throw std::runtime_error("positive LUH3 cropland has no crop-type target");
		int closure_st = -1;
		double largest = -1.0;
		long double assigned = 0.0L;
		for (size_t st = 0; st < result.size(); ++st) {
			if (stlist[st].landcover != CROPLAND) continue;
			const int crop = crop_target_index(stlist[st].name);
			if (crop < 0)
				throw std::runtime_error("configured crop stand has no LUH3 target mapping");
			result[st] = target_lc[CROPLAND] * row.crop_st[crop] / raw_sum;
			assigned += static_cast<long double>(result[st]);
			if (result[st] > largest) { largest = result[st]; closure_st = static_cast<int>(st); }
		}
		if (closure_st < 0)
			throw std::runtime_error("positive LUH3 cropland has no configured crop stand");
		result[closure_st] += static_cast<double>(
			static_cast<long double>(target_lc[CROPLAND]) - assigned);
	}
	close_fraction_vector(result);
	return result;
}

Stand& import_stand(Stand& source, Gridcell& destination, double fraction) {
	std::stringstream buffer;
	ArchiveOutStream output(buffer);
	source.serialize(output);
	Stand& imported = destination.create_stand(source.landcover);
	ArchiveInStream input(buffer);
	imported.serialize(input);
	imported.set_gridcell_fraction(fraction);
	return imported;
}

void move_incompatible_individual_to_litter(Individual& individual) {
	const double carbon = individual.ccont();
	const double nitrogen = individual.ncont();
	if (!std::isfinite(carbon) || !std::isfinite(nitrogen) ||
			carbon < -1.0e-12 || nitrogen < -1.0e-12)
		throw std::runtime_error("incompatible vegetation has invalid C/N pools");

	// The regular harvest helper operates on the legacy annual crop pools.  A
	// restart can instead contain the daily-growing-season crop pools used by
	// Individual::ccont(), so using that helper here can silently lose carbon.
	// This is an offline state migration, not a land-use event: move exactly the
	// pools counted by the restart to litter and create no harvest/fire flux.
	const double total_c = std::max(0.0, carbon);
	const double total_n = std::max(0.0, nitrogen);
	double root_c = 0.0;
	double wood_c = 0.0;
	if (individual.has_daily_turnover() && individual.cropindiv) {
		root_c = std::max(0.0, individual.cropindiv->grs_cmass_root);
	} else {
		root_c = std::max(0.0, individual.cmass_root);
		wood_c = std::max(0.0,
			individual.cmass_sap + individual.cmass_heart - individual.cmass_debt);
	}
	const double c_structural = root_c + wood_c;
	if (c_structural > total_c && c_structural > 0.0) {
		const double scale = total_c / c_structural;
		root_c *= scale;
		wood_c *= scale;
	}
	const double leaf_c = total_c - root_c - wood_c;

	double root_n = std::max(0.0, individual.nmass_root +
		individual.nstore_longterm + individual.nstore_labile);
	double wood_n = std::max(0.0, individual.nmass_sap + individual.nmass_heart);
	const double n_structural = root_n + wood_n;
	if (n_structural > total_n && n_structural > 0.0) {
		const double scale = total_n / n_structural;
		root_n *= scale;
		wood_n *= scale;
	}
	const double leaf_n = total_n - root_n - wood_n;

	Patchpft& litter = individual.vegetation.patch.pft[individual.pft.id];
	litter.cmass_litter_root += root_c;
	litter.cmass_litter_sap += wood_c;
	litter.cmass_litter_leaf += leaf_c;
	litter.nmass_litter_root += root_n;
	litter.nmass_litter_sap += wood_n;
	litter.nmass_litter_leaf += leaf_n;
}

void configure_retyped_stand(Stand& stand, int stid, landcovertype origin) {
	const double carbon_before = stand.ccont();
	const double nitrogen_before = stand.ncont();
	StandType& destination = stlist[stid];
	stand.landcover = destination.landcover;
	stand.stid = stid;
	stand.lc_origin = origin;
	stand.st_origin = stid;
	stand.current_rot = 0;
	stand.nyears_in_rotation = 0;
	stand.ndays_in_rotation = 0;
	stand.infallow = false;
	stand.isrotationday = false;
	stand.isirrigated = destination.get_management(0).hydrology == IRRIGATED;
	stand.hasgrassintercrop = false;
	stand.pftid = pftlist.getpftid(destination.get_management(0).pftname);

	const bool natural_vegetation = destination.naturalveg == "ALL";
	const bool natural_grass = natural_vegetation || destination.naturalveg == "GRASSONLY";
	for (unsigned int pft = 0; pft < pftlist.nobj; ++pft) {
		Pft& definition = pftlist[pft];
		Standpft& state = stand.pft[pft];
		const bool allowed = definition.landcover == destination.landcover ||
			(natural_vegetation && definition.landcover == NATURAL) ||
			(natural_grass && definition.landcover == NATURAL &&
			 definition.lifeform == GRASS);
		state.active = allowed;
		state.reestab = allowed;
		state.plant = allowed && definition.lifeform == TREE;
		state.irrigated = stand.isirrigated;
	}

	// Retain compatible vegetation. Incompatible living PFTs (for example a
	// peat-only PFT in a newly mineral NATURAL stand) are transferred to litter
	// without fire, harvest, or atmospheric flux, then removed. This preserves
	// C and N while preventing a semantically invalid stand from entering the run.
	for (unsigned int patch = 0; patch < stand.npatch(); ++patch) {
		Patch& patch_state = stand[patch];
		if (destination.landcover == NATURAL)
			patch_state.managed = false;
		patch_state.vegetation.firstobj();
		while (patch_state.vegetation.isobj) {
			Individual& individual = patch_state.vegetation.getobj();
			const Pft& definition = individual.pft;
			const bool allowed = definition.landcover == destination.landcover ||
				(natural_vegetation && definition.landcover == NATURAL) ||
				(natural_grass && definition.landcover == NATURAL &&
				 definition.lifeform == GRASS);
			if (allowed) {
				stand.pft[individual.pft.id].active = true;
				patch_state.vegetation.nextobj();
			} else {
				move_incompatible_individual_to_litter(individual);
				patch_state.vegetation.killobj();
			}
		}
	}
	const double carbon_tolerance = 1.0e-9 * std::max(1.0, std::fabs(carbon_before));
	const double nitrogen_tolerance = 1.0e-9 * std::max(1.0, std::fabs(nitrogen_before));
	const double carbon_after = stand.ccont();
	const double nitrogen_after = stand.ncont();
	if (std::fabs(carbon_after - carbon_before) > carbon_tolerance ||
			std::fabs(nitrogen_after - nitrogen_before) > nitrogen_tolerance) {
		std::fprintf(stderr,
			"RETYPE POOL DETAIL from_lc=%d to_lc=%d stid=%d "
			"C_before=%.17g C_after=%.17g C_delta=%.17g "
			"N_before=%.17g N_after=%.17g N_delta=%.17g\n",
			static_cast<int>(origin), static_cast<int>(destination.landcover), stid,
			carbon_before, carbon_after, carbon_after - carbon_before,
			nitrogen_before, nitrogen_after, nitrogen_after - nitrogen_before);
		throw std::runtime_error("vegetation-to-litter conversion did not conserve C/N");
	}
}

Stand& clone_to_type_without_pool_change(Stand& source, int stid, double fraction) {
	Gridcell& gridcell = source.get_gridcell();
	const landcovertype origin = source.landcover;
	Stand& created = import_stand(source, gridcell, fraction);
	const StandWater source_water = capture_water(created);
	configure_retyped_stand(created, stid, origin);
	apply_water_to_current_geometry(created, source_water);
	return created;
}

Stand* largest_stand(Gridcell& gridcell, int stid) {
	Stand* result = 0;
	for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s) {
		Stand& stand = gridcell[s];
		if (stand.stid == stid && stand.get_gridcell_fraction() > AREA_TOLERANCE &&
				(!result || stand.get_gridcell_fraction() > result->get_gridcell_fraction()))
			result = &stand;
	}
	return result;
}

double donor_st_total(Gridcell& gridcell, int stid) {
	double result = 0.0;
	for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s)
		if (gridcell[s].stid == stid && gridcell[s].get_gridcell_fraction() > 0.0)
			result += gridcell[s].get_gridcell_fraction();
	return result;
}

double import_historical_stands(Gridcell& destination, Gridcell& donor,
		int stid, double requested, Pools& recovered) {
	const double donor_total = donor_st_total(donor, stid);
	if (donor_total <= 0.0 || requested <= 0.0)
		return 0.0;
	double imported_total = 0.0;
	Stand* last_imported = 0;
	for (unsigned int s = 0; s < donor.nbr_stands(); ++s) {
		Stand& source = donor[s];
		if (source.stid != stid || source.get_gridcell_fraction() <= 0.0)
			continue;
		const double fraction = requested * source.get_gridcell_fraction() / donor_total;
		Stand& imported = import_stand(source, destination, fraction);
		last_imported = &imported;
		imported_total += fraction;
		add_weighted_pools(recovered, imported, fraction);
	}
	if (last_imported) {
		const double residual = requested - imported_total;
		last_imported->set_gridcell_fraction(
			last_imported->get_gridcell_fraction() + residual);
		if (residual != 0.0)
			add_weighted_pools(recovered, *last_imported, residual);
		imported_total += residual;
	}
	return imported_total;
}

double synthesize_missing_stand(Gridcell& gridcell, int stid, double requested,
		Pools& recovered) {
	if (requested <= 0.0)
		return 0.0;
	Stand* source = largest_stand(gridcell, stid);
	if (!source) {
		for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s) {
			if (gridcell[s].get_gridcell_fraction() > AREA_TOLERANCE &&
					gridcell[s].landcover == stlist[stid].landcover) {
				source = &gridcell[s];
				break;
			}
		}
	}
	if (!source) {
		for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s) {
			if (gridcell[s].get_gridcell_fraction() > AREA_TOLERANCE) {
				source = &gridcell[s];
				break;
			}
		}
	}
	if (!source)
		throw std::runtime_error("cannot synthesize a missing stand without any donor state");
	Stand& created = clone_to_type_without_pool_change(*source, stid, requested);
	add_weighted_pools(recovered, created, requested);
	return requested;
}

void recover_area(Gridcell& gridcell, Gridcell& donor1878, Gridcell& donor1854,
		int stid, double requested, MigrationAccounting& accounting) {
	if (requested <= AREA_TOLERANCE)
		return;
	double remaining = requested;
	double imported = import_historical_stands(
		gridcell, donor1878, stid, remaining, accounting.recovered);
	accounting.recovered_1878 += imported;
	remaining -= imported;
	if (remaining > AREA_TOLERANCE) {
		imported = import_historical_stands(
			gridcell, donor1854, stid, remaining, accounting.recovered);
		accounting.recovered_1854 += imported;
		remaining -= imported;
	}
	if (remaining > AREA_TOLERANCE) {
		const double synthesized = synthesize_missing_stand(
			gridcell, stid, remaining, accounting.recovered);
		accounting.recovered_synthetic += synthesized;
		remaining -= synthesized;
	}
	if (std::fabs(remaining) > AREA_TOLERANCE)
		throw std::runtime_error("failed to reconstruct all missing stand area");
}

void transfer_between_types(Gridcell& gridcell, int from, int to, double amount) {
	double remaining = amount;
	while (remaining > AREA_TOLERANCE) {
		Stand* source = largest_stand(gridcell, from);
		if (!source)
			throw std::runtime_error("stand-type transfer lacks physical donor area");
		const double moved = std::min(remaining, source->get_gridcell_fraction());
		Stand& created = clone_to_type_without_pool_change(*source, to, moved);
		(void)created;
		source->set_gridcell_fraction(source->get_gridcell_fraction() - moved);
		remaining -= moved;
	}
}

void remove_surplus(Gridcell& gridcell, int stid, double amount,
		MigrationAccounting& accounting) {
	double remaining = amount;
	while (remaining > AREA_TOLERANCE) {
		Stand* source = largest_stand(gridcell, stid);
		if (!source)
			throw std::runtime_error("cannot remove residual stand surplus");
		const double removed = std::min(remaining, source->get_gridcell_fraction());
		add_weighted_pools(accounting.removed, *source, removed);
		accounting.removed_area += removed;
		source->set_gridcell_fraction(source->get_gridcell_fraction() - removed);
		remaining -= removed;
	}
}

void move_to_target(Gridcell& gridcell, const std::vector<double>& target,
		MigrationAccounting* accounting) {
	std::vector<double> current = physical_st(gridcell);
	std::vector<double> surplus(target.size(), 0.0);
	std::vector<double> deficit(target.size(), 0.0);
	for (size_t st = 0; st < target.size(); ++st) {
		const double difference = current[st] - target[st];
		if (difference > AREA_TOLERANCE)
			surplus[st] = difference;
		else if (difference < -AREA_TOLERANCE)
			deficit[st] = -difference;
	}

	for (size_t from = 0; from < target.size(); ++from) {
		for (size_t to = 0; to < target.size() && surplus[from] > AREA_TOLERANCE; ++to) {
			if (deficit[to] <= AREA_TOLERANCE)
				continue;
			const double moved = std::min(surplus[from], deficit[to]);
			transfer_between_types(gridcell, static_cast<int>(from), static_cast<int>(to), moved);
			surplus[from] -= moved;
			deficit[to] -= moved;
		}
	}

	for (size_t st = 0; st < target.size(); ++st) {
		if (deficit[st] > AREA_TOLERANCE) {
			if (!accounting) {
				if (deficit[st] <= CONSERVATIVE_CLOSURE_TOLERANCE)
					continue;
				std::fprintf(stderr,
					"TARGET RESIDUAL deficit st=%zu area=%.17g current_total=%.17g "
					"target_total=%.17g\n", st, deficit[st],
					vector_sum(physical_st(gridcell)), vector_sum(target));
				throw std::runtime_error("conservative target transfer left an area deficit");
			}
			synthesize_missing_stand(gridcell, static_cast<int>(st),
				deficit[st], accounting->recovered);
			accounting->recovered_synthetic += deficit[st];
		}
		if (surplus[st] > AREA_TOLERANCE) {
			if (!accounting) {
				if (surplus[st] <= CONSERVATIVE_CLOSURE_TOLERANCE)
					continue;
				std::fprintf(stderr,
					"TARGET RESIDUAL surplus st=%zu area=%.17g current_total=%.17g "
					"target_total=%.17g\n", st, surplus[st],
					vector_sum(physical_st(gridcell)), vector_sum(target));
				throw std::runtime_error("conservative target transfer left an area surplus");
			}
			remove_surplus(gridcell, static_cast<int>(st), surplus[st], *accounting);
		}
	}
}

void reconstruct_restart_area(Gridcell& gridcell, Gridcell& donor1878,
		Gridcell& donor1854, const std::vector<double>& target,
		MigrationAccounting& accounting) {
	const std::vector<double> current = physical_st(gridcell);
	std::vector<double> deficit(target.size(), 0.0);
	double deficit_total = 0.0;
	for (size_t st = 0; st < target.size(); ++st) {
		deficit[st] = std::max(0.0, target[st] - current[st]);
		deficit_total += deficit[st];
	}
	const double missing_total = std::max(0.0, 1.0 - vector_sum(current));
	if (missing_total > AREA_TOLERANCE && deficit_total > 0.0) {
		double recovered_total = 0.0;
		int last_deficit = -1;
		for (size_t st = 0; st < target.size(); ++st)
			if (deficit[st] > AREA_TOLERANCE)
				last_deficit = static_cast<int>(st);
		for (size_t st = 0; st < target.size(); ++st) {
			if (deficit[st] <= AREA_TOLERANCE)
				continue;
			double share = missing_total * deficit[st] / deficit_total;
			if (static_cast<int>(st) == last_deficit)
				share = missing_total - recovered_total;
			share = std::min(share, deficit[st]);
			recover_area(gridcell, donor1878, donor1854,
				static_cast<int>(st), share, accounting);
			recovered_total += share;
		}
	}
	move_to_target(gridcell, target, &accounting);
}

void remove_zero_area_stands(Gridcell& gridcell) {
	Gridcell::iterator stand = gridcell.begin();
	while (stand != gridcell.end()) {
		if ((*stand).get_gridcell_fraction() <= AREA_TOLERANCE)
			stand = gridcell.delete_stand(stand);
		else
			++stand;
	}
}

void close_physical_peat(Gridcell& gridcell, double prescribed_peat) {
	Stand* closure = 0;
	long double other_peat = 0.0L;
	for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s) {
		Stand& stand = gridcell[s];
		if (stand.landcover == PEATLAND && (!closure ||
				stand.get_gridcell_fraction() > closure->get_gridcell_fraction()))
			closure = &stand;
	}
	if (!closure) {
		if (prescribed_peat > AREA_TOLERANCE)
			throw std::runtime_error("positive prescribed peat has no physical stand");
		return;
	}
	for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s)
		if (&gridcell[s] != closure && gridcell[s].landcover == PEATLAND)
			other_peat += static_cast<long double>(gridcell[s].get_gridcell_fraction());
	const double corrected =
		static_cast<double>(static_cast<long double>(prescribed_peat) - other_peat);
	if (corrected < -AREA_TOLERANCE)
		throw std::runtime_error("peat closure would make a peat stand negative");
	closure->set_gridcell_fraction(std::max(0.0, corrected));
}

void close_physical_total(Gridcell& gridcell, double prescribed_peat) {
	if (prescribed_peat >= 1.0 - AREA_TOLERANCE)
		return;
	Stand* closure = 0;
	long double others = 0.0L;
	for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s) {
		Stand& stand = gridcell[s];
		if (stand.landcover != PEATLAND && (!closure ||
				stand.get_gridcell_fraction() > closure->get_gridcell_fraction()))
			closure = &stand;
	}
	if (!closure)
		throw std::runtime_error("no non-peat physical stand available for closure");
	for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s)
		if (&gridcell[s] != closure)
			others += static_cast<long double>(gridcell[s].get_gridcell_fraction());
	const double corrected = static_cast<double>(1.0L - others);
	if (corrected < -AREA_TOLERANCE)
		throw std::runtime_error("physical closure would make a stand negative");
	closure->set_gridcell_fraction(std::max(0.0, corrected));
}

void synchronize_metadata(Gridcell& gridcell, double prescribed_peat) {
	double lc[NLANDCOVERTYPES] = {0.0};
	std::vector<double> st(gridcell.st.nobj, 0.0);
	std::vector<int> counts(gridcell.st.nobj, 0);
	for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s) {
		Stand& stand = gridcell[s];
		const double area = stand.get_gridcell_fraction();
		st[stand.stid] += area;
		lc[stand.landcover] += area;
		counts[stand.stid]++;
		stand.frac_old = area;
		stand.frac_temp = area;
		stand.frac_change = 0.0;
		stand.gross_frac_increase = 0.0;
		stand.gross_frac_decrease = 0.0;
		stand.cloned_fraction = 0.0;
		stand.cloned = false;
	}
	for (unsigned int i = 0; i < gridcell.st.nobj; ++i) {
		Gridcellst& gcst = gridcell.st[i];
		const double represented = stlist[i].landcover == PEATLAND ?
			prescribed_peat : st[i];
		gcst.frac = represented;
		gcst.frac_old = represented;
		gcst.frac_old_orig = represented;
		gcst.frac_change = 0.0;
		gcst.gross_frac_increase = 0.0;
		gcst.gross_frac_decrease = 0.0;
		gcst.nstands = counts[i];
	}
	for (int i = 0; i < NLANDCOVERTYPES; ++i) {
		const double represented = i == PEATLAND ? prescribed_peat : lc[i];
		gridcell.landcover.frac[i] = represented;
		gridcell.landcover.frac_old[i] = represented;
		gridcell.landcover.frac_change[i] = 0.0;
	}
}

void validate_final(Gridcell& gridcell,
		const double target_lc[NLANDCOVERTYPES]) {
	const std::vector<double> st = physical_st(gridcell);
	double lc[NLANDCOVERTYPES] = {0.0};
	double total = 0.0;
	for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s) {
		Stand& stand = gridcell[s];
		lc[stand.landcover] += stand.get_gridcell_fraction();
		total += stand.get_gridcell_fraction();
	}
	if (std::fabs(total - 1.0) > AREA_TOLERANCE)
		throw std::runtime_error("final physical stand fractions do not sum to one");
	if (std::fabs(lc[PEATLAND] - target_lc[PEATLAND]) > AREA_TOLERANCE) {
		std::fprintf(stderr,
			"PEAT DETAIL map=%.17g physical=%.17g delta=%.17g total=%.17g "
			"stands=%u\n", target_lc[PEATLAND], lc[PEATLAND],
			lc[PEATLAND] - target_lc[PEATLAND], total, gridcell.nbr_stands());
		throw std::runtime_error("final physical peat differs from the map");
	}
	for (int i = 0; i < NLANDCOVERTYPES; ++i) {
		if (std::fabs(lc[i] - target_lc[i]) > AREA_TOLERANCE) {
			std::fprintf(stderr,
				"LANDCOVER DETAIL lc=%d target=%.17g physical=%.17g delta=%.17g "
				"total=%.17g stands=%u lon=%.17g lat=%.17g\n",
				i, target_lc[i], lc[i], lc[i] - target_lc[i], total,
				gridcell.nbr_stands(), gridcell.get_lon(), gridcell.get_lat());
			throw std::runtime_error("final physical land cover differs from prescribed target");
		}
	}
	for (unsigned int i = 0; i < gridcell.st.nobj; ++i)
		if (std::fabs(gridcell.st[i].frac - st[i]) > AREA_TOLERANCE)
			throw std::runtime_error("final stand-type metadata differs from physical stands");
	for (int i = 0; i < NLANDCOVERTYPES; ++i)
		if (std::fabs(gridcell.landcover.frac[i] - lc[i]) > AREA_TOLERANCE)
			throw std::runtime_error("final land-cover metadata differs from physical stands");
}

bool pools_conserved(const Pools& before, const Pools& after) {
	const double ctol = 1.0e-9 * std::max(1.0, std::fabs(before.c));
	const double ntol = 1.0e-9 * std::max(1.0, std::fabs(before.n));
	const double wtol = 1.0e-9 * std::max(1.0, std::fabs(before.water));
	return std::fabs(after.c - before.c) <= ctol &&
		std::fabs(after.n - before.n) <= ntol &&
		std::fabs(after.water - before.water) <= wtol;
}

} // namespace

int main(int argc, char** argv) {
	if (argc != 11) {
		std::fprintf(stderr,
			"usage: %s INS SOURCE_STATE DONOR1_STATE DONOR2_STATE LUH3_PEAT_TARGET SOIL_CODES OUTPUT_DIR "
			"RANK RANKS REPORT\n", argv[0]);
		return 2;
	}
	const int rank = std::atoi(argv[8]);
	const int ranks = std::atoi(argv[9]);
	if (rank < 0 || ranks <= 0 || rank >= ranks)
		return 2;

	try {
		set_shell(new CommandLineShell(file_log));
		std::auto_ptr<InputModule> input_module(
			InputModuleRegistry::get_instance().create_input_module("ece"));
		GuessOutput::OutputModuleContainer output_modules;
		GuessOutput::OutputModuleRegistry::get_instance().create_all_modules(output_modules);
		read_instruction_file(argv[1]);
		if (iforganicsoilproperties)
			throw std::runtime_error(
				"offline migration does not support iforganicsoilproperties=1");

		LandcoverTargetMap target_map(argv[5]);
		SoilCodeMap soil_codes(argv[6]);
		GuessDeserializer source_state(argv[2]);
		GuessDeserializer donor1_state(argv[3]);
		GuessDeserializer donor2_state(argv[4]);
		const std::vector<std::pair<double, double> > coords =
			source_state.get_source_coordinates();
		if (donor1_state.get_source_coordinates() != coords ||
				donor2_state.get_source_coordinates() != coords)
			throw std::runtime_error("restart coordinate sets differ");

		FILE* report = std::fopen(argv[10], "w");
		if (!report)
			throw std::runtime_error("cannot create migration report");
		std::fprintf(report,
			"index\tlon\tlat\tluh3_valid\ttarget_urban\ttarget_crop\ttarget_pasture\t"
			"target_forest\ttarget_natural\ttarget_peat\ttarget_barren\t"
			"old_target_peat\tbefore_area\tbefore_peat\t"
			"before_c\tbefore_n\tbefore_water\tbaseline_c\tbaseline_n\tbaseline_water\t"
			"final_area\tfinal_urban\tfinal_crop\tfinal_pasture\tfinal_forest\t"
			"final_natural\tfinal_peat\tfinal_barren\tfinal_c\tfinal_n\tfinal_water\t"
			"recovered_donor1_area\trecovered_donor2_area\trecovered_synthetic_area\t"
			"recovered_c\trecovered_n\trecovered_water\tremoved_area\tremoved_c\t"
			"removed_n\tremoved_water\n");

		// A focused migration smoke test may select one source cell. Full
		// production runs leave this variable unset and process every cell.
		const char* test_cell_text = std::getenv("LPJG_TEST_CELL_INDEX");
		const size_t test_cell = test_cell_text ?
			static_cast<size_t>(std::strtoul(test_cell_text, 0, 10)) : coords.size();
		if (test_cell_text && test_cell >= coords.size())
			throw std::runtime_error("LPJG_TEST_CELL_INDEX is outside the source grid");
		GuessSerializer serializer(argv[7], rank, ranks);
		for (size_t index = static_cast<size_t>(rank); index < coords.size();
				index += static_cast<size_t>(ranks)) {
			if (test_cell_text && index != test_cell) continue;
			Gridcell gridcell;
			Gridcell old1878;
			Gridcell old1854;
			gridcell.set_coordinates(coords[index].first, coords[index].second);
			old1878.set_coordinates(coords[index].first, coords[index].second);
			old1854.set_coordinates(coords[index].first, coords[index].second);
			source_state.deserialize_gridcell(gridcell);
			donor1_state.deserialize_gridcell(old1878);
			donor2_state.deserialize_gridcell(old1854);
			const int soil_code = soil_codes.lookup(
				coords[index].first, coords[index].second);
			if (!ifs_soil_parameters(gridcell, soil_code) ||
					!ifs_soil_parameters(old1878, soil_code) ||
					!ifs_soil_parameters(old1854, soil_code))
				throw std::runtime_error("IFS soil-code initialization failed");
			normalize_gridcell_hydrology(gridcell);
			normalize_gridcell_hydrology(old1878);
			normalize_gridcell_hydrology(old1854);

			const TargetRow& target_row =
				target_map.lookup(coords[index].first, coords[index].second);
			const double map_peat = target_row.peat;
			const std::vector<double> old_target = restart_target(gridcell);
			double old_target_peat = 0.0;
			for (size_t st = 0; st < old_target.size(); ++st)
				if (stlist[st].landcover == PEATLAND)
					old_target_peat += old_target[st];

			const std::vector<double> before_st = physical_st(gridcell);
			const double before_area = vector_sum(before_st);
			double before_peat = 0.0;
			for (size_t st = 0; st < before_st.size(); ++st)
				if (stlist[st].landcover == PEATLAND)
					before_peat += before_st[st];
			const Pools before = gridcell_pools(gridcell);

			MigrationAccounting accounting;
			reconstruct_restart_area(gridcell, old1878, old1854,
				old_target, accounting);
			remove_zero_area_stands(gridcell);
			normalize_gridcell_hydrology(gridcell);
			const Pools baseline = gridcell_pools(gridcell);

			double prescribed_lc[NLANDCOVERTYPES] = {0.0};
			const std::vector<double> final_target =
				prescribed_landcover_target(old_target, target_row, prescribed_lc);
			move_to_target(gridcell, final_target, 0);
			remove_zero_area_stands(gridcell);
			normalize_gridcell_hydrology(gridcell);
			close_physical_peat(gridcell, map_peat);
			close_physical_total(gridcell, map_peat);
			// The source restart can be short or long by a few ulps. The first
			// conservative transfer then leaves a sub-1e-10 class deficit, and
			// total closure necessarily places it in one particular non-peat
			// class. Re-run the same physical stand transfer after total closure
			// so that the canonical total and every prescribed class agree. This
			// moves state between stands and therefore conserves C/N/water; it is
			// not a metadata-only repair or a relaxed comparison.
			move_to_target(gridcell, final_target, 0);
			remove_zero_area_stands(gridcell);
			normalize_gridcell_hydrology(gridcell);
			close_physical_peat(gridcell, map_peat);
			close_physical_total(gridcell, map_peat);
			synchronize_metadata(gridcell, map_peat);
			validate_final(gridcell, prescribed_lc);
			validate_fixed_peat_restart_soil_state(gridcell,
				"offline migration final state");
			const Pools final_pools = gridcell_pools(gridcell);
			if (!pools_conserved(baseline, final_pools)) {
				std::fprintf(stderr,
					"POOL DETAIL index=%zu baseline=(%.17g %.17g %.17g) "
					"final=(%.17g %.17g %.17g) delta=(%.17g %.17g %.17g)\n",
					index, baseline.c, baseline.n, baseline.water,
					final_pools.c, final_pools.n, final_pools.water,
					final_pools.c - baseline.c, final_pools.n - baseline.n,
					final_pools.water - baseline.water);
				throw std::runtime_error("fixed-peat reclassification changed C/N/water pools");
			}

			double final_lc[NLANDCOVERTYPES] = {0.0};
			double final_area = 0.0;
			for (unsigned int s = 0; s < gridcell.nbr_stands(); ++s) {
				final_lc[gridcell[s].landcover] += gridcell[s].get_gridcell_fraction();
				final_area += gridcell[s].get_gridcell_fraction();
			}

			serializer.serialize_gridcell(gridcell);
			std::fprintf(report,
				"%zu\t%.17g\t%.17g\t%d\t"
				"%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t"
				"%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t"
				"%.17g\t%.17g\t%.17g\t%.17g\t"
				"%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t"
				"%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t"
				"%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
				index, coords[index].first, coords[index].second,
				target_row.luh3_valid ? 1 : 0,
				prescribed_lc[URBAN], prescribed_lc[CROPLAND],
				prescribed_lc[PASTURE], prescribed_lc[FOREST],
				prescribed_lc[NATURAL], prescribed_lc[PEATLAND],
				prescribed_lc[BARREN], old_target_peat, before_area, before_peat,
				before.c, before.n, before.water,
				baseline.c, baseline.n, baseline.water, final_area,
				final_lc[URBAN], final_lc[CROPLAND], final_lc[PASTURE],
				final_lc[FOREST], final_lc[NATURAL], final_lc[PEATLAND],
				final_lc[BARREN],
				final_pools.c, final_pools.n, final_pools.water,
				accounting.recovered_1878, accounting.recovered_1854,
				accounting.recovered_synthetic, accounting.recovered.c,
				accounting.recovered.n, accounting.recovered.water,
				accounting.removed_area, accounting.removed.c,
				accounting.removed.n, accounting.removed.water);
		}
		std::fclose(report);
	}
	catch (const std::exception& error) {
		std::fprintf(stderr, "MIGRATION ERROR rank=%d: %s\n", rank, error.what());
		return 1;
	}
	return 0;
}
