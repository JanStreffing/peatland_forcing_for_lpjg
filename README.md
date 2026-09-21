# Peat forcing for LPJ-GUESS in AWI-ESM3

Prescribed peat fraction per OpenIFS grid cell, and the tools to make the maps and to
bring an existing LPJ-GUESS restart in line with a new one.

The maps themselves are not in git. They are reproducible from `TL255_peat_frac.txt`
and the LUH3 states in minutes, and two of them are above 100 MB.

## The maps

| file | what |
|---|---|
| `TL255_peat_frac.txt` | the source, at full coordinate precision. Tracked. |
| `<RES>_peat_frac.txt` | the source moved onto a grid by nearest neighbour. Uncapped. |
| `<RES>_peat_frac_Lund_corrected_0850-2024.txt` | the same, capped. AWI-ESM3 v3.5 uses these. |

Resolutions: TL159, TL255, TL511, TCO95, TCO159, TCO199, TCO319, TCO639, TCO1279, TCO2559.

### Why the cap

LPJ-GUESS takes the fixed peat area out of natural land first. Where a map holds more
peat than LUH3 has natural land, natural is zero and the rest is squeezed out of crops
and pasture. LUH3 still prescribes its full gross transfers out of natural, which then no
longer fit. The capped map is

    natural = min(primf + primn + secdf + secdn, 1 - crops - pasture - urban)
    peat    = min(peat, smallest natural over 0850-2024)

Cells LUH3 does not cover get zero peat. The rule is a270270's
(`/work/ab0995/a270270/peat_LU_fractions_EC_EARTH_METHOD`), who made the TCO95 and
TCO319 maps. `scripts/correct_peat_map_luh3.py` reproduces both byte for byte and also
handles the grids their script cannot read.

The cap holds for the historical period only. A scenario needs the minimum taken over
historical plus scenario states.

### Making a map

    # 1. a target grid at full precision. LPJ-GUESS finds a cell in the map to within
    #    5e-5 degrees, and the regridded LUH3 files carry coordinates rounded to four
    #    decimals, or for TL511 and TCO639 a latitude array that was never written.
    scripts/make_gaussian_target_grid.py <luh3_file_on_that_grid.nc> grid.nc

    # 2. the source onto that grid (a270270's tool)
    scripts/generate_peat_map.py --source TL255_peat_frac.txt --target-grid grid.nc \
        --resolution <RES> --output-dir stage --expected-count <cells>

    # 3. the cap
    scripts/correct_peat_map_luh3.py stage/<RES>_peat_frac.txt \
        <RES>_peat_frac_Lund_corrected_0850-2024.txt <luh3_states.nc> [...]

Step 3 takes a regridded LUH3 states file or the native 0.25 degree one. The regrids
were made with `cdo remapnn`, and a minimum over time commutes with nearest neighbour, so
the native file gives the same result apart from a few coastal cells, which get zero.
That is how TCO1279 and TCO2559 were made; TCO2559 has no states regrid. Given several
LUH3 files, the cap is the smallest of them: those two maps fit v3-1-1 and v3-1-2.

TCO95, TCO319 and TL255 are on the base maps that were already in use, since restarts
exist against them. They differ from a fresh step 2 in 38 of 40320 cells at TCO95, all
exactly midway between two source points.

## Bringing a restart in line: `restart_correction/`

LPJ-GUESS compares the peat map with the peat stands in a restart and stops on any
difference. With fixed land use it also requires every stand type to match the forcing
to 1e-6, because it has no January transfer to close a gap with. So a state written under
another map, or another crop split, has to be corrected offline.

The tools are a270270's, from `corrected_state_files/script/` of the package above, where
they corrected the four TCO319 CMIP7 production states. Changed here: two TCO319 cell
counts removed, and one generic driver in place of the experiment-specific ones.

    # once, against the LPJ-GUESS build that will run the state (esm_master build tree);
    # needs the compiler environment of comp-lpj_guess-*_script.sh
    restart_correction/BUILD_TOOLS.sh <model_dir>/lpj_guess

    # the target: LUH3 for the land-use year the run uses, crop stand types, peat
    restart_correction/build_luh3_peat_target.py --luh3 <states.nc> --management <mgmt.nc> \
        --peat <capped map> --year <year> --output target.tsv --summary target.json

    restart_correction/correct_state.sh SETUP_DIR STATE DONOR1 DONOR2 target.tsv \
        <capped map> <state year> OUTPUT_BASE

`<year>` is the `fixed_LU` year, or for a transient run the last completed year, one
less than the state's name. The result is usable only if the driver ends with PASS.
About three minutes on a login node for TCO95.

`TCO95_PICAL/` is the setup for the TCO95 CORE3 pre-industrial runs: instruction files,
and soil codes keyed on the coordinates as the restart stores them. The state probe
matches those exactly, so a table written from the grid does not work.

First use: `PICAL_momixoff` 1940 to LUH3 1850 and the capped map. 10157 cells, 821 with
changed peat, 11 of them all peat before with no crop or pasture stand at all. Validator
PASS, C, N and water conserved to 4e-11, 5e-13 and 2e-10. One coupled year from it ran
clean, and the state the model wrote at the end passes the startup probe again.

This needs LPJ-GUESS with commit "Use LUH3 crop fractions for fixed land use restarts".
Older code splits cropland equally over its ten stand types, five of them irrigated,
which is 50% irrigated cropland against 5.5% in LUH3 1850.

## `superseded_2026-01/`

The January 2026 maps, with two-decimal coordinates and, for TL511 and TCO639, scrambled
latitudes. See the note inside.
