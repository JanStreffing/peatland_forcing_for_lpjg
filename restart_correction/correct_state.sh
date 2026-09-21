#!/usr/bin/env bash
# Bring an LPJ-GUESS restart in line with a peat map and the LUH3 crop mix.
#
#   correct_state.sh SETUP_DIR STATE DONOR1 DONOR2 TARGET_TABLE PEAT_MAP YEAR OUTPUT_BASE
#
# SETUP_DIR     instructions/guess.ins and *_restart_exact_soil_types.txt for the
#               grid and configuration the state was written with
# STATE         the lpjg_state_<year> directory to correct. It is only read.
# DONOR1/2      older states of the same run, read only where a stand is missing
# TARGET_TABLE  from build_luh3_peat_target.py, for the land-use year the run uses:
#               its fixed_LU year, or the last completed year of a transient run
# PEAT_MAP      the peat map TARGET_TABLE was built with. The startup probe reads it
#               through guess.ins, as the model will.
# YEAR          the year of STATE, for the startup probe
#
# Runs on a login node in a few minutes for TCO95. The result is usable only if the
# last line says PASS and the startup probe differences are far below 1e-6.
set -euo pipefail
[[ $# -eq 8 ]] || { sed -n 2,18p "$0" >&2; exit 2; }
setup=$(realpath "$1"); state=$(realpath "$2"); donor1=$(realpath "$3"); donor2=$(realpath "$4")
target=$(realpath "$5"); peat_map=$(realpath "$6"); year=$7; out=$8
[[ -e $out ]] && { echo "Refusing to use an existing output base: $out" >&2; exit 2; }
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ins=$setup/instructions/guess.ins
soil=$(ls "$setup"/*_restart_exact_soil_types.txt)
cells=$(($(wc -l < "$soil") - 1))
export LPJG_STATE_SHARDS=$(ls "$state"/*.state | wc -l)
mkdir -p "$out"; out=$(realpath "$out")
new=$out/$(basename "$state")_corrected
# guess.ins names its inputs relative to the run directory.
mkdir -p "$setup/landuse"
ln -sfn "$peat_map" "$setup/landuse/peat_frac.txt"
cd "$setup/instructions"

bash "$here/RUN_ONE_RESTART_MIGRATION.sh" "$ins" "$state" "$donor1" "$donor2" "$target" "$soil" "$new" "$out/reports/migration"
bash "$here/RUN_FRACTION_READBACK.sh" "$ins" "$new" "$out/reports/fraction_readback"
bash "$here/RUN_STATE_READBACK.sh" "$ins" "$new" "$soil" "$out/reports/state_readback"

# What a fixed land-use restart checks at load, run offline with LPJ-GUESS's own input code.
mkdir -p "$out/reports/startup_probe"
for first in $(seq 0 16 $((LPJG_STATE_SHARDS - 1))); do
    for rank in $(seq "$first" $((first + 15 < LPJG_STATE_SHARDS ? first + 15 : LPJG_STATE_SHARDS - 1))); do
        "$here/bin/fixed_lu_restart_startup_probe" "$ins" "$new" "$year" "$rank" "$LPJG_STATE_SHARDS" \
            "$out/reports/startup_probe/s_$rank.tsv" > "$out/reports/startup_probe/log_$rank.txt" 2>&1 &
    done
    wait
done
awk -F'\t' 'FNR>1{c+=$3; if($4>m)m=$4; if($5>l)l=$5}
    END{printf "startup probe: %d cells, largest stand-type difference %.3g, land-cover %.3g (LPJ-GUESS allows 1e-6)\n", c, m, l
        if (m > 1e-9 || l > 1e-9) exit 1}' "$out"/reports/startup_probe/s_*.tsv

"${PYTHON:-python3}" "$here/validate_corrected_restart.py" --migration "$out/reports/migration" \
    --fractions "$out/reports/fraction_readback" --state "$out/reports/state_readback" \
    --target "$target" --experiment "$(basename "$state")" --source-state "$state" \
    --expected-cells "$cells" --output "$out/reports/validation.json"
echo "corrected state: $new"
