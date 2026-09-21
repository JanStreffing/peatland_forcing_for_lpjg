#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 8 ]]; then
    echo "Usage: $0 INS SOURCE DONOR1 DONOR2 TARGET_TABLE SOIL_CODES OUTPUT_STATE REPORT_DIR" >&2
    exit 2
fi

ins=$1
source_state=$2
donor1=$3
donor2=$4
target=$5
soil_codes=$6
output_state=$7
reports=$8
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
tool=$script_dir/bin/synchronize_restart_to_luh3_peat
state_shards=${LPJG_STATE_SHARDS:-64}
parallel_jobs=${LPJG_PARALLEL_JOBS:-16}

if ! [[ $state_shards =~ ^[1-9][0-9]*$ && $parallel_jobs =~ ^[1-9][0-9]*$ ]]; then
    echo "LPJG_STATE_SHARDS and LPJG_PARALLEL_JOBS must be positive integers" >&2
    exit 2
fi

for path in "$ins" "$source_state/meta.bin" "$donor1/meta.bin" "$donor2/meta.bin" \
            "$target" "$soil_codes" "$tool"; do
    [[ -s $path ]] || { echo "Missing or empty required input: $path" >&2; exit 2; }
done

for state_dir in "$source_state" "$donor1" "$donor2"; do
    count=$(find "$state_dir" -maxdepth 1 -type f -name '*.state' | wc -l)
    if [[ $count -ne $state_shards ]]; then
        echo "Expected $state_shards state shards in $state_dir; found $count" >&2
        exit 2
    fi
done
if [[ -e $output_state || -e $reports ]]; then
    echo "Refusing existing output: $output_state or $reports" >&2
    exit 2
fi
mkdir -p "$output_state" "$reports"

# Sixteen simultaneous readers/writers is the tested default for one Levante
# node. LPJG_PARALLEL_JOBS can lower the concurrency without changing output.
for ((first_rank = 0; first_rank < state_shards; first_rank += parallel_jobs)); do
    last_rank=$((first_rank + parallel_jobs - 1))
    if (( last_rank >= state_shards )); then
        last_rank=$((state_shards - 1))
    fi
    pids=()
    for rank in $(seq "$first_rank" "$last_rank"); do
        "$tool" "$ins" "$source_state" "$donor1" "$donor2" "$target" \
            "$soil_codes" "$output_state" "$rank" "$state_shards" \
            "$reports/rank_${rank}.tsv" \
            >"$reports/rank_${rank}.log" 2>&1 &
        pids+=("$!")
    done
    failed=0
    for pid in "${pids[@]}"; do
        wait "$pid" || failed=1
    done
    if (( failed )); then
        echo "Migration failed in rank batch $first_rank-$last_rank" >&2
        exit 1
    fi
done

state_count=$(find "$output_state" -maxdepth 1 -type f -name '*.state' | wc -l)
report_count=$(find "$reports" -maxdepth 1 -type f -name 'rank_*.tsv' | wc -l)
[[ $state_count -eq $state_shards && $report_count -eq $state_shards && -s $output_state/meta.bin ]] || {
    echo "Incomplete migration: states=$state_count reports=$report_count" >&2
    exit 1
}
echo "Migration complete: $output_state"
