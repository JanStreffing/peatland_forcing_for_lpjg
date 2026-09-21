#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 INS STATE_DIR OUTPUT_DIR" >&2
    exit 2
fi
ins=$1
state=$2
output=$3
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
probe=$script_dir/bin/restart_fraction_probe
readback_shards=${LPJG_READBACK_SHARDS:-16}
parallel_jobs=${LPJG_PARALLEL_JOBS:-16}
for path in "$ins" "$state/meta.bin" "$probe"; do
    [[ -s $path ]] || { echo "Missing or empty required input: $path" >&2; exit 2; }
done
if ! [[ $readback_shards =~ ^[1-9][0-9]*$ && $parallel_jobs =~ ^[1-9][0-9]*$ ]]; then
    echo "LPJG_READBACK_SHARDS and LPJG_PARALLEL_JOBS must be positive integers" >&2
    exit 2
fi
[[ ! -e $output ]] || { echo "Refusing existing output: $output" >&2; exit 2; }
mkdir -p "$output"
for ((first = 0; first < readback_shards; first += parallel_jobs)); do
    last=$((first + parallel_jobs - 1))
    if (( last >= readback_shards )); then last=$((readback_shards - 1)); fi
    pids=()
    for shard in $(seq "$first" "$last"); do
        "$probe" "$ins" "$state" "$shard" "$readback_shards" "$output/shard_${shard}.tsv" \
            >"$output/shard_${shard}.log" 2>&1 &
        pids+=("$!")
    done
    failed=0
    for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
    if (( failed )); then
        echo "Fraction read-back failed in shard batch $first-$last" >&2
        exit 1
    fi
done
count=$(find "$output" -maxdepth 1 -type f -name 'shard_*.tsv' | wc -l)
[[ $count -eq $readback_shards ]] || {
    echo "Incomplete fraction read-back: expected $readback_shards tables, found $count" >&2
    exit 1
}
echo "Fraction read-back complete: $output"
