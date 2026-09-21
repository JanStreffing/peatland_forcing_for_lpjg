#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 LPJ_GUESS_SOURCE_DIRECTORY" >&2
    exit 2
fi

source_dir=$1
build_dir=$source_dir/build
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
flags=$build_dir/CMakeFiles/guess.dir/flags.make
link_file=$build_dir/CMakeFiles/guess.dir/link.txt

if [[ ! -f $flags || ! -f $link_file ]]; then
    echo "Build LPJ-GUESS in $build_dir before building the migration tool." >&2
    exit 1
fi

# Read the three CMake make-variable values as data. They contain compiler
# flags, not executable shell code.
CXX_DEFINES=$(sed -n 's/^CXX_DEFINES = //p' "$flags")
CXX_INCLUDES=$(sed -n 's/^CXX_INCLUDES = //p' "$flags")
CXX_FLAGS=$(sed -n 's/^CXX_FLAGS = //p' "$flags")
compiler=$(which mpicxx)
read -r -a definitions <<< "$CXX_DEFINES"
read -r -a includes <<< "$CXX_INCLUDES"
read -r -a compile_flags <<< "$CXX_FLAGS"
main_object=CMakeFiles/guess.dir/command_line_version/main.cpp.o
mkdir -p "$script_dir/bin"

build_tool() {
    local name=$1
    local source=$script_dir/$name.cpp
    local output=$script_dir/bin/$name
    local object=$script_dir/bin/$name.o
    "$compiler" "${definitions[@]}" "${includes[@]}" "${compile_flags[@]}" \
        -c "$source" -o "$object"

    local link_command
    link_command=$(<"$link_file")
    if [[ $link_command != *"$main_object"* ]]; then
        echo "Could not find the LPJ-GUESS main object in $link_file" >&2
        exit 1
    fi
    link_command=${link_command/"$main_object"/"$object"}
    link_command=${link_command/' -o guess '/' -o '"$output"' '}
    (cd "$build_dir" && eval "$link_command")
    rm -f "$object"
    echo "Built $output"
}

build_tool synchronize_restart_to_luh3_peat
build_tool restart_state_probe
build_tool restart_fraction_probe
build_tool fixed_lu_restart_startup_probe

branch=$(git -C "$source_dir" branch --show-current 2>/dev/null || echo unknown)
commit=$(git -C "$source_dir" rev-parse HEAD 2>/dev/null || echo unknown)
status=$(git -C "$source_dir" status --porcelain 2>/dev/null || true)
{
    echo "source_directory=$(cd -- "$source_dir" && pwd)"
    echo "branch=$branch"
    echo "commit=$commit"
    echo "built_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    if [[ -n $status ]]; then
        echo "working_tree=modified"
        printf '%s\n' "$status"
    else
        echo "working_tree=clean"
    fi
} >"$script_dir/bin/build_provenance.txt"
echo "Wrote $script_dir/bin/build_provenance.txt"
