#!/usr/bin/env python3
"""Command-line entry point for complete peat-map generation."""

from __future__ import annotations

import argparse
import sys

from peat_map_tool import (
    DEFAULT_GRIDLIST_TOLERANCE,
    PeatMapError,
    format_report,
    generate_map,
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a complete resolution-wide LPJ-GUESS peat map with "
            "periodic nearest-neighbour interpolation."
        )
    )
    parser.add_argument("--source", required=True, help="source Lon/Lat/Peat_Frac table")
    parser.add_argument(
        "--target-grid", required=True, help="NetCDF containing the complete target grid"
    )
    parser.add_argument("--resolution", required=True, help="resolution identifier, e.g. TCO319")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="staging directory; output is <resolution>_peat_frac.txt",
    )
    parser.add_argument(
        "--expected-count",
        type=int,
        help="required number of complete target-grid rows",
    )
    parser.add_argument(
        "--gridlist",
        action="append",
        default=[],
        help="optional Lat/Lon gridlist that the complete output must cover (repeatable)",
    )
    parser.add_argument(
        "--gridlist-tolerance",
        type=float,
        default=DEFAULT_GRIDLIST_TOLERANCE,
        help=(
            "maximum absolute longitude and latitude difference in degrees for "
            "gridlist coverage (default 5e-5)"
        ),
    )
    parser.add_argument("--lon-var", default="lon", help="target longitude variable name")
    parser.add_argument("--lat-var", default="lat", help="target latitude variable name")
    parser.add_argument(
        "--coordinate-precision",
        type=int,
        default=6,
        help="decimal places for coordinates (minimum 6; default 6)",
    )
    parser.add_argument(
        "--fraction-precision",
        type=int,
        default=8,
        help="decimal places for peat fraction (minimum 8; default 8)",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=200_000,
        help="number of target coordinates queried per interpolation chunk",
    )
    parser.add_argument(
        "--replace",
        action="store_true",
        help="atomically replace an existing staged output",
    )
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    try:
        output, report = generate_map(
            source_path=arguments.source,
            target_grid_path=arguments.target_grid,
            resolution=arguments.resolution,
            output_directory=arguments.output_dir,
            expected_count=arguments.expected_count,
            gridlists=arguments.gridlist,
            coordinate_precision=arguments.coordinate_precision,
            fraction_precision=arguments.fraction_precision,
            lon_variable=arguments.lon_var,
            lat_variable=arguments.lat_var,
            chunk_size=arguments.chunk_size,
            gridlist_tolerance=arguments.gridlist_tolerance,
            replace=arguments.replace,
        )
    except (PeatMapError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(format_report(output, report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
