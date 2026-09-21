#!/usr/bin/env python3
"""Command-line entry point for complete peat-map validation."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from peat_map_tool import (
    DEFAULT_GRIDLIST_TOLERANCE,
    PeatMapError,
    format_report,
    read_target_grid,
    validate_map,
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate an LPJ-GUESS peat map against the complete target NetCDF grid."
        )
    )
    parser.add_argument("--map", required=True, help="Lon/Lat/Peat_Frac map to validate")
    parser.add_argument(
        "--target-grid", required=True, help="NetCDF containing the complete target grid"
    )
    parser.add_argument(
        "--resolution",
        required=True,
        help="resolution identifier; map basename must be <resolution>_peat_frac.txt",
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
        help="optional Lat/Lon gridlist that the map must cover (repeatable)",
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
        help="coordinate comparison precision (minimum 6; default 6)",
    )
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    map_path = Path(arguments.map)
    expected_name = f"{arguments.resolution}_peat_frac.txt"
    if map_path.name != expected_name:
        print(
            f"ERROR: map basename must be {expected_name!r}, got {map_path.name!r}",
            file=sys.stderr,
        )
        return 2
    try:
        target_lon, target_lat = read_target_grid(
            arguments.target_grid,
            lon_variable=arguments.lon_var,
            lat_variable=arguments.lat_var,
        )
        report = validate_map(
            map_path,
            target_lon,
            target_lat,
            arguments.coordinate_precision,
            expected_count=arguments.expected_count,
            gridlists=arguments.gridlist,
            gridlist_tolerance=arguments.gridlist_tolerance,
        )
    except (PeatMapError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(format_report(map_path, report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
