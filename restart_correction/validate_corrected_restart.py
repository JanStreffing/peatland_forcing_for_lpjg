#!/usr/bin/env python3
"""Independent validation of an offline-corrected LPJ-GUESS restart."""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
from collections import defaultdict


LAND_CLASSES = ("urban", "crop", "pasture", "forest", "natural", "peat", "barren")
CROP_TYPES = ("CC3ann", "CC3per", "CC3nfx", "CC4ann", "CC4per",
              "CC3anni", "CC3peri", "CC3nfxi", "CC4anni", "CC4peri")
AREA_TOL = 1.0e-10
POOL_REL_TOL = 1.0e-9
MAP_TOL = 1.0e-3


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_sharded_tables(directory: str, prefix: str) -> dict[int, dict[str, str]]:
    paths = sorted(glob.glob(os.path.join(directory, f"{prefix}_*.tsv")))
    require(paths, f"no {prefix}_*.tsv files in {directory}")
    rows: dict[int, dict[str, str]] = {}
    expected_header = None
    for path in paths:
        with open(path, newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            require(reader.fieldnames is not None, f"missing header in {path}")
            if expected_header is None:
                expected_header = reader.fieldnames
            require(reader.fieldnames == expected_header, f"inconsistent header in {path}")
            for row in reader:
                index = int(row["index"])
                require(index not in rows, f"duplicate grid-cell index {index} in {directory}")
                rows[index] = row
    return rows


def index_rows_by_coordinate(
    rows: dict[int, dict[str, str]], name: str
) -> dict[tuple[float, float], dict[str, str]]:
    """Index rows by the serialized grid-cell coordinates.

    Restart serialization may repartition/reorder cells across the 64 output
    shards.  Therefore, the integer report index is meaningful only within one
    report set and must not be used to join pre-write and post-write tables.
    """
    indexed: dict[tuple[float, float], dict[str, str]] = {}
    for row in rows.values():
        key = (float(row["lon"]), float(row["lat"]))
        require(key not in indexed, f"duplicate coordinate {key} in {name}")
        indexed[key] = row
    return indexed


def periodic_longitude(lon: float) -> float:
    return (lon + 180.0) % 360.0 - 180.0


def longitude_distance(a: float, b: float) -> float:
    difference = abs(periodic_longitude(a) - periodic_longitude(b))
    return min(difference, 360.0 - difference)


def bucket_key(lon: float, lat: float) -> tuple[int, int]:
    return (math.floor(periodic_longitude(lon) * 1000.0), math.floor(lat * 1000.0))


def read_target(path: str) -> tuple[list[dict[str, str]], dict[tuple[int, int], list[int]]]:
    rows: list[dict[str, str]] = []
    buckets: dict[tuple[int, int], list[int]] = defaultdict(list)
    expected = [
        "Lon", "Lat", "LUH3Valid", "RawNatural", "RawCropland",
        "RawPasture", "RawUrban", "CorrectedPeat",
    ] + [f"Crop_{name}" for name in CROP_TYPES]
    with open(path, encoding="utf-8") as stream:
        require(stream.readline().split() == expected, f"unexpected target header in {path}")
        for line_number, line in enumerate(stream, start=2):
            fields = line.split()
            require(len(fields) == len(expected), f"invalid target row {line_number} in {path}")
            row = dict(zip(expected, fields))
            index = len(rows)
            rows.append(row)
            buckets[bucket_key(float(row["Lon"]), float(row["Lat"]))].append(index)
    require(len(rows) > 0, "target has no rows")
    return rows, buckets


def target_lookup(
    lon: float,
    lat: float,
    target_rows: list[dict[str, str]],
    buckets: dict[tuple[int, int], list[int]],
) -> dict[str, str]:
    center = bucket_key(lon, lat)
    matches: list[dict[str, str]] = []
    for di in range(-2, 3):
        for dj in range(-2, 3):
            for index in buckets.get((center[0] + di, center[1] + dj), ()):
                row = target_rows[index]
                if (
                    longitude_distance(float(row["Lon"]), lon) <= MAP_TOL
                    and abs(float(row["Lat"]) - lat) <= MAP_TOL
                ):
                    matches.append(row)
    require(len(matches) == 1, f"target lookup at ({lon}, {lat}) found {len(matches)} rows")
    return matches[0]


def close(a: float, b: float, absolute: float = AREA_TOL, relative: float = 0.0) -> bool:
    return abs(a - b) <= max(absolute, relative * max(1.0, abs(a), abs(b)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--migration", required=True)
    parser.add_argument("--fractions", required=True)
    parser.add_argument("--state", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--experiment", required=True)
    parser.add_argument("--source-state", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--expected-cells", type=int, default=107132)
    args = parser.parse_args()

    migration = read_sharded_tables(args.migration, "rank")
    fractions = read_sharded_tables(args.fractions, "shard")
    state = read_sharded_tables(args.state, "shard")
    expected_indexes = set(range(args.expected_cells))
    for name, rows in (("migration", migration), ("fraction read-back", fractions), ("state read-back", state)):
        require(set(rows) == expected_indexes, f"{name} contains {len(rows)} cells, not indexes 0..{args.expected_cells - 1}")

    migration_by_coordinate = index_rows_by_coordinate(migration, "migration")
    fractions_by_coordinate = index_rows_by_coordinate(fractions, "fraction read-back")
    state_by_coordinate = index_rows_by_coordinate(state, "state read-back")
    migration_coordinates = set(migration_by_coordinate)
    require(
        set(fractions_by_coordinate) == migration_coordinates,
        "fraction read-back coordinates do not match migration coordinates",
    )
    require(
        set(state_by_coordinate) == migration_coordinates,
        "state read-back coordinates do not match migration coordinates",
    )

    target_rows, target_buckets = read_target(args.target)
    maxima = defaultdict(float)
    totals = defaultdict(float)
    changed_peat_cells = 0

    reordered_cells = 0
    for coordinate, written in migration_by_coordinate.items():
        reread_fraction = fractions_by_coordinate[coordinate]
        reread_state = state_by_coordinate[coordinate]
        source_index = int(written["index"])
        fraction_index = int(reread_fraction["index"])
        state_index = int(reread_state["index"])
        location = f"source index {source_index}, coordinate {coordinate}"
        require(fraction_index == state_index, f"read-back index mismatch at {location}")
        if fraction_index != source_index:
            reordered_cells += 1
        lon = float(written["lon"])
        lat = float(written["lat"])

        for reread_name, reread in (("fraction", reread_fraction), ("state", reread_state)):
            require(close(float(reread["lon"]), lon, 1.0e-12), f"{reread_name} longitude mismatch at {location}")
            require(close(float(reread["lat"]), lat, 1.0e-12), f"{reread_name} latitude mismatch at {location}")

        direct_target = target_lookup(lon, lat, target_rows, target_buckets)
        target_peat = float(written["target_peat"])
        direct_peat = float(direct_target["CorrectedPeat"])
        require(close(target_peat, direct_peat), f"corrected peat target mismatch at {location}")
        require(int(written["luh3_valid"]) == int(direct_target["LUH3Valid"]), f"LUH3 validity mismatch at {location}")

        target_sum = sum(float(written[f"target_{land}"]) for land in LAND_CLASSES)
        final_sum = sum(float(written[f"final_{land}"]) for land in LAND_CLASSES)
        physical_total = float(reread_fraction["physical_total"])
        metadata_total = float(reread_fraction["metadata_total"])
        state_total = float(reread_state["physical_total"])
        for label, value in (
            ("target", target_sum), ("writer final", final_sum),
            ("physical read-back", physical_total), ("metadata read-back", metadata_total),
            ("state read-back", state_total),
        ):
            maxima[f"{label}_area_error"] = max(maxima[f"{label}_area_error"], abs(value - 1.0))
            require(close(value, 1.0), f"{label} fractions do not sum to one at {location}: {value:.17g}")

        for land in LAND_CLASSES:
            target_value = float(written[f"target_{land}"])
            final_value = float(written[f"final_{land}"])
            physical_value = float(reread_fraction[f"physical_{land}"])
            metadata_value = float(reread_fraction[f"metadata_{land}"])
            require(min(target_value, final_value, physical_value, metadata_value) >= -AREA_TOL,
                    f"negative {land} fraction at {location}")
            maxima["writer_target_fraction_difference"] = max(
                maxima["writer_target_fraction_difference"], abs(final_value - target_value)
            )
            maxima["physical_target_fraction_difference"] = max(
                maxima["physical_target_fraction_difference"], abs(physical_value - target_value)
            )
            maxima["metadata_physical_fraction_difference"] = max(
                maxima["metadata_physical_fraction_difference"], abs(metadata_value - physical_value)
            )
            require(close(final_value, target_value), f"writer {land} differs from target at {location}")
            require(close(physical_value, target_value), f"read-back {land} differs from target at {location}")
            require(close(metadata_value, physical_value), f"metadata {land} differs from physical state at {location}")

        raw_crop_sum = sum(float(direct_target[f"Crop_{crop}"]) for crop in CROP_TYPES)
        prescribed_crop = float(written["target_crop"])
        require(raw_crop_sum > 0.0 or close(prescribed_crop, 0.0),
                f"positive prescribed cropland has no LUH3 crop types at {location}")
        for crop in CROP_TYPES:
            expected_crop = (prescribed_crop * float(direct_target[f"Crop_{crop}"]) / raw_crop_sum
                             if raw_crop_sum > 0.0 else 0.0)
            physical_crop = float(reread_fraction[f"physical_{crop}"])
            metadata_crop = float(reread_fraction[f"metadata_{crop}"])
            maxima["physical_crop_type_target_difference"] = max(
                maxima["physical_crop_type_target_difference"], abs(physical_crop - expected_crop)
            )
            maxima["metadata_physical_crop_type_difference"] = max(
                maxima["metadata_physical_crop_type_difference"], abs(metadata_crop - physical_crop)
            )
            require(close(physical_crop, expected_crop),
                    f"physical crop type {crop} differs from LUH3 at {location}")
            require(close(metadata_crop, physical_crop),
                    f"metadata crop type {crop} differs from physical state at {location}")

        require(close(float(reread_state["physical_peat"]), target_peat), f"state peat differs from target at {location}")
        if not close(float(written["before_peat"]), target_peat):
            changed_peat_cells += 1

        for pool, state_column in (("c", "c_total"), ("n", "n_total"), ("water", "water_total")):
            baseline = float(written[f"baseline_{pool}"])
            final = float(written[f"final_{pool}"])
            reread = float(reread_state[state_column])
            writer_delta = abs(final - baseline)
            readback_delta = abs(reread - final)
            maxima[f"writer_{pool}_conservation_error"] = max(maxima[f"writer_{pool}_conservation_error"], writer_delta)
            maxima[f"readback_{pool}_difference"] = max(maxima[f"readback_{pool}_difference"], readback_delta)
            require(close(final, baseline, 1.0e-9, POOL_REL_TOL), f"{pool} not conserved by writer at {location}")
            require(close(reread, final, 1.0e-9, POOL_REL_TOL), f"{pool} changed after serialization at {location}")

        for name in (
            "recovered_donor1_area", "recovered_donor2_area", "recovered_synthetic_area",
            "recovered_c", "recovered_n", "recovered_water", "removed_area",
            "removed_c", "removed_n", "removed_water",
        ):
            totals[name] += float(written[name])

    result = {
        "status": "PASS",
        "experiment": args.experiment,
        "source_state": os.path.realpath(args.source_state),
        "target_table": os.path.realpath(args.target),
        "validated_cells": args.expected_cells,
        "serialization_reordered_cells": reordered_cells,
        "peat_changed_cells": changed_peat_cells,
        "tolerances": {
            "area_absolute": AREA_TOL,
            "pool_relative": POOL_REL_TOL,
            "target_coordinate_degrees": MAP_TOL,
        },
        "maximum_errors": dict(sorted(maxima.items())),
        "source_repair_accounting_sums": dict(sorted(totals.items())),
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(f"PASS: independently validated {args.expected_cells} cells for {args.experiment}")
    print(f"Report: {args.output}")


if __name__ == "__main__":
    main()
