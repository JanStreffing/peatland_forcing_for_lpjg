#!/usr/bin/env python3
"""Generate and validate a complete, resolution-wide LPJ-GUESS peat map."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile
from typing import Iterable, Sequence

import netCDF4
import numpy as np
from scipy.spatial import cKDTree


class PeatMapError(RuntimeError):
    """Raised when an input or generated peat map violates its contract."""


DEFAULT_GRIDLIST_TOLERANCE = 5.0e-5


def canonical_longitude(longitude: np.ndarray) -> np.ndarray:
    """Return longitudes in the half-open interval [-180, 180)."""
    longitude = np.asarray(longitude, dtype=np.float64)
    return np.remainder(longitude + 180.0, 360.0) - 180.0


def _coordinate_keys(
    longitude: np.ndarray, latitude: np.ndarray, precision: int
) -> np.ndarray:
    """Return exact integer keys matching fixed-decimal text output."""
    scale = 10**precision
    lon = canonical_longitude(longitude)
    lon_key = np.rint(lon * scale).astype(np.int64)
    lat_key = np.rint(np.asarray(latitude, dtype=np.float64) * scale).astype(
        np.int64
    )

    # A value just below +180 can round to +180. Canonicalize that written key
    # back to -180 so output and validation have one dateline representation.
    lon_key[lon_key >= 180 * scale] -= 360 * scale

    keys = np.empty(lon_key.size, dtype=[("lon", "<i8"), ("lat", "<i8")])
    keys["lon"] = lon_key
    keys["lat"] = lat_key
    return keys


def _format_key(key: np.void, precision: int) -> str:
    scale = float(10**precision)
    return f"({key['lon'] / scale:.{precision}f}, {key['lat'] / scale:.{precision}f})"


def _check_coordinate_values(
    longitude: np.ndarray,
    latitude: np.ndarray,
    label: str,
) -> tuple[np.ndarray, np.ndarray]:
    longitude = np.asarray(longitude, dtype=np.float64).reshape(-1)
    latitude = np.asarray(latitude, dtype=np.float64).reshape(-1)
    if longitude.size == 0:
        raise PeatMapError(f"{label}: coordinate set is empty")
    if longitude.size != latitude.size:
        raise PeatMapError(
            f"{label}: longitude/latitude sizes differ "
            f"({longitude.size} != {latitude.size})"
        )
    if not np.all(np.isfinite(longitude)) or not np.all(np.isfinite(latitude)):
        raise PeatMapError(f"{label}: coordinates contain a non-finite value")
    # Both conventional longitude representations are accepted. Values outside
    # one complete turn are almost certainly corrupt coordinates, not wrapping.
    if np.any(longitude < -360.0) or np.any(longitude > 360.0):
        raise PeatMapError(f"{label}: longitude is outside [-360, 360]")
    if np.any(latitude < -90.0) or np.any(latitude > 90.0):
        raise PeatMapError(f"{label}: latitude is outside [-90, 90]")
    return canonical_longitude(longitude), latitude


def _check_unique_coordinates(
    longitude: np.ndarray,
    latitude: np.ndarray,
    precision: int,
    label: str,
) -> np.ndarray:
    keys = _coordinate_keys(longitude, latitude, precision)
    unique, counts = np.unique(keys, return_counts=True)
    duplicate = np.flatnonzero(counts != 1)
    if duplicate.size:
        examples = ", ".join(
            _format_key(unique[index], precision) for index in duplicate[:5]
        )
        raise PeatMapError(
            f"{label}: {duplicate.size} duplicate coordinate key(s) at "
            f"{precision} decimal places; examples: {examples}"
        )
    return keys


def read_peat_table(path: os.PathLike[str] | str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a strict Lon/Lat/Peat_Frac whitespace table."""
    path = Path(path)
    longitude: list[float] = []
    latitude: list[float] = []
    peat: list[float] = []

    try:
        stream = path.open("r", encoding="utf-8")
    except OSError as error:
        raise PeatMapError(f"cannot open peat table {path}: {error}") from error

    with stream:
        header = stream.readline().split()
        if header != ["Lon", "Lat", "Peat_Frac"]:
            raise PeatMapError(
                f"{path}: expected header 'Lon Lat Peat_Frac', got {header!r}"
            )
        for line_number, raw_line in enumerate(stream, start=2):
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) != 3:
                raise PeatMapError(
                    f"{path}:{line_number}: expected three columns, got {len(fields)}"
                )
            try:
                lon_value, lat_value, peat_value = map(float, fields)
            except ValueError as error:
                raise PeatMapError(
                    f"{path}:{line_number}: a field is not numeric"
                ) from error
            longitude.append(lon_value)
            latitude.append(lat_value)
            peat.append(peat_value)

    lon_array, lat_array = _check_coordinate_values(longitude, latitude, str(path))
    peat_array = np.asarray(peat, dtype=np.float64)
    if not np.all(np.isfinite(peat_array)):
        raise PeatMapError(f"{path}: peat fraction contains a non-finite value")
    if np.any(peat_array < 0.0) or np.any(peat_array > 1.0):
        minimum = float(np.min(peat_array))
        maximum = float(np.max(peat_array))
        raise PeatMapError(
            f"{path}: peat fraction is outside [0, 1] (range {minimum} to {maximum})"
        )
    return lon_array, lat_array, peat_array


def _netcdf_array(variable: netCDF4.Variable, label: str) -> np.ndarray:
    values = variable[:]
    if np.ma.isMaskedArray(values):
        mask = np.ma.getmaskarray(values)
        if np.any(mask):
            raise PeatMapError(f"{label}: coordinate variable contains missing values")
        values = values.data
    result = np.asarray(values, dtype=np.float64)
    if result.size == 0:
        raise PeatMapError(f"{label}: coordinate variable is empty")
    return result


def _find_coordinate_variable(
    dataset: netCDF4.Dataset, requested_name: str, standard_name: str
) -> netCDF4.Variable | None:
    if requested_name in dataset.variables:
        return dataset.variables[requested_name]
    candidates = [
        variable
        for variable in dataset.variables.values()
        if getattr(variable, "standard_name", "") == standard_name
    ]
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        names = ", ".join(variable.name for variable in candidates)
        raise PeatMapError(
            f"NetCDF has multiple variables with standard_name={standard_name!r}: {names}"
        )
    return None


def read_target_grid(
    path: os.PathLike[str] | str,
    lon_variable: str = "lon",
    lat_variable: str = "lat",
) -> tuple[np.ndarray, np.ndarray]:
    """Read paired, rectilinear, two-dimensional, or reduced-grid coordinates."""
    path = Path(path)
    try:
        dataset = netCDF4.Dataset(path, mode="r")
    except OSError as error:
        raise PeatMapError(f"cannot open target NetCDF {path}: {error}") from error

    try:
        lat_var = _find_coordinate_variable(dataset, lat_variable, "latitude")
        lon_var = _find_coordinate_variable(dataset, lon_variable, "longitude")
        if lat_var is None:
            raise PeatMapError(
                f"{path}: cannot find latitude variable {lat_variable!r} or a unique "
                "standard_name='latitude' variable"
            )

        latitude = _netcdf_array(lat_var, f"{path}:{lat_var.name}")
        if lon_var is None:
            if "reduced_points" not in dataset.variables or latitude.ndim != 1:
                raise PeatMapError(
                    f"{path}: cannot find longitude variable {lon_variable!r} or a "
                    "unique standard_name='longitude' variable"
                )
            reduced_points_raw = _netcdf_array(
                dataset.variables["reduced_points"], f"{path}:reduced_points"
            )
            if reduced_points_raw.ndim != 1 or reduced_points_raw.size != latitude.size:
                raise PeatMapError(
                    f"{path}: reduced_points must be one-dimensional and match latitude"
                )
            if not np.all(reduced_points_raw == np.rint(reduced_points_raw)):
                raise PeatMapError(f"{path}: reduced_points contains a non-integer value")
            reduced_points = reduced_points_raw.astype(np.int64)
            if np.any(reduced_points <= 0):
                raise PeatMapError(f"{path}: reduced_points contains a non-positive row")
            longitude_rows = [
                np.arange(count, dtype=np.float64) * (360.0 / float(count))
                for count in reduced_points
            ]
            latitude_rows = [
                np.full(count, latitude[index], dtype=np.float64)
                for index, count in enumerate(reduced_points)
            ]
            longitude = np.concatenate(longitude_rows)
            latitude = np.concatenate(latitude_rows)
        else:
            longitude = _netcdf_array(lon_var, f"{path}:{lon_var.name}")
            if longitude.shape == latitude.shape:
                longitude = longitude.reshape(-1)
                latitude = latitude.reshape(-1)
            elif longitude.ndim == 1 and latitude.ndim == 1:
                lon_mesh, lat_mesh = np.meshgrid(longitude, latitude, indexing="xy")
                longitude = lon_mesh.reshape(-1)
                latitude = lat_mesh.reshape(-1)
            else:
                try:
                    lon_broadcast, lat_broadcast = np.broadcast_arrays(
                        longitude, latitude
                    )
                except ValueError as error:
                    raise PeatMapError(
                        f"{path}: longitude shape {longitude.shape} and latitude shape "
                        f"{latitude.shape} are neither paired nor broadcastable"
                    ) from error
                longitude = lon_broadcast.reshape(-1)
                latitude = lat_broadcast.reshape(-1)
    finally:
        dataset.close()

    return _check_coordinate_values(longitude, latitude, str(path))


def read_gridlist(path: os.PathLike[str] | str) -> tuple[np.ndarray, np.ndarray]:
    """Read either supported LPJ-GUESS gridlist representation.

    Static two-column files use ``Lat Lon``. Coupler-created four-column files
    use ``Lon Lat soil_code cell_index``. The different column counts make this
    distinction unambiguous; mixed formats in one file are rejected.
    """
    path = Path(path)
    longitude: list[float] = []
    latitude: list[float] = []
    try:
        stream = path.open("r", encoding="utf-8")
    except OSError as error:
        raise PeatMapError(f"cannot open gridlist {path}: {error}") from error
    detected_columns: int | None = None
    with stream:
        for line_number, raw_line in enumerate(stream, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) not in (2, 4):
                raise PeatMapError(
                    f"{path}:{line_number}: expected two columns (Lat Lon) or four "
                    f"columns (Lon Lat soil_code cell_index), got {len(fields)}"
                )
            if detected_columns is None:
                detected_columns = len(fields)
            elif len(fields) != detected_columns:
                raise PeatMapError(
                    f"{path}:{line_number}: gridlist mixes {detected_columns}- and "
                    f"{len(fields)}-column rows"
                )
            try:
                if len(fields) == 2:
                    lat_value, lon_value = map(float, fields)
                else:
                    lon_value, lat_value = map(float, fields[:2])
                    # Parse metadata as numeric too so malformed rows cannot hide
                    # behind ignored columns.
                    float(fields[2])
                    float(fields[3])
            except ValueError as error:
                raise PeatMapError(
                    f"{path}:{line_number}: a field is not numeric"
                ) from error
            longitude.append(lon_value)
            latitude.append(lat_value)
    return _check_coordinate_values(longitude, latitude, str(path))


def interpolate_periodic_nearest(
    source_longitude: np.ndarray,
    source_latitude: np.ndarray,
    source_peat: np.ndarray,
    target_longitude: np.ndarray,
    target_latitude: np.ndarray,
    chunk_size: int = 200_000,
) -> np.ndarray:
    """Nearest-neighbour interpolation with periodic longitude.

    Distance remains Euclidean in longitude/latitude degrees, matching the
    previous scientific method. Three longitude copies make the dateline
    periodic. Exact distance ties are resolved by the earliest source row for
    deterministic output.
    """
    if chunk_size <= 0:
        raise PeatMapError("chunk_size must be positive")
    source_lon, source_lat = _check_coordinate_values(
        source_longitude, source_latitude, "source peat coordinates"
    )
    target_lon, target_lat = _check_coordinate_values(
        target_longitude, target_latitude, "target grid coordinates"
    )
    source_peat = np.asarray(source_peat, dtype=np.float64).reshape(-1)
    if source_peat.size != source_lon.size:
        raise PeatMapError("source peat size does not match source coordinate size")

    source_count = source_lon.size
    augmented_lon = np.concatenate((source_lon - 360.0, source_lon, source_lon + 360.0))
    augmented_lat = np.tile(source_lat, 3)
    augmented_to_source = np.tile(np.arange(source_count, dtype=np.int64), 3)
    tree = cKDTree(np.column_stack((augmented_lon, augmented_lat)))
    result = np.empty(target_lon.size, dtype=np.float64)

    for begin in range(0, target_lon.size, chunk_size):
        end = min(begin + chunk_size, target_lon.size)
        target = np.column_stack((target_lon[begin:end], target_lat[begin:end]))
        try:
            distances, indices = tree.query(target, k=2, workers=-1)
        except TypeError:  # scipy versions before the workers argument
            distances, indices = tree.query(target, k=2)

        chosen = augmented_to_source[indices[:, 0]].copy()
        tied = np.isclose(distances[:, 0], distances[:, 1], rtol=0.0, atol=1.0e-12)
        tied &= augmented_to_source[indices[:, 0]] != augmented_to_source[indices[:, 1]]
        for local_index in np.flatnonzero(tied):
            point = target[local_index]
            nearest_distance = distances[local_index, 0]
            candidate_indices = tree.query_ball_point(
                point, nearest_distance + 1.0e-12
            )
            candidates = np.asarray(candidate_indices, dtype=np.int64)
            candidate_points = tree.data[candidates]
            exact_distances = np.sqrt(np.sum((candidate_points - point) ** 2, axis=1))
            equally_near = candidates[
                np.isclose(
                    exact_distances, nearest_distance, rtol=0.0, atol=1.0e-12
                )
            ]
            chosen[local_index] = int(
                np.min(augmented_to_source[equally_near])
            )
        result[begin:end] = source_peat[chosen]

    if not np.all(np.isfinite(result)) or np.any(result < 0.0) or np.any(result > 1.0):
        raise PeatMapError("interpolation produced an invalid peat fraction")
    return result


def _compare_coordinate_sets(
    map_keys: np.ndarray,
    target_keys: np.ndarray,
    precision: int,
    map_label: str,
    target_label: str,
) -> None:
    sorted_map = np.sort(map_keys, order=("lon", "lat"))
    sorted_target = np.sort(target_keys, order=("lon", "lat"))
    if np.array_equal(sorted_map, sorted_target):
        return

    missing = np.setdiff1d(sorted_target, sorted_map, assume_unique=True)
    extra = np.setdiff1d(sorted_map, sorted_target, assume_unique=True)
    missing_examples = ", ".join(_format_key(key, precision) for key in missing[:5])
    extra_examples = ", ".join(_format_key(key, precision) for key in extra[:5])
    raise PeatMapError(
        f"coordinate mismatch between {map_label} and {target_label}: "
        f"{missing.size} missing, {extra.size} extra; "
        f"missing examples [{missing_examples}]; extra examples [{extra_examples}]"
    )


def validate_map(
    map_path: os.PathLike[str] | str,
    target_longitude: np.ndarray,
    target_latitude: np.ndarray,
    coordinate_precision: int,
    expected_count: int | None = None,
    gridlists: Iterable[os.PathLike[str] | str] = (),
    gridlist_tolerance: float = DEFAULT_GRIDLIST_TOLERANCE,
) -> dict[str, int | float]:
    """Validate values, full target equality, and optional gridlist coverage."""
    if coordinate_precision < 6:
        raise PeatMapError("coordinate_precision must be at least 6")
    if not np.isfinite(gridlist_tolerance) or gridlist_tolerance < 0.0:
        raise PeatMapError("gridlist_tolerance must be finite and non-negative")
    map_lon, map_lat, peat = read_peat_table(map_path)
    target_lon, target_lat = _check_coordinate_values(
        target_longitude, target_latitude, "target grid"
    )
    target_count = target_lon.size
    if expected_count is not None and target_count != expected_count:
        raise PeatMapError(
            f"target grid has {target_count} rows, expected {expected_count}"
        )
    if map_lon.size != target_count:
        raise PeatMapError(
            f"{map_path}: map has {map_lon.size} rows, target grid has {target_count}"
        )

    map_keys = _check_unique_coordinates(
        map_lon, map_lat, coordinate_precision, str(map_path)
    )
    target_keys = _check_unique_coordinates(
        target_lon, target_lat, coordinate_precision, "target grid"
    )
    _compare_coordinate_sets(
        map_keys,
        target_keys,
        coordinate_precision,
        str(map_path),
        "target grid",
    )

    map_points = np.column_stack((map_lon, map_lat))
    periodic_map_points = np.vstack(
        (
            np.column_stack((map_lon - 360.0, map_lat)),
            map_points,
            np.column_stack((map_lon + 360.0, map_lat)),
        )
    )
    periodic_to_map = np.tile(np.arange(map_lon.size, dtype=np.int64), 3)
    map_tree = cKDTree(periodic_map_points)
    gridlist_count = 0
    for gridlist in gridlists:
        grid_lon, grid_lat = read_gridlist(gridlist)
        _check_unique_coordinates(
            grid_lon, grid_lat, coordinate_precision, str(gridlist)
        )
        grid_points = np.column_stack((grid_lon, grid_lat))
        # Coverage is a per-coordinate contract, not a radial distance test:
        # |periodic delta longitude| and |delta latitude| must each be within
        # the tolerance. Chebyshev distance expresses exactly that condition.
        distances, indices = map_tree.query(grid_points, k=2, p=np.inf)
        covered = distances[:, 0] <= gridlist_tolerance
        if not np.all(covered):
            missing = _coordinate_keys(
                grid_lon[~covered], grid_lat[~covered], coordinate_precision
            )
            examples = ", ".join(
                _format_key(key, coordinate_precision) for key in missing[:5]
            )
            raise PeatMapError(
                f"{gridlist}: {missing.size} coordinates are absent from the peat map; "
                f"nearest distance exceeds {gridlist_tolerance:g} degree; "
                f"examples: {examples}"
            )
        ambiguous = (distances[:, 1] <= gridlist_tolerance) & (
            periodic_to_map[indices[:, 1]] != periodic_to_map[indices[:, 0]]
        )
        if np.any(ambiguous):
            raise PeatMapError(
                f"{gridlist}: {np.count_nonzero(ambiguous)} coordinate(s) match more "
                f"than one map row within {gridlist_tolerance:g} degree"
            )
        gridlist_count += 1

    return {
        "rows": int(map_lon.size),
        "positive_peat_rows": int(np.count_nonzero(peat > 0.0)),
        "minimum_peat": float(np.min(peat)),
        "maximum_peat": float(np.max(peat)),
        "gridlists_validated": gridlist_count,
    }


def _fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY
    directory_fd = os.open(path, flags)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def write_map_atomically(
    output_path: os.PathLike[str] | str,
    longitude: np.ndarray,
    latitude: np.ndarray,
    peat: np.ndarray,
    coordinate_precision: int,
    fraction_precision: int,
    expected_count: int | None,
    gridlists: Sequence[os.PathLike[str] | str],
    gridlist_tolerance: float,
    replace: bool,
) -> dict[str, int | float]:
    """Write, fsync, validate, and atomically publish a completed map."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists() and not replace:
        raise PeatMapError(
            f"refusing to replace existing output {output_path}; use --replace explicitly"
        )

    keys = _coordinate_keys(longitude, latitude, coordinate_precision)
    scale = float(10**coordinate_precision)
    written_lon = keys["lon"].astype(np.float64) / scale
    written_lat = keys["lat"].astype(np.float64) / scale

    descriptor, temporary_name = tempfile.mkstemp(
        dir=output_path.parent,
        prefix=f".{output_path.name}.",
        suffix=".tmp",
        text=True,
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as stream:
            # The input is intended to be reviewable and usable by collaborating
            # model users; mkstemp starts at 0600, so publish it read-only to
            # group/others while retaining owner write permission.
            os.fchmod(stream.fileno(), 0o644)
            stream.write("Lon\tLat\tPeat_Frac\n")
            for lon_value, lat_value, peat_value in zip(
                written_lon, written_lat, peat, strict=True
            ):
                stream.write(
                    f"{lon_value:.{coordinate_precision}f}\t"
                    f"{lat_value:.{coordinate_precision}f}\t"
                    f"{peat_value:.{fraction_precision}f}\n"
                )
            stream.flush()
            os.fsync(stream.fileno())

        report = validate_map(
            temporary_path,
            longitude,
            latitude,
            coordinate_precision,
            expected_count=expected_count,
            gridlists=gridlists,
            gridlist_tolerance=gridlist_tolerance,
        )
        if output_path.exists() and not replace:
            raise PeatMapError(
                f"output {output_path} appeared while generating; refusing to replace it"
            )
        os.replace(temporary_path, output_path)
        _fsync_directory(output_path.parent)
        return report
    except BaseException:
        try:
            temporary_path.unlink()
        except FileNotFoundError:
            pass
        raise


def generate_map(
    source_path: os.PathLike[str] | str,
    target_grid_path: os.PathLike[str] | str,
    resolution: str,
    output_directory: os.PathLike[str] | str,
    expected_count: int | None = None,
    gridlists: Sequence[os.PathLike[str] | str] = (),
    coordinate_precision: int = 6,
    fraction_precision: int = 8,
    lon_variable: str = "lon",
    lat_variable: str = "lat",
    chunk_size: int = 200_000,
    gridlist_tolerance: float = DEFAULT_GRIDLIST_TOLERANCE,
    replace: bool = False,
) -> tuple[Path, dict[str, int | float]]:
    """Generate one complete `<resolution>_peat_frac.txt` map."""
    if not resolution or any(character not in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-" for character in resolution):
        raise PeatMapError(
            "resolution must be a non-empty identifier containing only letters, "
            "digits, underscore, or hyphen"
        )
    if coordinate_precision < 6:
        raise PeatMapError("coordinate_precision must be at least 6")
    if fraction_precision < 8:
        raise PeatMapError("fraction_precision must be at least 8")

    source_lon, source_lat, source_peat = read_peat_table(source_path)
    _check_unique_coordinates(source_lon, source_lat, 12, str(source_path))
    target_lon, target_lat = read_target_grid(
        target_grid_path, lon_variable=lon_variable, lat_variable=lat_variable
    )
    _check_unique_coordinates(
        target_lon, target_lat, coordinate_precision, str(target_grid_path)
    )
    if expected_count is not None and target_lon.size != expected_count:
        raise PeatMapError(
            f"target grid has {target_lon.size} rows, expected {expected_count}"
        )

    peat = interpolate_periodic_nearest(
        source_lon,
        source_lat,
        source_peat,
        target_lon,
        target_lat,
        chunk_size=chunk_size,
    )
    output_path = Path(output_directory) / f"{resolution}_peat_frac.txt"
    report = write_map_atomically(
        output_path,
        target_lon,
        target_lat,
        peat,
        coordinate_precision,
        fraction_precision,
        expected_count,
        gridlists,
        gridlist_tolerance,
        replace,
    )
    return output_path, report


def format_report(path: os.PathLike[str] | str, report: dict[str, int | float]) -> str:
    return (
        f"VALID: {path}\n"
        f"rows={report['rows']}\n"
        f"positive_peat_rows={report['positive_peat_rows']}\n"
        f"peat_range=[{report['minimum_peat']}, {report['maximum_peat']}]\n"
        f"gridlists_validated={report['gridlists_validated']}"
    )
