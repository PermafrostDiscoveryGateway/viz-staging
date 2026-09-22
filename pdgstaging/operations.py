import os
from pathlib import Path

import geopandas as gpd
import pandas as pd

from .Deduplicator import deduplicate_by_footprint, deduplicate_neighbors
from .TilePathManager import TilePathManager
from .TileStager import TileStager


def _validate_source_key(value):
    """Validate that a source key is a safe single path component.

    Parameters
    ----------
    value : str
        The candidate source key; must be a non-empty string with no
        directory separators and not ``"."`` or ``".."``.

    Raises
    ------
    ValueError
        If the value is not a safe single path component.
    """
    if (
        not isinstance(value, str)
        or not value
        or value in (".", "..")
        or Path(value).name != value
    ):
        raise ValueError("source_key must be a non-empty path component")


def _valid_geopackage(path):
    """Check whether a path is a readable GeoPackage file.

    Parameters
    ----------
    path : pathlib.Path
        The file path to probe with ``geopandas.read_file``.

    Returns
    -------
    bool
        True when the path is a file that geopandas can read, False otherwise.
    """
    if not path.is_file():
        return False
    try:
        gpd.read_file(path)
    except Exception:
        return False
    return True


def _ensure_below(path, root):
    """Ensure a path stays inside the expected root directory.

    Parameters
    ----------
    path : str or pathlib.Path
        The candidate output path to check.
    root : str or pathlib.Path
        The directory the path must be equal to or contained in.

    Raises
    ------
    ValueError
        If the path escapes the root directory.
    """
    root = Path(root).absolute()
    path = Path(path).absolute()
    if root != path and root not in path.parents:
        raise ValueError(f"output escapes source directory: {path}")


def stage_source(
    input_path,
    shard_root,
    source_key,
    *,
    tms_id,
    z,
    path_structure,
    properties,
    input_crs=None,
    tolerance=None,
    overwrite=False,
):
    """Stage one source into tile shards; return the complete shard list.

    Parameters
    ----------
    input_path : str or pathlib.Path
        The source geospatial file to tile.
    shard_root : str or pathlib.Path
        The root directory under which per-source shard trees are written.
    source_key : str
        Safe single path component naming this source under ``shard_root``.
    tms_id : str
        Tile Matrix Set identifier used to build the tiling grid.
    z : int
        Maximum zoom level for tiling.
    path_structure : tuple
        Path template components passed to ``TilePathManager``.
    properties : dict
        Property-name mapping with ``"tile"`` and ``"centroid_tile"`` keys.
    input_crs : optional
        Coordinate reference system assigned to the input when needed.
    tolerance : float, optional
        Simplification tolerance passed to ``simplify_geoms``.
    overwrite : bool, default False
        When True, rewrite existing shards; otherwise reuse valid ones and
        raise on stale files.

    Returns
    -------
    list of pathlib.Path
        Sorted shard paths published for this source, empty when the input
        holds no Polygon or MultiPolygon geometry.

    Raises
    ------
    FileNotFoundError
        If the input file does not exist.
    ValueError
        If the source key is unsafe or the input cannot be loaded.
    FileExistsError
        If an existing shard is stale and ``overwrite`` is False.
    """
    _validate_source_key(source_key)
    source_root = Path(shard_root) / source_key
    input_path = Path(input_path)
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    tiles = TilePathManager(
        tms_id=tms_id,
        path_structure=path_structure,
        base_dirs={"staged": {"path": str(source_root), "ext": ".gpkg"}},
    )
    stager = TileStager(tiles=tiles, props=properties, max_z_level=z)
    gdf = stager.get_data(str(input_path))
    if gdf is None:
        raise ValueError(f"could not load input: {input_path}")
    gdf = gdf[gdf.geometry.type.isin(["Polygon", "MultiPolygon"])]
    if gdf.empty:
        return []
    gdf = stager.simplify_geoms(gdf, tolerance)
    gdf = stager.set_crs(gdf, input_crs)
    stager.grid = stager.make_tms_grid(gdf)
    gdf = stager.add_properties(gdf, str(input_path))
    published = []
    written = []
    try:
        for tile, data in gdf.groupby(properties["tile"]):
            data = data.copy()
            data[properties["tile"]] = data[properties["tile"]].astype(str)
            data[properties["centroid_tile"]] = data[
                properties["centroid_tile"]
            ].astype(str)
            output = Path(tiles.path_from_tile(tile, base_dir="staged"))
            _ensure_below(output, source_root)
            if output.exists():
                if not overwrite:
                    if _valid_geopackage(output):
                        published.append(output)
                        continue
                    raise FileExistsError(f"shard already exists: {output}")
            temporary = output.with_name(f"{output.stem}.tmp{output.suffix}")
            output.parent.mkdir(parents=True, exist_ok=True)
            try:
                data.to_file(temporary)
                if not _valid_geopackage(temporary):
                    raise ValueError(f"invalid shard: {temporary}")
                os.replace(temporary, output)
                published.append(output)
                written.append(output)
            finally:
                temporary.unlink(missing_ok=True)
        return sorted(published)
    except Exception:
        for output in written:
            output.unlink(missing_ok=True)
        raise


def _read_shard(path):
    """Read one staged shard into a GeoDataFrame.

    Parameters
    ----------
    path : str or pathlib.Path
        The shard file to read.

    Returns
    -------
    geopandas.GeoDataFrame
        The shard contents, guaranteed to include a ``"geometry"`` column.

    Raises
    ------
    FileNotFoundError
        If the shard file does not exist.
    ValueError
        If the shard is unreadable or has no ``"geometry"`` column.
    """
    shard = Path(path)
    if not shard.is_file():
        raise FileNotFoundError(f"missing shard: {shard}")
    try:
        data = gpd.read_file(shard)
    except Exception as error:
        raise ValueError(f"unreadable shard: {shard}") from error
    if "geometry" not in data.columns:
        raise ValueError(f"invalid shard: {shard}")
    return data


def _check_matching(frames):
    """Check that shard frames share CRS, schema, and tile identity.

    Parameters
    ----------
    frames : list of geopandas.GeoDataFrame
        The shard contents to compare; the first frame is the reference.

    Raises
    ------
    ValueError
        If the CRS values differ, the column sets differ, or the
        single-valued ``"Tile(...)"`` columns disagree across frames.
    """
    first, rest = frames[0], frames[1:]
    if any(frame.crs != first.crs for frame in rest):
        raise ValueError("shard CRS mismatch")
    columns = set(first.columns)
    if any(set(frame.columns) != columns for frame in rest):
        raise ValueError("shard schema mismatch")
    tile_columns = [
        column
        for column in first.columns
        if column != "geometry"
        and all(
            frame[column].astype(str).str.startswith("Tile(").all()
            and frame[column].nunique() == 1
            for frame in frames
        )
    ]
    if tile_columns and not any(
        all(frame[column].iloc[0] == first[column].iloc[0] for frame in rest)
        for column in tile_columns
    ):
        raise ValueError("shard tile mismatch")


def _deduplicated_data(shards, deduplication):
    """Concatenate shards and deduplicate the combined features.

    Parameters
    ----------
    shards : list of pathlib.Path
        The shard files to read and combine.
    deduplication : str or dict
        Either a method name (``"neighbors"`` or ``"footprints"``) or a
        mapping with a ``"method"`` key and optional ``"options"`` dict
        forwarded to the deduplicator.

    Returns
    -------
    geopandas.GeoDataFrame
        The combined, deduplicated features.

    Raises
    ------
    ValueError
        If the deduplication method is unknown.
    """
    data = gpd.GeoDataFrame(
        pd.concat([gpd.read_file(path) for path in shards], ignore_index=True)
    )
    method = deduplication if isinstance(deduplication, str) else deduplication.get("method")
    options = {} if isinstance(deduplication, str) else deduplication.get("options", {})
    if method == "neighbors":
        return deduplicate_neighbors(data, **options)
    if method == "footprints":
        return deduplicate_by_footprint(data, **options)
    raise ValueError(f"unknown deduplication method: {method}")


def merge_staged_tile(shard_paths, output_path, *, deduplication=None, overwrite=False):
    """Merge an explicit shard list into one tile; return the output path.

    Parameters
    ----------
    shard_paths : list of str or pathlib.Path
        The shard files to merge; duplicates are collapsed and entries sorted.
    output_path : str or pathlib.Path
        The destination tile file to write.
    deduplication : str or dict, optional
        Deduplication spec forwarded to ``_deduplicated_data``; when None the
        frames are concatenated without deduplication.
    overwrite : bool, default False
        When True, rewrite an existing output; otherwise reuse a valid output
        and raise on a stale one.

    Returns
    -------
    pathlib.Path
        The output tile path.

    Raises
    ------
    ValueError
        If no shards are given, a shard is missing or invalid, or the merged
        tile fails validation.
    FileExistsError
        If the output exists but is stale and ``overwrite`` is False.
    """
    shards = sorted({Path(path) for path in shard_paths})
    if not shards:
        raise ValueError("at least one shard is required")
    frames = [_read_shard(path) for path in shards]
    _check_matching(frames)
    output = Path(output_path)
    if output.exists():
        if not overwrite:
            if _valid_geopackage(output):
                return output
            raise FileExistsError(f"output already exists: {output}")
    temporary = output.with_name(f"{output.stem}.tmp{output.suffix}")
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        if deduplication is None:
            merged = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True))
        else:
            merged = _deduplicated_data(shards, deduplication)
        merged.to_file(temporary)
        if not _valid_geopackage(temporary):
            raise ValueError(f"invalid merged tile: {temporary}")
        os.replace(temporary, output)
        return output
    finally:
        temporary.unlink(missing_ok=True)
