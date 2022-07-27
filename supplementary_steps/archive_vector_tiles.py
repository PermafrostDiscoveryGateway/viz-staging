"""
This file prepares staged, vector tiles for archiving in the DataONE network.
It has methods that do the following:
    1. Find all the paths of vector tiles in the staged directory.
    2. Open each path in GeoPandas
    3. Remove the polygons whose centroids are not contained within the
       boundaries of the tile. (This avoids archiving the exact same polygons
       in two different tiles, in the case that the polygon intersects a tile's
       boundary.)
    4. Remove the "centroids_within_tile" and "tile" columns, which are no
       longer needed after the above step. The tile is still identified using
       by the "centroid_tile" property.
    5. Write the GeoPandas object to a file in a given archive directory, still
       maintaining the tile path structure (e.g.
       {TileMatrix}/{TileRow}/{TileCol})
"""

import os

import geopandas as gpd
import pandas as pd
from pdgstaging.TileStager import TileStager
import parsl
from parsl import python_app


def filter_gdf(gdf, filter_dict):
    """
    Filter a GeoPandas object by a dictionary of column names and values.

    Parameters
    ----------
    gdf : GeoPandas object filter_dict : dict
        A dictionary of column names and values to filter by.

    Returns
    -------
    gdf : GeoPandas object
    """
    return gdf.loc[(gdf[list(filter_dict)] ==
                   pd.Series(filter_dict)).all(axis=1)]


def remove_columns(gdf, columns):
    """
    Remove columns from a GeoPandas object.

    Parameters
    ----------
    gdf : GeoPandas object columns : list of str
        A list of column names to remove.

    Returns
    -------
    gdf : GeoPandas object
    """
    return gdf.drop(columns, axis=1)


def prepare_for_archiving(path, stager):
    """
    Prepare a vector tile for archiving. This involves:
        1. Removing the polygons whose centroids are not contained within the
           boundaries of the tile.
        2. Removing the "centroids_within_tile" and "tile" columns

    Parameters
    ----------
    path : str
        The path to the vector tile to prepare for archiving.
    stager : TileStager object or dict
        The TileStager object used to stage the vector tiles, or a dictionary
        of the TileStager configuration.

    Returns
    -------
    gdf : GeoPandas object
    """
    gdf = gpd.read_file(path)
    if isinstance(stager, (dict, str)):
        stager = TileStager(stager)
    config = stager.config
    prop_centroid_within_tile = config.polygon_prop('centroid_within_tile')
    prop_tile = stager.config.polygon_prop('tile')
    gdf = filter_gdf(gdf, {prop_centroid_within_tile: True})
    gdf = remove_columns(gdf, [prop_centroid_within_tile, prop_tile])
    return gdf


def archive_vector_tile(path, stager, archive_dir, ext=None):
    """
    Archive a single vector tile.

    Parameters
    ----------
    path : str
        The path to the vector tile to archive.
    stager : TileStager object or dict
        The TileStager object used to stage the vector tiles, or a dictionary
        of the TileStager configuration.
    archive_dir : str
        The path to the directory to archive the vector tile to.
    ext : str, optional
        The extension to use for the archive file. If not provided, the
        extension will be determined from the path.

    Returns
    -------
    archive_path : str
        The path to the archived file.

    """
    if isinstance(stager, (dict, str)):
        stager = TileStager(stager)
    if ext is None:
        ext = os.path.splitext(path)[1]
    path_manager = stager.tiles
    if 'archive' not in path_manager.base_dirs:
        path_manager.add_base_dir('archive', archive_dir, ext)
    out_path = path_manager.path_from_tile(tile=path, base_dir='archive')
    gdf = prepare_for_archiving(path, stager)
    # mk di
    stager.tiles.create_dirs(out_path)
    gdf.to_file(out_path)
    return out_path


def archive_vector_tiles(paths, stager, archive_dir, ext=None):
    """
    Archive a list of vector tiles.

    Parameters
    ----------
    path : str
        The path to the vector tile to archive.
    stager : TileStager object or dict
        The TileStager object used to stage the vector tiles, or a dictionary
        of the TileStager configuration.
    archive_dir : str
        The path to the directory to archive the vector tile to.
    ext : str, optional
        The extension to use for the archive file. If not provided, the
        extension will be determined from the path.

    Returns
    -------
    archive_path : list
        A list of paths to the archived files.

    """
    out_paths = []
    for path in paths:
        out_path = archive_vector_tile(path, stager, archive_dir, ext)
        out_paths.append(out_path)
    return out_paths


def archive_all_vector_tiles(stager, archive_dir, ext=None):
    """
    Archive all vector tiles in the staged directory.

    Parameters
    ----------
    stager : TileStager object or dict
        The TileStager object used to stage the vector tiles, or a dictionary
        of the TileStager configuration.
    archive_dir : str
        The path to the directory to archive the vector tile to.
    ext : str, optional
        The extension to use for the archive file. If not provided, the
        extension will be determined from the path.

    Returns
    -------
    archive_paths : list
        A list of paths to the archived files.

    """
    if isinstance(stager, (dict, str)):
        stager = TileStager(stager)
    path_manager = stager.tiles
    paths = path_manager.get_filenames_from_dir('staged')
    return archive_vector_tiles(paths, stager, archive_dir, ext)


@python_app
def archive_vector_tile_parsl(path, stager, archive_dir, ext=None):
    """
    Archive a single vector tile as a Parsl task.

    Parameters
    ----------
    path : str
        The path to the vector tile to archive.
    stager : dict
        A dictionary of the TileStager configuration.
    archive_dir : str
        The path to the directory to archive the vector tile to.
    ext : str, optional
        The extension to use for the archive file. If not provided, the
        extension will be determined from the path.

    Returns
    -------
    archive_path : parsl.app.futures.DataFuture
        A future that will contain the path to the archived file as the result.

    """
    if isinstance(stager, (dict, str)):
        stager = TileStager(stager)
    if ext is None:
        ext = os.path.splitext(path)[1]
    path_manager = stager.tiles
    if 'archive' not in path_manager.base_dirs:
        path_manager.add_base_dir('archive', archive_dir, ext)
    out_path = path_manager.path_from_tile(tile=path, base_dir='archive')
    gdf = prepare_for_archiving(path, stager)
    stager.tiles.create_dirs(out_path)
    gdf.to_file(out_path)

    # prepare for archiving
    gdf = gpd.read_file(path)
    config = stager.config
    prop_centroid_within_tile = config.polygon_prop('centroid_within_tile')
    prop_tile = stager.config.polygon_prop('tile')
    gdf = filter_gdf(gdf, {prop_centroid_within_tile: True})
    gdf = remove_columns(gdf, [prop_centroid_within_tile, prop_tile])

    gdf.to_file(out_path)
    return out_path


def archive_vector_tiles_parsl(paths, stager, archive_dir, ext=None):
    """
    Archive a list of vector tiles in parallel with Parsl.

    Parameters
    ----------
    path : str
        The path to the vector tile to archive.
    stager : dict
        A dictionary of the TileStager configuration.
    archive_dir : str
        The path to the directory to archive the vector tile to.
    ext : str, optional
        The extension to use for the archive file. If not provided, the
        extension will be determined from the path.

    Returns
    -------
    archive_paths : list of parsl.app.futures.DataFuture
        A list of futures that will contain the paths to the archived files as
        the result of each future.

    """
    data_futures = []
    for path in paths:
        data_future = archive_vector_tile_parsl(path, stager, archive_dir, ext)
        data_futures.append(data_future)
    return data_futures


def archive_all_vector_tiles_parsl(stager, archive_dir, ext=None):
    """
    Archive all vector tiles in the staged directory in parallel with Parsl.

    Parameters
    ----------
    path : str
        The path to the vector tile to archive.
    stager : dict
        A dictionary of the TileStager configuration.
    archive_dir : str
        The path to the directory to archive the vector tile to.
    ext : str, optional
        The extension to use for the archive file. If not provided, the
        extension will be determined from the path.

    Returns
    -------
    archive_paths : list of parsl.app.futures.DataFuture
        A list of futures that will contain the paths to the archived files as
        the result of each future.

    """
    if isinstance(stager, (dict, str)):
        stager = TileStager(stager)
    path_manager = stager.tiles
    paths = path_manager.get_filenames_from_dir('staged')
    config = stager.config.input_config
    return archive_vector_tiles_parsl(paths, config, archive_dir, ext)
