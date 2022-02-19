"""
This module comprises methods that are used in multiple steps of the tile
generation process.
"""

import os
from shapely.geometry import Polygon


def get_tile_path(
        prefix='',
        tms=None,
        tile=None,
        path_structure=['z', 'x', 'y'],
        ext='.png'):
    """
    Creates path to where tiles should be saved that follows a standardized
    directory structure

    Attributes
    -----------
    prefix : str, optional
        Any string to add as a prefix to the start of the path. e.g. a base
        directory path.
    tms : morecantile.TileMatrixSet, required
        Tile Matrix Set that the tiles use.
    tile :  morecantile.Tile, required
        The TMS tile that this path is for
    path_structure : list, optional
        A list of strings that represent the directory structure of last
        segment of the path that uses the x (TileCol), y (TileRow), and z
        (TileMatrix) integer values of the tile. By default, the path will be
        in the format of {TileMatrix}/{TileCol}/{TileRow}.ext, configured as
        ['z', 'x', 'y'].
    ext : str, optional
        The extension for the tile image. Defaults to '.png'.

    Returns
    -------
    path : str
        A string with a path to where the tile image should be saved
    """

    if (tile is None) or (tms is None):
        raise TypeError('Both a tile and tms is required to create a path')

    # Iterate through the path_structure array to create the path
    tile_path = [str(tile.__getattribute__(i)) for i in path_structure]
    tile_path = os.path.join(*tile_path)

    return os.path.join(
        prefix,
        tms.identifier,
        tile_path + ext
    )


def polygon_from_bb(north, east, south, west):
    """
        Create a Shapely polygon object from bounding box coordinates
        Parameters ---------- north : float
            North coordinate of bounding box
        east : float
            East coordinate of bounding box
        south : float
            South coordinate of bounding box
        west : float
            West coordinate of bounding box

        Returns
        -------
        polygon : Shapely Polygon
            Shapely polygon object
    """
    return Polygon([[east, north], [west, north],
                    [west, south], [east, south]])
