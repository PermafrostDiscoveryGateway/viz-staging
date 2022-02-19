import uuid
import os
from datetime import datetime
import logging

import numpy as np
import pandas as pd
import geopandas as gpd
import morecantile

from common import get_tile_path, polygon_from_bb

logger = logging.getLogger(__name__)


class TileStager():
    """
        Prepares vector files containing polygon geometries for the tiling
        process. Splits input data into smaller tile-specific files according
        to the given OGC TileMatrixSet ("TMS").

        The following steps are performed:
            1. Projects the input vector file to the CRS of the TMS. (Also sets
               the CRS initially if needed.)
            2. Simplifies the polygons using the Douglas-Peucker algorithm
                with the specified tolerance.
            3. Adds properties to each vector including centroid coordinates,
               area, and a UUID.
            4. Identifies which tile in the given TMS the polygon belongs to.
            5. Creates a new vector file that contain only the polygons that
                belong to the identified tile.
    """

    def __init__(
        self,
        input_dir=None,
        input_ext='.shp',
        output_dir=None,
        output_ext='.shp',
        input_crs=None,
        simplify_tolerance=0.000001,
        tms_identifier='WorldCRS84Quad',
        zoom_level=14,
        tile_path_structure=['z', 'x', 'y'],
        summary_path='staging-summary.csv',
    ):
        """
            Initialize the TileStager object.

            Parameters
            ----------
            input_dir : str
                The directory containing all of the input vector files to
                process. The directory is searched recursively for all files
                with the given input extension.
            input_ext : str
                The extension of the input vector files
            output_dir : str
                The directory to save the output vector files. The output files
                will be saved in a subdirectory named after the TMS.
            output_ext : str
                The extension of the output vector files (e.g. '.shp',
                '.geojson', '.gpkg')
            input_crs : str
                The CRS of the input vector files, if the CRS is not specified
                the file itself. Setting an input_crs does NOT re-project the
                data.
            simplify_tolerance : float
                The tolerance to use when simplifying the geometries
            tms_identifier : str
                The name of the TMS to use. Must be one that is supported by
                the morecantile package.
            zoom_level : int
                The zoom level (or 'TileMatrix index') to create tiles for.
                This should be the highest zoom level that will be used in the
                resulting raster and 3D tiles.
            tile_path_structure : list
                A list of strings that represent the format of last segment of
                the path that uses the x (TileCol), y (TileRow), and z
                (TileMatrix) integer values of the tile. By default, the path
                will be in the format of {TileMatrix}/{TileCol}/{TileRow}.ext,
                configured as ['z', 'x', 'y'].
            summary_path : str
                The path to save a CSV file to that contains a summary of the
                tiles created during the staging process. If a summary exists
                from a previous run, new summary rows will be appended to the
                existing file.
        """
        self.input_dir = input_dir
        self.input_ext = input_ext
        self.output_dir = output_dir
        self.output_ext = output_ext
        self.input_crs = input_crs
        self.simplify_tolerance = simplify_tolerance
        self.tms_identifier = tms_identifier
        self.zoom_level = zoom_level
        self.tile_path_structure = tile_path_structure
        self.summary_path = summary_path

        # Create the TileMatrixSet for the output tiles
        self.tms = morecantile.tms.get(tms_identifier)
        self.output_crs = self.tms.crs

        # Vectorize some functions
        self.which_tiles = np.vectorize(
            self.which_tile, otypes=[morecantile.commons.Tile]
        )
        self.get_all_tile_properties = np.vectorize(
            self.get_tile_properties, otypes=[dict]
        )

    def stage_all(self):
        """
            Create tiles using all of the available shapefiles
        """
        overall_start_time = datetime.now()
        # Create the output directory if it doesn't exist
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        self.get_input_paths()
        for path in self.input_paths:
            self.stage_file(path)
        # Calculate the total time to stage all files
        total_time = datetime.now() - overall_start_time
        total_files = len(self.input_paths)
        avg_time = total_time / total_files
        logger.info(
            f'Staged {total_files} files in {total_time} '
            f'({avg_time} per file on average)'
        )

    def stage_file(self, path):
        """
            Process and create tiles for a single vector file

            Parameters
            ----------
            path : str
                The path to the vector file to process and create tiles for.
        """
        gdf = self.get_data(path)
        if (gdf is not None) and (len(gdf) > 0):
            gdf = self.set_crs(gdf)
            gdf = self.simplify_geoms(gdf)
            gdf = self.add_properties(gdf, path)
            self.save_tiles(gdf)
        else:
            logger.warning(f'No features in {path}')

    def get_tile_path(self, tile):
        """
            Returns the path to save a tile to.

            Parameters
            ----------
            tile : morecantile.commons.Tile
                The tile to get the path for.

            Returns
            -------
            str
                The path to save the tile to.
        """
        return get_tile_path(
            prefix=self.output_dir,
            tms=self.tms,
            tile=tile,
            path_structure=self.tile_path_structure,
            ext=self.output_ext
        )

    def get_input_paths(self):
        """
            Get the paths for all of the input data files. Looks recursively
            through the input directory for all files with the given input
            extension.

            Returns
            -------
            list
                A list of paths to all of the files found in the input_dir with
                the given input_ext.
        """
        start_time = datetime.now()
        input_paths = []
        for root, dirs, files in os.walk(self.input_dir):
            for file in files:
                if file.endswith(self.input_ext):
                    path = os.path.join(root, file)
                    input_paths.append(path)
        logger.info(
            f'Found {len(input_paths)} vector files in directory: '
            f'{self.input_dir} in {datetime.now() - start_time}'
        )
        self.input_paths = input_paths
        return input_paths

    def get_data(self, input_path=None):
        """
            Read a vector file as a GeoDataFrame

            Parameters
            ----------
            input_path : str
                The path to the vector file to read

            Returns
            -------
            GeoDataFrame
                The GeoDataFrame representation of the vector file.
        """

        start_time = datetime.now()
        logger.info(f'Reading vector file: {input_path}')
        gdf = gpd.read_file(input_path)
        logger.info(
            f'Read in {input_path} in {(datetime.now() - start_time)}'
        )
        return gdf

    def set_crs(self, gdf):
        """
            Set the CRS of the GeoDataFrame to the input CRS, if there is one.
            Re-project to the CRS of the TMS.

            Parameters
            ----------
            gdf : GeoDataFrame
                The GeoDataFrame to set the CRS of

            Returns
            -------
            GeoDataFrame
                The GeoDataFrame with the updated CRS

        """
        start_time = datetime.now()
        # Set the input CRS if it hasn't been set
        if self.input_crs:
            gdf.set_crs(
                self.input_crs, inplace=True, allow_override=True
            )
        # Re-project the geoms
        if self.output_crs:
            gdf.to_crs(self.output_crs, inplace=True)
        logger.info(
            f'Re-projected {len(gdf.index)} polygons in '
            f'{datetime.now() - start_time}'
        )
        return gdf

    def simplify_geoms(self, gdf):
        """
            Simplify all geometries in a GeoDataFrame using the Douglas-Peucker
            algorithm

            Parameters
            ----------
            gdf : GeoDataFrame
                The GeoDataFrame to simplify

            Returns
            -------
            GeoDataFrame
                The GeoDataFrame with the simplified geometries
        """
        start_time = datetime.now()
        gdf['geometry'] = gdf['geometry'].simplify(self.simplify_tolerance)
        logger.info(
            f'Simplified {len(gdf.index)} polygons in '
            f'{datetime.now() - start_time}'
        )
        return gdf

    def add_properties(self, gdf, path=''):
        """
            Add area, centroid coordinates, filename, uuid to each polygon in
            the GeoDataFrame. Identify the tile that each polygon belongs to.

            Parameters
            ----------
            gdf : GeoDataFrame
                The GeoDataFrame to add the tile properties to
            path : str
                The path to the vector file that the GeoDataFrame originated
                from (used for the filename)

            Returns
            -------
            GeoDataFrame
                The GeoDataFrame with the tile properties added to each row
        """

        start_time = datetime.now()

        # Get the centroids and area, identify the tile for each polygon.
        # Calculating these properties will give a warning when the CRS is not
        # a projected CRS. Since the geometries are small, the resulting
        # numbers should be accurate enough.
        centroids = gdf.centroid
        gdf['centroid_tile'] = self.which_tiles(centroids.x, centroids.y)
        gdf['area'] = gdf.area
        gdf['centroid_x'] = centroids.x
        gdf['centroid_y'] = centroids.y
        # Add the original file name and an identifier for each polygon
        gdf['filename'] = path
        gdf['identifier'] = [
            str(uuid.uuid4()) for _ in range(len(gdf.index))
        ]
        gdf = self.assign_tile(gdf)
        gdf['centroid_within_tile'] = gdf['tile'] == gdf['centroid_tile']

        logger.info(
            f'Added properties for {len(gdf.index)} vectors in '
            f'{datetime.now() - start_time} from file {path}'
        )
        return gdf

    def assign_tile(self, gdf):
        """
            Assign each polygon in the GeoDataFrame to the tile that it
            intersects. When a polygon intersects multiple tiles, it will be
            duplicated.

            Parameters
            ----------
            gdf : GeoDataFrame
                The GeoDataFrame to assign tiles to
        """
        grid = self.make_tms_grid(gdf)
        gridded_gdf = gdf.sjoin(grid, how='left', predicate='intersects')
        gridded_gdf = gridded_gdf.rename(columns={'index_right': 'tile'})

        return gridded_gdf

    def make_tms_grid(self, gdf):
        """
            Create a GeoDataFrame of rectangles that creates a TMS grid that
            covers the given GeoDataFrame.

            Parameters
            ----------
            gdf : GeoDataFrame
                The GeoDataFrame from which a bounding box will be calculated
                and used for the extent of the TMS grid.
        """
        bbox = gdf.total_bounds
        west, south, east, north = bbox

        grid_cell_geoms = []
        tiles = []

        for tile in self.tms.tiles(west, south, east, north, self.zoom_level):
            tile_bb = self.tms.bounds(tile)
            poly = polygon_from_bb(
                tile_bb.top, tile_bb.right, tile_bb.bottom, tile_bb.left
            )
            grid_cell_geoms.append(poly)
            tiles.append(tile)

        grid = gpd.GeoDataFrame({'geometry': grid_cell_geoms, 'tile': tiles})
        grid = grid.set_index('tile')
        grid = grid.set_crs(self.tms.crs)

        return grid

    def which_tile(self, x, y):
        """
            Identify which tile a pair of x and y coordinates belongs to. x, y
            coordinates must be in the same CRS as the TMS

            Parameters
            ----------
            x : float
                The x coordinate to identify the tile for
            y : float
                The y coordinate to identify the tile for
        """
        return self.tms.tile(x, y, self.zoom_level)

    def save_tiles(self, gdf=None):
        """
            Given a processed GeoDataFrame, save vectors into tiled vector
            files.

            Parameters
            ----------
            gdf : GeoDataFrame
                The GeoDataFrame to save
        """

        if gdf is None:
            # TODO give warning
            return None

        for tile, data in gdf.groupby('tile'):

            # Get the tile path
            tile_path = self.get_tile_path(tile)

            # Track the start time, the tile, and the number of vectors
            start_time = datetime.now()
            logger.info(
                f'Saving {len(data.index)} vectors to tile {tile_path}'
            )

            tile_dir = os.path.dirname(tile_path)

            # Create the directory if it doesn't exist
            if not os.path.exists(tile_dir):
                os.makedirs(tile_dir, exist_ok=True)

            # Save a copy of the tile column for the summary
            tiles = data['tile'].copy()

            # Tile must be a string for saving as attribute
            data['tile'] = data['tile'].astype('str')
            data['centroid_tile'] = data['tile'].astype('str')

            # mode = 'a' will append data to the file if it already exists.
            mode = 'w'
            if os.path.isfile(tile_path):
                mode = 'a'
            data.to_file(tile_path, mode=mode)

            # Track the end time, the total time, and the number of vectors
            logger.info(
                f'Saved {tile_path} in {datetime.now() - start_time}'
            )

            # Record what was saved
            data['tile'] = tiles
            self.summarize(data)

    def get_tile_properties(self, tile):
        """
            For a given TMS tile, return certain properties of that tile as a
            dict for summary purposes

            Parameters
            ----------
            tile : morecantile.commons.Tile
                The TMS tile to get properties for
        """
        bounds = self.tms.bounds(tile)
        return {
            'tile_x': tile.x,
            'tile_y': tile.y,
            'tile_z': tile.z,
            'tile_top': bounds.top,
            'tile_right': bounds.right,
            'tile_bottom': bounds.bottom,
            'tile_left': bounds.left
        }

    def summarize(self, gdf=None):
        """
            For a given file, count how many vectors there are per tile, the
            area of vectors per tile, and get information about the tile itself
            (bounds, x, y, z values). Append these details to the summary
            dataframe.

            Parameters
            ----------
            gdf : GeoDataFrame
                The GeoDataFrame to summarize
        """

        if gdf is None:
            return None

        # Log the summary event, including how long it takes
        start_time = datetime.now()

        # Since summarized is called for each tile that is created, grouping
        # should create 1 row of data to add to the summary dataframe
        gdf_summary = gdf.groupby(['filename', 'tile'], as_index=False).agg(
            num_polygons=('geometry', 'count'),
            area_polygons=('area', 'sum')
        )
        # Add the date time that the tile was saved
        gdf_summary['datetime'] = datetime.now()

        tiles = gdf_summary['tile']
        tile_props = self.get_all_tile_properties(tiles)

        gdf_summary = pd.concat(
            [gdf_summary, pd.DataFrame(list(tile_props))],
            axis=1
        )
        gdf_summary['tms_identifier'] = self.tms_identifier
        gdf_summary = gdf_summary.drop(columns=['tile'])

        # Save the summary to a file
        header = False
        mode = 'a'
        if not os.path.isfile(self.summary_path):
            header = True
            mode = 'w'
        gdf_summary.to_csv(self.summary_path, mode=mode,
                           index=False, header=header)

        # Log the total time to create the summary
        logger.info(
            'Summarized Tile(z={}, x={}, y={}) in {}'
            .format(
                tile_props[0]['tile_z'],
                tile_props[0]['tile_x'],
                tile_props[0]['tile_y'],
                (datetime.now() - start_time)
            )
        )
