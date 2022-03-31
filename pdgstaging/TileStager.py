import uuid
import os
from datetime import datetime
import logging
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import box

from . import ConfigManager
from . import TilePathManager


logger = logging.getLogger(__name__)


class TileStager():
    """
        Prepares vector files containing polygon geometries for the tiling
        process. Splits input data into tile-specific files according
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
        config
    ):
        """
            Initialize the TileStager object.
        """

        self.config = ConfigManager(config)
        self.tiles = TilePathManager(**self.config.get_path_manager_config())

        # Configured names of properties that will be added to each polygon
        # during either staging or rasterization
        self.props = self.config.props

        # Create tiles for the maximum z-level configured
        self.z_level = self.config.get_max_z()

        # Vectorize some functions
        self.which_tiles = np.vectorize(
            self.which_tile, otypes=[self.tiles.tile_type]
        )
        self.get_all_tile_properties = np.vectorize(
            self.get_tile_properties, otypes=[dict]
        )

    def stage_all(self):
        """
            Create tiles using all of the available shapefiles
        """

        overall_start_time = datetime.now()
        input_paths = self.tiles.get_filenames_from_dir('input')
        num_paths = len(input_paths)

        if num_paths == 0:
            logger.error(
                'No vector files found for staging.'
            )
            return

        logger.info(
            f'Begin staging {num_paths} input vector files. '
        )

        for path in input_paths:
            self.stage(path)

        # Calculate the total time to stage all files
        total_time = datetime.now() - overall_start_time
        avg_time = total_time / num_paths
        logger.info(
            f'Staged {num_paths} files in {total_time} '
            f'({avg_time} per file on average)'
        )

    def stage(self, path):
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
        try:
            gdf = gpd.read_file(input_path)
        except FileNotFoundError:
            gdf = None
            logger.warning(f'{input_path} not found. It will be skipped.')
        logger.info(
            f'Read in {input_path} in {(datetime.now() - start_time)}'
        )

        # Check that none of the existing properties match the configured
        # properties that will be created. Finoa (used by geopandas) gives
        # error if properties are duplicated regardless of case.
        to_create = self.props.values()
        existing = gdf.columns.values
        duplicated = [b for b in to_create if b.lower() in (a.lower()
                                                            for a in existing)]

        if len(duplicated) > 0:
            error_msg = 'The following properties already exist in the '
            error_msg += f'input vector file: {duplicated}'
            error_msg += '\nPlease remove them or change the configured '
            error_msg += 'property names.'
            logger.error(error_msg)
            raise ValueError(error_msg)

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

        input_crs = self.config.get('input_crs')
        output_crs = self.tiles.crs

        # Set the input CRS if it hasn't been set
        if input_crs:
            gdf.set_crs(
                input_crs, inplace=True, allow_override=True
            )
        # Re-project the geoms
        if output_crs:
            gdf.to_crs(output_crs, inplace=True)

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
        tolerance = self.config.get('simplify_tolerance')
        if tolerance is not None:
            gdf['geometry'] = gdf['geometry'].simplify(tolerance)
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

        props = self.props

        # Get the centroids and area, identify the tile for each polygon.
        # Calculating these properties will give a warning when the CRS is not
        # a projected CRS. Since the geometries are small, the resulting
        # numbers should be accurate enough.

        # Catch the UserWarning that is raised if centroid is calculated
        # using the a non-projected coordinate system.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            centroids = gdf.centroid
            gdf[props['area']] = gdf.area

        gdf[props['centroid_tile']] = self.which_tiles(
            centroids.x, centroids.y)

        gdf[props['centroid_x']] = centroids.x
        gdf[props['centroid_y']] = centroids.y
        # Add the original file name and an identifier for each polygon
        gdf[props['filename']] = path
        gdf[props['identifier']] = [
            str(uuid.uuid4()) for _ in range(len(gdf.index))
        ]
        gdf = self.assign_tile(gdf)
        centroid_only = gdf[props['tile']] == gdf[props['centroid_tile']]
        gdf[props['centroid_within_tile']] = centroid_only

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
        gridded_gdf = gridded_gdf.rename(
            columns={'index_right': self.props['tile']})

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

        for tile in self.tiles.tms.tiles(
                west, south, east, north, self.z_level):
            tile_bb = self.tiles.get_bounding_box(tile)
            grid_cell_geoms.append(box(
                maxy=tile_bb['top'],
                miny=tile_bb['bottom'],
                maxx=tile_bb['right'],
                minx=tile_bb['left']))
            tiles.append(tile)

        grid = gpd.GeoDataFrame({'geometry': grid_cell_geoms, 'tile': tiles})
        grid = grid.set_index('tile')
        grid = grid.set_crs(self.tiles.crs)

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
        return self.tiles.tms.tile(x, y, self.z_level)

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

        for tile, data in gdf.groupby(self.props['tile']):

            # Get the tile path
            tile_path = self.tiles.path_from_tile(tile, base_dir='staged')
            self.tiles.create_dirs(tile_path)

            # Track the start time, the tile, and the number of vectors
            start_time = datetime.now()
            logger.info(
                f'Saving {len(data.index)} vectors to tile {tile_path}'
            )

            self.tiles.create_dirs(tile_path)

            # Save a copy of the tile column for the summary
            tiles = data[self.props['tile']].copy()

            # Tile must be a string for saving as attribute
            data[self.props['tile']] = data[self.props['tile']].astype('str')
            tile_strings = data[self.props['tile']].astype('str')
            data[self.props['centroid_tile']] = tile_strings

            # mode = 'a' will append data to the file if it already exists.
            mode = 'w'
            if os.path.isfile(tile_path):
                mode = 'a'
            # Ignore the FutureWarning that is raised because of geopandas
            # issue 2347
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', FutureWarning)
                data.to_file(tile_path, mode=mode)

            # Track the end time, the total time, and the number of vectors
            logger.info(
                f'Saved {tile_path} in {datetime.now() - start_time}'
            )

            # Record what was saved
            data[self.props['tile']] = tiles
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
        bounds = self.tiles.get_bounding_box(tile)
        return {
            'tile_x': tile.x,
            'tile_y': tile.y,
            'tile_z': tile.z,
            'tile_top': bounds['top'],
            'tile_right': bounds['right'],
            'tile_bottom': bounds['bottom'],
            'tile_left': bounds['left']
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

        prop_file = self.props['filename']
        prop_tile = self.props['tile']
        prop_area = self.props['area']

        # Since summarized is called for each tile that is created, grouping
        # should create 1 row of data to add to the summary dataframe
        gdf_grouped = gdf.groupby([prop_file, prop_tile], as_index=False)
        gdf_summary = gdf_grouped.agg(
            num_polygons=('geometry', 'count'),
            area_polygons=(prop_area, 'sum')
        )
        gdf_summary = gdf_summary.rename(
            columns={
                prop_file: 'filename',
                prop_tile: 'tile'
            })
        # Add the date time that the tile was saved
        gdf_summary['datetime'] = datetime.now()

        tiles = gdf_summary['tile']
        tile_props = self.get_all_tile_properties(tiles)

        gdf_summary = pd.concat(
            [gdf_summary, pd.DataFrame(list(tile_props))],
            axis=1
        )
        gdf_summary['tms_identifier'] = self.tiles.tms_id
        gdf_summary = gdf_summary.drop(columns=['tile'])

        # Save the summary to a file
        summary_path = self.config.get('filename_staging_summary')
        header = False
        mode = 'a'
        if not os.path.isfile(summary_path):
            header = True
            mode = 'w'
        gdf_summary.to_csv(summary_path, mode=mode,
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
