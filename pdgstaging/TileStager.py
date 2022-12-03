import logging
import os
import uuid
import warnings
from datetime import datetime

import geopandas as gpd
import numpy as np
import pandas as pd
from filelock import FileLock

from . import ConfigManager, TilePathManager, TMSGrid

logger = logging.getLogger(__name__)


class TileStager():
    """
        Prepares vector files containing polygon geometries for the tiling
        process. Splits input data into tile-specific files according
        to the given OGC TileMatrixSet ("TMS").

        The following steps are performed:
            1. Removes any non-polygon features.
            2. Projects the input vector file to the CRS of the TMS. (Also sets
               the CRS initially if needed.)
            3. Simplifies the polygons using the Douglas-Peucker algorithm
                with the specified tolerance.
            4. Adds properties to each vector including centroid coordinates,
               area, and a UUID.
            5. Identifies which tile in the given TMS the polygon belongs to.
            6. Creates a new vector file that contain only the polygons that
                belong to the identified tile.
    """

    def __init__(self, config=None, check_footprints=True):
        """
            Initialize the TileStager object.

            Parameters
            ----------
            config : dict
                The configuration dictionary to use. If not provided,
                the default configuration will be used.

            check_footprints : bool
                When True, upon initialization, checks that a footprint file
                can be found for each of the input files when the
                deduplicate_method is set to 'footprints'. If not, raises
                an exception. Default is True.
        """

        self.config = ConfigManager(config)
        self.tiles = TilePathManager(**self.config.get_path_manager_config())

        if(check_footprints and
           self.config.get('deduplicate_method') == 'footprints'):
            logger.info('Checking for footprint files...')
            missing = self.check_footprints()
            num_missing = len(missing)
            if 0 < num_missing < 20:
                logger.warning(
                    f'Missing footprint files for {num_missing} files: '
                    f'{missing}')
            elif num_missing > 20:
                logger.warning(
                    f'Missing footprint files for {num_missing} files: '
                    f'{len(missing)}')
                logger.warning(
                    f'Printing first 30 missing footprint: '
                    f'{missing[0:30]}')

        # Configured names of properties that will be added to each polygon
        # during either staging or rasterization
        self.props = self.config.props

        # Create tiles for the maximum z-level configured
        self.z_level = self.config.get_max_z()

        self.get_all_tile_properties = np.vectorize(
            self.get_tile_properties, otypes=[dict]
        )

    def stage_all(self):
        """
            Process and create tiles from all of the vector files in the
            configured input directory.
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
        # Remove any geometries that are not polygons
        gdf = gdf[gdf.geometry.type == 'Polygon']
        if (gdf is not None) and (len(gdf) > 0):
            gdf = self.simplify_geoms(gdf)
            gdf = self.set_crs(gdf)
            self.grid = self.make_tms_grid(gdf)
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
        num_polygons = len(gdf.index)

        # Get the centroids and area. Calculating these properties will give a
        # warning when the CRS is not a projected CRS. Since the geometries are
        # small, the resulting numbers should be accurate enough.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            centroids = gdf.centroid
            gdf[props['area']] = gdf.area
        gdf[props['centroid_x']] = centroids.x
        gdf[props['centroid_y']] = centroids.y

        # Add source file path, excluding the input dir and any leading slashes
        dir_input = self.config.get('dir_input')
        gdf[props['filename']] = path.removeprefix(dir_input).strip(os.sep)

        # Add identifier
        gdf[props['identifier']] = [
            str(uuid.uuid4()) for _ in range(num_polygons)]

        # Find the tile the polygon falls within
        gdf = self.assign_tile_by_centroid(gdf)
        # Find the tiles the polygon intersects
        gdf = self.assign_tile(gdf)
        centroid_only = gdf[props['tile']] == gdf[props['centroid_tile']]
        gdf[props['centroid_within_tile']] = centroid_only

        # Add the column to flag duplicate polygons. This will be set to True
        # later if duplicates are found.
        dedup_method = self.config.get_deduplication_method()
        if dedup_method is not None:
            # Mark all the polygons as not duplicated
            gdf[self.config.polygon_prop('duplicated')] = False

        logger.info(
            f'Added properties for {num_polygons} vectors in '
            f'{datetime.now() - start_time} from file {path}'
        )
        return gdf

    def assign_tile_by_centroid(self, gdf):
        """
            Assign the tile that the centroid of each polygon falls within.

            Parameters
            ----------
            gdf : GeoDataFrame
                The GeoDataFrame to assign the centroid tile to

            Returns
            -------
            GeoDataFrame
                The GeoDataFrame with the centroid tile assigned to each row
        """
        gdf[self.props['centroid_tile']] = self.grid.tiles_from_xy(
            gdf[self.props['centroid_x']], gdf[self.props['centroid_y']]
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

            Returns
            -------
            GeoDataFrame
                The GeoDataFrame with the tile assigned to each row

        """
        return self.grid.sjoin(
            gdf,
            how='left',
            predicate='intersects',
            as_tile=True)

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
        grid = TMSGrid(
            tms_id=self.tiles.tms_id,
            z=self.z_level,
            bounds=gdf.total_bounds
        )
        grid.TILE_NAME = self.props['tile']
        return grid

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

            # Lock the tile so that multiple processes don't try to write to
            lock = self.__lock_file(tile_path)

            # Track the start time, the tile, and the number of vectors
            start_time = datetime.now()
            logger.info(
                f'Saving {len(data.index)} vectors to tile {tile_path}'
            )

            # Tile must be a string for saving as attribute
            data[self.props['tile']] = data[self.props['tile']].astype('str')
            tile_strings = data[self.props['tile']].astype('str')
            data[self.props['centroid_tile']] = tile_strings

            # Open the file in write mode by default
            mode = 'w'

            if os.path.isfile(tile_path):
                # If the file exists and config is set to deduplicate, then
                # open the file, append the new data, and identify duplicates.
                # Remove the data if the config is set to remove duplicates
                # during staging. Overwrite existing file.
                dedup_method = self.config.get_deduplication_method()
                if dedup_method is not None:
                    data = self.combine_and_deduplicate(data, tile_path)
                # If the file exists and config is not set to deduplicate, then
                # just append the new data to the existing file.
                else:
                    mode = 'a'

            try:
                # Ignore the FutureWarning raised from geopandas issue 2347
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', FutureWarning)
                    data.to_file(tile_path, mode=mode)

                # convert each tile from string format to morecantile format 
                # so it can be added to summary
                # first create series of tiles in str format
                tiles_str = data[self.props['tile']].copy()

                tiles_morecantile = []
                for tile in tiles_str:
                    tile_morecantile = self.tiles.tile_from_str(tile)
                    tiles_morecantile.append(tile_morecantile)                

                # Record what was saved
                data[self.props['tile']] = tiles_morecantile
                self.summarize(data)
            finally:
                # Track the end time, the total time, and the number of vectors
                logger.info(
                    f'Saved {tile_path} in {datetime.now() - start_time}'
                )
                self.__release_file(lock)

    def combine_and_deduplicate(self, gdf, tile_path):
        """
            Combine existing data for a tile with the new data in a
            GeoDataFrame, and deduplicate the result.

            Parameters
            ----------
            gdf : GeoDataFrame
                The GeoDataFrame with new data to add to the existing tile
            tile_path : str
                The path to the existing data

            Returns
            -------
            gdf : GeoDataFrame
                The combined and deduplicated GeoDataFrame
        """

        dedup_start_time = datetime.now()

        dedup_method = self.config.get_deduplication_method()
        existing_gdf = gpd.read_file(tile_path)

        # Projection info can be lost during saving & reopening geopackage
        # files for the CRS used for some TMSs, see details:
        # https://github.com/PermafrostDiscoveryGateway/viz-staging/issues/11
        to_concat = [gdf, existing_gdf]
        num_unique_crs = len({g.crs for g in to_concat})
        if num_unique_crs != 1:
            existing_gdf.to_crs(gdf.crs, inplace=True)

        gdf = pd.concat(to_concat)
        dedup_config = self.config.get_deduplication_config(gdf)
        if dedup_method is None:
            return gdf

        logger.info(
            f'Starting deduplication in tile {tile_path} with {len(gdf)} '
            'polygons.'
        )
        gdf = dedup_method(gdf, **dedup_config)

        # drop duplicated polygons, if config is set to deduplicate here
        if self.config.deduplicate_at('staging'):
            prop_duplicated = self.config.polygon_prop('duplicated')
            if prop_duplicated in gdf.columns:
                gdf = gdf[~gdf[prop_duplicated]]

        logger.info(
            f'Finished deduplication in {datetime.now() - dedup_start_time}'
        )

        return gdf

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

    def check_footprints(self):
        """
        Check if all the footprint files exist.

        Returns
        -------
        list of str
            Returns a list of input filenames that are missing footprint files.
            Returns an empty list if all footprint files exist.
        """
        missing_footprints = []
        matching_footprints = []
        input_paths = self.tiles.get_filenames_from_dir('input')
        for path in input_paths:
            try:
                footprint = self.config.footprint_path_from_input(
                    path, check_exists=True)
            except FileNotFoundError:
                missing_footprints.append(path)
            else:
                matching_footprints.append(footprint)
        num_missing = len(missing_footprints)
        num_found = len(matching_footprints)
        logging.info(f'Found {num_found} matching footprints. '
                     f'{num_missing} missing.')
        return missing_footprints

    def __lock_file(self, path):
        """
            Lock a file for writing.

            Parameters
            ----------
            path : str
                The path to the file to lock

            Returns
            -------
            lock : FileLock
                The lock object
        """
        lock_path = path + '.lock'
        lock = FileLock(lock_path)
        lock.acquire()
        return lock

    def __release_file(self, lock):
        """
            Release a file lock. Remove the lock file.

            Parameters
            ----------
            lock : FileLock
                The lock to release
        """
        lock.release()
        if os.path.exists(lock.lock_file):
            os.remove(lock.lock_file)
