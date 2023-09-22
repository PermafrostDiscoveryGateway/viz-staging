import logging
from . import logging_config
import os
import uuid
import warnings
from datetime import datetime

import geopandas as gpd
import numpy as np
import pandas as pd
from filelock import FileLock

from . import ConfigManager, TilePathManager, TMSGrid
from .Deduplicator import clip_gdf

logger = logging_config.logger


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
        logger.info(f"Staging file {path}")
        # Remove any geometries that are not polygons
        gdf = gdf[gdf.geometry.type == 'Polygon']

        if (gdf is not None) and (len(gdf) > 0):
            gdf = self.simplify_geoms(gdf)
            # clip to footprint before CRS of IWP data is transformed
            # to EPSG:4326
            gdf = self.clip_to_footprint(gdf, path)
            gdf = self.set_crs(gdf)
            self.grid = self.make_tms_grid(gdf)
            gdf = self.add_properties(gdf, path)
            self.save_tiles(gdf)
        else:
            logger.warning(f'No features in {path}')
    
    def clip_to_footprint(self, gdf, path):
        """
            If the config is set to clip_to_footprint=True,
            and config is set to deduplicate at any step in the workflow,
            find the footprint file associated with the gdf, 
            determine which polygons fall outside the footprint,
            and label the polygons as True or False in a new column
            that defualts to 'staging_duplicated'.

            Parameters
            ----------
            gdf : GeoDataFrame
                The GDF undergoing staging.
            
            path: string
                Path to the input file. This is used to retrieve the footprint.

            Returns
            -------
            The GeodataFrame with polygons that fall outside
            the footprint labeled as duplicates.

        """

        # check if the config is set to clip to footprint
        clip_to_footprint = self.config.get('deduplicate_clip_to_footprint')
        logger.info(f"clip_to_footprint is {clip_to_footprint}")
        # check if the config is set to label duplicates
        dedup = self.config.get('deduplicate_method')
        # if the config is set to do so, clip to footprint
        if clip_to_footprint == True and dedup is not None:
            logger.info(f' Starting clipping_to_footprint() for file {path}.')
            # pull in footprint as a gdf called fp
            fp_path = self.config.footprint_path_from_input(path, check_exists=True)
            fp = self.get_data(fp_path)
            logger.info(f' Checking CRSs of polygons and footprint.')

            iwp_crs = gdf.crs
            fp_crs = fp.crs

            if iwp_crs == fp_crs:
                logger.info(f" CRSs match. They are both {iwp_crs}.")
            else:
                logger.info(f" CRSs do not match.\n IWP's CRS is {iwp_crs}."
                            f" Footprint's CRS is {fp_crs}.")
                # transform the footprint to the CRS of the polygon data
                fp.to_crs(iwp_crs, inplace = True)
                # check again
                fp_crs_transformed = fp.crs
                if iwp_crs == fp_crs_transformed:
                    logger.info("Footprint CRS has been transformed to CRS of polygons.")
                else:
                    logger.error("Failed to transform footprint CRS to CRS of polygons.")
                    return

            # determine if polygons fall within or outside the footprint
            # first retrieve the name of the column we will use to label duplicates
            prop_duplicated = self.config.polygon_prop('duplicated')
            gdf_with_labels = clip_gdf(
                gdf = gdf.copy(), # the gdf to clip
                boundary = fp.copy(), # the footprint
                method = 'intersects',
                prop_duplicated = prop_duplicated
            )

            return gdf_with_labels
        else:
            logger.info(f" Either clip_to_footprint was set to False, or config"
                        f" was not set to deduplicate at any step. Returning original GDF"
                        f" without clipping to footprint.")
            return gdf
    
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
                The GeoDataFrame to set the CRS of the Tile Matrix Set

            Returns
            -------
            GeoDataFrame
                The GeoDataFrame with the updated CRS

        """
        start_time = datetime.now()

        input_crs = self.config.get('input_crs')
        # assign the CRS of the TMS to `output_crs`
        output_crs = self.tiles.crs

        # log the CRS that already exists in the input data
        crs_in_input = gdf.crs

        if crs_in_input is not None:
            logger.info(f"CRS of input data is {crs_in_input}.\nIf input_crs is set in config," 
                        f" setting CRS to that.")
        else: 
            logger.info(f"No CRS set in input data. Setting to input_crs specified in config.")

        if input_crs:
            gdf.set_crs(
                input_crs, inplace=True, allow_override=True
            )
        # Re-project the geoms to the CRS of the TMS
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
            tile_strings = data[self.props['centroid_tile']].astype('str')
            data[self.props['centroid_tile']] = tile_strings

            dedup_method = self.config.get_deduplication_method()
            if os.path.isfile(tile_path):
                if dedup_method is not None:
                    # If the file exists and config is set to deduplicate
                    # at any step, then open the file, append the new data, 
                    # and identify duplicates. 
                    # (If deduplicating by footprint:
                    # duplicates where polygons fall outside
                    # the footprint were already labeled earlier).
                    # Then remove all the duplicated data if the config is 
                    # set to remove duplicates during staging.
                    logger.info(f"Tile exists and dedup is set to occur at some step,"
                                f" so executing `combine_and_deduplicate()`")
                    data = self.combine_and_deduplicate(data, tile_path)

                    mode = 'w'
                    # Overwrite existing file.
                    self.save_new_tile(data = data,
                                        tile_path = tile_path,
                                        mode = mode,
                                        start_time = start_time,
                                        lock = lock)
                else:
                    # If the file exists and config is not set to deduplicate
                    # at any step, simply append the polygons to the existing file.
                    # no polygons will be labeled as duplicates or not.
                    # If deduplicating by footprint: 
                    # neither file has been clipped to footprint
                    logger.info(f"Tile exists but dedup is not set to occur, so just "
                                f"appending the polygons.")
                    
                    # Append to existing tile
                    mode = 'a'
                    self.save_new_tile(data = data,
                                        tile_path = tile_path,
                                        mode = mode,
                                        start_time = start_time,
                                        lock = lock)
            else:
                if dedup_method is not None:
                    if self.config.deduplicate_at('staging'):
                        # If the file does not yet exist and the config is set to deduplicate
                        # at staging, just remove the polygons that were already labeled as 
                        # duplicates because they fell outside the footprint,
                        # and save the tile as a new tile.

                        # dedup_start_time = datetime.now()
                        # logger.info(f'Starting deduplication in tile {tile_path} with {len(gdf)}'
                        #             'polygons.')
                        # dedup_config = self.config.get_deduplication_config(gdf)
                        # gdf = dedup_method(gdf, **dedup_config)
                        logger.info(f"Tile does not yet exist and config is set to deduplicate "
                                    f"at staging, so removing polygons that fell outside the footprint "
                                    f"if deduplicating by footpint, and removing overlapping polygons\n"
                                    f"if deduplicating by neighbor.")
                        # retreive the name of the duplicated column from config
                        prop_duplicated = self.config.polygon_prop('duplicated')
                        if prop_duplicated in data.columns:
                            data = data[~data[prop_duplicated]]

                        mode = 'w'
                        self.save_new_tile(data = data,
                                           tile_path = tile_path,
                                           mode = mode,
                                           start_time = start_time,
                                           lock = lock)
                    else:
                        # If the file does not yet exist and the config is set to deduplicate
                        # at a step after staging, just save as a new tile.
                        # If deduplicating by footprint:
                        # The gdf already has polys labeled as duplicates if they fall outside 
                        # the footprint, and the next time a poly in this tile is produced, 
                        # it will be checked for duplicates where the two footprints overlap.
                        # If deduplicating by neighbor:
                        # the prop_duplicated col will be added, 
                        # with all values set to false so all staged files have same properties.
                        # This will be overwritten if this file overlaps with others later,
                        # with combine_and_deduplicate().
                        logger.info(f"Tile does not yet exist and config is set to deduplicate at a step "
                                    f"after staging, so just saving the new tile." 
                                    f"\nIf deduplicating by footprint:\n"
                                    f"Duplicates from `clip_gdf` were identified.")

                        dedup_method = self.config.get_deduplication_method()
                        prop_duplicated = self.config.polygon_prop('duplicated')
                        logger.info(f"Checking for presence of {prop_duplicated} col with dedup method set"
                                    f" to {dedup_method}")
                        if dedup_method is not None and prop_duplicated not in data.columns:
                            logger.info(f"Adding {prop_duplicated} column with False values to "
                                        f"{data['staging_tile']}\nbecause property did not "
                                        f"already exist.")
                            data[prop_duplicated] = False
                        else:
                            logger.info(f"File {data['staging_tile']} that did not yet exist did have "
                                        f" {prop_duplicated} col before saving.")
                        
                        mode = 'w'
                        self.save_new_tile(data = data,
                                           tile_path = tile_path,
                                           mode = mode,
                                           start_time = start_time,
                                           lock = lock)
                else:
                    # If the tile does not yet exist and the config is not set to deduplicate
                    # at any step, just save as a new tile. 
                    # If deduplicating by footprint:
                    # No duplicates were labeled earlier either,
                    # because the workflow did not check if any fell outside the footprint.
                    logger.info(f"Tile does not yet exist and config is not set to deduplicate, so just "
                                f"saving the new tile.\nIf deduplicating by footprint:\n"
                                f"Duplicates from `clip_gdf` were not identified.")
                    
                    mode = 'w'
                    self.save_new_tile(data = data,
                                       tile_path = tile_path,
                                       mode = mode,
                                       start_time = start_time,
                                       lock = lock)

    def save_new_tile(self, data, tile_path, mode, start_time, lock):
        """
            Save data as a new tile, if the tile does not already exist.
            Adds a row to staging_summary.csv.

            Parameters
            ----------
            data: GeoDataframe
                Subset of the original input GDF, which was split by self.props['tile']
                at the start of save_tiles().
            
            tile_path: str
                Path to tile file, considering the standardized directory
                structure set in the Path Manager.

            mode: str
                Specifies the way in which a file should be opened. "w" overwrites the
                contents of the file if it already exists, and creates a new file if it 
                doesn't exist. "a" appends data to the end of the file if it already exists,
                and creates a new file if it doesn't exist. 

            start_time: datetime
                The datetime when save_tiles() began (after the file was locked).

            lock: lock file object instance
                Ensures that only 1 process can access the file while it is being written.

            Returns
            -------
            Nothing. Saves new tile in `staged` directory.

        """

        try: 
            # Ignore the FutureWarning raised from geopandas issue 2347
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', FutureWarning)
                data.to_file(tile_path, mode=mode)

            # convert each tile from string format to morecantile format 
            # so it can be added to summary
            # first create series of tiles in str format
            tiles_str = data[self.props['tile']].copy()

            tiles_morecantile = [self.tiles.tile_from_str(tile) for tile in tiles_str]

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
            GeoDataFrame, label polygons as duplicates using one of 
            the deduplication methods, and remove tiles labeled as 
            duplicates if the config is set to do so at staging.

            Parameters
            ----------
            gdf : GeoDataFrame
                The GeoDataFrame with new data to add to the existing tile
            tile_path : str
                The path to the existing data

            Returns
            -------
            gdf : GeoDataFrame
                The combined GeoDataFrame without duplicates labeled or with 
                deduplicates labeled, depending on how the config is set.
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
            # if the CRS's don't match, transform the new data to the
            # CRS of the existing data
            # note: unlikely this is executed because both tiles have already
            # been input into set_crs()
            gdf.to_crs(existing_gdf.crs, inplace=True)

        gdf = pd.concat(to_concat, ignore_index=True)
        gdf.reset_index(drop = True, inplace = True)
        dedup_config = self.config.get_deduplication_config(gdf)

        logger.info(
            f'Starting deduplication in tile {tile_path} with {len(gdf)} '
            'polygons.'
        )

        # label duplicates depending on which deduplication type is set in config
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
        logger.info(f'Found {num_found} matching footprints. '
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
