import logging
import os
import uuid
import warnings
import gc
from datetime import datetime
from .Deduplicator import deduplicate_neighbors, deduplicate_by_footprint, clip_gdf

import geopandas as gpd
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from filelock import FileLock
from typing import Optional

from . import TilePathManager, TMSGrid


class TileStager:
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

    def __init__(
        self,
        tiles: Optional[TilePathManager] = None,
        props=[],
        max_z_level=13,
        tms_id="WGS1984Quad",
        path_structure=["style", "tms", "z", "x", "y"],
        base_dirs={},
    ):
        """
        Initialize the TileStager object.

        Parameters
        ----------
        tiles : TilePathManager, optional
            An instance of TilePathManager (defaults to None)

        props : list
            A list of all the properties associated with the tiling data

        max_z_level : int
            The max z level for tiling the input data

        tms_id : str
            The tile matrix set id for constructing the tile path manager object

        path_structure : list
            The path structure used for generating the tiling output

        base_dirs : dict
            Directory structure and extenstions for storing the output
        """
        self.logger = logging.getLogger(__name__)

        # Configured names of properties that will be added to each polygon
        # during either staging or rasterization
        self.props = props

        # Create tiles for the maximum z-level configured
        self.z_level = max_z_level

        if tiles is None:
            self.tiles = self.set_tiling_config(tms_id, path_structure, base_dirs)
        else:
            self.tiles = tiles

        staged_root = self.tiles.base_dirs["staged"]["path"]
        summary_root = os.path.dirname(staged_root)
        self.summary_path = os.path.join(summary_root, "staging_summary.csv")

        self.get_all_tile_properties = np.vectorize(
            self.get_tile_properties, otypes=[dict]
        )

    def stage_all(self):
        """
        Process and create tiles from all of the vector files in the
        configured input directory.
        """

        overall_start_time = datetime.now()
        input_paths = self.tiles.get_filenames_from_dir("input")
        num_paths = len(input_paths)

        if num_paths == 0:
            self.logger.error("No vector files found for staging.")
            return

        self.logger.info(f"Begin staging {num_paths} input vector files. ")

        for i, path in enumerate(input_paths):
            self.stage(path)
            # Periodic garbage collection to prevent memory buildup
            if (i + 1) % 50 == 0:
                gc.collect()
                self.logger.debug(f"Garbage collected after {i + 1} files")
        try:
            self.csv_to_parquet()
        except Exception:
            self.logger.exception("Failed to convert CSVs to Parquet")

        # Calculate the total time to stage all files
        total_time = datetime.now() - overall_start_time
        avg_time = total_time / num_paths
        self.logger.info(
            f"Staged {num_paths} files in {total_time} "
            f"({avg_time} per file on average)"
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
        self.logger.info(f"Staging file {path}")
        # Remove any geometries that are not polygons
        gdf = gdf[gdf.geometry.type.isin(["Polygon", "MultiPolygon"])]

        if (gdf is not None) and (len(gdf) > 0):
            gdf = self.simplify_geoms(gdf)
            # If clipping to footprint, do so before CRS is transformed
            gdf = self.clip_to_footprint(gdf, path)
            gdf = self.set_crs(gdf)
            self.grid = self.make_tms_grid(gdf)
            gdf = self.add_properties(gdf, path)
            self.save_tiles(gdf)
            # Clean up to free memory
            del gdf
            self.grid = None
        else:
            self.logger.warning(f"No features in {path}")
            # Clean up empty GeoDataFrame
            if gdf is not None:
                del gdf

    def clip_to_footprint(
        self,
        gdf,
        path,
        clip_to_footprint=False,
        dedup=None,
        fp_path="footprints",
        prop_duplicated="staging_duplicated",
    ):
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

        clip_to_footprint: boolean
            Flag to indicate if out of bounds geometries need to be clipped off.

        dedup: string
            The deduplication method for processing input geometries

        fp_path: str
            The path to the footprint file.

        prop_duplicated : str
            The name of the boolean property that indicates if a polygon
            was identified as a duplicate or not

        Returns
        -------
        The GeodataFrame with polygons that fall outside
        the footprint labeled as duplicates.

        """

        self.logger.info(f"clip_to_footprint is {clip_to_footprint}")
        # check if the config is set to label duplicates
        # if the config is set to do so, clip to footprint
        if clip_to_footprint == True and dedup is not None:
            self.logger.info(f" Starting clipping_to_footprint() for file {path}.")
            # pull in footprint as a gdf called fp
            fp = self.get_data(fp_path)
            self.logger.info(f" Checking CRSs of polygons and footprint.")

            data_crs = gdf.crs
            fp_crs = fp.crs

            if data_crs == fp_crs:
                self.logger.info(f" CRSs match. They are both {data_crs}.")
            else:
                self.logger.info(
                    f" CRSs do not match.\n Data's CRS is {data_crs}."
                    f" Footprint's CRS is {fp_crs}."
                )
                # transform the footprint to the CRS of the polygon data
                fp.to_crs(data_crs, inplace=True)
                # check again
                fp_crs_transformed = fp.crs
                if data_crs == fp_crs_transformed:
                    self.logger.info(
                        "Footprint CRS has been transformed to CRS of polygons."
                    )
                else:
                    self.logger.error(
                        "Failed to transform footprint CRS to CRS of polygons."
                    )
                    return

            # determine if polygons fall within or outside the footprint
            # first retrieve the name of the column we will use to label duplicates
            gdf_with_labels = clip_gdf(
                gdf=gdf.copy(),  # the gdf to clip
                boundary=fp.copy(),  # the footprint
                method="intersects",
                prop_duplicated=prop_duplicated,
            )
            # Clean up footprint GeoDataFrame to free memory
            del fp

            return gdf_with_labels
        else:
            self.logger.info(
                f" Either clip_to_footprint was set to False, or config"
                f" was not set to deduplicate at any step. Returning original GDF"
                f" without clipping to footprint."
            )
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
        self.logger.info(f"Reading vector file: {input_path}")
        try:
            gdf = gpd.read_file(input_path)
        except FileNotFoundError:
            gdf = None
            self.logger.warning(f"{input_path} not found. It will be skipped.")
            return None
        except Exception as e:
            self.logger.error(f"Error reading {input_path}: {e}")
            return None
        self.logger.info(f"Read in {input_path} in {(datetime.now() - start_time)}")

        # Check that none of the existing properties match the configured
        # properties that will be created. Finoa (used by geopandas) gives
        # error if properties are duplicated regardless of case.
        to_create = self.props.values()
        existing = gdf.columns.values
        duplicated = [
            b for b in to_create if b.lower() in (a.lower() for a in existing)
        ]

        if len(duplicated) > 0:
            error_msg = "The following properties already exist in the "
            error_msg += f"input vector file: {duplicated}"
            error_msg += "\nPlease remove them or change the configured "
            error_msg += "property names."
            self.logger.error(error_msg)
            raise ValueError(error_msg)

        return gdf

    def set_crs(self, gdf, input_crs=None):
        """
        Set the CRS of the GeoDataFrame to the input CRS, if there is one.
        Re-project to the CRS of the TMS.

        Parameters
        ----------
        gdf : GeoDataFrame
            The GeoDataFrame to set the CRS of the Tile Matrix Set

        input_crs : str
            If the input data is lacking CRS information, then the CRS of
            the input data. This will overwrite existing CRS data, if
            GeoPandas detects any. Input data will not be reprojected to
            this CRS.

        Returns
        -------
        GeoDataFrame
            The GeoDataFrame with the updated CRS

        """
        start_time = datetime.now()

        # assign the CRS of the TMS to `output_crs`
        output_crs = self.tiles.crs

        # log the CRS that already exists in the input data
        crs_in_input = gdf.crs

        if crs_in_input is not None:
            self.logger.info(
                f"CRS of input data is {crs_in_input}.\nIf input_crs is set in config,"
                f" setting CRS to that."
            )
        else:
            self.logger.info(
                f"No CRS set in input data. Setting to input_crs specified in config."
            )

        if input_crs:
            gdf.set_crs(input_crs, inplace=True, allow_override=True)
        # Re-project the geoms to the CRS of the TMS
        if output_crs:
            gdf.to_crs(output_crs, inplace=True)

        self.logger.info(
            f"Re-projected {len(gdf.index)} polygons in "
            f"{datetime.now() - start_time}"
        )
        return gdf

    def simplify_geoms(self, gdf, tolerance=None):
        """
        Simplify all geometries in a GeoDataFrame using the Douglas-Peucker
        algorithm

        Parameters
        ----------
        gdf : GeoDataFrame
            The GeoDataFrame to simplify

        tolerance: : float
            The tolerance to use when simplifying the input polygons.
            Defaults to 0.0001. Set to None to skip simplification.

        Returns
        -------
        GeoDataFrame
            The GeoDataFrame with the simplified geometries
        """
        start_time = datetime.now()
        if tolerance is not None:
            # Simplify in place and ensure old geometry data is released
            simplified = gdf["geometry"].simplify(tolerance)
            gdf["geometry"] = simplified
            del simplified
            self.logger.info(
                f"Simplified {len(gdf.index)} polygons in "
                f"{datetime.now() - start_time}"
            )
        return gdf

    def add_properties(self, gdf, path="", dir_input="input"):
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

        dir_input : str
            The directory to read input vector files from.

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
            warnings.simplefilter("ignore", UserWarning)
            centroids = gdf.centroid
            gdf[props["area"]] = gdf.area
        gdf[props["centroid_x"]] = centroids.x
        gdf[props["centroid_y"]] = centroids.y
        # Clean up centroids to free memory
        del centroids

        # Add source file path, excluding the input dir and any leading slashes
        gdf[props["filename"]] = path.removeprefix(dir_input).strip(os.sep)

        # Add identifier
        gdf[props["identifier"]] = [str(uuid.uuid4()) for _ in range(num_polygons)]

        # Find the tile the polygon falls within
        gdf = self.assign_tile_by_centroid(gdf)
        # Find the tiles the polygon intersects
        gdf = self.assign_tile(gdf)
        centroid_only = gdf[props["tile"]] == gdf[props["centroid_tile"]]
        gdf[props["centroid_within_tile"]] = centroid_only
        # Clean up intermediate boolean series
        del centroid_only

        self.logger.info(
            f"Added properties for {num_polygons} vectors in "
            f"{datetime.now() - start_time} from file {path}"
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
        gdf[self.props["centroid_tile"]] = self.grid.tiles_from_xy(
            gdf[self.props["centroid_x"]], gdf[self.props["centroid_y"]]
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
        return self.grid.sjoin(gdf, how="left", predicate="intersects", as_tile=True)

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
            tms_id=self.tiles.tms_id, z=self.z_level, bounds=gdf.total_bounds
        )
        grid.TILE_NAME = self.props["tile"]
        return grid

    def save_tiles(
        self,
        gdf=None,
        dedup=None,
        deduplicate_at="staging",
        prop_duplicated="staging_duplicated",
    ):
        """
        Given a processed GeoDataFrame, save vectors into tiled vector
        files.

        Parameters
        ----------
        gdf : GeoDataFrame
            The GeoDataFrame to save

        dedup: string
            The deduplication method for processing input geometries

        deduplicate_at : str
            Step at which de-duplication should occur

        prop_duplicated : str
            The name of the boolean property that indicates if a polygon
            was identified as a duplicate or not
        """

        if gdf is None:
            # TODO give warning
            return None

        if dedup is not None:
            if dedup == "neighbors":
                dedup_method = deduplicate_neighbors
            if dedup == "footprints":
                dedup_method = deduplicate_by_footprint
        else:
            dedup_method = None

        for tile, data in gdf.groupby(self.props["tile"]):
            # Create a copy to avoid modifying the grouped data
            data = data.copy()

            # Get the tile path
            tile_path = self.tiles.path_from_tile(tile, base_dir="staged")
            self.tiles.create_dirs(tile_path)

            # Lock the tile so that multiple processes don't try to write to
            lock = self.__lock_file(tile_path)

            # Track the start time, the tile, and the number of vectors
            start_time = datetime.now()
            self.logger.info(f"Saving {len(data.index)} vectors to tile {tile_path}")

            # Tile must be a string for saving as attribute
            data[self.props["tile"]] = data[self.props["tile"]].astype("str")
            tile_strings = data[self.props["centroid_tile"]].astype("str")
            data[self.props["centroid_tile"]] = tile_strings

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
                    self.logger.info(
                        f"Tile exists and dedup is set to occur at some step,"
                        f" so executing `combine_and_deduplicate()`"
                    )
                    data = self.combine_and_deduplicate(data, tile_path)

                    mode = "w"
                    # Overwrite existing file.
                    self.save_new_tile(
                        data=data,
                        tile_path=tile_path,
                        mode=mode,
                        start_time=start_time,
                        lock=lock,
                    )
                else:
                    # If the file exists and config is not set to deduplicate
                    # at any step, simply append the polygons to the existing file.
                    # no polygons will be labeled as duplicates or not.
                    # If deduplicating by footprint:
                    # neither file has been clipped to footprint
                    self.logger.info(
                        f"Tile exists but dedup is not set to occur, so"
                        f" appending polygons."
                    )

                    # Append to existing tile
                    mode = "a"
                    self.save_new_tile(
                        data=data,
                        tile_path=tile_path,
                        mode=mode,
                        start_time=start_time,
                        lock=lock,
                    )
            else:
                if dedup_method is not None:
                    if deduplicate_at == "staging":
                        # If the file does not yet exist and the config is set to deduplicate
                        # at staging:
                        # If deduplicatinf by footprint:
                        # remove the polygons that were already labeled as
                        # duplicates because they fell outside the footprint,
                        # and save the tile as a new tile.
                        # If deduplicating by neighbor:
                        # the prop_duplicated col has not been created yet,
                        # so create it and set all values to False
                        self.logger.info(
                            f"Tile does not yet exist and config is set to deduplicate "
                            f"at staging, so removing polygons that fell outside the footprint "
                            f"if deduplicating by footpint, and removing overlapping polygons\n"
                            f"if deduplicating by neighbor if column already existed.\n "
                            f"Creating column with all false values it it did not exist."
                        )
                        self.logger.info(
                            f"Checking for presence of {prop_duplicated} column."
                        )
                        if prop_duplicated in data.columns:
                            data = data[~data[prop_duplicated]]
                        else:
                            self.logger.info(
                                f"Adding {prop_duplicated} column because property did not "
                                f"already exist."
                            )
                            data[prop_duplicated] = False

                        mode = "w"
                        self.save_new_tile(
                            data=data,
                            tile_path=tile_path,
                            mode=mode,
                            start_time=start_time,
                            lock=lock,
                        )
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
                        self.logger.info(
                            f"Tile does not yet exist and config is set to deduplicate at a step "
                            f"after staging, so just saving the new tile."
                            f"\nIf deduplicating by footprint: "
                            f"Duplicates from `clip_gdf` were identified."
                        )

                        self.logger.info(
                            f"Checking for presence of {prop_duplicated} column."
                        )
                        if prop_duplicated not in data.columns:
                            self.logger.info(
                                f"Adding {prop_duplicated} column because property did not "
                                f"already exist."
                            )
                            data[prop_duplicated] = False
                        else:
                            self.logger.info(
                                f"Tile that did not yet exist did have "
                                f"{prop_duplicated} column before saving."
                            )

                        mode = "w"
                        self.save_new_tile(
                            data=data,
                            tile_path=tile_path,
                            mode=mode,
                            start_time=start_time,
                            lock=lock,
                        )
                else:
                    # If the tile does not yet exist and the config is not set to deduplicate
                    # at any step, just save as a new tile.
                    # If deduplicating by footprint:
                    # No duplicates were labeled earlier either,
                    # because the workflow did not check if any fell outside the footprint.
                    self.logger.info(
                        f"Tile does not yet exist and config is not set to deduplicate, so just "
                        f"saving the new tile.\nIf deduplicating by footprint:\n"
                        f"Duplicates from `clip_gdf` were not identified."
                    )

                    mode = "w"
                    self.save_new_tile(
                        data=data,
                        tile_path=tile_path,
                        mode=mode,
                        start_time=start_time,
                        lock=lock,
                    )

    def save_new_tile(self, data, tile_path, mode, start_time, lock):
        """
        Save data as a new tile, if the tile does not already exist.
        Adds a row to staging_summary.csv

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
                warnings.simplefilter("ignore", FutureWarning)
                data.to_file(tile_path, mode=mode)

            # convert each tile from string format to morecantile format
            # so it can be added to summary
            # first create series of tiles in str format
            tiles_str = data[self.props["tile"]].copy()

            tiles_morecantile = [self.tiles.tile_from_str(tile) for tile in tiles_str]

            # Clean up intermediate tiles_str to free memory
            del tiles_str

            # Record what was saved
            data[self.props["tile"]] = tiles_morecantile
            summary_csv_path = self.summary_path
            self.summarize(data, summary_csv_path)
        except Exception as e:
            self.logger.error(f"Error saving tile {tile_path}: {e}")
            raise
        finally:
            # Track the end time, the total time, and the number of vectors
            self.logger.info(f"Saved {tile_path} in {datetime.now() - start_time}")
            self.__release_file(lock)

    def combine_and_deduplicate(
        self,
        gdf,
        tile_path,
        dedup=None,
        deduplicate_at="staging",
        prop_duplicated="staging_duplicated",
        dedup_config={},
    ):
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

        dedup: str
            The deduplication method for processing input geometries

        deduplicate_at : str
            Step at which de-duplication should occur

        prop_duplicated : str
            The name of the boolean property that indicates if a polygon
            was identified as a duplicate or not

        dedup_config : str
            Config for constructing the deduplication object

        Returns
        -------
        gdf : GeoDataFrame
            The combined GeoDataFrame without duplicates labeled or with
            deduplicates labeled, depending on how the config is set.
        """

        dedup_start_time = datetime.now()

        if dedup is not None:
            if dedup == "neighbors":
                dedup_method = deduplicate_neighbors
            if dedup == "footprints":
                dedup_method = deduplicate_by_footprint
        else:
            dedup_method = None

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
        # Clean up intermediate objects to free memory
        del existing_gdf
        del to_concat
        gdf.reset_index(drop=True, inplace=True)

        self.logger.info(
            f"Starting deduplication in tile {tile_path} with {len(gdf)} " "polygons."
        )

        # label duplicates depending on which deduplication type is set in config
        gdf = dedup_method(gdf, **dedup_config)

        # drop duplicated polygons, if config is set to deduplicate here
        if deduplicate_at == "staging":
            if prop_duplicated in gdf.columns:
                before_dedup = len(gdf)
                gdf = gdf[~gdf[prop_duplicated]]
                # Force garbage collection after dropping rows
                if before_dedup - len(gdf) > 1000:
                    gc.collect()

        self.logger.info(
            f"Finished deduplication in {datetime.now() - dedup_start_time}"
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
            "tile_x": tile.x,
            "tile_y": tile.y,
            "tile_z": tile.z,
            "tile_top": bounds["top"],
            "tile_right": bounds["right"],
            "tile_bottom": bounds["bottom"],
            "tile_left": bounds["left"],
        }

    def summarize(self, gdf=None, summary_path="staging_summary.csv"):
        """
        For a given file, count how many vectors there are per tile, the
        area of vectors per tile, and get information about the tile itself
        (bounds, x, y, z values). Append these details to the summary
        dataframe.

        Parameters
        ----------
        gdf : GeoDataFrame
            The GeoDataFrame to summarize

        summary_path : str
            The path and filename to save a CSV file that summarizes the
            tiled files that were created during the staging process.
        """

        if gdf is None:
            return None

        # Log the summary event, including how long it takes
        start_time = datetime.now()

        prop_file = self.props["filename"]
        prop_tile = self.props["tile"]
        prop_area = self.props["area"]

        # Since summarized is called for each tile that is created, grouping
        # should create 1 row of data to add to the summary dataframe
        gdf_grouped = gdf.groupby([prop_file, prop_tile], as_index=False)
        gdf_summary = gdf_grouped.agg(
            num_polygons=("geometry", "count"), area_polygons=(prop_area, "sum")
        )
        gdf_summary = gdf_summary.rename(
            columns={prop_file: "filename", prop_tile: "tile"}
        )
        # Add the date time that the tile was saved
        gdf_summary["datetime"] = datetime.now()

        tiles = gdf_summary["tile"]
        tile_props = self.get_all_tile_properties(tiles)

        gdf_summary = pd.concat([gdf_summary, pd.DataFrame(list(tile_props))], axis=1)
        gdf_summary["tms_identifier"] = self.tiles.tms_id
        gdf_summary = gdf_summary.drop(columns=["tile"])

        # Save the summary to a file
        header = False
        mode = "a"
        if not os.path.isfile(summary_path):
            header = True
            mode = "w"
        gdf_summary.to_csv(summary_path, mode=mode, index=False, header=header)
        # Clean up DataFrame after writing
        del gdf_summary

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
        lock_path = path + ".lock"
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

    def get_default_base_dir(self):
        """
        Returns the default base dir for the TPM class

        Parameters
        ----------
        None

        Returns
        -------
        base_dir : dictionary
            Default base directory representing the TPM structure for the WGS1984Quad
        """
        tmp_base_dirs = {}

        # Set up input base dirs and output for the raster base dir
        tmp_base_dirs["input"] = {}
        tmp_base_dirs["input"]["path"] = "input"
        tmp_base_dirs["input"]["ext"] = ".shp"

        tmp_base_dirs["staged"] = {}
        tmp_base_dirs["staged"]["path"] = "staged"
        tmp_base_dirs["staged"]["ext"] = ".gpkg"

        return tmp_base_dirs

    def set_tiling_config(
        self,
        tms_id="WGS1984Quad",
        path_structure=["style", "tms", "z", "x", "y"],
        base_dirs={},
    ):
        """
        Updates the tiling config for the TilePathmanager class

        Parameters
        ----------
        tms_id : str
            The tile matrix set id for constructing the tile path manager object

        path_structure : list
            The path structure used for generating the tiling output

        base_dirs : dict
            Directory structure and extenstions for storing the output

        Returns
        -------
        None
        """
        self.tms_id = tms_id
        self.path_structure = path_structure

        if not base_dirs:
            self.base_dirs = self.get_default_base_dir()

        for key, value in base_dirs.items():
            if key in self.base_dirs:
                self.base_dirs[key].update(value)
            else:
                self.base_dirs[key] = value

        self.tiles = TilePathManager(
            tms_id=self.tms_id,
            path_structure=self.path_structure,
            base_dirs=self.base_dirs,
        )

        return None

    def csv_to_parquet(self):
        """
        This method is to converst the CSV files containing rasterization events and raster
        summaries to Parquet format, keeping the original CSVs.
        """
        csv_path = self.summary_path

        if not os.path.isfile(csv_path):
            self.logger.warning(f"CSV not found → {csv_path}")

        root, _ = os.path.splitext(csv_path)
        parquet_path = f"{root}.parquet"

        self.logger.warning(f"Converting {csv_path} to {parquet_path}.")

        df = pd.read_csv(csv_path)
        pq.write_table(
            pa.Table.from_pandas(df, preserve_index=False),
            parquet_path,
            compression="snappy",
        )
