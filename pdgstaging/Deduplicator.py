import operator
import geopandas as gpd
import pandas as pd
import uuid
import itertools
import warnings

import os
from datetime import datetime
import numpy as np
from filelock import FileLock
import logging
from . import logging_config
# NOTE: DO NOT IMPORT ConfigManager, TilePathManager, Grid
# because causes config import error for rasterization step 

logger = logging_config.logger


def keep_rules_to_sort_order(keep_rules):
    """
    Convert a list of keep rules to the format required for the pandas
    sort_values function.

    Parameters
    ----------
    keep_rules : list
        A list of tuples in the form (property, operator), as required for the
        deduplication methods.

    Returns
    -------
    tuple
        The first element in the tuple is the list of properties to sort by,
        and the second element is the list of booleans indicating whether the
        corresponding property should be sorted in ascending or descending
        order. To be used for the `by` and `ascending` arguments of the pandas
        sort_values function.
    """
    sort_props = [x[0] for x in keep_rules]
    sort_order = [x[1] == 'smaller' for x in keep_rules]
    return sort_props, sort_order

def clip_gdf(gdf = None, boundary = None, method = 'intersects', prop_duplicated = None):
    """
        Remove polygons from a GeoDataFrame that fall outside of some boundary.
        Determine if the polygons in the GeoDataFrame are within the boundary
        by using an sjoin operation, with a specified method.

        Parameters
        ----------
        gdf : GeoDataFrame
            The GeoDataFrame to clip.

        boundary : GeoDataFrame
            A GeoDataFrame that contains polygons representing the boundaries.

        method : str
            The predicate to use for the sjoin operation. Can be one of:
            'contains', 'contains_properly', 'covers', 'crosses', 'intersects',
            'overlaps', 'touches', 'within' (any option listed by
            geopandas.GeoDataFrame.sindex.valid_query_predicates)

        Returns
        -------
            A dictionary containing the clipped GeoDataFrame (the value of the
            'to_keep' key), and the GeoDataFrame of polygons that were removed
            (the value of the 'to_remove' key).

    """

    # Temporary property to use during the sjoin
    # assigns an identifier to 
    prop_in_fp_temp = 'WITHIN_BOUNDARY_' + uuid.uuid4().hex

    # define the boundary as the footprint gdf's geometry column
    boundary = boundary.copy().filter(['geometry'])

    boundary[prop_in_fp_temp] = True
    # use sjoin to determine if the polygons are within the footprint
    # then drop the column called `index_right`
    gdf = gdf.sjoin(boundary, how='left', predicate=method) \
             .drop(['index_right'], axis=1)

    # create a gdf called `within` that contains only the gdf rows where the
    # values in prop_in_fp_temp are not null (are within the footprint)
    within = gdf[gdf[prop_in_fp_temp].notnull()]
    # drop the column from `within` object that's used 
    # to identify polys that fall outside the footprint
    within = within.drop([prop_in_fp_temp], axis=1)

    # create a gdf called `outside` that contains only the gdf rows where the
    # values in prop_in_fp_temp are null (are not within the footprint)
    outside = gdf[~gdf[prop_in_fp_temp].notnull()]
    outside = outside.drop([prop_in_fp_temp], axis=1)

    outside[prop_duplicated] = True
    within[prop_duplicated] = False

    # stack the gdf's
    gdf_with_labels = pd.concat([outside, within], ignore_index = True)
    gdf_with_labels.reset_index(drop = True, inplace = True)

    return gdf_with_labels


def deduplicate_neighbors(
    gdf,
    split_by=None,
    prop_area=None,
    prop_centroid_x=None,
    prop_centroid_y=None,
    keep_rules=[],
    overlap_tolerance=0.5,
    overlap_both=True,
    centroid_tolerance=None,
    distance_crs='EPSG:3857',
    return_intersections=False,
    prop_duplicated='staging_duplicated'
     # defaults to these options only if aren't already specified in config
):
    """

    Identify duplicated polygons from two sources in a single GeoDataFrame, even
    when geometries are not identical.

    Two polygons are considered duplicates of one another when they:
        1. have a different value for a specific property (e.g. 'source_file'),
           AND
        2. are overlapping, AND
        3. have centroids within a given distance, OR
        4. at least a certain proportion of one or both polygons is overlapped
           by the other

   The deduplication can be set to ignore centroid distance or ignore the
   overlap proportion.

    When two polygons are identified as duplicates, only one is kept in the returned 'to_keep'
    dataframe. Which one is kept is determined by the 'keep_rules' parameter,
    which indicates the property to compare, and the rule to apply ('larger' or
    'smaller'). For example, to keep the newest polygon, set a keep_rule for a
    more recent date property, e.g. ('Date', 'larger').

    A dictionary is returned with the deduplicated GeoDataFrame ('to_keep'), the
    GeoDataFrame with the dropped duplicates ('to_remove'), and optionally, the
    intersection geometries of duplicated polygons ('intersections').

    Parameters
    ----------
    gdf : GeoDataFrame
        The GeoDataFrame to deduplicate.
    split_by : str
        The name of the property in the GeoDataFrame to use to identify
        duplicated polygons. Two polygons must have a different value for this
        property to be considered duplicated. e.g. 'source_file.' Every
        pairwise combination of unique values for this property will be
        compared, so ideally, this property should not have too many unique
        values, or the deduplication will be very slow.
    prop_area : str, optional
        If the area has already been calculated for the GeoDataFrame, give the
        name of the property containing the area. Area must be in the same unit
        as the CRS of the GeoDataFrame. If set to None, the area will be
        calculated.
    prop_centroid_x : str, optional
        If the centroid has already been calculated for the GeoDataFrame, give
        the name of the property containing the centroid x-coordinate. Centroid
        coordinates must be in the same CRS of the GeoDataFrame. If set to
        None, the centroid will be calculated.
    prop_centroid_y : str, optional
        If the centroid has already been calculated for the GeoDataFrame, give
        the name of the property containing the centroid y-coordinate. Centroid
        coordinates must be in the same CRS of the GeoDataFrame. If set to
        None, the centroid will be calculated.
    keep_rules : list, optional
        Rules that define which of the polygons to keep when two or more are
        duplicates. A list of tuples of the form (property, operator), where
        property is the name of a property in the GDF to use for the
        comparison. The operator is either 'larger' or 'smaller'. If the rule
        is 'larger', the polygons with the largest value for the property will
        be kept. If the rule is 'smaller', the polygon with the smallest value
        for the property will be kept. When two properties are equal, then the
        next property in the list will be checked.
    overlap_tolerance : float, optional
        The minimum proportion of a polygon's area that must be overlapped by
        another polygon to be considered a duplicate. Default is 0.5. Set to
        None to ignore overlap proportions when comparing polygons, and set a
        centroid threshold instead. Note that both an overlap_tolerance AND a
        centroid_tolerance can be used.
    overlap_both : bool, optional
        If True, then the overlap_tolerance proportion must be True for both of
        the intersecting polygons to be considered a duplicate. If False, then
        the overlap_tolerance proportion must be True for only one of the
        intersecting polygons to be considered a duplicate. Default is True.
    centroid_tolerance : float, optional
        The maximum distance between the centroids of two polygons to be
        considered a duplicate. Default is None. Set to None to ignore centroid
        distances when comparing polygons, and set an overlap threshold
        instead. Note that both an overlap_tolerance AND a centroid_tolerance
        can be used. The unit of the distance is the unit of the distance_crs
        property (e.g. meters for EPSG:3857), or the unit of the GeoDataFrame
        if distance_crs is None.
    distance_crs : str, optional
        The CRS to use for the centroid distance calculation. Default is
        EPSG:3857. Centroid points will be re-projected to this CRS before
        calculating the distance between them. centroid_tolerance will use the
        units of this CRS. Set to None to skip the re-projection and use the
        CRS of the GeoDataFrame.
    return_intersections : bool, optional
        If True, the GeoDataFrame with the intersection geometries is returned
        in addition to the GeoDataFrames with the deduplicated and removed
        geometries. Default is False.
    prop_duplicated : str, optional
        Defaults to "staging_duplicated". The column name / property to use to flag
        duplicates when label is True.

    Returns
    -------
    GeoDataFrame
        Returns the input GDF with the polygons flagged as duplicates. 
        Duplicates are marked as True in the prop_duplicated column.

        'intersections' represents the GeoDataFrame with the intersections. 
        It has not been integrated into the function again since the deduplication 
        approach changed from returning a dictionary to returning a labeled GDF.
    """

    if split_by is None:
        raise ValueError('An overlap property must be specified for '
                         'deduplication.')

    if (keep_rules is None) or (len(keep_rules) == 0):
        raise ValueError('A list of keep rules must be specified for '
                         'deduplication.')

    # Use uuid for the properties that will be temporarily created so they do
    # not conflict with existing properties
    rand_id = uuid.uuid4().hex
    prop_int_area = 'intersect_area_' + rand_id
    prop_int_id = 'intersect_id_' + rand_id
    prop_id = 'temp_id_' + rand_id

    gdf = gdf.copy()

    # Add an ID that we can refer to later, and make a copy of the original
    gdf[prop_id] = list(range(len(gdf)))

    gdf_original = gdf.copy()

    # Save the CRS
    crs = gdf.crs

    if not crs:
        # TODO: raise error?
        return

    # Calculate centroids if needed
    if (prop_centroid_x is None) or (prop_centroid_y is None):
        if centroid_tolerance is not None:
            prop_centroid_x = 'cent_x_' + rand_id
            prop_centroid_y = 'cent_y_' + rand_id
            centroids = gdf.centroid
            gdf[prop_centroid_x] = centroids.x
            gdf[prop_centroid_y] = centroids.y

    # Calculate area if needed
    if (prop_area is None):
        if overlap_tolerance is not None:
            prop_area = 'area_' + rand_id
            gdf[prop_area] = gdf.area

    # Organize keep rules for pandas sort_values
    sort_props, sort_order = keep_rules_to_sort_order(keep_rules)

    # Remove all unnecessary columns
    working_cols = [
        'geometry',
        prop_id,
        split_by,
        prop_area,
        prop_centroid_x,
        prop_centroid_y] + sort_props
    working_cols = [i for i in working_cols if i is not None]
    gdf.drop(gdf.columns.difference(working_cols), axis=1, inplace=True)

    # Split the gdf by the overlap property (e.g. source file)
    gdf_grouped = gdf.groupby(split_by)
    gdfs = [gdf_grouped.get_group(x) for x in gdf_grouped.groups]

    to_return = {
        'to_keep': None,
        'intersections': None,
        'to_remove': None
    }

    # If there is only one group, then there is nothing to deduplicate
    if len(gdfs) < 2:
        to_return['to_keep'] = gdf_original.drop(columns=[prop_id])

        to_return = label_duplicates(to_return, prop_duplicated)
        # end the deduplicate_neighbors() here, no need to continue
        return to_return

    # Set will hold all of the polygons IDs that we will remove
    ids_to_remove = set()

    # For every pairwise combination of two (subsetted) GDFs
    for gdf_pair in itertools.combinations(gdfs, 2):

        g1 = gdf_pair[0]
        g2 = gdf_pair[1]

        # Remove any of the ids_to_remove (so we don't process twice)
        g1 = g1[~g1[prop_id].isin(ids_to_remove)]
        g2 = g2[~g2[prop_id].isin(ids_to_remove)]

        if (len(g1) == 0) or (len(g2) == 0):
            continue

        # Polygons must intersect to be considered duplicated
        duplicates = g1.overlay(g2, how='intersection')
        # Add a unique ID for each new intersection geometry
        duplicates[prop_int_id] = list(range(0, len(duplicates)))

        # Identify a polygon as duplicated if the area of overlap is >
        # threshold proportion
        if overlap_tolerance is not None:
            # Filter out any intersections that do not meet the threshold %
            # overlap Get the area of the intersection
            # Catch the UserWarning that is raised if area is calculated
            # using the a non-projected coordinate system.
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', UserWarning)
                duplicates[prop_int_area] = duplicates.area
            duplicates[prop_int_area + '_1'] = duplicates[prop_int_area] / \
                duplicates[prop_area + '_1']
            duplicates[prop_int_area + '_2'] = duplicates[prop_int_area] / \
                duplicates[prop_area + '_2']
            # Must the overlap propotion be > threshold for both or just one
            # polygon
            op = operator.and_ if overlap_both else operator.or_
            # Remove polygons from list of duplicates that intersect but have
            # less than the threshold overlap proportion.
            duplicates = duplicates[op(
                (duplicates[prop_int_area + '_1'] > overlap_tolerance),
                (duplicates[prop_int_area + '_2'] > overlap_tolerance)
            )]
            # Remove columns no longer needed
            duplicates = duplicates.drop(
                columns=[
                    prop_int_area + '_1',
                    prop_int_area + '_2',
                    prop_int_area
                ])

        # Identify a polygon as duplicated when the distance between two
        # centroids is <= the threshold distance
        if centroid_tolerance is not None:
            # Calculate distance between centroids of each polygon, in unit of
            # measurement CRS
            points = []
            for i in [1, 2]:
                i = str(i)
                x = duplicates[prop_centroid_x + '_' + i]
                y = duplicates[prop_centroid_y + '_' + i]
                p = gpd.GeoSeries(gpd.points_from_xy(x, y, crs=crs))
                if(distance_crs):
                    p = p.to_crs(distance_crs)
                points.append(p)

            centroid_distances = points[0].distance(points[1], align=False)
            duplicated_by_centroid = centroid_distances <= centroid_tolerance
            duplicated_by_centroid = duplicated_by_centroid.tolist()
            # Remove polygons from list of duplicates that intersect but have
            # centroids further apart than the threshold distance
            duplicates = duplicates[duplicated_by_centroid]

        # If we want to return the intersections
        if return_intersections:
            to_return['intersections'] = duplicates.copy()
            to_return['intersections'].drop(
                [prop_int_id, prop_id + '_1', prop_id + '_2'],
                axis=1, inplace=True
            )

        # Can treat as a regular DF now
        duplicates = duplicates.drop(['geometry'], axis=1)

        # Reshape from wide to long
        reshape_dict = {}
        for c in working_cols:
            cols = [col for col in duplicates.columns if col.startswith(c)]
            if len(cols) > 1:
                reshape_dict[c] = cols

        duplicates = pd.lreshape(duplicates, reshape_dict)
        # Sort to determine which of the duplicates to keep, which should be
        # removed
        duplicates.sort_values(
            by=sort_props,
            ascending=sort_order,
            inplace=True)
        # 'remove' here is correct, should not be 'to_remove'
        duplicates['remove'] = duplicates.duplicated(prop_int_id, keep='first')

        # Get the IDs of only the polygons to remove, and add them to the list
        # 'remove' here is correct, should not be 'to_remove'
        ids_to_remove_pair = duplicates[duplicates['remove']][prop_id]
        ids_to_remove.update(ids_to_remove_pair)

    # Remove polygons from the original gdf when they have an ID that is in the
    # remove list
    remove = gdf_original[prop_id].isin(ids_to_remove)

    # Remove the prop_id column, no longer needed
    gdf_original.drop([prop_id], axis=1, inplace=True)

    to_return['to_remove'] = gdf_original[remove]
    to_return['to_keep'] = gdf_original[~remove]

    to_return = label_duplicates(to_return, prop_duplicated)

    if prop_duplicated in to_return.columns:
        if True in to_return[prop_duplicated].values:
            sum_true = (to_return[prop_duplicated] == True).value_counts()[True]
            logger.info(f"Sum of True values in the {prop_duplicated} col is: {sum_true}")
        else:
            sum_true = 0
            logger.info(f"Sum of True values in the {prop_duplicated} col is: {sum_true}")
    else:
        logger.info(f"{prop_duplicated} is not a column present after labeling.")

    return to_return


def deduplicate_by_footprint(
    gdf,
    split_by,
    footprints,
    keep_rules=[],
    return_intersections=False,
    prop_duplicated='staging_duplicated' 
    # defaults to these options only if aren't already specified in config
):
    """

    Label polygons to remove if they occur within an area where associated footprints
    overlap. This method can be used to remove polygons that originate from
    different overlapping files.

    Parameters
    ----------
    gdf : GeoDataFrame
        The GeoDataFrame to deduplicate
    split_by : str
        The name of the property in the GeoDataFrame that links the polygon to
        the footprint, e.g. 'source_file'. This property will be used to divide
        the GeoDataFrame into groups. The groups with the preferred polygons
        will be identified using the `rank` parameter. Polygons from less
        preferred groups will be removed when they occur within an area where
        footprints overlap. Every pairwise combination of unique values for
        this property will be compared, so ideally, this property should not
        have too many unique values, or the deduplication will be very slow.
    footprints : dict
        A dictionary that maps the unique values of `split_by` to a
        GeoDataFrame of or path (string) to the footprint of the associated
        file
    keep_rules : str
        One or more rules that define which footprint is preferred. Polygons
        that occur in an area where footprints overlap will be removed if they
        are not associated with the preferred footprint. Keep_rules are given
        as a list of tuples in the form (property, operator). The property is
        the name of a property in the footprint file. The operator is either
        'larger' or 'smaller'. If the rule is 'larger', then polygons from the
        footprint with the largest value for the property will be kept, and
        vice versa for smaller. When two properties are equal, then the next
        property in the list will be checked.
    return_intersections : bool, optional
        If true, the polygons that represent the intersections between
        footprints will be returned. Default is False. Not currently available
        in this release 0.9.0. return_intersections is to be integrated again 
        in future releases.
    prop_duplicated : str, optional
        Defaults to "staging_duplicated". The column name / property to use to 
        flag duplicates.

    Returns
    -------
    GeoDataFrame
        Returns input GDF with the polygons flagged as duplicates if they meet
        the criteria outlined in the docs. Duplicates are marked as True in the
        prop_duplicated column.

        `intersections` represents the polygon area where the footprints overlap. 
        It has not been integrated into the function again since the deduplication 
        approach changed from returning a dictionary to returning a labeled GDF.
        This will be integrated again in releases after 0.9.0.
    """

    logger.info(f"Executing deduplicate_by_footprint() for {gdf}")

    gdf = gdf.copy()

    # `to_remove` list will hold the polygons that fit either of the criteria:
    # 1. were previously labeled True for `duplicated` col because poly
    # fell outside footprint
    # 2. are labeled as True for `duplicated` col within this function
    # based on overlap of footprints
    to_remove = []

    # subset the gdf to just polys that were identified as dups because
    # they fell outside the footprint, earlier in the staging step
    known_dups = gdf[gdf[prop_duplicated] == True]
    logger.info(f"Length of known_dups is {len(known_dups)}.")

    # Will hold the polygons that defined the footprint intersections
    intersections = []

    # To make sure footprints and overlap GeoDataFrame are in the correct CRS
    crs = gdf.crs

    # Get the unique values of the split_by property. Divide the GeoDataFrame.
    # split_by is the filename, so data is grouped by input file 
    gdf_grouped = gdf.groupby(split_by)
    # create a dict with each key being an input file,
    # and the value of each key is the polygons that came from that file
    gdf_dict = {}
    for g in gdf_grouped.groups:
        gdf_dict[g] = gdf_grouped.get_group(g)
    names = list(gdf_grouped.groups)

    # If the values of footprints dict are strings, then assume the strings are
    # paths and load the footprints as individual GeoDataFrames
    if all([isinstance(v, str) for v in footprints.values()]):
        for name, path in footprints.items():
            try:
                footprints[name] = gpd.read_file(path)
            except Exception:
                footprints[name] = None
                warnings.warn(f'Footprint missing for {name}') 

    # Add a column to the footprint GeoDataFrame that contains the filename of the footprint
    prop_filename_temp = 'filename_' + uuid.uuid4().hex
    for name, fp_gdf in footprints.items():
        fp_gdf[prop_filename_temp] = name
        # Make sure the footprints are in the same CRS as the GeoDataFrame
        fp_gdf.to_crs(crs, inplace=True)
        footprints[name] = fp_gdf

    # Rank the footprints according to the keep_rules
    # to determine which file's polygons to keep where two
    # footprints overlap
    footprints_concat = gpd.GeoDataFrame(pd.concat(
        footprints.values(), ignore_index=True))

    sort_props, sort_order = keep_rules_to_sort_order(keep_rules)
    footprints_concat.sort_values(
        by=sort_props,
        ascending=sort_order,
        inplace=True)
    rank = footprints_concat[prop_filename_temp].tolist()

    # Find overlapping section of footprints
    # For every pairwise combination of two (subsetted) GDFs
    for pair in itertools.combinations(names, 2):
        name1 = pair[0]
        name2 = pair[1]
        footprint1 = footprints.get(name1)
        footprint2 = footprints.get(name2)
        if(footprint1 is None or footprint2 is None):
            continue

        # Get overlap between two footprints
        overlap = gpd.GeoDataFrame(
            geometry=footprint1.intersection(footprint2))
        overlap['overlap'] = True
        overlap.to_crs(crs, inplace=True)
        intersections.append(overlap)

        # Identify which GDF is not preferred in the pair. This is the GDF that
        # we will remove polygons from.
        least_preferred = name1
        if rank.index(name2) > rank.index(name1):
            least_preferred = name2
        to_reduce = gdf_dict.get(least_preferred)

        # Find polygons in the non-preferred GDF that intersect with overlap
        # and remove them.
        overlap_boolean = to_reduce.sjoin(overlap, how='left')['overlap']
        overlap_boolean = overlap_boolean.fillna(False)
        # subset the least_preferred GDF for only polygons that do not fall within
        # the region where the footprints overlap
        reduced = to_reduce[~overlap_boolean]

        # Update the dictionary with the reduced GDF
        gdf_dict[least_preferred] = reduced

        # Save the removed polygons to the list defined at start of function
        # that holds all dataframes with polys labeled as duplicates
        to_remove.append(to_reduce[overlap_boolean])

    # Recombine the GDFs from the dictionary:
    keep = pd.concat(gdf_dict.values(), ignore_index=True)
    keep.reset_index(drop = True, inplace = True)

    # combine the list of dataframes (called `to_remove`) into one dataframe
    # along rows (stack the dataframes)
    to_remove = pd.concat(to_remove, ignore_index = True)
    to_remove.reset_index(drop = True, inplace = True) 

    # the two duplicate dataframes have the same columns, so concatenate them (stack)
    to_remove_all = pd.concat([to_remove, known_dups], ignore_index = True)
    to_remove_all.reset_index(drop = True, inplace = True)

    # ensure that the `keep` gdf has the `prop_duplicated` col too,
    # and values in the `prop_duplicated` column = False
    if keep is not None:
        # if the column `prop_duplicated` is already present,
        # set values to True if already set to True
        # if the value is not already True, set it to False
        if prop_duplicated in keep.columns:
            keep[prop_duplicated] = \
                keep[prop_duplicated].apply(
                    lambda x: True if x is True else False)
        # if the column `prop_duplicated` is not already present, 
        # create it and set all values to False
        else: 
            keep[prop_duplicated] = False

    # ensure that the labels for all of `prop_duplicated` are True for to_remove_all
    if to_remove_all is not None:
        to_remove_all[prop_duplicated] = True

    # combine the keep and to_remove_all into 1 gdf, 
    # stack them because already have the same columns
    to_return = pd.concat([keep, to_remove_all], ignore_index = True)
    to_return.reset_index(drop = True, inplace = True)

    #if return_intersections:
    #    to_return['intersections'] = pd.concat(intersections, ignore_index=True)
    # need to re-integrate this optional step later, 
    # not currently available since removing dict step during cliping to FP

    if prop_duplicated in to_return.columns:
        if True in to_return[prop_duplicated].values:
            sum_true = (to_return[prop_duplicated] == True).value_counts()[True]
            logger.info(f"Sum of True values in the {prop_duplicated} col is: {sum_true}")
        else:
            sum_true = 0
            logger.info(f"Sum of True values in the {prop_duplicated} col is: {sum_true}")
    else:
        logger.info(f"{prop_duplicated} is not a column present after labeling.")

    return to_return


def label_duplicates(deduplicate_output, prop_duplicated):
    """Recombine the keep & removed GDFs and mark the removed as duplicates

    Parameters
    ----------
    deduplicate_output : dict
        The output of one of the deduplication methods

    prop_duplicated : str
        The column name to use to flag duplicate polygons

    Returns
    -------

    GeoDataFrame
    The combined GDF with duplicated polygons flagged.

    """

    not_duplicates = deduplicate_output['to_keep']
    duplicates = deduplicate_output['to_remove']

    if duplicates is not None:
        duplicates[prop_duplicated] = True

    if not_duplicates is not None:
        if prop_duplicated in not_duplicates.columns:
            # if the gdf of polygons to keep exists and
            # contains the column `prop_duplicated`
            # set values to True if already set to True
            # if the value is not already True, set it to False
            not_duplicates[prop_duplicated] = \
                not_duplicates[prop_duplicated].apply(
                    lambda x: True if x is True else False)
        else:
            # if the column `prop_duplicated` is not already present, 
            # create it and set all values to False
            not_duplicates[prop_duplicated] = False 

    combined = pd.concat([not_duplicates, duplicates], ignore_index=True)
    combined.reset_index(drop=True, inplace=True)

    if prop_duplicated not in combined.columns:
        logger.warning(f"Config set to deduplicate but {prop_duplicated} col not in gdf.")
    
    if combined[prop_duplicated].isna().any():
        logger.warning(f"{prop_duplicated} col has NA values.")

    return combined


# def plot_duplicates(deduplicate_output, split_by):
#     """
#     Plot the output of deduplication (useful for testing & choosing thresholds)
#     """
#     do = deduplicate_output
#     ax = do['to_keep'].plot(column=split_by, cmap='Dark2', alpha=0.6)
#     do['to_remove'].plot(
#         ax=ax,
#         facecolor='none',
#         edgecolor='k',
#         alpha=0.6,
#         linewidth=0.5)
#     do['intersections'].plot(
#         ax=ax,
#         facecolor='none',
#         edgecolor='red',
#         alpha=0.6,
#         linewidth=0.5)
#     return ax
