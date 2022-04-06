import operator
import geopandas as gpd
import pandas as pd
import uuid
import itertools


def deduplicate(
    gdf,
    prop_overlap=None,
    prop_area=None,
    prop_centroid_x=None,
    prop_centroid_y=None,
    keep_rules=[],
    overlap_tolerance=0.5,
    overlap_both=True,
    centroid_tolerance=None,
    distance_crs='EPSG:3857',
    return_intersections=False
):
    """

    Remove duplicated polygons from two sources in a single GeoDataFrame, even
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

    When two polygons are duplicated, only one is kept in the returned 'keep'
    dataframe. Which one is kept is determined by the 'keep_rules' parameter,
    which indicates the property to compare, and the rule to apply ('larger' or
    'smaller'). For example, to keep the newest polygon, set a keep_rule for a
    more recent date property, e.g. ('Date', 'larger').

    A dictionairy is returned with the deduplicated GeoDataFrame ('keep'), the
    GeoDataFrame with the dropped duplicates ('remove'), and optionally, the
    intersection geometries of duplicated polygons ('intersection').

    Parameters
    ----------
    gdf : GeoDataFrame
        The GeoDataFrame to deduplicate.
    prop_overlap : str
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
        duplicates. A list of tuples of the form (property, operator). The
        property is the property that exists in the geodataframe to use for the
        comparison. The operator is either 'larger' or 'smaller'. If the rule
        is 'larger', the polygon with the largest value for the property will
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

    Returns
    -------
    dict
        A dictionary with the following keys:

        - 'keep' : GeoDataFrame
            The deduplicated GeoDataFrame.
        - 'removed' : GeoDataFrame
            The GeoDataFrame with the removed features.
        - 'intersections' : GeoDataFrame
            The GeoDataFrame with the intersections.
    """

    if prop_overlap is None:
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
    sort_props = [x[0] for x in keep_rules]
    sort_order = [x[1] == 'smaller' for x in keep_rules]

    # Remove all unnecessary columns
    working_cols = [
        'geometry',
        prop_id,
        prop_overlap,
        prop_area,
        prop_centroid_x,
        prop_centroid_y] + sort_props
    working_cols = [i for i in working_cols if i is not None]
    gdf.drop(gdf.columns.difference(working_cols), axis=1, inplace=True)

    # Split the gdf by the overlap property (e.g. source file)
    gdf_grouped = gdf.groupby(prop_overlap)
    gdfs = [gdf_grouped.get_group(x) for x in gdf_grouped.groups]

    to_return = {
        'keep': None,
        'intersections': None,
        'remove': None
    }

    if len(gdfs) < 2:
        to_return['keep'] = gdf_original.drop(columns=[prop_id])
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
        duplicates['remove'] = duplicates.duplicated(prop_int_id, keep='first')

        # Get the IDs of polygons to remove, and add them to the list
        ids_to_remove_pair = duplicates[duplicates['remove']][prop_id]
        ids_to_remove.update(ids_to_remove_pair)

    # Remove polygons from the original gdf when they have an ID that is in the
    # to_remove list
    remove = gdf_original[prop_id].isin(ids_to_remove)

    # Remove the prop_id column, no longer needed
    gdf_original.drop([prop_id], axis=1, inplace=True)

    to_return['remove'] = gdf_original[remove]
    to_return['keep'] = gdf_original[~remove]

    return to_return


def plot_duplicates(deduplicate_output, prop_overlap):
    """
    Plot the output of deduplication (useful for testing & choosing thresholds)
    """
    do = deduplicate_output
    ax = do['keep'].plot(column=prop_overlap, cmap='Dark2', alpha=0.6)
    do['remove'].plot(
        ax=ax,
        facecolor='none',
        edgecolor='k',
        alpha=0.6,
        linewidth=0.5)
    do['intersections'].plot(
        ax=ax,
        facecolor='none',
        edgecolor='red',
        alpha=0.6,
        linewidth=0.5)
    return ax
