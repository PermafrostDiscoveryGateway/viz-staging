# Deduplication

If input files contain overlapping regions with the same features, then the
resulting tiles will contain duplicate polygons. When polygons are duplicated,
then the resulting raster & 3d tile products will tend to over-count polygons
in the areas of overlap.

Tiles may be deduplicated using one of two methods:

1. Deduplication by source file footprint overlap ("footprints")
2. Deduplication by neighbor proximity ("neighbor")

Neither method is perfect: The "footprints" method will fairly evenly sample
all areas, but will tend to under-count and underestimate the size of polygons
at file boundaries. The "neighbor" method may over-sample areas of overlap, but
will more accurately reflect the polygons in the files.

### Example

Here are two overlapping files with thick lines representing the boundaries of those files:

![input files with boundary lines](images/tiles_input-with-footprint.png)

Without any deduplication, here is what the resulting tiles will look like:

![tiles without any deduplication](images/tiles_not_deduped.png)

Deduplicating with the **footprints** method removes all polygons from one file where the two files overlap:

![tiles deduplicated by the footprints method](images/tiles_deduped_footprint.png)

Deduplicating with the **neighbor** method leaves polygons from both files in the area where two files overlap,
but removes polygons from one file where they overlap or are nearby a polygon from another file:

![tiles deduplicated by the neighbor method](images/tiles_deduped_neighbor.png)

## Common options

Regardless of which method is selected, the following options should be set in the configuration:

* `deduplicate_at` : (list of str or None)
    When to deduplicate the input data. Options are 'staging',
    'raster', '3dtiles', or None to skip deduplication.

* `deduplicate_method` : Set to `neighbor` or `footprints`

* `deduplicate_keep_rules` : (list of tuple: []) Rules that determine which
    polygons should be kept and which should be removed when they are
    considered duplicates. The keep rules are provided as a list of tuples in
    the format `(property, operator)`, where each tuple is a rule. `property` is a
    property that exists in the footprint file (in the case of 'footprint'
    method) or in the input data (in the case of the 'neighbor' method) that
    can be used to compare polygons/footprints. `operator` can be set to
    'larger' or 'smaller'. If the rule is 'larger', then polygon with the
    largest value for the property will be kept, and vice versa for smaller.
    When two properties are equal, then the next property in the list will be
    checked.

# Footprint method

The footprint method removes polygons that originate from one source file in
areas where two source files overlap. To use this method, each input vector
file must have a matching footprint vector file with the same name, located in
the directory specified by `dir_footprints`. Each footprint file should contain
one or more polygons that specify the boundary of the associate input file. The
footprint polygons must also have a property that can be used to rank the files
in order of preference, such that only polygons from the preferred file are
kept in areas of overlap. (This property should be specified in `keep_rules`).

For more details, see [footprints for deduplication](footprints.md)

# Neighbor method

The neighbor method removes nearby or overlapping polygons when they come from
a different source file.

The following configuration options are available for this method:

* `deduplicate_overlap_tolerance` : (float, optional)
    The minimum proportion of a polygon's area that must be
    overlapped by another polygon to be considered a duplicate.
    Default is 0.5. Set to None to ignore overlap proportions
    when comparing polygons, and set a centroid threshold
    instead. Note that both an overlap_tolerance AND a
    centroid_tolerance can be used.

* `deduplicate_overlap_both` : (bool, optional)
    If True, then the overlap_tolerance proportion must be True
    for both of the intersecting polygons to be considered a
    duplicate. If False, then the overlap_tolerance proportion
    must be True for only one of the intersecting polygons to
    be considered a duplicate. Default is True.

* `deduplicate_centroid_tolerance` : (float, optional)
    The maximum distance between the centroids of two polygons
    to be considered a duplicate. Default is None. Set to None
    to ignore centroid distances when comparing polygons, and
    set an overlap threshold instead. Note that both an
    overlap_tolerance AND a centroid_tolerance can be used. The
    unit of the distance is the unit of the distance_crs
    property (e.g. meters for EPSG:3857), or the unit of the
    GeoDataFrame if distance_crs is None.

* `deduplicate_distance_crs` : (str, optional)
    The CRS to use for the centroid distance calculation.
    Default is EPSG:3857. Centroid points will be re-projected
    to this CRS before calculating the distance between them.
    centroid_tolerance will use the units of this CRS. Set to
    None to skip the re-projection and use the CRS of the
    GeoDataFrame.

## Example 

```python
overlap_tolerance=0.5
overlap_both=True
centroid_tolerance=5
```

### Before

![Before example](images/iwp_before.png)

### After

![After example](images/iwp_deduplicated.png)


## How it works

![two-polygons-two-files](images/1_two-polygons-two-files.png)

### Configuration options

![centroid-tolerance](images/2_centroid-tolerance.png)

![overlap-tolerance](images/3_overlap-tolerance.png)

### Minimal examples

![same-source-file](images/4_same-source-file.png)

![centroid-tolerance-example](images/5_centroid-tolerance-example.png)

![overlap-tolerance-example](images/6_overlap-tolerance-example.png)

![centroid-and-overlap-example](images/7_centroid-and-overlap-example.png)

![overlap-both-example](images/7_overlap-both-example.png)

## Testing

```python
import geopandas as gpd
import pandas as pd
from pdgstaging.Deduplicator import deduplicate, plot_duplicates

# Import two files that overlap
file1 = 'path/to/first/vector/file.shp'
file2 = 'path/to/second/vector/file.shp'

gdf1 = gpd.read_file(file1)
gdf2 = gpd.read_file(file2)

# Add properties (if they don't already exist)
gdf1['source_file'] = file1
gdf2['source_file'] = file2

gdf1['Date'] = '2019-01-01'
gdf2['Date'] = '2020-01-01'

# Combine the two files
gdf = pd.concat([gdf1, gdf2])

gdf['Area'] = gdf.area
gdf['Centroid_x'] = gdf1.centroid.x
gdf['Centroid_y'] = gdf1.centroid.y

# Deduplicate
output = deduplicate_neighbors(
    gdf,
    split_by='source_file',
    prop_area='Area',
    prop_centroid_x='Centroid_x',
    prop_centroid_y='Centroid_y',
    keep_rules=[['Date', 'larger']],
    overlap_tolerance=0.5,
    overlap_both=True,
    centroid_tolerance=5,
    distance_crs='EPSG:3413',
    return_intersections=True
)

# See the output
plot_duplicates(output, 'source_file')

# Save the deduplicated file
output['keep'].to_file('path/to/output/file.shp')

```
