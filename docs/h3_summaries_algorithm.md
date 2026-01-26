# H3 Grid Summary Generator

Goal: Convert feature geometries into H3 cell indices, then aggregate per-cell metrics.

1. For any input geometry file: compute "count" per h3 cell (how many features intersect with the cell).
2. For polygon features: additionally compute "area_km2" per h3 cell (how much polygon area falls inside the cell).
3. Optionally compute land-only statistics: "land_area_km2", "land_fraction", and "land_coverage_fraction".

### Inputs

* `input_path` : Path to input vector dataset (.shp or .gpkg)
* `output_path` : Path to output GeoPackage
* `h3_res` : H3 resolution (1–15)
* `land_polygons_path` (optional): Land/coastline polygon dataset path
* `area_epsg` : Equal-area EPSG code for area computations (default 6933)

### Inputs

* `out_gdf` : H3 polygon layer with aggregated attributes written to output_path (GeoPackage)

### Function 1: build_h3_summary(input_path, output_path, h3_res, land_polygons_path, area_epsg=6933)

1. Load and normalize CRS
    * Read gdf = `gpd.read_file(input_path)``.
    * Ensure geometry is in WGS84 (EPSG:4326) because H3 indexing expects lon/lat.

2. Decide which metrics are to be computed
    * Always compute "count" (feature-to-cell intersection count).
    * Compute "area_km2" only if input is set of polygonal features: create `gdf_area = gdf.to_crs(epsg=area_epsg)`` for reliable area computations.

3. Generate each feature's statistics corresponding to each cell into a `records` table
    * Initialize `records = []`.
    * For each feature (row) with geometry geom:
        - Compute the set of H3 cells that the feature intersects: 
            • `h3_indices = feature_to_h3_indices(geom, h3_res)``
            • (Logic: point → one cell; polygon/line → polyfill; if polyfill empty → fallback to representative point cell)
        -  For each `h3_index = h in h3_indices`, append a per-cell statistic:
            • Always append `{"h3_index": h, "_count": 1}`
            • If geom is polygonal, compute how much of that polygon lies inside the H3 cell:
                - `cell_poly = h3_to_polygon(h)` (WGS84)
                - Reproject `cell_poly` to area_epsg
                - `inter = intersection(geom_area, cell_poly_area)`
                - Append `area_km2 = area(inter)/1e6` (or 0.0 if empty)

    Note: `records` is a table where each feature can generate multiple rows (one per H3 cell).

4. Aggregate contributions per H3 cell