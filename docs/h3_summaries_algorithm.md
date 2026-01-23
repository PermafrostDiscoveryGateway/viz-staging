# H3 Grid Summary Generator

Goal: Convert feature geometries into H3 cell indices, then aggregate per-cell metrics.
	•	For any input geometry file: compute "_count" per h3 cell (how many features intersect with the cell).
	•	For polygon features: additionally compute "area_km2" per h3 cell (how much polygon area falls inside the cell).
	•	Optionally compute land-only denominators: "land_area_km2", "land_fraction", and "land_coverage_fraction".

Inputs

	•	input_path: Path to input vector dataset (.shp or .gpkg)
	•	output_path: Path to output GeoPackage
	•	h3_res: H3 resolution (1–15)
	•	land_polygons_path (optional): Land/coastline polygon dataset path
	•	area_epsg: Equal-area EPSG code for area computations (default 6933)

Output

	•	out_gdf: H3 grid (polygons) with per-cell summary attributes written to output_path (GeoPackage)

### Function 1: `build_h3_summary(input_path, output_path, h3_res, land_polygons_path, area_epsg=6933)`

	1.	Read input
		* gdf ← gpd.read_file(input_path)
		* If gdf.crs is missing: raise error
		* If gdf.crs.to_epsg() != 4326: gdf ← gdf.to_crs(epsg=4326)

	2. Initialize
		* records ← []
		* gdf_area ← None

	3. Pre-project for polygon area calculations (if needed)
		3.1. If any geometry type in gdf.geom_type is Polygon or MultiPolygon:
			 * gdf_area ← gdf.to_crs(epsg=area_epsg)

	
