H3 Grid Summary Generator

Inputs

	•	input_path: Path to input vector dataset (.shp or .gpkg)
	•	output_path: Path to output GeoPackage
	•	h3_res: H3 resolution (1–15)
	•	land_polygons_path (optional): Land/coastline polygon dataset path
	•	area_epsg: Equal-area EPSG code for area computations (default 6933)

Output

	•	out_gdf: H3 grid (polygons) with per-cell summary attributes written to output_path (GeoPackage)

Algorithm 1: build_h3_summary(input_path, output_path, h3_res, land_polygons_path=None, area_epsg=6933)

	1.	Read input
	1.	gdf ← gpd.read_file(input_path)
	2.	If gdf.crs is missing: raise error
	3.	If gdf.crs.to_epsg() != 4326: gdf ← gdf.to_crs(epsg=4326)
