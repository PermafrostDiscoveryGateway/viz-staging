# viz-staging

## Vector file staging for the PDG tiling pipeline

This repository contains code that prepares vector data (e.g. shapefiles) for subsequent steps in the [PDG](https://permafrost.arcticdata.io/) tiling pipeline (such as [viz-3dtiles](https://github.com/PermafrostDiscoveryGateway/viz-3dtiles) and [viz-raster](https://github.com/PermafrostDiscoveryGateway/viz-raster)). The staging step creates output vector files that conform to a specified [OGC Two Dimensional Tile Matrix Set](http://docs.opengeospatial.org/is/17-083r2/17-083r2.html) ("TMS"). Specifically, for each input file, the staging process:

1. Simplifies polygons and re-projects them to the Coordinate Reference System ("CRS") used by the desired TMS.
2. Assigns area, centroid, and other properties to each polygon.
3. Saves polygons to one file for each tile in the specified level of the TMS.

Polygons are assigned to a tile file if the polygon is within the tile or if it intersects with the bounding box of the tile (i.e. if it is at least *partially* within that tile). This means that polygons that fall within two or more tiles will be duplicated in the output. (This allows subsequent rasterization steps to measure the area of polygons that are only partially within the tile - otherwise some area is lost.)

However, polygon-tile relationships are also identified using the centroid of each polygon: The `centroid_tile` property assigned to polygons identifies the tile within which the polygon's centroid falls. (In the rare event that a polygon's centroid falls exactly on a tile boundary, the polygon will be added to the southern/eastern tile.)

The `centroid_within_tile` property is `True` when the polygon's centroid is within the same tile as the output file. To avoid using duplicated polygons in subsequent tiling steps (e.g. when generating 3D tiles), first filter out all polygons where `centroid_within_tile` is `False`.

The area and centroid of each polygon are calculated in the CRS of the TMS. When this is a geographic coordinate rather than a projected coordinate, the resulting values may vary slightly from the real values. 

## Polygon properties

After being run through this staging process, each polygon will have the following new properties:
- **tile** (string) - The tile that this polygon has been assigned to (matches the file). Indicated as a string in the format `Tile(x=6, y=10, z=4)`.
- **centroid_tile** (string) - The tile within which the centroid of this polygon falls (may differ from file). Indicated as a string in the format `Tile(x=6, y=10, z=4)`.
- **area** (number) - The area of the tile in the units of the CRS
- **centroid_x**, **centroid_y** (number) - The x and y coordinate for the centroid of the tile (in the units of the CRS)
- **filename** (string) - The path to the file from which this polygon originated
- **identifier** (string) - A unique identifier for the polygon
- **centroid_within_tile** (boolean) - True when the `centroid_tile` property matches the `tile`, i.e. when the centroid of the polygon is within the same tile as the file it is saved in

## Summary fields

The staging process will also output a summary CSV file with one row for each tile created from each file. The fields in the CSV are:
- **filename**: The path to the input vector file that was processed
- **num_polygons**: The number of polygons from that file that are in given tile
- **area_polygons**: The total area of the polygons from the file that are in the given tile
- **datetime**: The time at which polygons from this file were saved to the output tile
- **tile_x, tile_y, tile_z**: The x, y, and z indices for the tile
- **tile_top, tile_right, tile_bottom, tile_left**: The coordinates (bounding box) for the tile
- **tms_identifier**: The TMS that this tile belongs to

## Assumptions

- It is assumed that incoming vector data comprises only valid polygons. This has not yet been tested with multi-polygons, points, lines, or other geometries.
- It's also assumed that each incoming vector file is staged only once. If a file passes through the staging step twice, then all polygons from that file will be duplicated in the output (but with a different identifier). This is due to the fact that when a tile file already exists, additional polygons that belong to this tile will be appended to the file.

## Running

Uses Python version `3.9` and packages listed under `requirements.txt`

1. Make sure you have Python version 3.9 installed and that you know the path to that version of Python (try `which python3.9`).
2. Clone this repository.
3. From within the newly created directory, create a virtual environment: `/usr/local/bin/python3.9 -m venv .venv` (where `/usr/local/bin/python3.9` is your path to version 3.9 of Python).
4. Activate the environment: `source .venv/bin/activate`.
5. Install dependencies: `pip install -r requirements.txt`.
6. Run: `python main.py -i /path/to/input/shapefiles -o /path/to/output/directory`.
