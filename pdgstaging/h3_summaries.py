#!/usr/bin/env python3

"""
H3 Grid Summary Generator
-------------------------

This script takes an input vector dataset (e.g., .shp, .gpkg),
computes H3 cell coverage at a specified resolution, aggregates
feature-level attributes into H3 cells, and writes the resulting H3
grid as a vector .gpkg layer.

Notes
----------------

- Input geometries are reprojected to WGS84 (EPSG:4326) for H3 indexing.
- H3 resolution must be an integer between 1 and 15 (validated via
  argparse).
- For polygons that contain no H3 cell centroids at coarse resolutions,
  a fallback assigns the H3 cell containing the polygon’s representative
  interior point to ensure complete coverage.
- Area calculations are done in an equal-area CRS (default EPSG:6933).
  For Arctic-datasets it maybe preferable to switch to EPSG:3575.
- If the input geometry is Polygon/MultiPolygon, computes the **area of
  feature overlap within each H3 cell** (km²) and aggregates it per cell.
- Optionally computes **terrestrial (land) area** of each H3 cell using a
  coastline/land polygon dataset (e.g., Natural Earth land polygons), and
  adds:
    - land_area_km2 (total area of input polygons overlapping an H3 cell)
    - land_fraction
    - count_per_land_km2
    - area_per_land_km2 (if polygon overlap area was computed)

H3 Cell Summaries
--------------------
The following summary statistics are computed:

- _count - The number of input features intersecting an H3 cell.
- sum_area_m2 - The total area (in square meters) of polygon features intersecting the H3 cell.
                Only applies when input geometry is polygonal.
- area_km2 - The same value as sum_area_m2, but converted to square kilometers.
- cell_area_km2 - The total surface area of the H3 cell itself, in km².
- land_area_km2 - The portion of the H3 cell that is land, excluding water.
- land_fraction - The fraction of the H3 cell that is land (land_fraction = land_area_km2 / cell_area_km2).
- area_per_land_km2 - Fraction of land in a H3 cell that is covered by features.


Author
------
Shirly Stephen
"""


import argparse
from pathlib import Path

import geopandas as gpd
from shapely.geometry import Polygon
from shapely.ops import unary_union
import h3
import pandas as pd
import numpy as np


def feature_to_h3_indices(geom, res: int) -> set[str]:
    """Return a set of H3 indices for the geometry at resolution res.

    - Point: single H3 cell using lat/lon.
    - Polygon/MultiPolygon: polyfill coverage (geo_to_cells).
    - Other: attempts geo_to_cells coverage; falls back to representative point.
    """
    if geom is None or geom.is_empty:
        return set()

    if geom.geom_type == "Point":
        return {h3.latlng_to_cell(geom.y, geom.x, res)}

    if geom.geom_type in ("Polygon", "MultiPolygon"):
        return polygon_to_h3_cells(geom, res)

    # LineString, MultiLineString, GeometryCollection, etc.
    try:
        geo = geom.__geo_interface__
    except AttributeError:
        return set()

    cells = set(h3.geo_to_cells(geo, res))

    # Optional: keep a fallback for weird edge-cases
    if not cells:
        rp = geom.representative_point()
        cells = {h3.latlng_to_cell(rp.y, rp.x, res)} 

    return cells

def polygon_to_h3_cells(geom, res: int) -> set[str]:
    """Return H3 cells covering a Polygon/MultiPolygon at the given resolution.

    Uses H3 polyfill (geo_to_cells). If polyfill returns empty at very coarse
    resolutions, falls back to the H3 cell containing a representative point.
    """
    if geom is None or geom.is_empty:
        return set()

    try:
        geo = geom.__geo_interface__  # GeoJSON mapping (lon/lat order)
    except AttributeError:
        return set()

    cells = set(h3.geo_to_cells(geo, res))

    if not cells:
        rp = geom.representative_point()
        cells = {h3.latlng_to_cell(rp.y, rp.x, res)}
    return cells

def h3_to_polygon(h: str) -> Polygon:
    """Convert an H3 cell ID to a Shapely Polygon in lon/lat (EPSG:4326)."""
    boundary = h3.cell_to_boundary(h)  # list[(lat, lon)]
    boundary_xy = [(lon, lat) for lat, lon in boundary]
    return Polygon(boundary_xy)

def add_land_metrics(
    h3_gdf: gpd.GeoDataFrame,
    land_polygons_path: str | Path,
    area_epsg: int = 6933,
) -> gpd.GeoDataFrame:
    """Compute land_area_km2 and land_fraction for each H3 cell.

    Option A: polygon intersection between H3 cell polygons and a land polygon dataset.
    """
    land = gpd.read_file(land_polygons_path)
    print("Finished reading land polygons");
    if land.crs is None:
        raise ValueError("Land polygon dataset has no CRS. Please define it before use.")
    if h3_gdf.crs is None:
        raise ValueError("H3 GeoDataFrame has no CRS. Expected EPSG:4326.")

    # Project to equal-area CRS for area computation
    h3_ea = h3_gdf.to_crs(epsg=area_epsg)
    land_ea = land.to_crs(epsg=area_epsg)

    # Compute full cell area (km²)
    h3_ea["cell_area_km2"] = h3_ea.geometry.area / 1e6

    # Clip land polygons to the H3 extent to speed up overlay
    extent_geom = unary_union(h3_ea.geometry)
    land_clip = gpd.overlay(
        land_ea[["geometry"]],
        gpd.GeoDataFrame(geometry=[extent_geom], crs=h3_ea.crs),
        how="intersection",
    )

    # Intersect each H3 cell with land polygons
    inter = gpd.overlay(
        h3_ea[["h3_index", "geometry"]],
        land_clip[["geometry"]],
        how="intersection",
    )
    if len(inter) == 0:
        # No land overlap found (e.g., all ocean) – fill zeros
        h3_ea["land_area_km2"] = 0.0
        h3_ea["land_fraction"] = 0.0
    else:
        inter["land_area_km2"] = inter.geometry.area / 1e6
        land_area = inter.groupby("h3_index")["land_area_km2"].sum()

        h3_ea["land_area_km2"] = h3_ea["h3_index"].map(land_area).fillna(0.0)
        h3_ea["land_fraction"] = (
            h3_ea["land_area_km2"] / h3_ea["cell_area_km2"]
        ).clip(0, 1)
    

##    # Add normalized count metric (avoid div-by-zero)
##    denom = h3_ea["land_area_km2"].replace({0.0: pd.NA})
##    if "_count" in h3_ea.columns:
##        h3_ea["count_per_land_km2"] = (h3_ea["_count"] / denom).astype("float")

    # If polygon overlap area exists, normalize it too
##    if "area_km2" in h3_ea.columns:
##        h3_ea["area_per_land_km2"] = (h3_ea["area_km2"] / denom).astype("float")
    
    if {"area_km2", "land_area_km2"}.issubset(h3_ea.columns):
        denom = h3_ea["land_area_km2"].replace(0.0, np.nan)

        h3_ea["land_coverage_fraction"] = h3_ea["area_km2"] / denom

        # Optional but recommended: cap to [0, 1] to behave like a fraction
        h3_ea["land_coverage_fraction"] = h3_ea["land_coverage_fraction"].clip(lower=0, upper=1)


    # Return to EPSG:4326 for output consistency
    return h3_ea.to_crs(h3_gdf.crs)


def build_h3_summary(
    input_path: str | Path,
    output_path: str | Path,
    h3_res: int,
##    attr_to_sum=None,
##    attr_to_mean=None,
    land_polygons_path: str | Path | None = None,
    area_epsg: int = 6933,
):
    """Build an H3 summary grid for an input vector dataset."""

    gdf = gpd.read_file(input_path)

    # Ensure CRS is WGS84 for H3 indexing
    if gdf.crs is None:
        raise ValueError("Input dataset has no CRS; please define or reproject to EPSG:4326.")
    if gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(epsg=4326)

##    attr_to_sum = attr_to_sum or []
##    attr_to_mean = attr_to_mean or []

    records = []

    # Pre-project to equal-area once for polygon area calculations (if needed)
    gdf_area = None
    if any(gt in ("Polygon", "MultiPolygon") for gt in gdf.geom_type.unique()):
        gdf_area = gdf.to_crs(epsg=area_epsg)

    for idx, row in gdf.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue

        h3_indices = feature_to_h3_indices(geom, h3_res)
        if not h3_indices:
            continue

        # If polygon, use area-projected geometry for intersections
        geom_area = None
        if gdf_area is not None and geom.geom_type in ("Polygon", "MultiPolygon"):
            geom_area = gdf_area.loc[idx].geometry

        for h in h3_indices:
            rec = {"h3_index": h, "_count": 1}

            # Attribute aggregation passthrough
##            for col in attr_to_sum:
##                rec[f"sum_{col}"] = row[col]
##            for col in attr_to_mean:
##                rec[f"mean_{col}"] = row[col]

            # Polygon overlap area per cell (km²)
            if geom_area is not None:
                cell_poly = h3_to_polygon(h)
                cell_poly_area = gpd.GeoSeries([cell_poly], crs="EPSG:4326").to_crs(epsg=area_epsg).iloc[0]
                inter = geom_area.intersection(cell_poly_area)
                if not inter.is_empty:
                    rec["area_km2"] = inter.area / 1e6
                else:
                    rec["area_km2"] = 0.0

            records.append(rec)

    if not records:
        raise RuntimeError("No H3 coverage generated. Check geometries / resolution.")

    df = pd.DataFrame(records)

    # Aggregate by H3 index
    agg_dict = {"_count": "sum"}
    for col in attr_to_sum:
        agg_dict[f"sum_{col}"] = "sum"
    for col in attr_to_mean:
        agg_dict[f"mean_{col}"] = "mean"
    if "area_km2" in df.columns:
        agg_dict["area_km2"] = "sum"

    grouped = df.groupby("h3_index", as_index=False).agg(agg_dict)

    # Build H3 polygons
    grouped["geometry"] = grouped["h3_index"].apply(h3_to_polygon)
    out_gdf = gpd.GeoDataFrame(grouped, geometry="geometry", crs="EPSG:4326")

    # Optional: compute land normalization metrics
    if land_polygons_path is not None:
        out_gdf = add_land_metrics(out_gdf, land_polygons_path, area_epsg=area_epsg)

    # Write to file
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_gdf.to_file(output_path, driver="GPKG")
    print(f"Finished writing H3 summary grid to {output_path}")

            

def valid_h3_resolution(value: str) -> int:
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("H3 resolution must be an integer.")

    if ivalue < 1 or ivalue > 15:
        raise argparse.ArgumentTypeError(
            "H3 resolution must be between 1 and 15."
        )

    return ivalue

def main():
    parser = argparse.ArgumentParser(
        description="Generate H3 vector summary grids for input vector data."
    )
    parser.add_argument("input", help="Input vector file (.shp or .gpkg)")
    parser.add_argument("output", help="Output vector file (.gpkg)")
    parser.add_argument(
        "res",
        type=valid_h3_resolution,
        help="H3 resolution (1–15). See https://h3geo.org/docs/core-library/restable/"
    )
    parser.add_argument(
        "--sum-cols",
        nargs="*",
        default=[],
        help="Attribute columns to sum by H3 cell",
    )
    parser.add_argument(
        "--mean-cols",
        nargs="*",
        default=[],
        help="Attribute columns to average by H3 cell",
    )
    parser.add_argument(
        "--land-polygons",
        default=None,
        help="Optional: path to a land/coastline polygon dataset (e.g., Natural Earth land polygons) "
             "used to compute land_area_km2 and land_fraction per H3 cell.",
    )
    parser.add_argument(
        "--area-epsg",
        type=int,
        default=6933,
        help="EPSG code for equal-area projection used for area calculations (default: 6933). "
             "For pan-Arctic work consider 3575.",
    )


    args = parser.parse_args()

    build_h3_summary(
        input_path=args.input,
        output_path=args.output,
        h3_res=args.res,
##        attr_to_sum=args.sum_cols,
##        attr_to_mean=args.mean_cols,
        land_polygons_path=args.land_polygons,
        area_epsg=args.area_epsg,
    )

if __name__ == "__main__":
    TEST_MODE = True  # switch to True when testing
    print("Control here");

    if TEST_MODE:
        build_h3_summary(
            input_path="data/DARTS_NitzeEtAl_v1-2_features_2018-2023_level1.gpkg",
            output_path="data/h3/darts_2001-2023_level1_h3level9_alternate.gpkg",
            h3_res=9, # h3 cell counts and cell areas (https://h3geo.org/docs/core-library/restable/)
##            attr_to_sum=["area_m2"],
##            attr_to_mean=[],
            land_polygons_path="data/land/",
            area_epsg=6933,
            
        )
    else:
        main()
