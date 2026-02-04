# H3GridSummaryGenerator.py
#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

import geopandas as gpd
from shapely.geometry import Polygon
from shapely.ops import unary_union
import h3
import pandas as pd
import numpy as np

from . import TilePathManager


class H3GridSummaryGenerator:
    def __init__(
        self,
        area_epsg: int = 6933,
        land_polygons_path: str | Path | None = None,
        logger: logging.Logger | None = None,
        tiles: TilePathManager | None = None,
        out_base_dir: str = "h3",
        attr_to_sum=None,
        attr_to_mean=None,
    ):
        self.area_epsg = area_epsg
        self.land_polygons_path = land_polygons_path
        self.logger = logger or logging.getLogger(__name__)
        self.tiles = tiles
        self.out_base_dir = out_base_dir
        self.attr_to_sum = attr_to_sum or []
        self.attr_to_mean = attr_to_mean or []

    def polygon_to_h3_cells(self, geom, res: int) -> set[str]:
        if geom is None or geom.is_empty:
            return set()

        try:
            geo = geom.__geo_interface__
        except AttributeError:
            return set()

        cells = set(h3.geo_to_cells(geo, res))
        if not cells:
            rp = geom.representative_point()
            cells = {h3.latlng_to_cell(rp.y, rp.x, res)}
        return cells

    def feature_to_h3_indices(self, geom, res: int) -> set[str]:
        if geom is None or geom.is_empty:
            return set()

        if geom.geom_type == "Point":
            return {h3.latlng_to_cell(geom.y, geom.x, res)}

        if geom.geom_type in ("Polygon", "MultiPolygon"):
            return self.polygon_to_h3_cells(geom, res)

        try:
            geo = geom.__geo_interface__
        except AttributeError:
            return set()

        cells = set(h3.geo_to_cells(geo, res))
        if not cells:
            rp = geom.representative_point()
            cells = {h3.latlng_to_cell(rp.y, rp.x, res)}
        return cells

    def h3_to_polygon(self, h: str) -> Polygon:
        boundary = h3.cell_to_boundary(h)
        boundary_xy = [(lon, lat) for lat, lon in boundary]
        return Polygon(boundary_xy)

    def add_land_metrics(
        self,
        h3_gdf: gpd.GeoDataFrame,
        land_polygons_path: str | Path,
        area_epsg: int = 6933,
    ) -> gpd.GeoDataFrame:
        land = gpd.read_file(land_polygons_path)
        self.logger.info("Finished reading land polygons")

        if land.crs is None:
            raise ValueError("Land polygon dataset has no CRS. Please define it before use.")
        if h3_gdf.crs is None:
            raise ValueError("H3 GeoDataFrame has no CRS. Expected EPSG:4326.")

        h3_ea = h3_gdf.to_crs(epsg=area_epsg)
        land_ea = land.to_crs(epsg=area_epsg)

        h3_ea["cell_area_km2"] = h3_ea.geometry.area / 1e6

        extent_geom = unary_union(h3_ea.geometry)
        land_clip = gpd.overlay(
            land_ea[["geometry"]],
            gpd.GeoDataFrame(geometry=[extent_geom], crs=h3_ea.crs),
            how="intersection",
        )

        inter = gpd.overlay(
            h3_ea[["h3_index", "geometry"]],
            land_clip[["geometry"]],
            how="intersection",
        )

        if len(inter) == 0:
            h3_ea["land_area_km2"] = 0.0
            h3_ea["land_fraction"] = 0.0
        else:
            inter["land_area_km2"] = inter.geometry.area / 1e6
            land_area = inter.groupby("h3_index")["land_area_km2"].sum()
            h3_ea["land_area_km2"] = h3_ea["h3_index"].map(land_area).fillna(0.0)
            h3_ea["land_fraction"] = (h3_ea["land_area_km2"] / h3_ea["cell_area_km2"]).clip(0, 1)

        if {"area_km2", "land_area_km2"}.issubset(h3_ea.columns):
            denom = h3_ea["land_area_km2"].replace(0.0, np.nan)
            h3_ea["land_coverage_fraction"] = (h3_ea["area_km2"] / denom).clip(lower=0, upper=1)

        return h3_ea.to_crs(h3_gdf.crs)

    def build_h3_summary(
        self,
        input_path: str | Path,
        output_path: str | Path,
        h3_res: int,
        land_polygons_path: str | Path | None = None,
        area_epsg: int | None = None,
        attr_to_sum=None,
        attr_to_mean=None,
    ):
        if area_epsg is None:
            area_epsg = self.area_epsg
        if land_polygons_path is None:
            land_polygons_path = self.land_polygons_path

        if attr_to_sum is None:
            attr_to_sum = self.attr_to_sum
        if attr_to_mean is None:
            attr_to_mean = self.attr_to_mean

        gdf = gpd.read_file(input_path)

        if gdf.crs is None:
            raise ValueError("Input dataset has no CRS; please define or reproject to EPSG:4326.")
        if gdf.crs.to_epsg() != 4326:
            gdf = gdf.to_crs(epsg=4326)

        records = []

        gdf_area = None
        if any(gt in ("Polygon", "MultiPolygon") for gt in gdf.geom_type.unique()):
            gdf_area = gdf.to_crs(epsg=area_epsg)

        for idx, row in gdf.iterrows():
            geom = row.geometry
            if geom is None or geom.is_empty:
                continue

            h3_indices = self.feature_to_h3_indices(geom, h3_res)
            if not h3_indices:
                continue

            geom_area = None
            if gdf_area is not None and geom.geom_type in ("Polygon", "MultiPolygon"):
                geom_area = gdf_area.loc[idx].geometry

            for h in h3_indices:
                rec = {"h3_index": h, "_count": 1}

                # optional attribute passthrough
                for col in attr_to_sum:
                    rec[f"sum_{col}"] = row[col]
                for col in attr_to_mean:
                    rec[f"mean_{col}"] = row[col]

                if geom_area is not None:
                    cell_poly = self.h3_to_polygon(h)
                    cell_poly_area = (
                        gpd.GeoSeries([cell_poly], crs="EPSG:4326")
                        .to_crs(epsg=area_epsg)
                        .iloc[0]
                    )
                    inter = geom_area.intersection(cell_poly_area)
                    rec["area_km2"] = (inter.area / 1e6) if not inter.is_empty else 0.0

                records.append(rec)

        if not records:
            raise RuntimeError("No H3 coverage generated. Check geometries / resolution.")

        df = pd.DataFrame(records)

        agg_dict = {"_count": "sum"}
        for col in attr_to_sum:
            agg_dict[f"sum_{col}"] = "sum"
        for col in attr_to_mean:
            agg_dict[f"mean_{col}"] = "mean"
        if "area_km2" in df.columns:
            agg_dict["area_km2"] = "sum"

        grouped = df.groupby("h3_index", as_index=False).agg(agg_dict)

        grouped["geometry"] = grouped["h3_index"].apply(self.h3_to_polygon)
        out_gdf = gpd.GeoDataFrame(grouped, geometry="geometry", crs="EPSG:4326")

        if land_polygons_path is not None:
            out_gdf = self.add_land_metrics(out_gdf, land_polygons_path, area_epsg=area_epsg)

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        out_gdf.to_file(output_path, driver="GPKG")
        self.logger.info(f"Finished writing H3 summary grid to {output_path}")

    def valid_h3_resolution(self, value: str) -> int:
        try:
            ivalue = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError("H3 resolution must be an integer.")
        if ivalue < 1 or ivalue > 15:
            raise argparse.ArgumentTypeError("H3 resolution must be between 1 and 15.")
        return ivalue

    def main(self):
        parser = argparse.ArgumentParser(description="Generate H3 vector summary grids for input vector data.")
        parser.add_argument("input", help="Input vector file (.shp or .gpkg)")
        parser.add_argument("output", help="Output vector file (.gpkg)")
        parser.add_argument("res", type=self.valid_h3_resolution, help="H3 resolution (1–15).")
        parser.add_argument("--sum-cols", nargs="*", default=[], help="Attribute columns to sum by H3 cell")
        parser.add_argument("--mean-cols", nargs="*", default=[], help="Attribute columns to average by H3 cell")
        parser.add_argument("--land-polygons", default=None, help="Optional: path to land/coastline polygon dataset")
        parser.add_argument(
            "--area-epsg",
            type=int,
            default=self.area_epsg,
            help="EPSG code for equal-area projection used for area calculations (default: 6933).",
        )
        args = parser.parse_args()

        self.build_h3_summary(
            input_path=args.input,
            output_path=args.output,
            h3_res=args.res,
            attr_to_sum=args.sum_cols,
            attr_to_mean=args.mean_cols,
            land_polygons_path=args.land_polygons,
            area_epsg=args.area_epsg,
        )
