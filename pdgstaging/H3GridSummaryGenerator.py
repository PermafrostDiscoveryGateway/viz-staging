#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path
from typing import Optional, Union, Set, List

import geopandas as gpd
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.validation import make_valid
from shapely.errors import GEOSException
from shapely.ops import transform
import h3
import pandas as pd
import numpy as np

from . import TilePathManager


PathLike = Union[str, Path]


class H3GridSummaryGenerator:
    def __init__(
        self,
        config,
        area_epsg: int = 6933,
        land_polygons_path: Optional[PathLike] = None,
        logger: Optional[logging.Logger] = None,
        tiles: Optional[TilePathManager] = None,
        out_base_dir: str = "h3",
        attr_to_sum: Optional[List[str]] = None,
        attr_to_mean: Optional[List[str]] = None,
    ):
        self.area_epsg = area_epsg
        self.land_polygons_path = land_polygons_path
        self.logger = logger or logging.getLogger(__name__)
        self.tiles = tiles
        self.config = config
        self.out_base_dir = out_base_dir
        self.attr_to_sum = attr_to_sum or []
        self.attr_to_mean = attr_to_mean or []

    def feature_to_h3_index(self, geom, res: int) -> Optional[str]:
        """
        Routes a given geometry to a single H3 cell based on its center point.
        
        Points use their exact coordinates. Polygons, MultiPolygons, and other
        geometries use their computed centroid.
        """
        if geom is None or geom.is_empty:
            return None

        if geom.geom_type == "Point":
            return h3.latlng_to_cell(geom.y, geom.x, res)

        # for Polygons, MultiPolygons, and all other geometries:
        # calculate the centroid and assign the entire feature to that single cell
        centroid = geom.centroid
        return h3.latlng_to_cell(centroid.y, centroid.x, res)

    def h3_to_polygon(self, h: str) -> Polygon:
        """
        Assigns geometry to an H3 cell, needed to get the cells back on a map.
        """
        boundary = h3.cell_to_boundary(h)  # list[(lat, lon)]
        boundary_xy = [(lon, lat) for lat, lon in boundary]
        return Polygon(boundary_xy)

    def add_land_metrics(
        self,
        h3_gdf: gpd.GeoDataFrame,
        land_polygons_path: PathLike,
        area_epsg: int = 6933,
    ) -> gpd.GeoDataFrame:
        """
        Overlays the generated H3 grid with a land polygon dataset to calculate 
        land-based area metrics for each cell.
        
        The method projects both datasets to an equal-area coordinate reference 
        system to ensure accurate area calculations. Calculates the total land area in 
        square kilometers (`land_area_km2`) and the percentage of the cell covered by 
        land (`land_fraction`). If feature area data is present, it also calculates the
        feature's coverage relative strictly to the land area.
        """
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
        input_path: PathLike,
        output_paths: dict,
        h3_res: List[int],
        land_polygons_path: Optional[PathLike] = None, # Not used until combination step, but kept for signature
        area_epsg: Optional[int] = None,
        attr_to_sum: Optional[List[str]] = None,
        attr_to_mean: Optional[List[str]] = None,
    ) -> None:
        if area_epsg is None:
            area_epsg = self.area_epsg
        if attr_to_sum is None:
            attr_to_sum = self.attr_to_sum
        if attr_to_mean is None:
            attr_to_mean = self.attr_to_mean

        input_path = Path(input_path)
        gdf = gpd.read_file(input_path)

        if gdf.crs is None:
            raise ValueError("Input dataset has no CRS; please define or reproject to EPSG:4326.")
        if gdf.crs.to_epsg() != 4326:
            gdf = gdf.to_crs(epsg=4326)

        records_by_res = {res: [] for res in h3_res}

        # filter out polygons duplicated across tiles from staging step
        if 'staging_centroid_within_tile' in gdf.columns:
            gdf = gdf[gdf['staging_centroid_within_tile'] == True]
        if gdf.empty:
            self.logger.info("All features filtered out (centroids in other tiles). Skipping file.")
            return

        has_polygons = any(gt in ("Polygon", "MultiPolygon") for gt in gdf.geom_type.unique())
        if has_polygons:
            gdf["area_km2"] = gdf.to_crs(epsg=area_epsg).geometry.area / 1e6
        else:
            gdf["area_km2"] = 0.0

        for row in gdf.itertuples(index=True):
            idx = row.Index
            geom = row.geometry
            if geom is None or geom.is_empty:
                continue

            # calculate metrics for all resolutions requested while the feature is loaded
            for res in h3_res:
                h3_index = self.feature_to_h3_index(geom, res)
                if not h3_index:
                    continue

                rec = {"h3_index": h3_index, "_count": 1}

                for col in attr_to_sum:
                    rec[f"sum_{col}"] = getattr(row, col)
                for col in attr_to_mean:
                    # TODO: get math right here
                    rec[f"mean_{col}"] = getattr(row, col)

                if has_polygons:
                    rec["area_km2"] = row.area_km2

                records_by_res[res].append(rec)

        # process each resolution into its own PARQUET chunk
        for res, records in records_by_res.items():
            if not records:
                self.logger.warning("No H3 coverage generated for resolution %s. Skipping.", res)
                continue

            df = pd.DataFrame(records)

            agg_dict = {"_count": "sum"}
            for col in attr_to_sum:
                agg_dict[f"sum_{col}"] = "sum"
            for col in attr_to_mean:
                agg_dict[f"mean_{col}"] = "mean"
            if "area_km2" in df.columns:
                agg_dict["area_km2"] = "sum"

            grouped = df.groupby("h3_index", as_index=False).agg(agg_dict)

            # Define output chunks directory (e.g., adjacent to where the final gpkg will live)
            output_path = Path(output_paths[res])
            chunk_dir = output_path.parent / "chunks"
            chunk_dir.mkdir(parents=True, exist_ok=True)
            
            # Save directly to parquet (NO geometries generated here!)
            chunk_path = chunk_dir / f"res_{res}_chunk_{input_path.stem}.parquet"
            grouped.to_parquet(chunk_path)

    def combine_h3_summaries(
        self,
        output_paths: dict,
        h3_res: List[int],
        land_polygons_path: Optional[PathLike] = None,
        area_epsg: Optional[int] = None,
        attr_to_sum: Optional[List[str]] = None,
        attr_to_mean: Optional[List[str]] = None,
    ) -> None:
        """Runs once after all chunks are generated to produce the final GPKGs."""
        if area_epsg is None:
            area_epsg = self.area_epsg
        if land_polygons_path is None:
            land_polygons_path = self.land_polygons_path
        if attr_to_sum is None:
            attr_to_sum = self.attr_to_sum
        if attr_to_mean is None:
            attr_to_mean = self.attr_to_mean

        agg_dict = {"_count": "sum"}
        for col in attr_to_sum:
            agg_dict[f"sum_{col}"] = "sum"
        for col in attr_to_mean:
            agg_dict[f"mean_{col}"] = "mean" # Reminder: true mean requires dividing final sum by final count

        for res in h3_res:
            output_path = Path(output_paths[res])
            chunk_dir = output_path.parent / "chunks"
            
            if not chunk_dir.exists():
                self.logger.warning("No chunks directory found for res %s. Skipping combination.", res)
                continue
                
            chunk_files = list(chunk_dir.glob(f"res_{res}_chunk_*.parquet"))
            if not chunk_files:
                continue
            
            self.logger.info("Combining %d chunks for resolution %s...", len(chunk_files), res)
            
            combined_df = pd.concat([pd.read_parquet(f) for f in chunk_files], ignore_index=True)
            
            if "area_km2" in combined_df.columns:
                agg_dict["area_km2"] = "sum"

            final_grouped = combined_df.groupby("h3_index", as_index=False).agg(agg_dict)
            
            final_grouped["geometry"] = final_grouped["h3_index"].apply(self.h3_to_polygon)
            out_gdf = gpd.GeoDataFrame(final_grouped, geometry="geometry", crs="EPSG:4326")

            if land_polygons_path is not None:
                out_gdf = self.add_land_metrics(out_gdf, land_polygons_path, area_epsg=area_epsg)
            
            out_gdf.to_file(output_path, driver="GPKG", mode="w")
            self.logger.info("Finished writing final GeoPackage for resolution %s.", res)


    
    def valid_h3_resolution(self, value: str) -> int:
        try:
            ivalue = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError("H3 resolution must be an integer.")
        if ivalue < 1 or ivalue > 15:
            raise argparse.ArgumentTypeError("H3 resolution must be between 1 and 15.")
        return ivalue

    
