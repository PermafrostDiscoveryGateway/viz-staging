#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path
from typing import Optional, Union, List

import geopandas as gpd
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely import transform
import pyproj
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

    def feature_to_h3_index(self, geom, res: int, feature_split: bool) -> Optional[Union[str, set[str]]]:
        """Routes a given geometry to H3 cells based on the split configuration.
        
        Points use their exact coordinates. For polygons and other geometries:
        - If feature_split is False, routes to a single H3 cell based on the centroid.
        - If feature_split is True, routes to all intersecting H3 cells.

        Args:
            geom: The Shapely geometry object to process (e.g., Point, Polygon).
            res (int): The H3 resolution level to use for the cell index.
            feature_split (bool): Whether to return a set of all intersecting cells 
                or just a single cell based on the centroid.

        Returns:
            Optional[Union[str, set[str]]]: A single H3 string if feature_split=False, 
                a set of H3 strings if feature_split=True, or None if geometry is empty.
        """
        if geom is None or geom.is_empty:
            return None

        if geom.geom_type == "Point":
            cell = h3.latlng_to_cell(geom.y, geom.x, res)
            return {cell} if feature_split else cell

        # for Polygons, MultiPolygons, and all other geometries:
        # if no feature split, calculate the centroid and assign the entire feature to that single cell
        centroid = geom.centroid
        if not feature_split:
            return h3.latlng_to_cell(centroid.y, centroid.x, res)
        # else, get all cells associated with it
        else:
            try:
                # convert Shapely (x, y) to H3Shape (lat, lng) for experimental fn
                if geom.geom_type == 'Polygon':
                    outer = [(y, x) for x, y in geom.exterior.coords]
                    holes = [[(y, x) for x, y in interior.coords] for interior in geom.interiors]
                    h3_geom = h3.LatLngPoly(outer, *holes)
                    
                elif geom.geom_type == 'MultiPolygon':
                    h3_polys = []
                    for p in geom.geoms:
                        outer = [(y, x) for x, y in p.exterior.coords]
                        holes = [[(y, x) for x, y in interior.coords] for interior in p.interiors]
                        h3_polys.append(h3.LatLngPoly(outer, *holes))
                    h3_geom = h3.LatLngMultiPoly(*h3_polys)
                    
                else:
                    h3_geom = geom  # fallback for unexpected types

                # catch all overlapping cells using the native H3 shape
                # this ensures we get all H3 cells that touch the polygon, not just ones that touch the 
                # center of the H3 cell
                cells = set(h3.h3shape_to_cells_experimental(h3_geom, res, contain='overlap'))
                
            except (AttributeError, ValueError) as e:
                # fallback to standard fill if experimental fails or isn't supported
                cells = set(h3.geo_to_cells(geom, res))
            if not cells:
                cells = {h3.latlng_to_cell(centroid.y, centroid.x, res)}
            return cells


    def h3_to_polygon(self, h: str) -> Polygon:
        """Assigns geometry to an H3 cell, needed to get the cells back on a map.

        Converts the boundary of an H3 cell index from latitude/longitude 
        coordinates into a projected Shapely Polygon.

        Args:
            h (str): The H3 cell index string to convert.

        Returns:
            Polygon: A Shapely Polygon representing the geographic boundary 
                of the H3 cell.
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
        """Overlays the generated H3 grid with a land polygon dataset to calculate 
        land-based area metrics for each cell.
        
        The method projects both datasets to an equal-area coordinate reference 
        system to ensure accurate area calculations. Calculates the total land area in 
        square kilometers (`land_area_km2`) and the percentage of the cell covered by 
        land (`land_fraction`). If feature area data is present, it also calculates the
        feature's coverage relative strictly to the land area.

        Args:
            h3_gdf (gpd.GeoDataFrame): The GeoDataFrame containing H3 grid cells.
            land_polygons_path (PathLike): The file path to the vector dataset 
                representing land boundaries.
            area_epsg (int, optional): The EPSG code for the equal-area coordinate 
                reference system used for accurate area calculations. Defaults to 6933.

        Returns:
            gpd.GeoDataFrame: The updated GeoDataFrame with appended land metric 
                columns (`land_area_km2`, `land_fraction`, and potentially others).
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
        input_path: Union[str, Path],
        output_paths: dict,
        h3_res: List[int],
        feature_split: bool = False,
        land_polygons_path: Optional[Union[str, Path]] = None,
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

        if 'staging_centroid_within_tile' in gdf.columns:
            gdf = gdf[gdf["staging_centroid_within_tile"].astype(bool)]
        if gdf.empty:
            self.logger.info("All features filtered out (centroids in other tiles). Skipping file.")
            return

        if 'staging_duplicated' in gdf.columns:
            gdf = gdf[gdf["staging_duplicated"].astype(bool)]
        if gdf.empty:
            self.logger.info("All features filtered out (duplicated features). Skipping file.")
            return

        has_polygons = any(gt in ("Polygon", "MultiPolygon") for gt in gdf.geom_type.unique())
        if has_polygons:
            gdf["area_km2"] = gdf.to_crs(epsg=area_epsg).geometry.area / 1e6
        else:
            gdf["area_km2"] = 0.0

        # Pre-compile the reprojection transformer for fast sliver area calculations
        project_to_area = None
        if feature_split and has_polygons:
            project_to_area = pyproj.Transformer.from_crs(
                pyproj.CRS("EPSG:4326"), 
                pyproj.CRS(f"EPSG:{area_epsg}"), 
                always_xy=True
            ).transform

        for row in gdf.itertuples(index=True):
            geom = row.geometry
            if geom is None or geom.is_empty:
                continue

            for res in h3_res:
                h3_indices = self.feature_to_h3_index(geom, res, feature_split)
                if not h3_indices:
                    continue
                
                # Normalize to an iterable so we can handle both single strings and sets cleanly
                if isinstance(h3_indices, str):
                    h3_indices = {h3_indices}
                
                needs_split = (
                    feature_split 
                    and has_polygons 
                    and geom.geom_type not in ("Point", "MultiPoint") 
                    and len(h3_indices) > 1
                )
                count_cell = None
                # get the cell to assign the feature count to
                if needs_split:
                    centroid = geom.centroid
                    centroid_cell = h3.latlng_to_cell(centroid.y, centroid.x, res)
                    if centroid_cell in h3_indices:
                        count_cell = centroid_cell
                    else:
                        count_cell = next(iter(h3_indices))
                else:
                    count_cell = next(iter(h3_indices))

                for h3_cell in h3_indices:
                    piece_area = row.area_km2 if has_polygons else 0.0
                    ratio = 1.0

                    # do the splitting if required
                    if needs_split:
                        cell_boundary = h3.cell_to_boundary(h3_cell)
                        h3_poly = Polygon([(lng, lat) for lat, lng in cell_boundary])
                        
                        try:
                            intersection = geom.intersection(h3_poly)
                        except Exception as e:
                            self.logger.warning(f"Intersection geometry error: {e}")
                            continue
                            
                        if intersection.is_empty:
                            continue
                            
                        # reproject the sliver to get it's true area
                        # have to do this because intersection is done on the non-equal area
                        # projections
                        intersection_proj = transform(intersection, project_to_area, interleaved=False)
                        piece_area = intersection_proj.area / 1e6
                        
                        # calculate what percentage of the original feature this sliver represents
                        if row.area_km2 > 0:
                            ratio = piece_area / row.area_km2
                        else:
                            ratio = 0.0

                    
                    count_val = 1 if h3_cell == count_cell else 0
                    rec = {"h3_index": h3_cell, "_count": count_val}

                    for col in attr_to_sum:
                        val = getattr(row, col)
                        # Prorate extensive variables (like population or crop yield) based on area split
                        rec[f"sum_{col}"] = val * ratio if feature_split else val
                        
                    for col in attr_to_mean:
                        # Intensive variables (like temperature) do not scale with area
                        rec[f"mean_{col}"] = getattr(row, col)

                    if has_polygons:
                        rec["area_km2"] = piece_area

                    records_by_res[res].append(rec)

        for res, records in records_by_res.items():
            if not records:
                self.logger.warning("No H3 coverage generated for resolution %s. Skipping.", res)
                continue

            df = pd.DataFrame(records)

            agg_dict = {"_count": "sum"}
            for col in attr_to_sum:
                agg_dict[f"sum_{col}"] = "sum"
            for col in attr_to_mean:
                agg_dict[f"mean_{col}"] = "sum"
            if "area_km2" in df.columns:
                agg_dict["area_km2"] = "sum"

            grouped = df.groupby("h3_index", as_index=False).agg(agg_dict)

            output_path = Path(output_paths[res])
            chunk_dir = output_path.parent / "chunks"
            chunk_dir.mkdir(parents=True, exist_ok=True)

            tile_id = f"{input_path.parent.parent.name}_{input_path.parent.name}_{input_path.stem}"
            chunk_filename = f"res_{res}_chunk_{tile_id}.parquet"
            chunk_path = chunk_dir / chunk_filename
            
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
        """Aggregates intermediate Parquet chunks into final GeoPackage datasets.

        Runs once after all chunks are generated. This method reads the intermediate 
        Parquet files for each resolution, aggregates the counts and attribute sums, 
        calculates final attribute means (by dividing the aggregated sum by the total 
        count), and generates spatial geometries for the H3 cells. Optionally calculates 
        land coverage metrics before writing the final output to GeoPackage files.

        Args:
            output_paths (dict): A dictionary mapping H3 resolutions (int) to their 
                target GeoPackage output file paths (PathLike).
            h3_res (List[int]): A list of H3 resolution levels to combine and output.
            land_polygons_path (Optional[PathLike], optional): File path to a vector 
                dataset of land boundaries used to calculate land coverage metrics. 
                Defaults to the instance's land_polygons_path.
            area_epsg (Optional[int], optional): The EPSG code for the equal-area 
                coordinate reference system used in area calculations. Defaults to 
                the instance's area_epsg.
            attr_to_sum (Optional[List[str]], optional): List of column names that 
                were aggregated by summation. Defaults to the instance's attr_to_sum.
            attr_to_mean (Optional[List[str]], optional): List of column names that 
                require a final mean calculation. Defaults to the instance's attr_to_mean.
        """
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
            agg_dict[f"mean_{col}"] = "sum"

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

            # finalize means
            for col in attr_to_mean:
                final_grouped[f"mean_{col}"] = final_grouped[f"mean_{col}"] / final_grouped["_count"]
            
            # generate geometries for h3 cells
            final_grouped["geometry"] = final_grouped["h3_index"].apply(self.h3_to_polygon)
            out_gdf = gpd.GeoDataFrame(final_grouped, geometry="geometry", crs="EPSG:4326")

            # calculate % cover
            if "area_km2" in out_gdf.columns and area_epsg is not None:
                # calculate cell area using the same projection as the features
                cell_area_km2 = out_gdf.to_crs(epsg=area_epsg).geometry.area / 1e6
                
                # calculate cover fraction
                out_gdf["percent_cover"] = (out_gdf["area_km2"] / cell_area_km2) * 100
                
                # Clip to 100 to handle microscopic floating-point projection artifacts
                #out_gdf["percent_cover"] = out_gdf["percent_cover"].clip(upper=100)

            if land_polygons_path is not None:
                out_gdf = self.add_land_metrics(out_gdf, land_polygons_path, area_epsg=area_epsg)
            
            out_gdf.to_file(output_path, driver="GPKG", mode="w")
            self.logger.info("Finished writing final GeoPackage for resolution %s.", res)


    
    def valid_h3_resolution(self, value: str) -> int:
        """Validates and parses an H3 resolution value.

        Args:
            value (str): The string representation of the H3 resolution to validate.

        Raises:
            argparse.ArgumentTypeError: If the value is not a valid integer or 
                if it falls outside the range of 1 to 15.

        Returns:
            int: The validated H3 resolution as an integer.
        """
        try:
            ivalue = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError("H3 resolution must be an integer.")
        if ivalue < 1 or ivalue > 15:
            raise argparse.ArgumentTypeError("H3 resolution must be between 1 and 15.")
        return ivalue

    
