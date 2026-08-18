import pytest
from pathlib import Path
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
from pdgstaging import H3GridSummaryGenerator 


@pytest.fixture
def sample_staging_data(tmp_path: Path):
    """
    Sets up test fixtures by creating two overlapping GeoPackages 
    in a temporary directory. Yields the list of file paths.
    """
    # File 1: 2 overlapping features
    p1 = Polygon([(-150, 65), (-150, 66), (-149, 66), (-149, 65)])
    p2 = Polygon([(-149.5, 65.5), (-149.5, 66.5), (-148.5, 66.5), (-148.5, 65.5)])
    
    df1 = pd.DataFrame({
        "carbon_stock": [100, 200],  
        "temperature": [-5.0, -4.0], 
        "staging_centroid_within_tile": [True, True]
    })
    gdf1 = gpd.GeoDataFrame(df1, geometry=[p1, p2], crs="EPSG:4326")
    file1 = tmp_path / "staging_part1.gpkg"
    gdf1.to_file(file1, driver="GPKG")

    # File 2: 2 more features right in the same area
    p3 = Polygon([(-150.2, 65.2), (-150.2, 66.2), (-149.2, 66.2), (-149.2, 65.2)])
    p4 = Polygon([(-149.8, 65.8), (-149.8, 66.8), (-148.8, 66.8), (-148.8, 65.8)])
    
    df2 = pd.DataFrame({
        "carbon_stock": [300, 400],
        "temperature": [-6.0, -7.0],
        "staging_centroid_within_tile": [True, True]
    })
    gdf2 = gpd.GeoDataFrame(df2, geometry=[p3, p4], crs="EPSG:4326")
    file2 = tmp_path / "staging_part2.gpkg"
    gdf2.to_file(file2, driver="GPKG")

    return [file1, file2]


@pytest.fixture
def h3_stager():
    """
    Initializes and returns the H3 generator instance.
    """
    return H3GridSummaryGenerator(
        config=None,
        attr_to_sum=["carbon_stock"],
        attr_to_mean=["temperature"]
    )


def test_h3_feature_count_single_resolution(tmp_path: Path, sample_staging_data: list, h3_stager: H3GridSummaryGenerator):
    """
    Tests that the final aggregated H3 GeoPackage correctly accounts 
    for all features across multiple input files for a single zoom level.
    """

    h3_res_list = [3]
    out_paths = {3: tmp_path / "final_h3_res3.gpkg"}

    for file_path in sample_staging_data:
        h3_stager.build_h3_summary(
            input_path=file_path, 
            output_paths=out_paths, 
            h3_res=h3_res_list
        )

    h3_stager.combine_h3_summaries(
        output_paths=out_paths, 
        h3_res=h3_res_list
    )

    final_file = out_paths[3]
    assert final_file.exists(), "Final GeoPackage was not created by the combiner."

    final_gdf = gpd.read_file(final_file)
    
    total_features_counted = final_gdf["_count"].sum()
    
    assert total_features_counted == 4, f"Expected exactly 4 total features counted, but got {total_features_counted}"

def test_h3_feature_count_multi_resolution(tmp_path: Path, sample_staging_data: list, h3_stager: H3GridSummaryGenerator):
    """
    Tests that the final aggregated H3 GeoPackage correctly accounts 
    for all features across multiple input files for multiple zoom levels.
    """

    h3_res_list = [3,4]
    out_paths = {3: tmp_path / "final_h3_res3.gpkg",
                 4: tmp_path / "final_h3_res4.gpkg"}

    for file_path in sample_staging_data:
        h3_stager.build_h3_summary(
            input_path=file_path, 
            output_paths=out_paths, 
            h3_res=h3_res_list
        )

    h3_stager.combine_h3_summaries(
        output_paths=out_paths, 
        h3_res=h3_res_list
    )

    for res, final_file in out_paths.items():
        assert final_file.exists(), f"Final GeoPackage for resolution {res} was not created."

        final_gdf = gpd.read_file(final_file)
        assert not final_gdf.empty, f"GeoDataFrame for resolution {res} is empty."

        total_features_counted = final_gdf["_count"].sum()

        # At higher resolutions (like res 4), polygons can span multiple H3 cells,
        # so total `_count` across all cells will be >= 4.
        assert total_features_counted >= 4, (
            f"Expected 4 total feature counts for resolution {res}, "
            f"but got {total_features_counted}"
        )

def test_h3_area_calculation_multi_resolution(
    tmp_path: Path, sample_staging_data: list, h3_stager: H3GridSummaryGenerator
):
    """
    Tests that the total intersected area (area_km2) is calculated accurately
    and conserved across multiple H3 resolutions using the fixture data.
    """

    expected_area_km2 = 0.0
    
    for file_path in sample_staging_data:
        gdf = gpd.read_file(file_path)
        gdf_ea = gdf.to_crs(epsg=6933)
        expected_area_km2 += (gdf_ea.geometry.area / 1e6).sum()

    h3_res_list = [3, 4]
    out_paths = {
        3: tmp_path / "final_h3_area_res3.gpkg",
        4: tmp_path / "final_h3_area_res4.gpkg",
    }

    for file_path in sample_staging_data:
        h3_stager.build_h3_summary(
            input_path=file_path, 
            output_paths=out_paths, 
            h3_res=h3_res_list
        )

    h3_stager.combine_h3_summaries(
        output_paths=out_paths,
        h3_res=h3_res_list
    )

    for res, final_file in out_paths.items():
        assert final_file.exists(), f"Final GeoPackage for resolution {res} was not created."

        final_gdf = gpd.read_file(final_file)
        assert "area_km2" in final_gdf.columns, f"Res {res} output is missing 'area_km2' column."

        total_calculated_area_km2 = final_gdf["area_km2"].sum()

        assert total_calculated_area_km2 == pytest.approx(expected_area_km2, rel=1e-3), (
            f"Resolution {res} calculated total area {total_calculated_area_km2} km², "
            f"expected {expected_area_km2} km²"
        )