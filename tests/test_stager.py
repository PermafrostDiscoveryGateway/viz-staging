import geopandas as gpd
import geopandas.testing as gpd_test
from shapely import Polygon
from pdgstaging.TileStager import clip_to_footprint
# from pdgstaging.TileStager import df_intersects_antimeridian, clip_to_footprint

# def test_intersects_antimeridian():
#     tests = {
#         "test1": {
#             "gdf": gpd.GeoDataFrame({
#                 "geometry": [
#                     Polygon([(178, 80), (-178, 80), (-178, 65), (178, 65)])
#                 ]
#             }, crs="EPSG:4326"),
#             "want": [True],
#         }
#     }
#     for test in tests.values(): 
#         res = df_intersects_antimeridian(test["gdf"])
#         print(res)
#         for idx, val in enumerate(res): 
#             assert val == test["want"][idx] 

def test_clip_to_footprint():        
    tests = {
        "test1": {
            "gdf": gpd.GeoDataFrame({
                "geometry": [
                    # rough bounds of the city of Albuquerque
                    Polygon([(-106.667597, 35.217839), (-106.482254, 35.198681), (-106.485891, 35.042242), (-106.746253, 35.043637)]),
                    # a random triangle in West Texas
                    Polygon([(-104.137710, 31.140494), (-104.076552, 31.042003), (-104.203313, 31.028436)]),

                ],

            }, crs="EPSG:4326"),
            "fp": gpd.GeoDataFrame({
                # very rough bounds of New Mexico
                "geometry": [Polygon([(-109.030037, 36.993093), (-103.060255, 37.009544), (-103.084202, 32.067447), (-109.036903, 31.334431)])]
            }, crs="EPSG:4326"),
            "want": gpd.GeoDataFrame({
                "geometry": [
                    # a random triangle in West Texas
                    Polygon([(-104.137710, 31.140494), (-104.076552, 31.042003), (-104.203313, 31.028436)]),
                    # rough bounds of the city of Albuquerque
                    Polygon([(-106.667597, 35.217839), (-106.482254, 35.198681), (-106.485891, 35.042242), (-106.746253, 35.043637)]),

                ],
                "dup_label": [True, False],

            }, crs="EPSG:4326")
        }

    }
    for test in tests.values(): 
        test["gdf"].plot()
        clipped = clip_to_footprint(test["gdf"], test["fp"], "dup_label")
        gpd_test.assert_geodataframe_equal(clipped, test["want"])