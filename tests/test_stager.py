import geopandas as gpd
import geopandas.testing as gpd_test
from shapely import Polygon, MultiPolygon, Point, GeometryCollection
from shapely.testing import assert_geometries_equal
from pdgstaging.TileStager import * 

def test_clip_to_footprint():        
    tests = {
        "clip one polygon": {
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
        },
        "clip to multipolygon footprint": {
            "gdf": gpd.GeoDataFrame({
                "geometry": [
                    # rough bounds of the city of Albuquerque
                    Polygon([(-106.667597, 35.217839), (-106.482254, 35.198681), (-106.485891, 35.042242), (-106.746253, 35.043637)]),
                    # a random triangle in West Texas
                    Polygon([(-104.137710, 31.140494), (-104.076552, 31.042003), (-104.203313, 31.028436)]),
                    # Garden City, KS
                    Polygon([(-100.936106, 38.001829), (-100.820050, 38.002620), (-100.820404, 37.944815), (-100.941119, 37.932927)])

                ],

            }, crs="EPSG:4326"),
            "fp": gpd.GeoDataFrame({
                "geometry": [MultiPolygon([
                    # very rough bounds of New Mexico
                    Polygon([(-109.030037, 36.993093), (-103.060255, 37.009544), (-103.084202, 32.067447), (-109.036903, 31.334431)]),
                    # very rough bounds of Kansas
                    Polygon([(-101.994710, 39.965017), (-95.011623, 39.829632), (-94.666142, 37.041031), (-101.991577, 37.017151)])])]
            }, crs="EPSG:4326"),
            "want": gpd.GeoDataFrame({
                "geometry": [
                    # a random triangle in West Texas, should be marked duplicate
                    Polygon([(-104.137710, 31.140494), (-104.076552, 31.042003), (-104.203313, 31.028436)]),
                    # rough bounds of the city of Albuquerque, should not be marked duplicate
                    Polygon([(-106.667597, 35.217839), (-106.482254, 35.198681), (-106.485891, 35.042242), (-106.746253, 35.043637)]),
                    # Garden City, KS, should not be marked duplicate
                    Polygon([(-100.936106, 38.001829), (-100.820050, 38.002620), (-100.820404, 37.944815), (-100.941119, 37.932927)])

                ],
                "dup_label": [True, False, False],

            }, crs="EPSG:4326")
        }

    }
    for name, test in tests.items(): 
        print(f'running test: {name}')
        test["gdf"].plot()
        clipped = clip_to_footprint(test["gdf"], test["fp"], "dup_label")
        gpd_test.assert_geodataframe_equal(clipped, test["want"])

def test_label_am_crossings(): 
    gdf = gpd.GeoDataFrame(
        {"name": ["cross1", "nocross1"],
            "geometry": [
            Polygon([(178, 80), (-178, 80), (-178, 65), (178, 65)]),
            Polygon([(170, 80), (178, 80), (178, 65), (170, 65)]),
            ]}, crs="EPSG:4326",
    )
    label = "crosses_am"
    want = gpd.GeoDataFrame(
        {"name": ["cross1", "nocross1"],
            "geometry": [
            Polygon([(178, 80), (-178, 80), (-178, 65), (178, 65)]),
            Polygon([(170, 80), (178, 80), (178, 65), (170, 65)]),
            ],
            "crosses_am": [True, False]},
            crs="EPSG:4326",
    )

    res = label_am_crossings(gdf, label)
    gpd_test.assert_geodataframe_equal(res, want)

def test_split_am_polygon():
    tests = {
        "test1": {
            "polygon": Polygon([(178, 80), (-178, 80), (-178, 65), (178, 65)]),
            "want_polygons": GeometryCollection([
                Polygon([(178, 80), (180, 80), (180, 65), (178, 65)]),
                Polygon([(180, 80), (182, 80), (182, 65), (180, 65)]),
            ])
        }, 
        "test2": {
            "polygon": Polygon([(-178, 80), (178, 80), (178, 65), (-178, 65)]),
            "want_polygons": GeometryCollection([
                Polygon([(-180, 80), (-178, 80), (-178, 65), (-180, 65)]),
                Polygon([(-180, 65), (-182, 65),(-182, 80),(-180, 80) ]),
            ])
        }
    }
    for _, test in tests.items(): 
        got = split_am_polygon(test["polygon"])
        print(got)
        assert_geometries_equal(got, test["want_polygons"]) 

def test_split_am_crossing_polygons(): 
    gdf = gpd.GeoDataFrame(
        {"name": ["cross1", "nocross1", "cross2"],
            "geometry": [
            Polygon([(178, 80), (-178, 80), (-178, 65), (178, 65)]),
            Polygon([(170, 80), (178, 80), (178, 65), (170, 65)]),
Polygon([(-178, 80), (178, 80), (178, 65), (-178, 65)]),
            ]}, crs="EPSG:4326",
    )
    label = "crosses_am"
    want = gpd.GeoDataFrame(
        {"name": ["nocross1", "cross1", "cross1", "cross2", "cross2"],
            "geometry": [
                Polygon([(170, 80), (178, 80), (178, 65), (170, 65)]),
                Polygon([(178, 80), (180, 80), (180, 65), (178, 65)]),
                Polygon([(-180, 80), (-178, 80), (-178, 65), (-180, 65)]),
                Polygon([(-180, 80), (-178, 80), (-178, 65), (-180, 65)]),
                Polygon([(180, 65), (178, 65),(178, 80),(180, 80) ]),
            ]},
            crs="EPSG:4326",
    )

    res = split_am_crossing_polygons(gdf, label)
    print(res.difference(want))
    gpd_test.assert_geodataframe_equal(res, want)