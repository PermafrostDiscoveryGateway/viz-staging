import geopandas as gpd
import shapely
from pdgstaging.ConfigManager import ConfigManager


def test_validate_dedup_rules():
    tests = {
        "valid_rule1": {
            "data": gpd.GeoDataFrame({
                "col1": ["a", "b"],
                "geometry": [shapely.Point(1, 1), shapely.Point(1,2)]
            }),
            "rules": [["col1", "smaller"]]
        }
    }
    for test in tests.values():
        ConfigManager.validate_dedup_rules(test["data"], test["rules"])
