import geopandas as gpd
from shapely.geometry import box

from pdgstaging.Deduplicator import (
    keep_rules_to_sort_order,
    clip_gdf,
    deduplicate_neighbors,
    label_duplicates,
)

def _gdf(polys, crs="EPSG:3857", **props):
    # Helper: build a GeoDataFrame with given geometries and arbitrary columns.
    gdf = gpd.GeoDataFrame(props, geometry=polys, crs=crs)
    return gdf

def test_keep_rules_to_sort_order_basic():
    # ("area","smaller")  -> ascending (True)
    # ("date","larger")   -> descending (False)
    rules = [("area", "smaller"), ("date", "larger")]
    props, order = keep_rules_to_sort_order(rules)
    assert props == ["area", "date"]
    assert order == [True, False]

def test_clip_gdf_intersects_labels_inside_outside():
    # clip_gdf(..., method="intersects") should mark features outside boundary
    # as duplicates=True and inside as False.

    # boundary: 0..10 square
    boundary = _gdf([box(0,0,10,10)])
    inside  = box(1,1,2,2)
    outside = box(20,20,21,21)
    # features: 'a' inside, 'b' outside
    gdf = _gdf([inside, outside], prop=["a","b"])
    out = clip_gdf(gdf=gdf, boundary=boundary, method="intersects", prop_duplicated="dup")
    # Expect: only 'b' is flagged dup=True
    assert set(out.columns) >= {"prop", "geometry", "dup"}
    assert out["dup"].sum() == 1
    assert out.loc[out["prop"]=="b","dup"].item() is True
    assert out.loc[out["prop"]=="a","dup"].item() is False

def test_label_duplicates_combines_and_marks():
    gdf_keep = _gdf([box(0,0,1,1)], a=[1])
    gdf_drop = _gdf([box(1,1,2,2)], a=[2])
    combined = label_duplicates({"to_keep": gdf_keep.copy(), "to_remove": gdf_drop.copy()}, "dup")
    assert set(combined["dup"].unique()) == {False, True}
    assert len(combined) == 2

def test_deduplicate_neighbors_overlap_by_area_keeps_larger():
    g1 = box(0,0,2,2)        
    g2 = box(1,1,2.5,2.5)     
    gdf = _gdf([g1,g2], source=["A","B"])
    out = deduplicate_neighbors(
        gdf,
        split_by="source",
        prop_area=None,
        prop_centroid_x=None,
        prop_centroid_y=None,
        keep_rules=[("area_", "larger")],  
        overlap_tolerance=0.2,             # they overlap > 20%
        overlap_both=False,
        centroid_tolerance=None,
        return_intersections=False,
        prop_duplicated="dup",
    )
    assert "dup" in out.columns
    # smaller polygon should be flagged True
    dup_sources = out.loc[out["dup"], "source"].tolist()
    assert dup_sources in (["B"], ["B"])  # B should be removed
