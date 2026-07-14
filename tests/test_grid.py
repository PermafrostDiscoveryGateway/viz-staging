import geopandas as gpd
from shapely.geometry import Point, box
from pdgstaging.Grid import Grid, TMSGrid

def test_grid_sjoin_points_fastpath():
    g = Grid(bounds=[0,0,10,10], crs=3857, nrows=2, ncols=2)
    pts = gpd.GeoDataFrame(geometry=[Point(1,1), Point(9,9)], crs=3857)
    out = g.sjoin(pts, how="left", predicate="intersects", as_index=False)
    assert {"grid_row_index","grid_column_index"}.issubset(out.columns)
    assert len(out) == 2

def test_tmsgrid_tiles_from_xy_roundtrip():
    tg = TMSGrid("WebMercatorQuad", 4, bounds=[-1e6,-1e6,1e6,1e6])
    # pick a point within bounds
    pts = gpd.GeoDataFrame(geometry=[Point(0,0)], crs=tg.crs)
    out = tg.sjoin(pts, as_tile=True)
    assert tg.TILE_NAME in out.columns
    tile = out[tg.TILE_NAME].iloc[0]
    # tile object has x,y,z
    assert tile.z == 4



