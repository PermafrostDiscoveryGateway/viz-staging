import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from pyproj import Transformer
from shapely.geometry import MultiPolygon, Polygon, box

from pdgstaging.antimeridian import (
    crosses_antimeridian,
    split_at_antimeridian,
    unwrap_antimeridian,
    wrap_longitude,
)

PROPERTIES = {
    "area": "area", "centroid_x": "centroid_x", "centroid_y": "centroid_y",
    "filename": "filename", "identifier": "identifier", "tile": "tile",
    "centroid_tile": "centroid_tile", "centroid_within_tile": "centroid_within_tile",
}
STAGE_KWARGS = dict(
    tms_id="WGS1984Quad", z=3, path_structure=("tms", "z", "x", "y"),
    properties=PROPERTIES,
)


def _wrapped_box(west=170.0, east=190.0, south=65.0, north=70.0):
    """A lon/lat box that straddles 180 degrees, written as a ring whose longitudes wrap."""
    ring = [(west, south), (180.0 - 1e-6, south), (-180.0 + 1e-6, south), (east - 360, south),
            (east - 360, north), (-180.0 + 1e-6, north), (180.0 - 1e-6, north), (west, north)]
    return Polygon(ring)


def _polar_source(path, west=170.0, east=190.0, south=65.0, north=70.0):
    """Write the box as a polygon in EPSG:3413 so reprojection to lon/lat makes the ring wrap."""
    t = Transformer.from_crs(4326, 3413, always_xy=True)
    lons = np.concatenate([np.linspace(west, east, 40), np.full(20, east),
                           np.linspace(east, west, 40), np.full(20, west)])
    lats = np.concatenate([np.full(40, south), np.linspace(south, north, 20),
                           np.full(40, north), np.linspace(north, south, 20)])
    x, y = t.transform(lons, lats)
    gpd.GeoDataFrame(geometry=[Polygon(zip(x, y))], crs=3413).to_file(path)
    return path


def test_crossing_ring_is_detected_and_unwrapped():
    wrapped = _wrapped_box()
    assert crosses_antimeridian(wrapped)

    unwrapped = unwrap_antimeridian(wrapped)

    assert unwrapped is not wrapped and unwrapped.is_valid
    assert unwrapped.bounds[0] == pytest.approx(170.0)
    assert unwrapped.bounds[2] == pytest.approx(190.0, abs=1e-3)
    assert unwrapped.area == pytest.approx(100.0, rel=1e-3)


def test_non_crossing_geometry_is_returned_unchanged():
    geom = box(-10, 60, 10, 70)
    assert not crosses_antimeridian(geom)
    assert unwrap_antimeridian(geom) is geom
    assert split_at_antimeridian(geom) is geom


def test_split_folds_east_part_into_range_and_keeps_area():
    unwrapped = unwrap_antimeridian(_wrapped_box())

    split = split_at_antimeridian(unwrapped)

    assert isinstance(split, MultiPolygon) and len(split.geoms) == 2
    assert split.bounds[0] >= -180.0 and split.bounds[2] <= 180.0
    assert split.area == pytest.approx(unwrapped.area)
    assert split.is_valid


def test_wrap_longitude():
    assert wrap_longitude([185.0, -185.0, 180.0, -180.0, 10.0]).tolist() == [
        -175.0, 175.0, -180.0, -180.0, 10.0,
    ]


def test_stage_source_keeps_polygon_across_antimeridian_whole(tmp_path):
    from pdgstaging.operations import stage_source

    source = _polar_source(tmp_path / "source.gpkg")
    shards = stage_source(source, tmp_path / "shards", "one", **STAGE_KWARGS)

    frames = [gpd.read_file(s) for s in shards]
    staged = gpd.GeoDataFrame(pd.concat(frames), crs=frames[0].crs)
    # Every stored geometry stays inside the TMS extent and nothing spans the globe.
    assert staged.total_bounds[0] >= -180.0 and staged.total_bounds[2] <= 180.0
    parts = [p for g in staged.geometry for p in g.geoms]
    assert all(p.bounds[2] - p.bounds[0] < 180 for p in parts)
    assert staged.geometry.is_valid.all()
    # Area and centroid describe the real ~20 x 5 degree polygon, not a 360 degree smear.
    owner = staged[staged[PROPERTIES["centroid_within_tile"]]]
    assert len(owner) == 1
    assert owner[PROPERTIES["area"]].iloc[0] == pytest.approx(100.0, rel=0.1)
    assert abs(owner[PROPERTIES["centroid_x"]].iloc[0]) > 170.0
    assert 65.0 < owner[PROPERTIES["centroid_y"]].iloc[0] < 70.0
    # No temporary helper columns leak into the staged tiles.
    assert not any(c.startswith("_antimeridian") for c in staged.columns)
