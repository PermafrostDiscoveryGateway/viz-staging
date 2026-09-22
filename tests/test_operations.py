from pathlib import Path

import geopandas as gpd
import pytest
from shapely.geometry import Point, box

from pdgstaging.TileStager import TileStager

PROPERTIES = {
    "area": "area", "centroid_x": "centroid_x", "centroid_y": "centroid_y",
    "filename": "filename", "identifier": "identifier", "tile": "tile",
    "centroid_tile": "centroid_tile", "centroid_within_tile": "centroid_within_tile",
}

STAGE_KWARGS = dict(
    tms_id="WGS1984Quad", z=0, path_structure=("tms", "z", "x", "y"),
    properties=PROPERTIES,
)


def _write_source(path, geometries, crs="EPSG:4326"):
    gpd.GeoDataFrame(geometry=list(geometries), crs=crs).to_file(path)
    return path


def test_tile_stager_defaults_are_not_shared():
    first = TileStager()
    second = TileStager()

    assert first.props is not second.props


def test_stage_source_requires_an_existing_input(tmp_path):
    from pdgstaging.operations import stage_source

    with pytest.raises(FileNotFoundError):
        stage_source(
            tmp_path / "missing", tmp_path / "shards", "source", **STAGE_KWARGS
        )


def test_source_key_must_be_a_safe_component(tmp_path):
    from pdgstaging.operations import stage_source

    source = _write_source(tmp_path / "source.gpkg", [box(-1, -1, 1, 1)])
    for bad in ("", ".", "..", "a/b", "/abs"):
        with pytest.raises(ValueError):
            stage_source(source, tmp_path / "shards", bad, **STAGE_KWARGS)


def test_overlapping_sources_stay_separate(tmp_path):
    from pdgstaging.operations import stage_source

    source = _write_source(tmp_path / "source.gpkg", [box(-1, -1, 1, 1)])
    root = tmp_path / "shards"
    first = stage_source(source, root, "one", **STAGE_KWARGS)
    second = stage_source(source, root, "two", **STAGE_KWARGS)

    assert first and second and len(first) == len(second)
    assert set(first).isdisjoint(second)
    assert all(path.is_relative_to(root / "one") for path in first)
    assert all(path.is_relative_to(root / "two") for path in second)


def test_stage_source_empty_input_returns_empty_list(tmp_path):
    from pdgstaging.operations import stage_source

    source = _write_source(tmp_path / "source.gpkg", [Point(0, 0)])
    assert stage_source(source, tmp_path / "shards", "one", **STAGE_KWARGS) == []


def test_stage_source_overwrite_reuses_and_regenerates(tmp_path):
    from pdgstaging.operations import stage_source

    source = _write_source(tmp_path / "source.gpkg", [box(-1, -1, 1, 1)])
    first = stage_source(source, tmp_path / "shards", "one", **STAGE_KWARGS)
    assert stage_source(source, tmp_path / "shards", "one", **STAGE_KWARGS) == first

    first[0].write_text("not a geopackage")
    with pytest.raises((FileExistsError, ValueError)):
        stage_source(source, tmp_path / "shards", "one", **STAGE_KWARGS)
    regenerated = stage_source(
        source, tmp_path / "shards", "one", overwrite=True, **STAGE_KWARGS
    )
    assert regenerated == first and len(gpd.read_file(regenerated[0])) == 1


def test_merge_uses_only_explicit_shards(tmp_path):
    from pdgstaging.operations import merge_staged_tile, stage_source

    source = _write_source(tmp_path / "source.gpkg", [box(0, 0, 1, 1)])
    root = tmp_path / "shards"
    first = stage_source(source, root, "one", **STAGE_KWARGS)
    second = stage_source(source, root, "two", **STAGE_KWARGS)
    third = stage_source(source, root, "three", **STAGE_KWARGS)

    merged = merge_staged_tile(
        [second[0], first[0]], tmp_path / "tile.gpkg",
    )
    assert merged == tmp_path / "tile.gpkg"
    assert len(gpd.read_file(merged)) == 2
    assert third[0].is_file() and len(gpd.read_file(merged)) == 2


def test_merge_rejects_missing_and_malformed_shards(tmp_path):
    from pdgstaging.operations import merge_staged_tile, stage_source

    source = _write_source(tmp_path / "source.gpkg", [box(0, 0, 1, 1)])
    (shard,) = stage_source(source, tmp_path / "shards", "one", **STAGE_KWARGS)
    malformed = tmp_path / "bad.gpkg"
    malformed.write_text("not a geopackage")

    with pytest.raises(FileNotFoundError):
        merge_staged_tile([tmp_path / "absent.gpkg"], tmp_path / "out.gpkg")
    with pytest.raises(ValueError):
        merge_staged_tile([malformed], tmp_path / "out.gpkg")
    with pytest.raises(ValueError):
        merge_staged_tile([shard, malformed], tmp_path / "out.gpkg")
    with pytest.raises(ValueError):
        merge_staged_tile([], tmp_path / "out.gpkg")


def test_merge_rejects_mismatched_shards(tmp_path):
    from pdgstaging.operations import merge_staged_tile, stage_source

    source = _write_source(tmp_path / "source.gpkg", [box(0, 0, 1, 1)])
    (shard,) = stage_source(source, tmp_path / "shards", "one", **STAGE_KWARGS)
    frame = gpd.read_file(shard)

    other_crs = tmp_path / "other-crs.gpkg"
    frame.set_crs(frame.crs, allow_override=True).to_crs("EPSG:3857").to_file(other_crs)
    with pytest.raises(ValueError, match="[Cc][Rr][Ss]"):
        merge_staged_tile([shard, other_crs], tmp_path / "crs.gpkg")

    other_schema = tmp_path / "other-schema.gpkg"
    frame.drop(columns=[PROPERTIES["identifier"]]).to_file(other_schema)
    with pytest.raises(ValueError, match="[Ss]chema"):
        merge_staged_tile([shard, other_schema], tmp_path / "schema.gpkg")

    other_tile = tmp_path / "other-tile.gpkg"
    relabeled = frame.copy()
    for column in relabeled.columns:
        if column != "geometry" and (
            relabeled[column].astype(str).str.startswith("Tile(").all()
        ):
            relabeled[column] = "Tile(x=999, y=999, z=0)"
    relabeled.to_file(other_tile)
    with pytest.raises(ValueError, match="[Tt]ile"):
        merge_staged_tile([shard, other_tile], tmp_path / "tile.gpkg")


def test_merge_overwrite_reuses_valid_output(tmp_path):
    from pdgstaging.operations import merge_staged_tile, stage_source

    source = _write_source(tmp_path / "source.gpkg", [box(0, 0, 1, 1)])
    (shard,) = stage_source(source, tmp_path / "shards", "one", **STAGE_KWARGS)
    output = tmp_path / "tile.gpkg"

    assert merge_staged_tile([shard], output) == output
    assert merge_staged_tile([shard], output) == output

    output.write_text("not a geopackage")
    with pytest.raises((FileExistsError, ValueError)):
        merge_staged_tile([shard], output)
    assert merge_staged_tile([shard], output, overwrite=True) == output
    assert len(gpd.read_file(output)) == 1


def test_merge_overwrite_failed_write_leaves_no_output(tmp_path, monkeypatch):
    from pdgstaging.operations import merge_staged_tile, stage_source

    source = _write_source(tmp_path / "source.gpkg", [box(0, 0, 1, 1)])
    (shard,) = stage_source(source, tmp_path / "shards", "one", **STAGE_KWARGS)
    output = tmp_path / "tile.gpkg"

    def _fail(*args, **kwargs):
        raise RuntimeError("write failed")

    monkeypatch.setattr(gpd.GeoDataFrame, "to_file", _fail)
    with pytest.raises(RuntimeError):
        merge_staged_tile([shard], output, overwrite=True)

    assert not output.exists()
    assert list(tmp_path.glob("*.tmp*")) == []
