import json
import subprocess
import sys

import geopandas as gpd
from shapely.geometry import box
from test_operations import PROPERTIES


def _run(*args, cwd):
    return subprocess.run(
        [sys.executable, "-m", "pdgstaging", *args],
        capture_output=True, text=True, cwd=cwd,
    )


def test_cli_stages_a_source(tmp_path):
    source = tmp_path / "source.gpkg"
    gpd.GeoDataFrame(geometry=[box(0, 0, 1, 1)], crs="EPSG:4326").to_file(source)

    staged = _run(
        "stage-source", str(source), str(tmp_path / "shards"), "one",
        "--z", "0", "--path-structure", '["tms", "z", "x", "y"]',
        "--properties", json.dumps(PROPERTIES), cwd=tmp_path,
    )
    assert staged.returncode == 0, staged.stderr
    shards = json.loads(staged.stdout)
    assert len(shards) == 1 and shards[0].endswith(".gpkg")


def test_cli_failure_exits_nonzero_with_a_message(tmp_path):
    result = _run(
        "stage-source", "missing.gpkg", str(tmp_path / "shards"), "one",
        "--z", "0", "--properties", json.dumps(PROPERTIES), cwd=tmp_path,
    )
    assert result.returncode == 2
    assert "missing.gpkg" in result.stderr


def test_cli_merges_explicit_shards(tmp_path):
    source = tmp_path / "source.gpkg"
    gpd.GeoDataFrame(geometry=[box(0, 0, 1, 1)], crs="EPSG:4326").to_file(source)
    staged = _run(
        "stage-source", str(source), str(tmp_path / "shards"), "one",
        "--z", "0", "--path-structure", '["tms", "z", "x", "y"]',
        "--properties", json.dumps(PROPERTIES), cwd=tmp_path,
    )
    output = tmp_path / "tile.gpkg"
    merged = _run("merge-tile", str(output), *json.loads(staged.stdout), cwd=tmp_path)
    assert merged.returncode == 0, merged.stderr
    assert json.loads(merged.stdout) == str(output)
    assert len(gpd.read_file(output)) == 1
    missing = _run("merge-tile", str(tmp_path / "x.gpkg"), "nope.gpkg", cwd=tmp_path)
    assert missing.returncode == 2 and "missing shard" in missing.stderr
