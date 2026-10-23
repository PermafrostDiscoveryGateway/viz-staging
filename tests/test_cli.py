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
