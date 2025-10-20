import os
import pytest
from pdgstaging.TilePathManager import TilePathManager

def test_validate_tms_and_path_structure_ok():
    tpm = TilePathManager(tms_id="WebMercatorQuad", path_structure=("style","tms","z","x","y"))
    assert tpm.tms_id == "WebMercatorQuad"

def test_path_roundtrip_from_dict_and_back(tmp_path):
    tpm = TilePathManager(tms_id="WebMercatorQuad", path_structure=("style","tms","z","x","y"))
    base = tmp_path.as_posix()
    d = {"style":"vec","tms":tpm.tms_id,"z":5,"x":10,"y":12,"base_dir":base,"ext":".gpkg"}
    p = tpm.path_from_dict(d)
    got = tpm.dict_from_path(p)
    for k in ("style","tms","z","x","y","ext"):
        assert got[k] == d[k]
    tile = tpm.tile_from_path(p)
    assert (tile.x, tile.y, tile.z) == (10,12,5)
