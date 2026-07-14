import pytest

@pytest.fixture(autouse=True)
def _silence_geopandas_warnings(recwarn):
    pass
