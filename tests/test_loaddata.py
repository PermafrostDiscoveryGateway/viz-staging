from pandas import DataFrame
from pdgstaging import ConfigManager

def test_init():
    """Initialize tests and show they are working.
    """
    assert 1 == 1

def test_load_data():
    """Load example testing data for staging tests.
    """
    df = DataFrame( dict( x=[1,2,3], y=[4,5,6] ) )
    assert df is not None
