import numpy as np
import xarray as xr
import pytest

from wrfvis import cfg, grid

@pytest.fixture
def sample_dataset():
    lon = np.linspace(-180, 180, 10)
    lat = np.linspace(-90, 90, 5)
    XLONG, XLAT = np.meshgrid(lon, lat)
    
    # Creating a more standard WRF-like dataset with coordinates
    ds = xr.Dataset({
        'XLONG': (['south_north', 'west_east'], XLONG),
        'XLAT': (['south_north', 'west_east'], XLAT)
    })

    # Add coordinates to the dataset
    ds = ds.assign_coords({
        'south_north': np.arange(5),
        'west_east': np.arange(10)
    })

    return ds

# sample_dataset = sample_dataset()

def test_haversine():
    c = grid.haversine(34, 42, 35, 42)
    np.testing.assert_allclose(c, 82633.46475287154)

    c = grid.haversine(34, 42, [35, 36], [42, 42])
    np.testing.assert_allclose(c, np.array([82633.46475287, 165264.11172113]))


def test_find_nearest_gridcell():

    # test dataset
    hgt = xr.open_dataarray(cfg.test_hgt)

    ind, dist = grid.find_nearest_gridcell(hgt.XLONG, hgt.XLAT, 11, 45)
    assert type(ind) == tuple
    assert ind == (229, 303)
    np.testing.assert_allclose(dist, 3024.5848211250755)
    

def test_find_grid_cells_in_radius(sample_dataset):
    
    # Assuming you have a function find_nearest_gridcell to get the nearest grid cell
    ngcind, ngcdist = grid.find_nearest_gridcell(sample_dataset.XLONG, sample_dataset.XLAT, 0, 0)

    # Test with a small radius
    with pytest.raises(ValueError, match="Radius has to be larger than distance between nearest grid cell and target location"):
        grid.find_grid_cells_in_radius(ngcind, ngcdist, rad=ngcdist - 1,
                                       lon=0, lat=0, ds=sample_dataset)
