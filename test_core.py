import os
import numpy as np

from wrfvis import core, cfg


def test_get_ts():

    df, hgt, col_names = core.get_wrf_timeseries('T', 11, 45, 300, rad = 10000)

    assert df.attrs['variable_name'] == 'T'
    assert df.attrs['variable_units'] == 'K'
    assert df.attrs['grid_point_elevation_time0'] < 400
    assert df.attrs['grid_point_elevation_time0'] > 10
    np.testing.assert_allclose(df['T'].attrs['distance_to_grid_point'].XLAT, df['T'].attrs['lat_grid_point'])
    np.testing.assert_allclose(df['T'].attrs['distance_to_grid_point'].XLONG, df['T'].attrs['lon_grid_point'])

    # dimensions of hgt
    assert hgt.dims == ('south_north', 'west_east')
    
    # dimesions of df
    assert df.shape == (36,len(col_names))


def test_mkdir(tmpdir):

    dir = str(tmpdir.join('html_dir'))
    core.mkdir(dir)
    assert os.path.isdir(dir)
