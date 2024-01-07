""" Configuration module containing settings and constants. """

import os
import sys

#wrfout = '/home/c707201/temp/WRF_output_project.nc'
wrfout = r'C:\Users\andre\Downloads\WRF_output_project.nc'
output_directory = r'C:\Users\mhde\2_Wissprog\Term_Project\Plots_and_HTML'

if os.path.isfile(wrfout):
    # location of data directory
    pkgdir = os.path.dirname(__file__)
    html_template = os.path.join(pkgdir, 'data', 'template.html')
    html_skewT_template = os.path.join(pkgdir, 'data', 'skewT_template.html')  
    test_ts_df = os.path.join(pkgdir, 'data', 'test_df_timeseries.pkl')
    test_hgt = os.path.join(pkgdir, 'data', 'test_hgt.nc')

    # minimum and maximum elevations for topography plot
    topo_min = 0
    topo_max = 3200
else:
    print('The specified WRF output file does not exist. Please set a valid path in cfg.py.')
    sys.exit(1)
