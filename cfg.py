"""
Configuration module containing settings and constants for a WRF (Weather Research and Forecasting) output analysis script.

Constants
---------
- `wrfout`: Path to the WRF output file.
- `output_directory`: Directory where generated plots and HTML files will be saved.

Files
-----
- `html_template`: Path to the general HTML template file.
- `html_skewT_template`: Path to the skewT HTML template file.
- `html_weather_template`: Path to the weather table HTML template file.
- `html_geographic_template`: Path to the geographic table HTML template file.
- `html_twotable_template`: Path to the two-table HTML template file.
- `html_combined_template`: Path to the combined HTML template file.
- `test_ts_df`: Path to the test timeseries DataFrame file.
- `test_hgt`: Path to the test topography file.

Parameters
----------
- `topo_min`: Minimum elevation for the topography plot.
- `topo_max`: Maximum elevation for the topography plot.

Note
----
If the specified WRF output file does not exist, an error message will be displayed, and the program will exit.
"""

import os
import sys

wrfout = r'C:\Users\andre\Downloads\WRF_output_project.nc'
output_directory = r'C:\Users\andre\OneDrive\Dokumente\UNI\Master_Semester_1\ScientificProgramming\FinalProject\FinalVersion\wrvis_main\wrfvis\data\Plots_and_HTML'

if os.path.isfile(wrfout):
    
    # location of data directory
    pkgdir = os.path.dirname(__file__)
    html_template = os.path.join(pkgdir, 'data', 'template.html') 
    html_skewT_template = os.path.join(pkgdir, 'data', 'skewT_template.html')
    html_weather_template = os.path.join(pkgdir, 'data', 'weather_table_template.html')
    html_geographic_template = os.path.join(pkgdir, 'data', 'geogr_table_template.html')
    html_twotable_template = os.path.join(pkgdir, 'data', 'twotable_template.html')
    html_combined_template = os.path.join(pkgdir, 'data', 'combined_template.html')
    test_ts_df = os.path.join(pkgdir, 'data', 'test_df_timeseries.pkl')
    test_hgt = os.path.join(pkgdir, 'data', 'test_hgt.nc')

    # minimum and maximum elevations for topography plot
    topo_min = 0
    topo_max = 3200
else:
    print('The specified WRF output file does not exist.'
          ' Please set a valid path in cfg.py.')
    sys.exit(1)
