"""
Configuration module containing settings and constants.

Constants
---------
- `wrfout`: Path to the WRF output file.
- `output_directory`: Directory where generated plots and HTML files will be saved.

Files
-----
- `html_template`: Path to the HTML template file.
- `html_topo_template`: Path to the topo HTML template file.
- `html_timeseries_template`: Path to the timeseries HTML template file.
- `html_skewT_template`: Path to the skewT HTML template file.
- `html_snowcheck_template`: Path to the snowcheck HTML template file.
- `test_ts_df`: Path to the test timeseries DataFrame file.
- `test_hgt`: Path to the test topography file.

Parameters
----------
- `topo_min`: Minimum elevation for topography plot.
- `topo_max`: Maximum elevation for topography plot.

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
    html_topo_template = os.path.join(pkgdir, 'data', 'topo_template.html')
    html_timeseries_template = os.path.join(pkgdir, 'data', 'timeseries_template.html')
    html_skewT_template = os.path.join(pkgdir, 'data', 'skewT_template.html')
    html_snowcheck_template = os.path.join(pkgdir, 'data', 'snowcheck_template.html')
    test_ts_df = os.path.join(pkgdir, 'data', 'test_df_timeseries.pkl')
    test_hgt = os.path.join(pkgdir, 'data', 'test_hgt.nc')

    # minimum and maximum elevations for topography plot
    topo_min = 0
    topo_max = 3200
else:
    print('The specified WRF output file does not exist.'
          ' Please set a valid path in cfg.py.')
    sys.exit(1)
