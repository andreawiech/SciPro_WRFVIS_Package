# This is a hard coded version string.
# Real packages use more sophisticated methods to make sure that this
# string is synchronised with `setup.py`, but for our purposes this is OK
__version__ = '0.0.2'

from wrfvis.core import write_html_multiple_gridcell, write_html_skewT, write_html_snowcheck, generate_combined_html
from wrfvis.core import get_wrf_timeseries
from wrfvis.grid import haversine
from wrfvis.grid import find_nearest_gridcell
from wrfvis.grid import find_nearest_vlevel

## makes the package a regular package instead of only a namespace package
