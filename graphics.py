"""
Module for visualizing Weather Research and Forecasting (WRF) model output data.
It contains functions to create plots related to WRF topography and timeseries.

Functions
---------
1. `plot_topo(topo, lonlat, filepath=None)`: 
    Plots WRF topography at a specified location and saves the plot if a filepath is provided.

2. `plot_ts(df, col_names, filepath=None)`: 
    Plots timeseries of WRF data and saves the plot if a filepath is provided.

Author
------
Manuela Lehner
Andrea Wiech

Dependencies
------------
- numpy as np
- matplotlib.pyplot as plt
- matplotlib.dates
- wrfvis.cfg

Usage
-----
1. Import the module: `import wrfvis.graphics`.
2. Call the desired functions for WRF data visualization.
3. Ensure that the required dependencies are installed before using the module.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import dates
from wrfvis import cfg 

#plt.ioff()

def plot_topo(topo, lonlat, filepath=None):
    ''' 
    Plot topography

    Parameters
    ----------
    topo: xarray DataArray
        WRF topography

    lonlat: tuple
        longitude, latitude of WRF grid cell
    '''
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_position([0.1, 0.1, 0.75, 0.85])
    ax.set_xlabel('Longitude ($^{\circ}$)')
    ax.set_ylabel('Latitude ($^{\circ}$)')

    clevels = np.arange(cfg.topo_min, cfg.topo_max, 200)
    hc = ax.contourf(topo.XLONG, topo.XLAT, topo.data, levels=clevels,
                     cmap='terrain', vmin=cfg.topo_min, vmax=cfg.topo_max,
                     extend='both')
    ax.scatter(*lonlat, s=30, c='black', marker='s')

    # colorbar
    cbax = fig.add_axes([0.88, 0.1, 0.02, 0.85])
    plt.axis('off')
    cb = plt.colorbar(hc, ax=cbax, fraction=1, format='%.0f')
    cb.ax.set_ylabel('$z$ (MSL)')

    plt.savefig(filepath, dpi=150)
    plt.close()

    return fig


def plot_ts(df, col_names, filepath=None):
    ''' 
    Plot timeseries

    Parameters
    ----------
    df: pandas dataframe
        timeseries of df.variable_name
    '''
    
    fig, ax = plt.subplots(figsize=(10, 4))
    if df.shape == (36,): 
        ax.plot(df, color='black')
        ax.set_ylabel(f"{df.attrs['variable_name']} ({df.attrs['variable_units']})")
         # format the datetime tick mark labels
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H%M'))
        ax.set_xlabel('Time (UTC)')
        
        # title contains information about lon, lat, z agl, and time
        try:
            title = ('WRF time series at location {:.2f}$^{{\circ}}$E/'
                     '{:.2f}$^{{\circ}}$N, grid point elevation at time 0: {:.2f}'
                     ' m a.g.l\nModel initialization time: {:%d %b %Y, %H%M} UTC')
            plt.title(title.format(df.attrs['lon_grid_point'], df.attrs['lat_grid_point'],
                                   df.attrs['grid_point_elevation_time0'],
                                   df.index[0]), loc='left')
        
            plt.savefig(filepath, dpi=150)
            
        except:   
            title = ('WRF time series at location {:.2f}$^{{\circ}}$E/'
                     '{:.2f}$^{{\circ}}$N\nModel initialization time:'
                     ' {:%d %b %Y, %H%M} UTC')
            plt.title(title.format(df.attrs['lon_grid_point'], df.attrs['lat_grid_point'],
                                   df.index[0]), loc='left')
            
            plt.savefig(filepath, dpi=150)
        
    else:
        for i in col_names:
            ax.plot(df[i], label=f"Lon: {df[i].attrs['lon_grid_point']:.2f}, "
                    f"Lat: {df[i].attrs['lat_grid_point']:.2f} ")
            ax.legend(loc='upper right')
        ax.set_ylabel(f"{df.attrs['variable_name']} ({df.attrs['variable_units']})")
        
        # format the datetime tick mark labels
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H%M'))
        ax.set_xlabel('Time (UTC)')
        # title contains information about lon, lat, z agl, and time
        try:
            title = ('WRF time series at location {:.2f}$^{{\circ}}$E/'
                     '{:.2f}$^{{\circ}}$N, grid point elevation at time 0: {:.2f}'
                     ' m a.g.l\nModel initialization time: {:%d %b %Y, %H%M} UTC')
            plt.title(title.format(float(df[df.attrs['variable_name']].attrs['lon_grid_point']),
                                   float(df[df.attrs['variable_name']].attrs['lat_grid_point']),
                                   df.attrs['grid_point_elevation_time0'],
                                   df.index[0]), loc='left')
          
            plt.savefig(filepath, dpi = 150)
        except:
            title = ('WRF time series at location {lonstr:.2f}$^{{\circ}}$E/'
                     '{latstr:.2f}$^{{\circ}}$N\nModel initialization time: {index:%d %b %Y, %H%M} UTC').format(
                        lonstr=df[df.attrs['variable_name']].attrs['lon_grid_point'],
                        latstr=df[df.attrs['variable_name']].attrs['lat_grid_point'],
                        index=df.index[0])
            plt.title(title)
            
            plt.savefig(filepath, dpi=150)

    return fig
