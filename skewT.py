"""
Module to create SkewT plots from WRF output data.

This module provides functions to read WRF (Weather Research and Forecasting) output files, extract vertical profiles, and generate SkewT plots. Additionally, it includes functionality to plot Moist Static Energy (MSE) diagrams.

Functions
---------
1. `skewt_and_mseplot_dataframe(lon, lat, time_index)`: 
    Reads the time series from the WRF output file and returns a DataFrame with vertical profiles and related attributes.

2. `skewt_and_mseplot(df_skewt, pressure, temperature, dewpoint, uwind, vwind, lcl_pressure, lcl_temperature, lfc_pressure, lfc_temperature, water_vapor, zlev, prof, filepath=None)`: 
    Plots SkewT and Moist Static Energy diagrams based on the provided data and saves the plots if a filepath is specified.

Doc Author
------
Matilda Achaab

Dependencies
------------
- sys
- pandas as pd
- xarray as xr
- numpy as np
- matplotlib.pyplot as plt
- metpy.calc as mpcalc
- metpy.units as units
- metpy.plots as SkewT
- MSEplots.plots as mpt
- wrfvis.cfg, grid

Usage
-----
1. Import the module: `import skewT `.
2. Call the desired functions:
    - For extracting vertical profiles: `skewT.skewt_and_mseplot_dataframe(lon, lat, time_index)`.
    - For plotting SkewT and MSE diagrams: `skewT.skewt_and_mseplot(df_skewt, pressure, temperature, dewpoint, uwind, vwind, lcl_pressure, lcl_temperature, lfc_pressure, lfc_temperature, water_vapor, zlev, prof, filepath=None)`.

Note
----
Ensure that the required dependencies are installed before using the module.
"""

# Contains functions to make SkewT plots
import sys
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
from metpy.plots import SkewT
from MSEplots import plots as mpt
from wrfvis import cfg, grid

plt.ioff()

def skewt_and_mseplot_dataframe(lon, lat, time_index):
    """
    @author: Matilda Achaab
    
    Read the time series from the WRF output file.
    
    Parameters
    ----------
    lon : float
        the longitude
    lat : float
        the latitude
    time_index: int
        the time for extraction of vertical profile

    Returns
    -------
    df_skewT: pd.DataFrame 
        timeseries variables of the vertical profile with additional attributes (grid cell lon, lat, dist, ...)
    lcl_pressure : float
        Pressure (hPa) at the Lifted Condensation Level (LCL).
    lcl_temperature : float
        Temperature (°C) at the Lifted Condensation Level (LCL).
    lfc_pressure : float
        Pressure (hPa) at the Level of Free Convection (LFC).
    lfc_temperature : float
        Temperature (°C) at the Level of Free Convection (LFC).
    zlev : array
        Altitude values (in meters).
    """
    with xr.open_dataset(cfg.wrfout) as ds:
        ngcind, ngcdist = grid.find_nearest_gridcell(
                          ds.XLONG[0,:,:], ds.XLAT[0,:,:], lon, lat)
        
        
        # convert binary times to datetime
        wrf_time = pd.to_datetime(
                   [bytes.decode(time) for time in ds.Times.data], 
                   format='%Y-%m-%d_%H:%M:%S') 
        # replace time coordinate (1-len(time)) with datetime times
        ds = ds.assign_coords({'Time': wrf_time})
        
        # Dynamically create variable names
        variable_names = ['P', 'T', 'QVAPOR', 'U', 'V']
        df_skewT= pd.DataFrame()
        try:
            for var_name in variable_names:
                # extratxting variables with the vertical profile for the skewT
                if var_name == 'P':
                    var = (ds.variables[var_name][time_index, :,ngcind[1],ngcind[0]] + ds.variables['PB'][0, :,ngcind[1],ngcind[0]]) * 0.01* units.hPa
                elif var_name == 'T':
                    var = ds.variables[var_name][time_index,:, ngcind[1],ngcind[0]] + 300
                else:
                    var = ds.variables[var_name][time_index, :, ngcind[1],ngcind[0]]
                
                    
                # Create a DataFrame for the current grid point
                df_skewT[var_name] = pd.DataFrame(var.values)
                
                # add information about the variable
                df_skewT[var_name].attrs['variable_name'] = var_name
                df_skewT.attrs['variable_units'] = ds[var_name].units
    
                # add information about the location
                df_skewT.attrs['distance_to_grid_point'] = ngcdist
                df_skewT.attrs['Xtime'] = ds.XTIME
                # Convert numpy.datetime64 to a datetime.datetime object
                time_datetime = pd.to_datetime(df_skewT.attrs['Xtime'].values)
                # Format the datetime as a string with the desired format
                df_skewT.attrs['time']=  time_datetime.tz_localize('UTC')
                df_skewT.attrs['lon_grid_point'] = ds.XLONG.to_numpy()[0, ngcind[0], ngcind[1]]
                df_skewT.attrs['lat_grid_point'] = ds.XLAT.to_numpy()[0, ngcind[0], ngcind[1]]
        except:
              print('choose a time index from 0-35')
              sys.exit()
                   
        #extracting the hght for the MSE plot
        Zlev=(ds['PHB'][:,:,ngcind[0],ngcind[1]] + ds['PH'][:,:,ngcind[0],ngcind[1]]) / 9.81
        
        #adding calculations for the LCL,LFC,Temperature_in_degC,dewpoint
        
        # Calculating the lifting condensation level
        df_skewT['Pressure'] = np.round(df_skewT['P'].values*units('hPa'),2)
        
        #calculating and converting actual temperature from potential temperature to celsuis
        actual_temp = mpcalc.temperature_from_potential_temperature(df_skewT['P'].values*units('hPa'), 
                                                                    df_skewT['T'].values*units('K'))
        df_skewT['Temperature_in_degC'] = actual_temp.to('degC')
        
        # Calculate dewpoint using specific humidity
        df_skewT['dewpoint'] = mpcalc.dewpoint_from_specific_humidity(df_skewT['P'].values*units('hPa'),
                                                                      df_skewT['Temperature_in_degC'].values*units('degC'),
                                                                      df_skewT['QVAPOR'].values * units('kg/kg'))
          
        # Calculate full parcel profile and add to plot as black line
        prof = mpcalc.parcel_profile(df_skewT['P'].values*units('hPa'),
                                                    df_skewT['Temperature_in_degC'][0]*units('degC'),
                                                    df_skewT['dewpoint'][0]*units('degC'))
        
        df_skewT['profile'] = prof.to('degC')    
        # Calculating the lifting condensation level
        lcl_pressure, lcl_temperature = mpcalc.lcl(df_skewT['P'][0]*units('hPa'),
                                                   df_skewT['Temperature_in_degC'][0]*units('degC'),
                                                   df_skewT['dewpoint'][0]*units('degC'))
        
        # Calculating the Level of free convection
        lfc_pressure, lfc_temperature = mpcalc.lfc(df_skewT['P'].values*units('hPa'),
                                                   df_skewT['Temperature_in_degC'].values*units('degC'),
                                                   df_skewT['dewpoint'].values*units('degC'))
        
    return df_skewT, Zlev,  lcl_pressure, lcl_temperature, lfc_pressure, lfc_temperature

   


def skewt_and_mseplot(df_skewT, pressure, temperature, dewpoint, uwind, vwind, lcl_pressure, lcl_temperature, lfc_pressure, lfc_temperature, water_vapor, zlev, prof, filepath=None):
    ''' 
    @authors: Matilda Achaab
    
    Plot SkewT and Moist static energy 

    Parameters
    ----------
    df_skewT: pandas dataframe
        attrs of the df.variable_names
    pressure : numpy.ndarray
        Atmospheric pressure values (in hPa) corresponding to the temperature data.
    temperature : pandas dataframe
        Atmospheric temperature values (in °C).
    dewpoint : pandas dataframe
        Dewpoint temperature values (in °C).
    uwind : numpy.ndarray
        Zonal wind component values.
    vwind : pandas dataframe
        Meridional wind component values.
    lcl_pressure : float
        Pressure (hPa) at the Lifted Condensation Level (LCL).
    lcl_temperature : float
        Temperature (°C) at the Lifted Condensation Level (LCL).
    lfc_pressure : float
        Pressure (hPa) at the Level of Free Convection (LFC).
    lfc_temperature : float
        Temperature (°C) at the Level of Free Convection (LFC).
    water_vapor : pandas dataframe
        Water vapor values.
    zlev : array
        Altitude values (in meters).
    filepath : str, optional
        Filepath to save the plot. If not provided, the plot will be displayed but not saved.

    Returns
    -------
    matplotlib.figure.Figure
        The generated Matplotlib Figure containing the SkewT and Moist Static Energy Diagrams.
    '''
    fig = plt.figure(figsize=(8, 6))
    skew = SkewT(fig, rotation=45)

    title = ('WRF time series at location {:.2f}$^{{\circ}}$E/{:.2f}$^{{\circ}}$N,'
             + '\nModel initialization time: {:%d %b %Y, %H%M} UTC')

    plt.title(title.format(df_skewT.attrs['lon_grid_point'], df_skewT.attrs['lat_grid_point'],
                           df_skewT.attrs['time'][0], loc='left'))

    # Customize labels
    skew.ax.set_ylabel('Pressure (hPa)')
    skew.ax.set_xlabel('Temperature (°C)')

    # Parcel Profile
    skew.plot(pressure.values, prof, 'k', linewidth=2)

    # Plot the temperature
    skew.plot(pressure, temperature, 'r', label='Temperature')
    skew.plot(pressure, dewpoint, 'g', label='dewpoint')
    skew.plot_barbs(pressure, uwind, vwind)
    skew.plot(lcl_pressure, lcl_temperature, 'ko', label='LCL')
    skew.plot(lfc_pressure, lfc_temperature, 'bo', label='LFC')
    # Additional Skew-T features
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()

    # Shade areas of CAPE and CIN
    skew.shade_cin(pressure.values * units('hPa'), temperature.values * units('degC'),
                   prof.values * units('degC'), dewpoint.values * units('degC'), label='CIN')

    skew.shade_cape(pressure.values * units('hPa'), temperature.values * units('degC'),
                    prof.values * units('degC'), label='CAPE')

   #plot the skewT
   if 'skewT.png' in filepath: 
        plt.savefig(filepath, dpi=150)
        plt.show()
        plt.close()
        print(f"Skew-T plot saved as: {filepath}")
        
    elif 'MSE.png' in filepath:    
        # Plot the MSE
        print('plotting MSE')
        fig,ax= mpt.msed_plots(pressure.values, temperature.values, water_vapor.values,
                            zlev.values * units.m, h0_std=2000, ensemble_size=20,
                            ent_rate=np.arange(0, 2, 0.05), entrain=False)
    
        # Save the MSE plot if filepath is provided
        plt.savefig(filepath, dpi=150)
        plt.show()
        plt.close()
        print(f"MSE plot saved as: {filepath}")
    return fig
