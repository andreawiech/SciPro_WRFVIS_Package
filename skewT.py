"""contains functions to make skewT plots"""

#import os
#from tempfile import mkdtemp
#import shutil
#import sys
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import dates

import metpy.calc as mpcalc
from metpy.units import units
import metpy
#import metpy.plots as plots
from metpy.plots import SkewT
from MSEplots import plots as mpt

from wrfvis import cfg, grid
#from datetime import timezone

plt.ioff()

def skewT_dataframe(lon,lat,time_index):
    """Read the time series from the WRF output file.
    
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
    LCL: pressure and temperature value
    LFC : preesure and Temperature value
    Zlev : the height 
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
                    var = (ds.variables[var_name][time_index, :,ngcind[0],ngcind[1]] + ds.variables['PB'][0, :,ngcind[0],ngcind[1]]) * 0.01* units.hPa
                elif var_name == 'T':
                    var = ds.variables[var_name][time_index,:, ngcind[0],ngcind[1]] + 300
                else:
                    var = ds.variables[var_name][time_index, :, ngcind[0],ngcind[1]]
                
                    
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
              print('choose a time index from 0-36')
                   
        #extracting the hght for the MSE plot
        Zlev=(ds['PHB'][:, :, ngcind[0], ngcind[1]] + ds['PH'][:, :, ngcind[0], ngcind[1]]) / 9.81
        
        #adding calculations for the LCL,LFC,Temperature_in_degC,dewpoint
        
        #calculating and converting actual temperature from potential temperature to celsuis
        actual_temp = mpcalc.temperature_from_potential_temperature(df_skewT['P'].values*units('hPa'), df_skewT['T'].values*units('K'))
        df_skewT['Temperature_in_degC'] = actual_temp - 273.15*units.K
        
        # Calculate dewpoint using specific humidity
        df_skewT['dewpoint'] = mpcalc.dewpoint_from_specific_humidity(df_skewT['P'].values*units('hPa'), actual_temp, df_skewT['QVAPOR'].values * units('kg/kg'))
        
        # Calculating the lifting condensation level

        lcl_pressure, lcl_temperature = mpcalc.lcl(df_skewT['P'][0]*units('hPa'), actual_temp[0], df_skewT['dewpoint'][0]*units('degC'))
        
        # Calculating the Level of free convection
        lfc_pressure, lfc_temperature =metpy.calc.lfc(df_skewT['P'].values*units('hPa'), actual_temp, df_skewT['dewpoint'].values*units('degC'))
        
    return df_skewT, Zlev,  lcl_pressure, lcl_temperature, lfc_pressure, lfc_temperature

   

def skewT_plot(df_skewT, pressure, temperature, dewpoint, uwind, vwind, lcl_pressure, lcl_temperature, lfc_pressure, lfc_temperature, water_vapor, zlev,filepath=None): 
    ''' plot SkewT and Moist static energy 

    Parameters
    ----------
    df_SkewT: pandas dataframe
        SkewT_and_MSE plot of df_SkewT.vertical_profiles_variables
    '''
    # Plot the Skew-T diagram
    fig = plt.figure(figsize=(12, 12))
    skew = SkewT(fig, rotation=45)
    
    title = ('WRF time series at location {:.2f}$^{{\circ}}$E/{:.2f}$^{{\circ}}$N,'
             + '\nModel initialization time: {:%d %b %Y, %H%M} UTC')

    plt.title(title.format(df_skewT.attrs['lon_grid_point'], df_skewT.attrs['lat_grid_point'], df_skewT.attrs['time'][0] , loc='left'))
    
    # Customize labels
    skew.ax.set_ylabel('Pressure (hPa)')
    skew.ax.set_xlabel('Temperature (°C)')

    
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
    
    #MSE plots
    #ax = mpt.msed_plots(pressure.values, temperature.values, water_vapor.values , zlev.values*units.m, h0_std=2000, ensemble_size=20, ent_rate=np.arange(0,2,0.05), entrain=False)
    
    # Show legend
    skew.ax.legend()
    
    plt.show()

    if filepath is not None:
        plt.savefig(filepath, dpi=150)
        plt.close()
    # Show the plot
    plt.show()
    
    print(f"Skew-T plot saved as: {filepath}")  
    return fig

   # Example usage:
   # skewT_plot(pressure_data, temperature_data, dewpoint_data,
