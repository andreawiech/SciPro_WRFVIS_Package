"""
Module to create SkewT plots from WRF output data.

This module provides functions to read WRF (Weather Research and Forecasting) output files,
extract vertical profiles, and generate SkewT plots. Additionally, it includes functionality 
to plot Moist Static Energy (MSE) diagrams.

Functions
---------
1. `extration_vertial_profile(param, lon, lat)`: 
    Reads the WRF output file and extracts the vertical profile of the specified parameter at the nearest grid cell.

2. `extration_skewT_variables(time_index, lon, lat)`: 
    Extracts the vertical profiles of multiple variables for a specific time, longitude, and latitude.

3. `convert_var_to_actual_values(df)`: 
    Converts perturbation variables to actual values and adds units.

4. `calculation_dewpoint(temp, pressure, mixing_ratio)`: 
    Calculates dewpoint temperature using specific humidity.

5. `parcel_profie(temp, pressure, dewpoint)`: 
    Calculates the parcel profile using MetPy.

6. `attris_of_skewT(pressure, temperature, dewpoint)`: 
    Calculates attributes for SkewT plots, including the Lifted Condensation Level (LCL) and Level of Free Convection (LFC).

7. `severe_weather_par(pressure, temperature, dewpoint, height)`: 
    Calculates severe weather parameters, including K-Index and Total Totals Index.

8. `mixed_layer_properties(pressure, temperature, prof)`: 
    Calculates mixed-layer parcel properties, including Mixed-Layer Convective Available Potential 
    Energy (MLCAPE) and Convective Inhibition (MLCIN).

9. `unstable_parcel_properties(pressure, temperature, dewpoint, zlev)`: 
    Calculates properties of the most unstable parcel, including Most Unstable Convective Available Potential Energy 
    (MUCAPE) and Convective Inhibition (MUCIN).

10. `surface_based_cape_cin(pressure, temperature, dewpoint)`: 
    Calculates Surface-based Convective Available Potential Energy (SBCAPE) and Convective Inhibition (SBCIN).

Author
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
1. Import the module: `import extraction_and_calculation as ec.
2. Call the desired functions:
    - For extracting vertical profiles: `ec.extration_skewT_variables(time_index, lon, lat)`.
    
    - For convert_var_to_actual_values: `ec.convert_var_to_actual_values(df)`.
    
    - For calculating dewpoint: `ec.calculation_dewpoint(temp, pressure, mixing_ratio)`.
    
    - For calculating parcel profile: `ec.parcel_profie(temp, pressure, dewpoint)`.
    
    - For calculating SkewT attributes: `ec.attris_of_skewT(pressure, temperature, dewpoint)`.
    
    - For calculating severe weather parameters: `ec.severe_weather_par(pressure, temperature,
                                                                        dewpoint, height)`.
    
    - For calculating mixed-layer properties: `ec.mixed_layer_properties(pressure, temperature,
                                                                         prof)`.
    - For calculating unstable parcel properties: `ec.unstable_parcel_properties(pressure, 
                                                                temperature, dewpoint, zlev)`.
    
    - For calculating surface-based CAPE and CIN: `ec.surface_based_cape_cin(pressure, 
                                                                temperature, dewpoint)`.


Note
----
Ensure that the required dependencies are installed before using the module.
"""


import sys
import pandas as pd
import xarray as xr
import numpy as np
import metpy.calc as mpcalc
from metpy.units import units
import cfg, grid



def extration_vertial_profile(param,lon,lat): 
    """
    Read the time series from the WRF output file.
    
    Parameters
    ----------
    param: str
        WRF output variable (only 3D variables implemented so far)
    lon : float
        the longitude
    lat : float
        the latitude

    Returns
    -------
    var_array: xarray Dataarray 
        datasets at naarest gridcell with attributes with additional attributes (grid cell lon, lat, dist, ...)
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

        # extratxting variables with the vertical profile for the skewT
        if param == 'T':
            var_array = ds[param][:,:,ngcind[0],ngcind[1]] + 300
        else:
            var_array = ds[param][:,:,ngcind[0],ngcind[1]]
        # add information about the location
        var_array.attrs['distance_to_grid_point'] = ngcdist
        var_array.attrs['lon_grid_point'] = ds.XLONG.to_numpy()[0, ngcind[0], ngcind[1]]
        var_array.attrs['lat_grid_point'] = ds.XLAT.to_numpy()[0, ngcind[0], ngcind[1]]
            #extracting the hght for the MSE plot
              
    return var_array
    
def extration_skewT_variables(time_index,lon,lat): 
    """
    Extract the vertical profile of specific variables.
    
    Parameters
    ----------
    time_index: integer
        the time for extracting the vertical profile
    lon : float
        the longitude
    lat : float
        the latitude

    Returns
    -------
    df: pd.DataFrame 
        vertical extraction of param with additional attributes ( variable name, units, time, ...)
    """
    
    df = pd.DataFrame()
    # Dynamically create variable names
    variable_names = ['P', 'T', 'QVAPOR', 'U', 'V','PB','PHB','PH']
    try:
        for var_name in variable_names:
            extract_var=extration_vertial_profile(var_name,lon,lat)
            var=extract_var[time_index,:]
            
            # Create a DataFrame for the current grid point
            df[var_name] = pd.DataFrame(var.values)
            
            # add information about the variable
            df[var_name].attrs['variable_name'] = var_name
            df[var_name].attrs['units'] = extract_var.attrs.get('units')
        
            df[var_name].attrs['Xtime'] = extract_var.XTIME
            # Convert numpy.datetime64 to a datetime.datetime object
            time_datetime = pd.to_datetime(df[var_name].attrs['Xtime'].values)
            # Format the datetime as a string with the desired format
            df.attrs['time']=  time_datetime.tz_localize('UTC')
            df.attrs['lon_grid_point'] = extract_var.lon_grid_point
            df.attrs['lat_grid_point'] = extract_var.lat_grid_point
    except:
          print('choose a time index from 0-35')
          sys.exit()
    return df


def convert_var_to_actual_values(df):
    """
    converting pertubations to actual values and adding units.
    
    Parameters
    ----------
    df: pd.Datafram
       The vertical profile of perturbation variables

    Returns
    -------
    df: pd.DataFrame 
        total pressure, temperature, mixing ratio and geo-potential height dataframe
    """
    
    #coverting pressure to total pressure
    df['Pressure'] = (df['P'] +df['PB']) *0.01
    
    # rounding to 2.d place
    df['Pressure'] = np.round(df['Pressure'].values*units('hPa'),2)
    
    #calculating geopotential height
    df['geopot_hgt'] = (df['PHB'] +df['PH'])/9.81
   
    #calculating and converting actual temperature from potential temperature to celsuis
    actual_temp= mpcalc.temperature_from_potential_temperature(df['Pressure'].values*units('hPa'), 
                                                               df['T'].values*units('K'))
    
    df['Temperature_to_deg'] = np.round(actual_temp.to('degC'),2)
     
    #assigning units to mixing ratio
    df['QVAPOR'] = df['QVAPOR'].values * units('kg/kg')

    return df['Temperature_to_deg'],df['Pressure'],df['QVAPOR'],df['geopot_hgt']

    
       
def calculation_dewpoint(temp,pressure,mixing_ratio):
    """
    calculating dewpoint temperature using metpy package.
    
    Parameters
    ----------
    temp: float
        temperature in celsuis
    pressure: float
         pressure in hPa
    mixing ratio: float
         mixing ratio dimensionless
    
    Returns
    -------
    dewpoint: pd.DataFrame 
        dewpoint temperature of atmosphere
    """
    
    # Calculate dewpoint using specific humidity
    dewpoint= mpcalc.dewpoint_from_specific_humidity(pressure.values*units('hPa'),
                                                     temp.values*units('degC'),
                                                      mixing_ratio.values*units('kg/kg'))
    #rounding to 2.d place
    dewpoint = np.round(dewpoint,2)
    return dewpoint

def parcel_profie(temp,pressure,dewpoint):
    """
    calculating the parcel profile using metpy.
    
    Parameters
    ----------
    temp: float
        temperature in celsuis
    pressure: float
         pressure in hPa
    dewpoint: float
         mdewpoint temperature in degC
    Returns
    -------
    prof: pd.DataFrame 
        parcel profile in degC
    """
    
    prof = mpcalc.parcel_profile(pressure.values*units('hPa'),
                                 temp[0]*units('degC'),dewpoint[0])
    
    #convert profile temperature to right remperature
    prof = prof.to('degC') 
    
    return prof

def attris_of_skewT(pressure,temperature,dewpoint): 
    """
    calculating lifting condesating level and level of free convection using metpy package.
    
    Parameters
    ----------
    temp:erature float
        temperature in celsuis
    pressure: float
         pressure in hPa
    dewpoint: float
         dewpoint temperature in degC
    
    Returns
    -------
    LCL: float
        pressure and temperature values at lifting condensation level
    LFC: float
        Pressure and tempearure values at level of free convection
    """
    
    
    # Calculating the lifting condensation level
    lcl_pressure, lcl_temperature = mpcalc.lcl(pressure[0]*units('hPa'),
                                               temperature[0]*units('degC'),
                                               dewpoint[0])
     
    # Calculating the Level of free convection
    lfc_pressure, lfc_temperature = mpcalc.lfc(pressure.values*units('hPa'),
                                               temperature.values*units('degC'),
                                               dewpoint)
    
    return lcl_pressure, lcl_temperature,lfc_pressure, lfc_temperature



def severe_weather_par(pressure,temperature,dewpoint,height):
    """
    ccalculating lsevere weather conditions using metpy package.
    
    Parameters
    ----------
    temperature : dataframe
        temperature in celsuis
    pressure: dataframe
         pressure in hPa
    dewpoint: dataframe
         dewpoint temperature in degC
   height : fdataframe
         geo-potential height in meters
    
    Returns
    -------
    kindex : float
        K-Index, a measure of thunderstorm potential based on atmospheric instability.
    total_totals : float
        Total Totals Index, an index combining various parameters to assess severe 
        weather risk.
    """
    
    # Here are some classic severe parameters!
    kindex = mpcalc.k_index(pressure.values*units('hPa'), 
                            temperature.values*units('degC'), 
                            dewpoint.magnitude*units('degC'))
    
    totals_index = mpcalc.total_totals_index(pressure.values*units('hPa'), 
                                             temperature.values*units('degC'),
                                             dewpoint.magnitude*units('degC'))
    return kindex, totals_index


def mixed_layer_properties(pressure, temperature, prof):
    
    """
    Calculate mixed-layer parcel properties using MetPy package.

    Parameters
    ----------
    pressure : dataframe
        Atmospheric pressure values in hPa.
    temperature : dataframe
        Atmospheric temperature values in CdegC.
    dewpoint : dataframe
        Dewpoint temperature values in degC.
    prof : ataframe
        Atmospheric profile temperature in degC.

    Returns
    -------
    mlcape : float
        Mixed-Layer Convective Available Potential Energy (CAPE), 
        a measure of atmospheric instability.
    mlcin : float
        Mixed-Layer Convective Inhibition (CIN), 
        an indicator of the inhibition of convection in the mixed layer.
    """
    
    # Mixed-layer parcel properties!
    mlcape, mlcin = mpcalc.mixed_layer_cape_cin(pressure.values * units('hPa'), 
                                                temperature.values * units('degC'),
                                                prof.magnitude * units('degC'), 
                                                depth=50 * units.hPa)
    
    return mlcape, mlcin




def unstable_parcel_properties(pressure, temperature, dewpoint, zlev):
    """
    Calculate properties of the most unstable parcel using MetPy package.

    Parameters
    ----------
    pressure : dataframe
        Atmospheric pressure values in hPa.
    temperature : dataframe
        Atmospheric temperature values in degC.
    dewpoint : dataframe
        Dewpoint temperature values in degC.
    zlev : dataframe
        Geo-potential height values in meters.

    Returns
    -------
    mucape : float
        Most Unstable Convective Available Potential Energy (CAPE), 
        a measure of atmospheric instability.
    mucin : float
        Most Unstable Convective Inhibition (CIN), an indicator 
        of the inhibition of convection in the most unstable parcel.
    """
    
    # Most unstable parcel properties!
    mu_p, mu_t, mu_td, _ = mpcalc.most_unstable_parcel(pressure.values * units('hPa'),
                                                       temperature.values * units('degC'), 
                                                       dewpoint.magnitude * units('degC'), 
                                                       depth=50 * units.hPa)
    
    mucape, mucin = mpcalc.most_unstable_cape_cin(pressure.values * units('hPa'), 
                                                  temperature.values * units('degC'), 
                                                  dewpoint.magnitude * units('degC'), 
                                                  depth=50 * units.hPa)
    
    return mucape, mucin

    
    

def surface_based_cape_cin(pressure, temperature, dewpoint):
    """
    Calculate Surface-based Convective Available Potential Energy (CAPE) 
    and Convective Inhibition (CIN) using MetPy package.

    Parameters
    ----------
    pressure : dataframe
        Atmospheric pressure values in hPa.
    temperature : dataframe
        Atmospheric temperature values in CdegC.
    dewpoint : dataframe
        Dewpoint temperature values in degC.

    Returns
    -------
    sbcape : float
        Surface-based Convective Available Potential Energy (CAPE),
        a measure of atmospheric instability.
    sbcin : float
        Surface-based Convective Inhibition (CIN), an indicator
        of the inhibition of convection at the surface.
    """
    
    # Compute Surface-based CAPE and CIN
    sbcape, sbcin = mpcalc.surface_based_cape_cin(pressure.values * units('hPa'),
                                                  temperature.values * units('degC'), 
                                                  dewpoint.magnitude * units('degC'))
    
    return sbcape, sbcin






     

    
    
    
    
    
    
    
