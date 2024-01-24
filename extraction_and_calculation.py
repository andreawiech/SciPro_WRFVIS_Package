"""
wrfvis: Module for SkewT and Moist Static Energy (MSE) Plots from WRF Output Data

This module provides functions to analyze WRF (Weather Research and Forecasting) output files, 
extract vertical profiles, and generate SkewT plots as well as Moist Static Energy (MSE) diagrams.

Functions:
----------
1. `skewT_and_MSE_dataframe(time_index, lon, lat)`: 
    Extracts vertical profile data from WRF output files and returns a DataFrame with relevant attributes.
    - Parameters:
        - `time_index`: Integer, the time index for extracting the vertical profile.
        - `lon`: Float, the longitude.
        - `lat`: Float, the latitude.
    - Returns:
        - A tuple containing DataFrame and profiles of pressure, temperature, 
        mixing ratio, geopotential height, dewpoint, and more.

2. `summary_severe_weather_par(pressure, temp, dewpoint, zlev, prof)`: 
    Calculates and summarizes various severe weather parameters using MetPy functions.
    - Parameters:
        - `pressure`: pd.Dataframe, atmospheric pressure values in hPa.
        - `temp`: pd.Dataframe,, atmospheric temperature values in degC.
        - `dewpoint`: pd.Dataframe,, dewpoint temperature values in degc.
        - `zlev`: pd.Dataframe, geopotential height in meters.
        - `prof`: pd.Dataframe, atmospheric profile.
    - Returns:
        - Various severe weather parameters such as CAPE, CIN, K-Index, Total Totals, etc.

3. `skewT_plot(df, pressure, temperature, dewpoint, uwind, vwind, lcl_pressure, lcl_temperature, 
               lfc_pressure, mlcape, lfc_temperature, prof, sbcape, sbcin, mucape, mucin, mlcin, 
               kindex, total_totals, filepath=None)`: 
    
    Plots SkewT diagrams based on the provided data and saves the plots if a filepath is specified.
    - Parameters:
        - `df`: DataFrame, attributes of the df.variable_names.
        - `pressure`: DataFrame, atmospheric pressure values (hPa).
        - `temperature`: DataFrame, atmospheric temperature values (°C).
        - `dewpoint`: DataFrame, dewpoint temperature values (°C).
        - ... (other atmospheric variables)
        - `filepath`: Str, optional, filepath to save the plot.

4. `MSE_plot(df, pressure, temperature, water_vapor, zlev, filepath=None)`: 
    Plots Moist Static Energy (MSE) diagrams based on the provided data and saves 
    the plots if a filepath is specified.
    
    - Parameters:
        - `df`: DataFrame, attributes of the df.variable_names.
        - `pressure`: DataFrame, atmospheric pressure values (hPa).
        - `temperature`: DataFrame, atmospheric temperature values (°C).
        - `water_vapor`: DataFrame, water vapor values.
        - `zlev`: DataFrame, altitude values (meters).
        - `filepath`: Str, optional, filepath to save the plot.

Author:
-------
Matilda Achaab

Dependencies:
-------------
- numpy as np
- matplotlib.pyplot as plt
- metpy.units as units
- metpy.plots as SkewT
- MSEplots.plots as mpt
- extration_and_calculation as ec

Usage:
------
1. Import the module: `import skewT_and_MSE as sm`.
2. Call the desired functions:
    - For extracting vertical profiles: `sm.skewT_and_MSE_dataframe(time_index, lon, lat)`.
    
    - Additional Severe Weather Summary: sm.summary_severe_weather_par(pressure, temperature,
                                                                     dewpoint, geo_hght, prof)

    - For Plotting SkewT Diagrams with Severe Weather Parameters:`skewT.skewT_plot(df, pressure, 
                                                                temperature,dewpoint, uwind, 
                                                                vwind, ...)`.
    
    - For Plotting Moist Static Energy (MSE) Diagrams: sm.skewT.MSE_plot(df, pressure, temperature,
                                                        mixing_ratio, geo_hght, filepath='mse_plot.png')

Note:
-----
Ensure that the required dependencies are installed before using the module.
"""

import numpy as np
from metpy.units import units
import matplotlib.pyplot as plt
from metpy.plots import SkewT
from MSEplots import plots as mpt
from wrfvis import extraction_and_calculation as ec


def skewT_and_Mse_dataframe(time_index,lon,lat):
    """
    Extraction of all the dataset required for the skewT.
    
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
        timeseries variables of the vertical profile with additional attributes (grid cell lon, lat, dist, ...)
    pressure: pd.Dtaframe
            pressure profile.
    temperature: pd.Dataframe
           temperature profile  
    mixing ratio: pd.Dataframe
            water_vapor profile 
    dewpoint: pd.Dataframe
            dewpoint temperature profile
    lcl_pressure : float
            Pressure (hPa) at the Lifted Condensation Level (LCL).
    lcl_temperature : float
            Temperature ('°'C) at the Lifted Condensation Level (LCL).
    lfc_pressure : float
            Pressure (hPa) at the Level of Free Convection (LFC).
    lfc_temperature : float
            Temperature ('°'C) at the Level of Free Convection (LFC).
    geopotential : pd.Dataframe
        Altitude values (in meters).
    """

    df = ec.extration_skewT_variables(time_index,lon,lat)
    temp,pressure,mixing_ratio,geo_hght =  ec.convert_var_to_actual_values(df)
    dewpoint = ec.calculation_dewpoint(temp,pressure,mixing_ratio)
    prof = ec.parcel_profie(temp,pressure,dewpoint)
    lcl_pressure, lcl_temperature,lfc_pressure,lfc_temperature = ec.attris_of_skewT(pressure,temp,dewpoint)

    return (df,pressure,temp,mixing_ratio,geo_hght,dewpoint,prof,
        lcl_pressure,lcl_temperature,lfc_pressure,lfc_temperature)


def summary_severe_weather_par(pressure,temp,dewpoint,geo_hght,prof):
    """
    Calculate and summarize various severe weather parameters using MetPy functions.

    Parameters
    ----------
    pressure : array-like
        Atmospheric pressure values in hPa.
    temp : array-like
        Atmospheric temperature values in Celsius.
    dewpoint : array-like
        Dewpoint temperature values in degC.
    prof : array-like
        Atmospheric profile.
    zlev : float
        Geo-potential height in meters.

    Returns
    -------
    sbcape : float
        Surface-based Convective Available Potential Energy (CAPE).
    sbcin : float
        Surface-based Convective Inhibition (CIN).
    mucape : float
        Most Unstable Convective Available Potential Energy (CAPE).
    mucin : float
        Most Unstable Convective Inhibition (CIN).
    mlcape : float
        Mixed Layer Convective Available Potential Energy (CAPE).
    mlcin : float
        Mixed Layer Convective Inhibition (CIN).
    kindex : float
        K-Index, a measure of thunderstorm potential.
    total_totals : float
        Total Totals Index, indicating severe weather potential.

    """
    
    kindex, totals_index = ec.severe_weather_par(pressure, temp, dewpoint,geo_hght)
    mlcape, mlcin = ec.mixed_layer_properties(pressure, temp,prof)
    mucape, mucin = ec.unstable_parcel_properties(pressure, temp, dewpoint,geo_hght)
    sbcape, sbcin = ec.surface_based_cape_cin(pressure, temp, dewpoint)

    return sbcape, sbcin, mucape, mucin, mlcape, mlcin, kindex, totals_index



def skewT_plot(df, pressure, temperature, dewpoint, uwind, vwind,
               lcl_pressure, lcl_temperature, lfc_pressure,mlcape, 
               lfc_temperature, prof,sbcape, sbcin, mucape, mucin, 
               mlcin, kindex, totals_index, filepath=None):
    
    """
    @authors: Matilda Achaab
    
    Plot SkewT

    Parameters
    ----------
    df : pd.DataFrame
        Attributes of the df.variable_names.
    pressure : pd.DataFrame
        Atmospheric pressure values (in hPa) corresponding to the temperature data.
    temperature : pd.DataFrame
        Atmospheric temperature values degC.
    dewpoint : pd.DataFrame
        Dewpoint temperature values (degC.
    uwind : d.DataFrame
        Zonal wind component values in m/s
    vwind : pd.DataFrame
        Meridional wind component values in m/s
    lcl_pressure : float
        Pressure (hPa) at the Lifted Condensation Level (LCL).
    lcl_temperature : float
        Temperature (°C) at the Lifted Condensation Level (LCL).
    lfc_pressure : float
        Pressure (hPa) at the Level of Free Convection (LFC).
    lfc_temperature : float
        Temperature (°C) at the Level of Free Convection (LFC).
    prof : pd.DataFrame
        Atmospheric profile.
    sbcape : float
        Surface-based Convective Available Potential Energy (CAPE).
    sbcin : float
        Surface-based Convective Inhibition (CIN).
    mucape : float
        Most Unstable Convective Available Potential Energy (CAPE).
    mucin : float
        Most Unstable Convective Inhibition (CIN).
    mlcape : float
        Mixed Layer Convective Available Potential Energy (CAPE).
    mlcin : float
        Mixed Layer Convective Inhibition (CIN).
    kindex : float
        K-Index, a measure of thunderstorm potential.
    total_totals : float
        Total Totals Index, indicating severe weather potential.
    filepath : str, optional
        Filepath to save the plot. 

    Returns
    -------
    matplotlib.figure.Figure
        The generated Matplotlib Figure containing the SkewT.
    """
    fig = plt.figure(figsize=(15, 8))
    skew = SkewT(fig, rotation=45, rect=(0.05, 0.05, 0.50, 0.90))

    title = ('WRF time series at location {:.2f}$^{{\circ}}$E/{:.2f}$^{{\circ}}$N,'
            + '\nModel initialization time: {:%d %b %Y, %H%M} UTC')

    plt.title(title.format(df.attrs['lon_grid_point'], df.attrs['lat_grid_point'],
                           df.attrs['time'][0], loc='left'))

    # Customize labels
    skew.ax.set_ylabel('Pressure (hPa)')
    skew.ax.set_xlabel('Temperature (°C)')

    # Parcel Profile
    skew.plot(pressure.values, prof, 'k', linewidth=2)

    # Plot the temperature
    skew.plot(pressure, temperature, 'r', label='Temperature')
    skew.plot(pressure, dewpoint, 'g', label='Dewpoint')
    skew.plot_barbs(pressure, uwind, vwind)
    skew.plot(lcl_pressure, lcl_temperature, 'ko', label='LCL')
    skew.plot(lfc_pressure, lfc_temperature, 'bo', label='LFC')
    # Additional Skew-T features
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()

    # Shade areas of CAPE and CIN
    skew.shade_cin(pressure.values * units('hPa'), temperature.values * units('degC'),
                   prof.magnitude * units('degC'), dewpoint.magnitude * units('degC'), label='CIN')

    skew.shade_cape(pressure.values * units('hPa'), temperature.values * units('degC'),
                    prof.magnitude * units('degC'), label='CAPE')
    
    # Show legend
    skew.ax.legend(loc='upper left')
    
    # Add a white rectangle as background for parameter values
    fig.patches.extend([plt.Rectangle((0.59, 0.30), 0.25, 0.50,
                                      edgecolor='black', facecolor='white',
                                      linewidth=1, alpha=1, transform=fig.transFigure,
                                      figure=fig)])
    
    # Display severe weather parameters
    plt.figtext(0.60, 0.75, f'SBCAPE: {sbcape:.0f~P}', weight='bold', fontsize=15,
                color='red')
    
    plt.figtext(0.60, 0.70, f'SBCIN: {sbcin:.0f~P}', weight='bold',
                fontsize=15, color='blue')
    
    plt.figtext(0.60, 0.65, f'MLCAPE: {mlcape:.0f~P}', weight='bold', fontsize=15,
                color='red')

    
    plt.figtext(0.60, 0.60, f'MLCIN: {mlcin:.0f~P}', weight='bold', fontsize=15,
                color='blue')
   
    
    plt.figtext(0.60, 0.55, f'MUCAPE: {mucape:.0f~P} ', weight='bold', fontsize=15,
                color='red')
    
    plt.figtext(0.60, 0.50, f'MUCIN: {mucin:.0f~P}', weight='bold', fontsize=15,
                color='blue', ha='left')
    
    plt.figtext(0.60, 0.45, f'TT-INDEX: {totals_index:.0f~P}', weight='bold', fontsize=15,
                color='black', ha='left')

    plt.figtext(0.60, 0.40, f'K-INDEX: {kindex:.0f~P}', weight='bold', fontsize=15,
                color='black', ha='left')
    
    #plot the MSE
    if filepath is None: 
        plt.savefig(filepath, dpi=150)
        print(f"Skew-T plot saved as: {filepath}")
    plt.savefig('SkewT.png')
        
    return fig



def mse_plot(df,pressure,temperature,water_vapor,zlev,filepath =None):
    """
    @authors: Matilda Achaab
    
    Plot Moist static energy 

    Parameters
    ----------
    df: pandas dataframe
        attrs of the df.variable_names
    pressure : numpy.ndarray
        Atmospheric pressure values (in hPa) corresponding to the temperature data.
    temperature : pandas dataframe
        Atmospheric temperature values (in °C).
    water_vapor : pandas dataframe
        Water vapor values.
    zlev : pd.Dataframe
        Altitude values (in meters).
    filepath : str, optional
        Filepath to save the plot. If not provided, the plot will be displayed but not saved.

    Returns
    -------
    matplotlib.figure.Figure
        The generated Matplotlib Figure containing Moist Static Energy Diagrams.
    """
    
    title = ('WRF time series at location {:.2f}$^{{\circ}}$E/{:.2f}$^{{\circ}}$N,'
             + '\nModel initialization time: {:%d %b %Y, %H%M} UTC')

    plt.title(title.format(df.attrs['lon_grid_point'], df.attrs['lat_grid_point'],
                           df.attrs['time'][0], loc='left'))

    fig = plt.figure(figsize=(8, 6))
    fig, ax= mpt.msed_plots(pressure.values, temperature.values,water_vapor.values, 
                       zlev.values, h0_std=2000, ensemble_size=20,
                       ent_rate=np.arange(0, 2, 0.05), entrain=False)
    
    #plot the Mse
    if filepath is None: 
        plt.savefig(filepath, dpi=150)
        plt.close()
        print(f"MSE plot saved as: {filepath}")
    plt.savefig('mse.png')
    return fig
