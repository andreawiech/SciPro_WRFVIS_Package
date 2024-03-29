"""
Plenty of useful functions doing useful things.

Module for visualizing Weather Research and Forecasting (WRF) model output data. 
It provides functions to extract time series from WRF output files, generate geographical and weather tables, 
and create HTML reports with topography maps, timeseries plots, SkewT diagrams, and Moist Static Energy (MSE) plots.

Functions
---------
1. `get_wrf_timeseries(param, lon, lat, zagl, rad=0)`: 
    Reads WRF output data and extracts time series for a given variable at a specified location. 
    Optionally, a radius can be defined to consider grid cells within a certain distance.

2. `mkdir(path, reset=False)`: 
    Checks if a directory exists and creates it if not. 
    If the directory already exists, it can be reset by erasing its content.

3. `write_html(param, lon, lat, zagl, rad=0, directory=None)`: 
    Generates an HTML report with a topography map and a timeseries plot for a given WRF variable at a specified location.

4. `skewT_html(lon, lat, time_index, directory=None)`: 
    Generates an HTML report with SkewT and Moist Static Energy (MSE) plots for a given location and time index.

5. `generate_combined_html(param, lon, lat, time_index, zagl, rad=0, directory=None)`: 
    Generates a comprehensive HTML report combining geographical and weather tables, model topography, timeseries plot, 
    SkewT diagram, and MSE diagram for a specified WRF variable, location, and time index.

Author
------
Manuela Lehner
Andrea Wiech
Malte Hildebrandt

Dependencies
------------
- os
- shutil
- numpy as np
- pandas as pd
- xarray as xr
- wrfvis.cfg, grid, graphics, tables
- skewT_and_mse as sm

Usage
-----
1. Import the module: `import wrfvis.core`.
2. Call the desired functions based on the visualization needs.
3. Ensure that the required dependencies are installed before using the module.
"""

import os
import shutil

import numpy as np
import pandas as pd
import xarray as xr

from wrfvis import cfg, grid, graphics,tables
from wrfvis import skewT_and_mse as sm

def get_wrf_timeseries(param, lon, lat, zagl, rad = 0):
    """
    Read the time series from the WRF output file.
    
    Parameters
    ----------
    param: str
        WRF output variable
    lon : float
        the longitude
    lat : float
        the latitude
    zagl : float
        height above ground level
    rad : float
        radius in which the gridcells of interest lie

    Returns
    -------
    df: pd.DataFrame 
        timeseries of param with additional attributes (grid cell lon, lat, dist, ...)
    wrf_hgt: xarray DataArray
        WRF topography
    col_names: list
        contains column names of all gridcells in the target radius
        empty if radius smaller then nearest gridcell distance
    """

    with xr.open_dataset(cfg.wrfout) as ds:
        # find nearest grid cell
        ngcind, ngcdist = grid.find_nearest_gridcell(
                          ds.XLONG[0,:,:], ds.XLAT[0,:,:], lon, lat)
        # convert binary times to datetime
        wrf_time = pd.to_datetime(
                   [bytes.decode(time) for time in ds.Times.data], 
                   format='%Y-%m-%d_%H:%M:%S') 
        # replace time coordinate (1-len(time)) with datetime times
        ds = ds.assign_coords({'Time': wrf_time})
        
        # If no radius is given or radius is smaller than nearest gridcell distance
        if rad <= ngcdist:
            # check if variable is 3-dim
            # variable is 3-dim if dimensions of dataframe variable >= 4
            if len(ds[param].dims)>=4:
                
                # find nearest vertical level
                nlind, nlhgt = grid.find_nearest_vlevel(
                               ds[['PHB', 'PH', 'HGT', param]], ngcind,
                               param, zagl)
                # extract time series
                if param == 'T':
                    # WRF output is perturbation potential temperature
                    vararray = ds[param][np.arange(len(ds.Time)), nlind,
                                         ngcind[0], ngcind[1]] + 300
                else:
                    vararray = ds[param][np.arange(len(ds.Time)),
                                         nlind, ngcind[0], ngcind[1]]
                df = vararray[:,0].to_dataframe()[param]
    
                # add information about the variable
                df.attrs['variable_name'] = param
                df.attrs['variable_units'] = ds[param].units
    
                # add information about the location
                df.attrs['distance_to_grid_point'] = ngcdist
                df.attrs['lon_grid_point'] = ds.XLONG.to_numpy()[0, ngcind[0],
                                                                 ngcind[1]]
                df.attrs['lat_grid_point'] = ds.XLAT.to_numpy()[0, ngcind[0], 
                                                                ngcind[1]]
                df.attrs['grid_point_elevation_time0'] = nlhgt[0]
    
                # terrain elevation
                wrf_hgt = ds.HGT[0,:,:]
                
                # col_names have to be returned that the number of returned values stays same
                col_names = []
                
                return df, wrf_hgt, col_names
        
            # extract time series
            vararray = ds[param][np.arange(len(ds.Time)), ngcind[0], ngcind[1]]
            df = vararray[:].to_dataframe()[param]

            # add information about the variable
            df.attrs['variable_name'] = param
            df.attrs['variable_units'] = ds[param].units

            # add information about the location
            df.attrs['distance_to_grid_point'] = ngcdist
            df.attrs['lon_grid_point'] = ds.XLONG.to_numpy()[0, ngcind[0],
                                                             ngcind[1]]
            df.attrs['lat_grid_point'] = ds.XLAT.to_numpy()[0, ngcind[0],
                                                            ngcind[1]]

            # terrain elevation
            wrf_hgt = ds.HGT[0,:,:]
            
            # col_names have to be returned that the number of returned values stays same
            col_names = []
            
            return df, wrf_hgt, col_names
        
        # if a radius is given
        else:        
            # find grid cells that lie inside specified radius
            gcind_inrad = grid.find_grid_cells_in_radius(ngcind, ngcdist, rad,
                                                         lon, lat, ds)
            
            # check if variable is 3-dim
            # variable is 3-dim if dimensions of dataframe variable >= 4
            if len(ds[param].dims)>=4:
                
                # find nearest vertical level
                nlind, nlhgt = grid.find_nearest_vlevel(
                           ds[['PHB', 'PH', 'HGT', param]], ngcind, param, zagl)
    
                # extract time series
                if param == 'T':
                    # for 3-dim varibale T (perturbation potential Temperature) 
                    # add 300 to get potential Temperature
                    # create an array with values equal to zero 
                    # with shape number of grids inside radius x 36 (time) x 36
                    # (bottom_top)
                    vararray = np.zeros((len(gcind_inrad),36,36))
                    n = 0
                    # Loop to write values of each grid in an array
                    for i in gcind_inrad:
                        vararray[n] = ds[param][np.arange(len(ds.Time)), nlind,
                                                i[0], i[1]] + 300
                        n += 1
                else:
                    # for all other 4-dim variables 
                    vararray = np.zeros((len(gcind_inrad),36,36))
                    n = 0
    
                    for i in gcind_inrad:
                        vararray[n] = ds[param][np.arange(len(ds.Time)), nlind,
                                                i[0], i[1]]
                        n += 1
                
                # create column lables needed for DataFrame
                col_names = [param]
                for i in range(1, len(gcind_inrad)):
                    col_names = col_names + [f'{param}{i}'] 
                    
                # write for each Grid the variables at hight level 0 to one DataFrame
                df = pd.DataFrame()
                
                for i in range(len(gcind_inrad)):
                    dfx = pd.DataFrame(data=vararray[i,:,0],
                                       index=wrf_time, 
                                       columns=[col_names[i]], dtype=None, copy=None)
                    df = pd.concat([df,dfx],axis=1)
    
                # add information about the variable
                df.attrs['variable_name'] = param
                df.attrs['variable_units'] = ds[param].units
    
                # add information about the location
                # information about location that changes for every gridcell
                n = 0
                for i in col_names:
                    df[i].attrs['distance_to_grid_point'] = grid.haversine(lon,lat,
                                                                           ds.XLONG[0,gcind_inrad[n,0],gcind_inrad[n,1]],
                                                                           ds.XLAT[0,gcind_inrad[n,0],gcind_inrad[n,1]])
                    df[i].attrs['lon_grid_point'] = ds.XLONG.to_numpy()[0,gcind_inrad[n,0],gcind_inrad[n,1]]
                    df[i].attrs['lat_grid_point'] = ds.XLAT.to_numpy()[0, gcind_inrad[n,0], gcind_inrad[n,1]]
                    n += 1
                # information about location that stays the same for all gridcells
                df.attrs['grid_point_elevation_time0'] = nlhgt[0]
                
                # terrain elevation
                wrf_hgt = ds.HGT[0,:,:]
                return df, wrf_hgt, col_names
                
            # if variable is 2D
            # extract time series
            vararray = np.zeros((len(gcind_inrad),36,36))
            n = 0
    
            for i in gcind_inrad:
                vararray[n] = ds[param][np.arange(len(ds.Time)), i[0], i[1]]
                n += 1
             
            # create column lables needed for DataFrame
            col_names = [param]
            for i in range(1, len(gcind_inrad)):
                col_names = col_names + [f'{param}{i}']   
             
            # write for each Grid the variables at hight level 0 to one DataFrame
            df = pd.DataFrame()
            
            for i in range(len(gcind_inrad)):
                dfx = pd.DataFrame(data=vararray[i,0,:],
                                   index=wrf_time, 
                                   columns=[col_names[i]], dtype=None, copy=None)
                df = pd.concat([df,dfx],axis=1)
                
            # add information about the variable
            df.attrs['variable_name'] = param
            df.attrs['variable_units'] = ds[param].units
    
            # add information about the location
            # information about location that changes for every gridcell
            n = 0
            for i in col_names:
                df[i].attrs['distance_to_grid_point'] = grid.haversine(lon,lat,
                                                                       ds.XLONG[0,gcind_inrad[n,0],gcind_inrad[n,1]],
                                                                       ds.XLAT[0,gcind_inrad[n,0],gcind_inrad[n,1]])
                df[i].attrs['lon_grid_point'] = ds.XLONG.to_numpy()[0,gcind_inrad[n,0],gcind_inrad[n,1]]
                df[i].attrs['lat_grid_point'] = ds.XLAT.to_numpy()[0, gcind_inrad[n,0], gcind_inrad[n,1]]
                n += 1
            
            # terrain elevation
            wrf_hgt = ds.HGT[0,:,:]

    return df, wrf_hgt, col_names


def mkdir(path, reset=False):
    """
    Check if directory exists and if not, create one.
        
    Parameters
    ----------
    path: str
        path to directory
    reset: bool 
        erase the content of the directory if it exists

    Returns
    -------
    path: str
        path to directory
    """
    
    if reset and os.path.exists(path):
        shutil.rmtree(path)
    try:
        os.makedirs(path)
    except FileExistsError:
        pass
    return path


def write_html_multiple_gridcell(param, lon, lat, zagl, rad=0, directory=None):
    """ 
    @ authors: Manuela Lehner, Malte Hildebrandt
    
    Create HTML with WRF topo and timeseries plot
    
    Returns
    -------
    html_content: str
        HTML content as a string
    outpath: str
        path to HTML file
    """

    # create directory for the plot
    if directory is None:
        directory = cfg.output_directory

    # extract timeseries from WRF output
    print('Extracting timeseries at grid cells')
    df, hgt, col_names = get_wrf_timeseries(param, lon, lat, zagl, rad)
    print('Plotting data')

    # plot a topography map
    png_topo = os.path.join(directory, 'topography.png')
    if df.shape == (36,):
        graphics.plot_topo(hgt, (df.attrs['lon_grid_point'],
                                 df.attrs['lat_grid_point']), filepath=png_topo)

    else:
        graphics.plot_topo(hgt, (df[param].attrs['lon_grid_point'], 
                                 df[param].attrs['lat_grid_point']), filepath=png_topo)
   
    # plot the timeseries
    png_ts = os.path.join(directory, 'timeseries.png')
    graphics.plot_ts(df,col_names, filepath=png_ts)


    # create HTML from template
    outpath = os.path.join(directory, 'topo_ts.html')
    
    
    with open(cfg.html_template, 'r') as infile:
        lines = infile.readlines()
        out = []
        for txt in lines:
            txt = txt.replace('[PLOTTYPE]', 'Timeseries')
            txt = txt.replace('[PLOTVAR]', param)
            txt = txt.replace('[IMGTYPE]', 'timeseries')
            out.append(txt)
    with open(outpath, 'w') as outfile:
            outfile.writelines(out)

    with open(outpath, 'r') as html_file:
        html_content = html_file.read()
        
    print(f"HTML file with topo and timeseries is saved at: {outpath}")

    return html_content, outpath

def write_html_snowcheck(lon, lat, directory=None):
    """
    @ author: Malte Hildebrandt
    Generates HTML tables containing weather and geographical information for skiing conditions.

    Parameters
    ----------
    lon : float
        User input longitude value in degrees East.
    lat : float
        User input latitude value in degrees North.
    directory : str, optional
        Directory where the HTML file and associated plots will be saved. If not provided, a temporary directory will be created.

    Returns
    -------
    html_content : str
        HTML content as a string containing weather and geographical tables.
    outpath : str
        Path to the saved HTML file.
    """

    # Create directory for the plot if not provided
    if directory is None:
        directory = cfg.output_directory

    # Open the WRF dataset
    with xr.open_dataset(cfg.wrfout) as ds:
        # Generate tables
        print("Generating geographical table...")
        table_geographical = tables.geographical_table(lon, lat, ds)
        print("Generating weather table...")
        table_weather = tables.weather_table(lon, lat, ds)

    # Use html_page function to combine tables into a single HTML content
    html_content = tables.html_page(table_weather, table_geographical, directory=directory)

    # Create HTML file path
    outpath = os.path.join(directory, 'snowcheck.html')

    print(f"HTML file with weather and geographic table is saved at: {outpath}")

    return html_content, outpath


def write_html_skewT(lon, lat, time_index, directory=None):
    """ 
    @ authors: Malte Hildebrandt, Matilda Achaab
    
    Create HTML with SkewT and MSE plots 
    
    Returns
    -------
    html_content: str
        HTML content as a string
    outpath: str
        path to HTML file
    """
    
    # create directory for the plot
    if directory is None:
        directory = cfg.output_directory

    # extract vertical profile for SkewT and MSE 
    print('Extracting vertical profile for SkewT and MSE')
    df,pressure,temp,mixing_ratio,geo_hght,dewpoint,prof, lcl_pressure,lcl_temperature,\
    lfc_pressure,lfc_temperature = sm.skewT_and_Mse_dataframe(time_index,lon,lat)
    
    print('Extracting severe weather parameters')
    sbcape, sbcin, mucape, mucin, mlcape, mlcin, kindex, totals_index \
        = sm.summary_severe_weather_par(pressure,temp,dewpoint,geo_hght,prof)

    # plot the SkewT 
    print('plotting SkewT')
    png_skewT = os.path.join(directory, 'skewt.png')
    sm.skewT_plot(df, pressure, temp, dewpoint, df['U'], df['V'],
                  lcl_pressure, lcl_temperature, lfc_pressure,mlcape, 
                  lfc_temperature, prof,sbcape, sbcin, mucape, mucin, 
                  mlcin, kindex, totals_index,filepath=png_skewT)
    
    #plot the MSE
    print('plotting MSE')
    png_mse = os.path.join(directory, 'MSE.png')
    sm.mse_plot(df,pressure,temp,mixing_ratio,geo_hght,filepath=png_mse)

    # create HTML from template
    outpath = os.path.join(directory, 'skewT_MSE.html')
    with open(cfg.html_skewT_template, 'r') as infile:
        lines = infile.readlines()
        out = []
        for txt in lines:
            if '[SKEWT_IMGPATH]' in txt and '[MSE_IMGPATH]' in txt:
                # Add section for SkewT plot
                txt = txt.replace('[SKEWT_IMGPATH]', os.path.relpath(png_skewT,
                                                                     directory))
                
                # Add section for MSE plot
                txt_mse = txt.replace('[MSE_IMGPATH]',
                                      os.path.relpath(png_mse, directory))
                
                out.append(txt)
                out.append(txt_mse)
            else:
                out.append(txt)

        with open(outpath, 'w') as outfile:
            outfile.writelines(out)

    with open(outpath, 'r') as html_file:
        html_content = html_file.read()

    print(f"HTML file with skewT and MSE plot is saved at: {outpath}")
    
    return html_content, outpath


def generate_combined_html(param, lon, lat, time_index, zagl, rad=0,
                           directory=None):
    """
    @ author: Malte Hildebrandt
    
    Generate an HTML file combining tables and plots for WRF model output visualization.
    
    Parameters:
        param (str): WRF variable name for time series.    
        lon (float): Longitude of the WRF grid cell.
        lat (float): Latitude of the WRF grid cell.
        time_index (int): Time index for the time series.
        zagl (float): Height above ground level for skewT plot.
        rad (float, optional): Radius for skewT plot. Default is 0.
        directory (str, optional): Output directory path. If None, the current working directory is used.
    
    Returns:
        str: Filepath to the generated HTML file.
    """
    print("Opening WRF dataset...")
    # Open WRF dataset
    with xr.open_dataset(cfg.wrfout) as ds:
        # Create the output directory if it doesn't exist
        if directory is None:
            directory = cfg.output_directory

        # Generate tables
        print("Generating geographical table...")
        table_geographical = tables.geographical_table(lon, lat, ds)
        print("Generating weather table...")
        table_weather = tables.weather_table(lon, lat, ds)

        # Generate plots
        print("Generating WRF timeseries plot...")
        html_timeseries, path_timeseries = write_html_multiple_gridcell(param, lon, lat, zagl,
                                                      rad=rad,
                                                      directory=directory)
        print("Generating SkewT and MSE plots...")
        html_skewT_mse, path_skewT_mse = write_html_skewT(lon, lat, time_index,
                                                    directory=directory)

        # Read the template
        with open(cfg.html_combined_template, 'r') as template_file:
            template_content = template_file.read()

        # Replace placeholders in the template
        template_content = template_content.replace("[GEOGRAPHICAL_TABLE_CONTENT]", table_geographical)
        template_content = template_content.replace("[WEATHER_TABLE_CONTENT]", table_weather)
        template_content = template_content.replace("[TOPOGRAPHY_IMG_PATH]", os.path.relpath("topography.png", directory))
        template_content = template_content.replace("[TIMESERIES_IMG_PATH]", os.path.relpath(path_timeseries, directory))
        template_content = template_content.replace("[SKEWT_IMG_PATH]", os.path.relpath(path_skewT_mse, directory))
        template_content = template_content.replace("[MSE_IMG_PATH]", os.path.relpath(path_skewT_mse, directory))

        # Save HTML file
        html_filepath = os.path.join(directory, 'Visualization_of_WRF_model.html')
        with open(html_filepath, 'w') as html_file:
            html_file.write(template_content)

        print(f"HTML file saved at: {html_filepath}")
        return html_filepath
