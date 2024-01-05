"""Plenty of useful functions doing useful things.  """

import os
from tempfile import mkdtemp
import shutil

import numpy as np
import pandas as pd
import xarray as xr

from wrfvis import cfg, grid, graphics, skewT

def get_wrf_timeseries(param, lon, lat, zagl, rad = 0):
    """Read the time series from the WRF output file.
    
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
                               ds[['PHB', 'PH', 'HGT', param]], ngcind, param, zagl)
                # extract time series
                if param == 'T':
                    # WRF output is perturbation potential temperature
                    vararray = ds[param][np.arange(len(ds.Time)), nlind, ngcind[0], ngcind[1]] + 300
                else:
                    vararray = ds[param][np.arange(len(ds.Time)), nlind, ngcind[0], ngcind[1]]
                df = vararray[:,0].to_dataframe()
    
                # add information about the variable
                df.attrs['variable_name'] = param
                df.attrs['variable_units'] = ds[param].units
    
                # add information about the location
                df.attrs['distance_to_grid_point'] = ngcdist
                df.attrs['lon_grid_point'] = ds.XLONG.to_numpy()[0, ngcind[0], ngcind[1]]
                df.attrs['lat_grid_point'] = ds.XLAT.to_numpy()[0, ngcind[0], ngcind[1]]
                df.attrs['grid_point_elevation_time0'] = nlhgt[0]
    
                # terrain elevation
                wrf_hgt = ds.HGT[0,:,:]
                
                # col_names have to be returned that the number of returned values stays same
                col_names = []
                
                return df, wrf_hgt, col_names
        
            # extract time series
            vararray = ds[param][np.arange(len(ds.Time)), ngcind[0], ngcind[1]]
            df = vararray[:].to_dataframe()

            # add information about the variable
            df.attrs['variable_name'] = param
            df.attrs['variable_units'] = ds[param].units

            # add information about the location
            df.attrs['distance_to_grid_point'] = ngcdist
            df.attrs['lon_grid_point'] = ds.XLONG.to_numpy()[0, ngcind[0], ngcind[1]]
            df.attrs['lat_grid_point'] = ds.XLAT.to_numpy()[0, ngcind[0], ngcind[1]]

            # terrain elevation
            wrf_hgt = ds.HGT[0,:,:]
            
            # col_names have to be returned that the number of returned values stays same
            col_names = []
            
            return df, wrf_hgt, col_names
        
        # if a radius is given
        else:        
            # find grid cells that lie inside specified radius
            gcind_inrad = grid.find_grid_cells_in_radius(ngcind, ngcdist, rad, lon, lat, ds)
            
            # check if variable is 3-dim
            # variable is 3-dim if dimensions of dataframe variable >= 4
            if len(ds[param].dims)>=4:
                
                # find nearest vertical level
                nlind, nlhgt = grid.find_nearest_vlevel(
                           ds[['PHB', 'PH', 'HGT', param]], ngcind, param, zagl)
    
                # extract time series
                if param == 'T':
                    # for 3-dim varibale T (perturbation potential Temperature) add 300 to get potential Temperature
                    # create an array with values equal to zero 
                    # with shape number of grids inside radius x 36 (time) x 36 (bottom_top)
                    vararray = np.zeros((len(gcind_inrad),36,36))
                    n = 0
                    # Loop to write values of each grid in an array
                    for i in gcind_inrad:
                        vararray[n] = ds[param][np.arange(len(ds.Time)), nlind, i[0], i[1]] + 300
                        n = n + 1
                else:
                    # for all other 4-dim variables 
                    vararray = np.zeros((len(gcind_inrad),36,36))
                    n = 0
    
                    for i in gcind_inrad:
                        vararray[n] = ds[param][np.arange(len(ds.Time)), nlind, i[0], i[1]]
                        n = n + 1
                
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
                    n = n + 1
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
                n = n + 1
             
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
                n = n + 1
            
            # terrain elevation
            wrf_hgt = ds.HGT[0,:,:]

    return df, wrf_hgt, col_names


def mkdir(path, reset=False):
    """Check if directory exists and if not, create one.
        
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


def write_html(param, lon, lat, zagl, rad = 0, directory=None):
    """ Create HTML with WRF plot 
    
    Returns
    -------
    outpath: str
        path to HTML file
    """

    # create directory for the plot
    if directory is None:
        directory = mkdtemp()
    mkdir(directory)

    # extract timeseries from WRF output
    print('Extracting timeseries at nearest grid cell')
    df, hgt, col_names = get_wrf_timeseries(param, lon, lat, zagl, rad)
   
    print('Plotting data')
    # plot the timeseries
    png = os.path.join(directory, 'timeseries.png')
    graphics.plot_ts(df,col_names, filepath=png)

    # plot a topography map
    png = os.path.join(directory, 'topography.png')
    for i in col_names:
        graphics.plot_topo(hgt, (df[param].attrs['lon_grid_point'], 
                       df[i].attrs['lat_grid_point']), filepath=png)
    #graphics.plot_topo(hgt, (df.attrs['lon_grid_point'], 
     #                  df.attrs['lat_grid_point']), filepath=png)

    # create HTML from template
    outpath = os.path.join(directory, 'index.html')
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

    return outpath

def skewT_html(lon, lat, time_index, directory=None):
    """ Create HTML with skewT plot 
    
    Returns
    -------
    outpath: str
        path to HTML file
    """
    
    #extract vertical profile for skewT 
    print('Extracting  vertical profile for skewT')
    df_skewT, zlev, lcl_pressure, lcl_temperature,lfc_pressure, lfc_temperature = skewT.skewT_and_MSEplot_dataframe(lon,lat,time_index)
    
    # plot the skewT
    print('plotting skewT')
    
    png = os.path.join('skewT.png')
    skewT.skewT_and_MSED_plot(df_skewT, df_skewT['P'], df_skewT['Temperature_in_degC'], df_skewT['dewpoint'], df_skewT['U'], df_skewT['V'], lcl_pressure, lcl_temperature, lfc_pressure, lfc_temperature, df_skewT['QVAPOR'], zlev, filepath ='skewT.png')
    
    #create HTML from template
    outpath = os.path.join(directory, 'index.html')
    with open(cfg.html_template, 'r') as infile:
        lines = infile.readlines()
        out = []
        for txt in lines:
            txt = txt.replace('[PLOTTYPE]', 'SkewT-Diagramm')
            #txt = txt.replace('[PLOTVAR]', param)
            txt = txt.replace('[IMGTYPE]', 'vertical profile')
            out.append(txt)
            
    
    with open(outpath, 'w') as outfile:
        outfile.writelines(out)
    
    
    return outpath
