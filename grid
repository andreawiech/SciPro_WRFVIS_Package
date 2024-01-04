"""contains functions related to the WRF grid """

import numpy as np


def haversine(lon1, lat1, lon2, lat2):
    """Great circle distance between two (or more) points on Earth

    Author: Fabien Maussion

    Parameters
    ----------
    lon1 : float
       scalar or array of point(s) longitude
    lat1 : float
       scalar or array of point(s) longitude
    lon2 : float
       scalar or array of point(s) longitude
    lat2 : float
       scalar or array of point(s) longitude

    Returns
    -------
    the distances

    Examples:
    ---------
    >>> haversine(34, 42, 35, 42)
    82633.46475287154
    >>> haversine(34, 42, [35, 36], [42, 42])
    array([ 82633.46475287, 165264.11172113])
    """

    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
        
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    return c * 6371000  # Radius of earth in meters


def find_nearest_gridcell(xlong, xlat, lon, lat):
    """ find nearest grid cell
    
    Parameters
    ----------
    xlong: xarray DataArray
        WRF grid longitude
    xlat: xarray DataArray
        WRF grid latitude
    lon: float
        target longitude
    lat: float
        target langitude

    Returns
    -------
    ngcind: tuple
        indeces of nearest grid cell
    ngcdist: float
        distance between nearest grid cell and target location
    """

    # find grid point with the smallest horizontal distance to target location
    horizdist = haversine(lon, lat, xlong, xlat)
    ngcdist = np.min(horizdist.data)
    ngcind = np.unravel_index(np.argmin(horizdist.data), horizdist.shape)

    return ngcind, ngcdist 


def find_nearest_vlevel(ds, gcind, param, ztarget):
    """ find nearest vertical level
    
    Parameters
    ----------
    ds: xarray dataset
        needs to contain the variables PHB, PH, HGT, and param
    gcind: tuple
        indeces of WRF grid cell
    param: str
        variable for which to determine the nearest vertical level
        needed to determine whether to unstagger the grid
    ztarget: float
        target height AGL

    Returns
    -------
    nlind: numpy array, integer
        index of nearest vertical level at each time
    nlhgt: numpy array, float
        height of nearest vertical level at each time
    """

    # geopotential height
    geopot_hgt = (ds['PHB'][:,:,gcind[0],gcind[1]] + ds['PH'][:,:,gcind[0],gcind[1]]) / 9.81

    # unstagger the geopotential height to half levels
    # for our purposes it is good enough to assume that the half levels are exactly 
    # in the middle of two full levels
    if 'bottom_top' in ds[param].dims:
        geopot_hgt = 0.5 * (geopot_hgt[:, :-1] + geopot_hgt[:, 1:])

    # subtract topo height to get height above ground
    wrf_zagl = geopot_hgt - ds['HGT'][0,gcind[0],gcind[1]]

    # find nearest vertical level
    nlind = np.argmin(np.abs(wrf_zagl - ztarget).to_numpy(), axis=1)
    nlhgt = wrf_zagl.to_numpy()[np.arange(len(nlind)), nlind]

    return nlind, nlhgt


def find_direct_neighbors(gcind, rad, lon, lat, ds):
    """ select direct neighboring grid cells that lie inside the radius of interest around a WRF grid cell
    
    Parameters
    ----------
    gcind : tuple
        indeces of WRF grid cell
    rad : float
        radius in which the grid cells of interest lie
    lon : float
        target longitude
    lat : float
        target langitude
    ds : xarray dataset
        needs to contain the variables XLONG and XLAT

    Returns
    -------
    insiderad : numpy array, integer
        array with indices of the direct neighboring grid cells that are inside the radius
    """
    # indices of all neighboring gridcells
    neighbor_gcind = np.array([(gcind[0],gcind[1]), (gcind[0]+1, gcind[1]),
                                (gcind[0]-1, gcind[1]), (gcind[0], gcind[1]+1), 
                                (gcind[0], gcind[1]-1), (gcind[0]-1, gcind[1]+1), 
                                (gcind[0]-1, gcind[1]-1), (gcind[0]+1, gcind[1]-1),
                                (gcind[0]+1,gcind[1]+1)])
    
    horizdist = np.zeros(0) # empty array
    
    for gcind in neighbor_gcind:
        # calculate distance between starting point and center of each the gridcell
        horizdist = np.append(horizdist, haversine(lon, lat, ds.XLONG[0, gcind[0], gcind[1]].data,
                                                   ds.XLAT[0, gcind[0], gcind[1]].data))
    
    # select gridcells that are inside the radius
    insiderad = neighbor_gcind[horizdist <= rad] 
    
    return insiderad


def find_grid_cells_in_radius(ngcind, ngcdist, rad, lon, lat, ds):
    """ select all grid cells in a radius around the target location
    

    Parameters
    ----------
    ngcind : tuple
        indeces of WRF grid cell
    ngcdist : float
        distance between nearest grid cell and target location
    rad : float
        radius in which the grid cells of interest lie
    lon : float
        target longitude
    lat : float
        target langitude
    ds : xarray dataset
        needs to contain the variables XLONG and XLAT

    Raises
    ------
    ValueError: 
        if the radius is smaller than the distance between nearest grid cell and target location
       
    Returns
    -------
    neighboring_ngcind: 2-dim numpy array, integer
        array with indices of all the neighboring grid cells inside the radius

    """
    
    if rad < ngcdist:
        raise ValueError('Radius has to be larger than distance between nearest grid cell and target location')
        
    else:
        neighboring_ngcind = find_direct_neighbors(ngcind, rad, lon, lat, ds)
        
        if len(neighboring_ngcind) != len([ngcind]):
            neighboring_ngcind_old = [ngcind]
            
            while len(neighboring_ngcind) != len(neighboring_ngcind_old):
                neighboring_ngcind_old = neighboring_ngcind
                for i in neighboring_ngcind_old:
                    neighboring_ngcind = np.append(neighboring_ngcind_old, find_direct_neighbors(i, rad, lon, lat, ds), axis=0)
                    neighboring_ngcind = np.unique(neighboring_ngcind, axis=0)
                
            return neighboring_ngcind
        
        else:
            return neighboring_ngcind
