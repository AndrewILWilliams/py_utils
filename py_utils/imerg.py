import xarray as xr 
import glob 

def barbados_subset(ds):
    '''subsample barbados region for relevant var'''
    # N.B. slice just returns everything in that range
    return ds['precipitationCal'].sel(lat=slice(10, 15), lon=slice(-60, -55))

def load_imerg_subset(years):
    """
    Loading Matt C's IMERG satellite data from the 
    Earth Observation shared storage.
    
    Data is in HDF5 format, so need to specify 
    'Grid' group for xarray to use. 
    """
    
    if type(years)!=list:
        raise Exception('"years" input needs to be a list,'
                        +f' is actually type {type(years)}')
        
    basepath='/gws/nopw/j04/eo_shared_data_vol1/satellite/imerg/Data/'
    files=[]
    for year in years:
        files.append(glob.glob(basepath+f'{year}/*/*/*'))
    
    da = xr.open_mfdataset(
        files, engine='h5netcdf', group='Grid',
        concat_dim="time", data_vars='minimal', 
        coords='minimal', compat='override',
        preprocess=barbados_subset, parallel=True
    )
    
    return da

def load_imerg_data(years, var='precipitationCal'):
    """
    Loading Matt C's IMERG satellite data from the 
    Earth Observation shared storage.
    
    Data is in HDF5 format, so need to specify 
    'Grid' group for xarray to use. 
    """
    
    if type(years)!=list:
        raise Exception('"years" input needs to be a list,'
                        +f' is actually type {type(years)}')
        
    basepath='/gws/nopw/j04/eo_shared_data_vol1/satellite/imerg/Data/'
    files=[]
    for year in years:
        files.append(glob.glob(basepath+f'{year}/*/*/*'))
    
    ds = xr.open_mfdataset(
        files, engine='h5netcdf', group='Grid',
        concat_dim="time", data_vars='minimal', 
        coords='minimal', compat='override', 
        parallel=True
    )
    
    return ds[var]