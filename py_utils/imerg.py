import xarray as xr 

def barbados_subset(ds):
    '''subsample barbados region for relevant var'''
    # N.B. slice just returns everything in that range
    return ds['precipitationCal'].sel(lat=slice(10, 15), lon=slice(-60, -55))

def get_data(in_file, var):
    """Load-and-chunk function for IMERG HDF files
    Data is in HDF5 format, so need to specify 
    'Grid' group for xarray to use. """
    data = xr.open_dataset(in_file, engine='h5netcdf', group='Grid', chunks={'time':1, 'lon':100, 'lat': 100})
    return data[var]

def load_imerg_data(years, var='precipitationCal'):
    """ 
    Loading Matt C's IMERG satellite data from the 
    Earth Observation shared storage. *Lazy*

    Args:
    * years (list):
        list of years (str) to load data in for
    * var (str):
        which variable to extract, default precipCal
    
    Output:
    * da (xr.DataArray):
        DataArray (time, lon, lat) with the global, 
        half-hourly data for each year requested
    """
    try:
        import multiprocessing as mp
        import glob 
    except:
        raise Exception("This function requires"
                       +"'glob' and 'multiprocessing'!")

    if type(years)!=list:
        raise Exception('"years" input needs to be a list,'
                        +f' is actually type {type(years)}')
        
    # Generate filelist
    basepath='/gws/nopw/j04/eo_shared_data_vol1/satellite/imerg/Data/'
    files=[]
    for year in years:
        files.append(glob.glob(basepath+f'{year}/*/*/*'))
        
    # In parallel, lazily open each dataset, extract relevant field
    # and append to list
    pool = mp.Pool(processes=20)
    p_list = []
    f_data = []
    
    # Files has shape (n_years, ~17520)
    for _ in range(len(years)):
        for f in files[_]:
            if f.split('/')[11].split('.')[0]=='3B-HHR':
                p_list.append(pool.apply_async(get_data, args=(f,var)))

    for p in p_list:
        t_res = p.get()
        f_data.append(t_res)

    # Finally, concatenate the many many dataarrays along time axis
    da = xr.concat(f_data, dim='time')
    return da
    


def load_imerg_data_OLD(years, var='precipitationCal'):
    """
    Loading Matt C's IMERG satellite data from the 
    Earth Observation shared storage.
    
    Data is in HDF5 format, so need to specify 
    'Grid' group for xarray to use. 
    """
    try:
        import glob 
    except:
        raise Exception("This function requires `glob`!")
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