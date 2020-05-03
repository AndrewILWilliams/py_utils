import numpy as np
import glob
import xarray as xr

def global_mean(ds):
    # calculate global avg temp, for example
    ds = unify_latlon(ds)
    # Calculate area weights
    weights      = np.cos(np.deg2rad(ds.lat))
    ds_weighted  = ds.weighted(weights)
    ds_weightavg = ds_weighted.mean(('lon', 'lat'))
    return ds_weightavg

def unify_latlon(ds):
    # Check dimensions are consistently named...
    for dim in ds.dims:
        if dim=='latitude' or dim=='longitude':
            ds=ds.rename({'longitude': 'lon','latitude': 'lat'})
    return ds

def prune_coords(ds):
    # Remove extraneous coords like 'height: 2m'
    for coord in ds.coords:
        if coord not in ['lat','lon','time']:
            ds = ds.drop(coord)
    return ds

def load_landmask(model):
    """
    The CMIP5 landmasks are usually stored here:
    /badc/cmip5/data/cmip5/
    + output1/MODELCENTRE/MODELNAME/historical/fx/atmos/fx/r0i0p0/latest/sftlf
    
    But occassionally they're stored in piControl or other
    experiments. 
    This function should organise all of these
    different paths so I don't have to worry about it. 
    """
    
    modelname_to_path = {
    'bcc-csm1-1':     'BCC/bcc-csm1-1/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'bcc-csm1-1-m':   'BCC/bcc-csm1-1-m/piControl/fx/atmos/fx/r0i0p0/latest/sftlf',
    'BNU-ESM':        'BNU/BNU-ESM/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'CanESM2':        'CCCma/CanESM2/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'CMCC-CESM':      'CMCC/CMCC-CESM/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'CNRM-CM5':       'CNRM-CERFACS/CNRM-CM5/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'CNRM-CM5-2':     'CNRM-CERFACS/CNRM-CM5-2/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'ACCESS1-0':      'CSIRO-BOM/ACCESS1-0/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'ACCESS1-3':      'CSIRO-BOM/ACCESS1-3/piControl/fx/atmos/fx/r0i0p0/latest/sftlf',
    'CSIRO-Mk3-6-0':  'CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'EC-EARTH':       'ICHEC/EC-EARTH/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'inmcm4':         'INM/inmcm4/amip/fx/atmos/fx/r0i0p0/latest/sftlf',
    'IPSL-CM5A-LR':   'IPSL/IPSL-CM5A-LR/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'IPSL-CM5A-MR':   'IPSL/IPSL-CM5A-MR/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'IPSL-CM5B-LR':   'IPSL/IPSL-CM5B-LR/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'FGOALS-g2':      'LASG-CESS/FGOALS-g2/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'HadCM3':         'MOHC/HadCM3/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'HadGEM2-CC':     'MOHC/HadGEM2-CC/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'HadGEM2-ES':     'MOHC/HadGEM2-ES/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'MPI-ESM-LR':     'MPI-M/MPI-ESM-LR/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'MPI-ESM-MR':     'MPI-M/MPI-ESM-MR/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'MPI-ESM-P':      'MPI-M/MPI-ESM-P/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'GISS-E2-H':      'NASA-GISS/GISS-E2-H/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'GISS-E2-H-CC':   'NASA-GISS/GISS-E2-H-CC/piControl/fx/atmos/fx/r0i0p0/latest/sftlf',
    'GISS-E2-R':      'NASA-GISS/GISS-E2-R/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'GISS-E2-R-C':    'NASA-GISS/GISS-E2-R-CC/piControl/fx/atmos/fx/r0i0p0/latest/sftlf',
    'CCSM4':          'NCAR/CCSM4/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'NorESM1-M':      'NCC/NorESM1-M/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'NorESM1-ME':     'NCC/NorESM1-ME/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'GFDL-CM2p1':     'NOAA-GFDL/GFDL-CM2p1/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'GFDL-CM3':       'NOAA-GFDL/GFDL-CM3/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'GFDL-ESM2G':     'NOAA-GFDL/GFDL-ESM2G/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'GFDL-ESM2M':     'NOAA-GFDL/GFDL-ESM2M/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'CESM1-BGC':      'NSF-DOE-NCAR/CESM1-BGC/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'CESM1-CAM5':     'NSF-DOE-NCAR/CESM1-CAM5/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'CESM1-FASTCHEM': 'NSF-DOE-NCAR/CESM1-FASTCHEM/historical/fx/atmos/fx/r0i0p0/latest/sftlf',
    'CESM1-WACCM':    'NSF-DOE-NCAR/CESM1-WACCM/historical/fx/atmos/fx/r0i0p0/latest/sftlf'
    }
    print(modelname_to_path[model])
    if model not in modelname_to_path.keys():
        raise Exception(f"Oops! Don't have landfrac information for {model}.")
        
    da = xr.open_mfdataset('/badc/cmip5/data/cmip5/output1/'
                           +modelname_to_path[model]+'/*.nc')['sftlf']
        
    landmask = (da==1)
    
    return landmask