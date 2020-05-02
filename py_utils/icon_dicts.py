"""
=================================================================
March 2020 --- Andrew Williams
=================================================================
Dictionaries to aid the loading in of variables for 
ICON-LAM simulations. 
Contains conversions from standard name-->long name
and long name --> standard name.
This isn't all of the variables 
but it's all the ones I'm using currently.
==================================================================
# Code snippet to generate dicts
ICON_3D_strd_to_long = {}
ICON_3D_long_to_strd = {}
for std_name in test3D.data_vars:
    if 'long_name' not in test3D[f'{std_name}'].attrs:
        continue
    else:
        lng_name = test3D[f'{std_name}'].attrs['long_name']
        ICON_3D_strd_to_long.update({f'{std_name}':f'{lng_name}'})
        ICON_3D_long_to_strd.update({f'{lng_name}':f'{std_name}'})
==================================================================
"""


ICON_2D_strd_to_long={'albdif':        'Shortwave albedo for diffuse radiation',
                     'albvisdif':      'UV visible albedo for diffuse radiation',
                     'albnirdif':      'Near IR albedo for diffuse radiation',
                     'sod_t':          'downward shortwave flux at TOA',
                     'sob_s':          'shortwave net flux at surface',
                     'sob_t':          'shortwave net flux at TOA',
                     'sou_t':          'shortwave upward flux at TOA',
                     'sou_s':          'shortwave upward flux at surface',
                     'swflx_par_sfc':  'downward photosynthetically active flux at surface',
                     'sodifd_s':       'shortwave diffuse downward flux at surface',
                     'thb_s':          'longwave net flux at surface',
                     'thu_s':          'longwave upward flux at surface',
                     'athb_t':         'TOA net thermal radiation mean since model start',
                     'asob_t':         'TOA net solar radiation mean since model start',
                     'athb_s':         'surface net thermal radiation mean since model start',
                     'asob_s':         'Surface net solar radiation mean since model start',
                     'asod_t':         'Top down solar radiation mean since model start',
                     'asou_t':         'Top up solar radiation mean since model start',
                     'athd_s':         'Surface down thermal radiationmean since model start',
                     'athu_s':         'Surface up thermal radiation mean since model start',
                     'asodird_s':      'Surface down solar direct rad.mean since model start',
                     'asodifd_s':      'Surface down solar diff. rad. mean since model start',
                     'asodifu_s':      'Surface up solar diff. rad. mean since model start',
                     'aswflx_par_sfc': 'Downward PAR flux mean since model start',
                     'tqv_dia':        'total column integrated water vapour (diagnostic)',
                     'tqc_dia':        'total column integrated cloud water (diagnostic)',
                     'tqi_dia':        'total column integrated cloud ice (diagnostic)',
                     't_2m':           'temperature in 2m',
                     'rh_2m':          'relative humidity in 2m',
                     'sp_10m':         'wind speed in 10m',
                     'u_10m':          'zonal wind in 10m',
                     'v_10m':          'meridional wind in 10m',
                     'shfl_s':         'surface sensible heat flux',
                     'lhfl_s':         'surface latent heat flux',
                     'tot_prec':       'total precip',
                     'clct':           'total cloud cover',
                     'pres_msl':       'mean sea level pressure',
                     'cape_ml':        'cape of mean surface layer parcel',
                     'cin_ml':         'convective inhibition of mean surface layer parcel'}

ICON_2D_long_to_strd={'Shortwave albedo for diffuse radiation': 'albdif',
                     'UV visible albedo for diffuse radiation': 'albvisdif',
                     'Near IR albedo for diffuse radiation': 'albnirdif',
                     'downward shortwave flux at TOA': 'sod_t',
                     'shortwave net flux at surface': 'sob_s',
                     'shortwave net flux at TOA': 'sob_t',
                     'shortwave upward flux at TOA': 'sou_t',
                     'shortwave upward flux at surface': 'sou_s',
                     'downward photosynthetically active flux at surface': 'swflx_par_sfc',
                     'shortwave diffuse downward flux at surface': 'sodifd_s',
                     'longwave net flux at surface': 'thb_s',
                     'longwave upward flux at surface': 'thu_s',
                     'TOA net thermal radiation mean since model start': 'athb_t',
                     'TOA net solar radiation mean since model start': 'asob_t',
                     'surface net thermal radiation mean since model start': 'athb_s',
                     'Surface net solar radiation mean since model start': 'asob_s',
                     'Top down solar radiation mean since model start': 'asod_t',
                     'Top up solar radiation mean since model start': 'asou_t',
                     'Surface down thermal radiationmean since model start': 'athd_s',
                     'Surface up thermal radiation mean since model start': 'athu_s',
                     'Surface down solar direct rad.mean since model start': 'asodird_s',
                     'Surface down solar diff. rad. mean since model start': 'asodifd_s',
                     'Surface up solar diff. rad. mean since model start': 'asodifu_s',
                     'Downward PAR flux mean since model start': 'aswflx_par_sfc',
                     'total column integrated water vapour (diagnostic)': 'tqv_dia',
                     'total column integrated cloud water (diagnostic)': 'tqc_dia',
                     'total column integrated cloud ice (diagnostic)': 'tqi_dia',
                     'temperature in 2m': 't_2m',
                     'relative humidity in 2m': 'rh_2m',
                     'wind speed in 10m': 'sp_10m',
                     'zonal wind in 10m': 'u_10m',
                     'meridional wind in 10m': 'v_10m',
                     'surface sensible heat flux': 'shfl_s',
                     'surface latent heat flux': 'lhfl_s',
                     'total precip': 'tot_prec',
                     'total cloud cover': 'clct',
                     'mean sea level pressure': 'pres_msl',
                     'cape of mean surface layer parcel': 'cape_ml',
                     'convective inhibition of mean surface layer parcel': 'cin_ml'}

ICON_3D_long_to_strd={'Zonal wind': 'u',
                     'Meridional wind': 'v',
                     'Vertical velocity': 'w',
                     'relative humidity': 'rh',
                     'Temperature': 'temp',
                     'cloud cover': 'clc',
                     'geopotential at full level cell centre': 'geopot',
                     'Specific humidity': 'qv',
                     'specific_cloud_water_content': 'qc',
                     'rain_mixing_ratio': 'qr',
                     'specific_cloud_ice_content': 'qi',
                     'Divergence': 'div'}

ICON_3D_strd_to_long={'u':    'Zonal wind',
                     'v':     'Meridional wind',
                     'w':     'Vertical velocity',
                     'rh':    'relative humidity',
                     'temp':  'Temperature',
                     'clc':   'cloud cover',
                     'geopot':'geopotential at full level cell centre',
                     'qv':    'Specific humidity',
                     'qc':    'specific_cloud_water_content',
                     'qr':    'rain_mixing_ratio',
                     'qi':    'specific_cloud_ice_content',
                     'div':   'Divergence'}