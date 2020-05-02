"""
=================================================================
March/April 2020 --- Andrew Williams
=================================================================
Utility functions for loading in data from the ensemble ICON-LAM 
simulations Guy conducted for doing correlation analyses.
Loads in variable, resamples to hourly freq, takes domain avg 
and then returns an array with that 
==================================================================
"""

import numpy as np
import xarray as xr

"""
===============================================================================
Functions to generate an hourly-resolution shallow/deep/high boolean cloud 
mask, based off of the maximum of the vertical qt profile.
Either for domain-mean data or for each pixel in the domain. 
===============================================================================
"""

def load_cloudmask(campaign, cdnc, cloud_type, freq='hrly'):
    """
    ================================================
    Loading hourly resolution shallow cloud mask 
    for domain averaged data...or calculate it.
    
    For SHALLOW clouds, we require:
    1) max in qt[P] profile >=700hPa
    2) sum(qt | P<600hPa) < 10**-6 kg/kg 
       ------------
       Note: the % cut-off can't be too harsh
            because sometimes cirrus clouds in
            a two-layer scenario can be ~10% 
            of the total!
       ------------
        
    N.B. This method doesn't classify high clouds.
         This is *not* a perfect method...!
    
    Inputs:
    campaign: 'NARVAL1' or 'NARVAL2'
    cdnc: '20' or '200'
    cloud_type: 'shal', 'deep'
    freq: 'hrly', '3hry', '6hrly', '12hrly', 'daily'
    ================================================
    """
    
    import os
    import time

    if os.path.exists('data/cloud_masks/shalcld_mask_{}_{}_{}CDNC.npy'.format(freq, campaign, cdnc)) and cloud_type=='shal':
        _=np.load('data/cloud_masks/shalcld_mask_{}_{}_{}CDNC.npy'.format(freq, campaign, cdnc))
        return _
    
    if os.path.exists('data/cloud_masks/deepcld_mask_{}_{}_{}CDNC.npy'.format(freq, campaign, cdnc)) and cloud_type=='deep':
        _=np.load('data/cloud_masks/deepcld_mask_{}_{}_{}CDNC.npy'.format(freq, campaign, cdnc))
        return _
        
    else:
        print(80*"=")
        print('NEWEST/FIXED VERSION -- using qt max and penalizing high values of sum(qt[P<600hPa])')
        print("Cloud masks not yet generated for that type/time-resolution: calculating now!")
        print("Once this is over, try loading again!")
        print(80*"=")
        
        basepath='/gws/nopw/j04/aopp/GuyDagan/barbedos_3deg_2/barbedos_3deg'
        
        if campaign=='NARVAL1': # Missing day 24
            missing_days=['24']
            fpath=basepath+'/NARVAL1/barbedos_3deg_12{}_{}CDNC_2/*3D*.nc'
        
        if campaign=='NARVAL2': 
            missing_days=[]
            fpath=basepath+'/barbedos_3deg_{}08_{}CDNC_2/*3D*.nc'
    
        days = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', 
                '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', 
                '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', 
                '31']
        
        freq_dict = {'hrly':'1H', '3hrly':'3H', '6hrly':'6H',
                     '12hrly':'12H', 'daily': '1D'}

        print("Processing {} cloud mask for campaign={}, {}CDNC".format(freq, 
                                                                        campaign, 
                                                                        cdnc)
             )
    
        # Initialize empty outp arrays and start clock
        pmax_vals = np.array([])
        frac_vals = np.array([])
        abs_vals  = np.array([])
        start_time = time.time()

        for idx, day in enumerate(days):

            if day in missing_days:
                print('Oops! Missing day {} in {}.'.format(day,campaign)
                      + '...filling with nans...')
                
                freq_dict2 = {'hrly': 1, '3hrly': 3, 
                              '6hrly': 6, '12hrly': 12,
                              'daily': 24}
                
                for dummy_idx in range(int(24/freq_dict2[freq])):
                    print(dummy_idx)
                    pmax_vals = np.append(pmax_vals, np.nan)
                    frac_vals = np.append(frac_vals, np.nan)
                    abs_vals  = np.append(abs_vals, np.nan)
                continue

            """Open the 3D files, get hrly means and take domain mean"""
            """and select only first 24 hrs"""
            _ = xr.open_mfdataset(fpath.format(day, cdnc), 
                                  combine='by_coords'
                                 )
            # Find maximum in the qt(=qc+qi+qr) profile for each period 
            q_profiles={}
            for var in ['qc', 'qi', 'qr']:
                arr = _[var].resample(time=freq_dict[freq]).mean().mean('ncells')
                arr = arr.transpose('time', ...)[:-1,...] # rmv extraneous timeidx
                q_profiles[var]=arr
                
            qt = q_profiles['qc']+q_profiles['qi']+q_profiles['qr']

            """(1) Where is the max in the qt profile?"""
            qt_max   = qt.argmax(dim='plev')
            pmax_tmp = qt.plev[qt_max].values
            
            pmax_vals = np.append(pmax_vals, pmax_tmp)
            #print(qt)
            """(2) What is the ratio between qt<600hPa and total qt?"""
            """Alternatively, just penalize absolute sum of qt<600hPa!!"""
            s  = xr.where(qt.plev<60000., qt, np.nan).sum('plev').values
            st = qt.sum('plev').values
            
            #frac_vals = np.append(frac_vals, s/st)
            abs_vals  = np.append(abs_vals, s)
            
            """Remove/unbind loaded vars"""
            del arr, _, q_profiles, s, st

            """Progress bar"""
            t_elapsed=time.time()-start_time
            print("Time elapsed = {} seconds.".format(np.round(t_elapsed, 1)))
            print('{:03d}%'.format(int( ((idx+1)/ (len(days))) *100)), 
                  "|"
                  +int( ( 100/ (len(days))) * (idx+1) )*"="
                  +int((100/(len(days))))*(len(days)-idx-1)*" "
                  +"|")
        
        """Data is loaded, now calculate the masks!"""
        print(80*"=")
        print("Now calculating masks for shal/deep clouds!")
        print("If there were missing days, there might be warnings here...IGNORE THEM.")
        print(80*"=")
        shalcld = (pmax_vals>=70000.) & (abs_vals <= 10**-6)
        deepcld = ~shalcld
        
        np.save('shalcld_mask_{}_{}_{}CDNC'.format(freq, campaign, cdnc), shalcld)
        np.save('deepcld_mask_{}_{}_{}CDNC'.format(freq, campaign, cdnc), deepcld)

        
        
def load_domavg_cloudmask_DEPRECATED(campaign, cdnc, cloud_type, freq='hrly', use_anyways=False):
    """
    ===============================================
    Loading hourly resolution cloud mask for domain
    averaged data.
    
    Inputs:
    campaign: 'NARVAL1' or 'NARVAL2'
    cdnc: '20' or '200'
    cloud_type: 'shal', 'deep', 'high'
    freq: 'hrly', '3hry', '6hrly', '12hrly'
    ===============================================
    """
    
    if use_anyways==False:
        raise Exception(
            + "\n"
            + "WARNING! This function is deprecated." 
            + "\n" 
            + "The `max_qt` condition is not robust enough to differentiate between shallow and deep convection."
            + "\n"
            + "If you really still want to use this function, set `use_anyways=True`.")
    
    import os
    import time

    if os.path.exists('data/cloud_masks/shalcld_mask_{}_{}_{}CDNC.npy'.format(freq, campaign, cdnc)) and cloud_type=='shal':
        _=np.load('data/cloud_masks/shalcld_mask_{}_{}_{}CDNC.npy'.format(freq, campaign, cdnc))
        return _
    
    if os.path.exists('data/cloud_masks/deepcld_mask_{}_{}_{}CDNC.npy'.format(freq, campaign, cdnc)) and cloud_type=='deep':
        _=np.load('data/cloud_masks/deepcld_mask_{}_{}_{}CDNC.npy'.format(freq, campaign, cdnc))
        return _
    
    if os.path.exists('data/cloud_masks/highcld_mask_{}_{}_{}CDNC.npy'.format(freq, campaign, cdnc)) and cloud_type=='high':
        _=np.load('data/cloud_masks/highcld_mask_{}_{}_{}CDNC.npy'.format(freq, campaign, cdnc))
        return _
        
    else:
        print(80*"=")
        print('NEWEST VERSION')
        print("Cloud masks not yet generated for that type/time-resolution: calculating now!")
        print("Once this is over, try loading again!")
        print(80*"=")
        
        basepath='/gws/nopw/j04/aopp/GuyDagan/barbedos_3deg_2/barbedos_3deg'
        
        if campaign=='NARVAL1': # Missing day 24
            missing_days=['24']
            fpath=basepath+'/NARVAL1/barbedos_3deg_12{}_{}CDNC_2/*3D*.nc'
        
        if campaign=='NARVAL2': 
            missing_days=[]
            fpath=basepath+'/barbedos_3deg_{}08_{}CDNC_2/*3D*.nc'
    
        days = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', 
                '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', 
                '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', 
                '31']
        
        freq_dict = {'hrly':'1H', '3hrly':'3H', '6hrly':'6H', '12hrly':'12H'}

        print("Processing {} cloud mask for campaign={}, {}CDNC".format(freq, 
                                                                        campaign, 
                                                                        cdnc)
             )
    
        # Initialize empty outp array and start clock
        pmax_vals = np.array([])
        start_time = time.time()

        for idx, day in enumerate(days):

            if day in missing_days:
                print('Oops! Missing day {} in {}, filling with nans.'.format(day, 
                                                                              campaign)
                     )
                freq_dict2 = {'hrly': 1, '3hrly': 3, 
                              '6hrly': 6, '12hrly': 12}
                
                for dummy_idx in range(int(24/freq_dict2[freq])):
                    print(dummy_idx)
                    pmax_vals = np.append(pmax_vals, np.nan)
                continue

            """Open the 3D files, get hrly means and take domain mean"""
            """and select only first 24 hrs"""
            _ = xr.open_mfdataset(fpath.format(day, cdnc), 
                                  combine='by_coords'
                                 )
            # Find maximum in the qt(=qc+qi+qr) profile for each period 
            q_profiles={}
            for var in ['qc', 'qi', 'qr']:
                arr = _[var].resample(time=freq_dict[freq]).mean().mean('ncells')
                arr = arr.transpose('time', ...)[:-1,...] # sample first 24hrs
                q_profiles[var]=arr
                
            qt = q_profiles['qc']+q_profiles['qi']+q_profiles['qr']

            qt_max = qt.argmax(dim='plev')

            pmax_tmp=qt.plev[qt_max].values
            
            pmax_vals = np.append(pmax_vals, pmax_tmp)
            
            """Remove/unbind loaded vars"""
            del arr, _, q_profiles

            """Progress bar"""
            t_elapsed=time.time()-start_time
            print("Time elapsed = {} seconds.".format(np.round(t_elapsed, 1)))
            print('{:03d}%'.format(int( ((idx+1)/ (len(days))) *100)), 
                  "|"
                  +int( ( 100/ (len(days))) * (idx+1) )*"="
                  +int((100/(len(days))))*(len(days)-idx-1)*" "
                  +"|")
        
        """Data is loaded, now calculate the masks!"""
        print(80*"=")
        print("Now calculating masks for shal/deep/high clouds!")
        print("If there were missing days, there might be warnings here...IGNORE THEM.")
        print(80*"=")
        shalcld = pmax_vals>=80000.
        deepcld = (pmax_vals < 80000.) & (pmax_vals > 20000.)
        highcld = pmax_vals<=20000.
        
        np.save('shalcld_mask_{}_{}_{}CDNC'.format(freq, campaign, cdnc), shalcld)
        np.save('deepcld_mask_{}_{}_{}CDNC'.format(freq, campaign, cdnc), deepcld)
        np.save('highcld_mask_{}_{}_{}CDNC'.format(freq, campaign, cdnc), highcld)
        
        
        
        
        
"""
===============================================================================
FUNCTIONS TO LOAD DATA VALUES INTO NUMPY ARRAY -- DROPPING ALL METADATA!
===============================================================================
"""

def load_hrly_3D_values(varname, campaign, cdnc, preslev):
    """
    ===============================================
    Load hourly, domain-mean data for the 3D vars
    
    NEED TO SPECIFY which plev to take!
    
    varname, campaign, cdnc -- string
    preslev -- list of ints (also string?) eg. [70000., 100000.] [Pa]
    ===============================================
    """
    import time
    
    basepath='/gws/nopw/j04/aopp/GuyDagan/barbedos_3deg_2/barbedos_3deg'
    
    if campaign=='NARVAL1': # Missing day 24
        missing_days=['24']
        fpath=basepath+'/NARVAL1/barbedos_3deg_12{}_{}CDNC_2/*3D*.nc'
        
    if campaign=='NARVAL2': 
        missing_days=[]
        fpath=basepath+'/barbedos_3deg_{}08_{}CDNC_2/*3D*.nc'
    
    days = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', 
            '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
            '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', 
            '31']

    print("Processing hourly data for var={}, pressure={}hPa, campaign={}, {}CDNC".format(varname, 
                                                                                          preslev/100, 
                                                                                          campaign, 
                                                                                          cdnc))
    
    # Initialize outp array and start clock
    var_outp=np.array([])
    start_time = time.time()
    
    for idx, day in enumerate(days):
        
        if day in missing_days:
            print('Oops! Missing day {} in {}, filling with nans.'.format(day, 
                                                                          campaign))
            for dummy_idx in range(24):
                var_outp=np.append(var_outp, np.nan)
            continue
        
        """Open the 3D files, get hrly means and take domain mean"""
        """and select plev"""
        _ = xr.open_mfdataset(fpath.format(day, cdnc), combine='by_coords')
        arr = _[varname].resample(time='1H').mean().mean('ncells').sel(plev=preslev).values
        
        """Save to array"""
        # Output format means the one day = 25 hrs, so discard last hour
        # OK because that hour will be picked up by the first hour of the next simulation day
        # Also, for accumulated vars like precip, differencing means we only get 24 hrs anyways
        # So, to be consistent we should only use 24 hrs of everything else too!
        var_outp = np.append(var_outp, arr[:-1])
        
        """Remove/unbind loaded vars"""
        del arr, _
        
        """Progress bar"""
        if idx==0:
            dt=time.time()-start_time

        t_elapsed=time.time()-start_time
        print("Time elapsed = {} seconds.".format(np.round(t_elapsed, 1)))
        print('{:03d}%'.format(int( ((idx+1)/ (len(days))) *100)), 
              "|"
              +int( ( 100/ (len(days))) * (idx+1) )*"="
              +int((100/(len(days))))*(len(days)-idx-1)*" "
              +"|")
    
    return var_outp



def load_hrly_2D_values(varname, campaign, cdnc):
    """
    ===============================================
    Load hourly, domain-mean data for the 2D vars
    
    Specify CDNC, or default is 20cm^-3
    ===============================================
    """
    import time
    
    basepath='/gws/nopw/j04/aopp/GuyDagan/barbedos_3deg_2/barbedos_3deg'
    
    if campaign=='NARVAL1': # Missing day 24
        missing_days=['24']
        fpath=basepath+'/NARVAL1/barbedos_3deg_12{}_{}CDNC_2/*2D*.nc'
        
    if campaign=='NARVAL2': 
        missing_days=[]
        fpath=basepath+'/barbedos_3deg_{}08_{}CDNC_2/*2D*.nc'
    
    days = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', 
            '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', 
            '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', 
            '31']
    
    print("Processing hourly data for var={}, campaign={}, CDNC={}.".format(varname, campaign, cdnc))
    
    # Initialize outp array and start clock
    var_outp=np.array([])
    start_time = time.time()
    
    for idx, day in enumerate(days):
        
        if day in missing_days:
            print('Oops! Missing day {} in {}, filling with nans.'.format(day, campaign))
            for dummy_idx in range(24):
                var_outp=np.append(var_outp, np.nan)
            continue
        
        """Open the 2D files, get hrly means and take domain mean"""
        _ = xr.open_mfdataset(fpath.format(day, cdnc), combine='by_coords')
        arr = _[varname].resample(time='1H').mean().mean('ncells').values
        
        """Save to array"""
        # Only have cumulative precip, so need to use time-differencing
        if varname=='tot_prec':
            arr=np.diff(arr)
            var_outp = np.append(var_outp, arr)
        else:
            # Output format means the one day = 25 hrs, so discard last hour
            # It's ok because that hour will be picked up by the first hour of the next simulation day
            # Also, for accumulated vars like precip, differencing means we only get 24 hrs anyways
            # So, to be consistent we should only use 24 hrs of everything else too!
            var_outp = np.append(var_outp, arr[:-1])
        
        """Remove/unbind loaded vars"""
        del arr, _
        
        """Progress bar"""
        if idx==0:
            dt=time.time()-start_time

        t_elapsed=time.time()-start_time
        print("Time elapsed = {} seconds.".format(np.round(t_elapsed, 1)))
        print('{:03d}%'.format(int( ((idx+1)/ (len(days))) *100)), 
              "|"
              +int( ( 100/ (len(days))) * (idx+1) )*"="
              +int((100/(len(days))))*(len(days)-idx-1)*" "
              +"|")
        
    return var_outp



def load_hrly_BC_values(varname, campaign):
    """
    ==============================================================================================
    Load hourly, domain-mean data for the BCs --- MOSTLY JUST SST !!
    
    NB Because the BCs are actually daily data, we just loop 24 times
    ==============================================================================================
    """
    
    print('Processing hourly BC data for {}, {}'.format(varname, campaign))
    
    basepath='/gws/nopw/j04/aopp/GuyDagan/barbedos_3deg_2/barbedos_3deg'
    
    if campaign=='NARVAL1': # Missing day 24
        missing_days=['24']
        fpath=basepath+'/NARVAL1/barbedos_3deg_12{}_20CDNC/initial_and_buondery/dei4_NARVALI_201312{}00_fg_DOM01_ML_0012.nc'
        
    if campaign=='NARVAL2': # Missing day 1
        missing_days=['01']
        fpath=basepath+'/barbedos_3deg_{}08_20CDNC/initial_and_buondery/dei4_NARVALII_201608{}00_fg_DOM01_ML_0012.nc'
        
    
    var_outp=np.array([])
    
    for day in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', 
                '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', 
                '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', 
                '31']:
        
        if day in missing_days:
            print('Missing day {}, filling with nans instead!'.format(day))
            for dummy_idx in range(24):
                var_outp=np.append(var_outp, np.nan)
            continue
        
        else:
            _ = xr.open_dataset(fpath.format(day, day))
            
            for dummy_idx in range(24):
                var_outp = np.append(var_outp, _[varname].mean('ncells').values)
            
            """Remove/unbind loaded vars"""
            del _
            
    return var_outp

"""
====================================================================================================================================
FUNCTIONS TO LOAD RELEVANT VARIABLES AS `xr.DataArray`s -- retaining metadata!
Need to decide if these are useful? 
Also, should I make functions to just load in single days? Or load in all the days? <-- If I make a single day one, can just loop??
====================================================================================================================================
"""

def load_hrly_3D_dataarray(varname, campaign, day, cdnc, preslev=None, hrs=24):
    """
    ===============================================
    Load hourly, domain-mean data for the 3D vars
    
    NEED TO SPECIFY which plev to take!
    
    varname, day, campaign, cdnc -- string
    preslev -- int eg. 70000. [Pa] 
        ^--> If None, return profile.
    ===============================================
    """
    import time
    
    if campaign=='NARVAL1' and day=='24': # Missing day 24
        raise Exception("Shit, sorry! We are missing that data at the moment. Try a different day?")
    
    days = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', 
            '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31']
    
    if day not in days:
        raise Exception("Oops! That's not a legit day! Try another, eg. '03', '23'.")

    if campaign=='NARVAL1': 
        fpath='/gws/nopw/j04/aopp/GuyDagan/barbedos_3deg_2/barbedos_3deg/NARVAL1/barbedos_3deg_12{}_{}CDNC_2/*3D*.nc'
        
    if campaign=='NARVAL2': 
        fpath='/gws/nopw/j04/aopp/GuyDagan/barbedos_3deg_2/barbedos_3deg/barbedos_3deg_{}08_{}CDNC_2/*3D*.nc'
    
    # Start clock
    start_time = time.time()
    print("Processing hourly data for var={}, pressure={}, campaign={}, {}CDNC".format(varname, 
                                                                                       'full_profile' if preslev==None else '{}hPa'.format(preslev/100), 
                                                                                       campaign, 
                                                                                       cdnc))

    """Open the 3D files, get hrly means and take domain mean"""
    """and select plev"""
    _ = xr.open_mfdataset(fpath.format(day, cdnc), combine='by_coords')
    
    if preslev==None:
        da = _[varname].resample(time='1H').mean().mean('ncells')
        
    else:
        da = _[varname].resample(time='1H').mean().mean('ncells').sel(plev=preslev)

    t_elapsed=time.time()-start_time
    print("Done! Time elapsed = {} seconds.".format(np.round(t_elapsed, 1)))

    if hrs==24:
        da=da.transpose('time', ...)[:-1,...]
    
    return da



def load_hrly_2D_dataarray(varname, campaign, day, cdnc=20, accumulated=False):
    """
    ===============================================
    Load hourly, domain-mean data for the 2D vars
    
    Specify CDNC, or default is 20cm^-3
    ===============================================
    """
    import time
    
    if campaign=='NARVAL1' and day=='24': # Missing day 24
        raise Exception("Shit, sorry! We are missing that data at the moment. Try a different day?")
    
    days = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', 
            '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31']
    
    if day not in days:
        raise Exception("Oops! That's not a legit day! Try another, eg. '03', '23'.")

    if campaign=='NARVAL1': 
        fpath='/gws/nopw/j04/aopp/GuyDagan/barbedos_3deg_2/barbedos_3deg/NARVAL1/barbedos_3deg_12{}_{}CDNC_2/*2D*.nc'
        
    if campaign=='NARVAL2': 
        fpath='/gws/nopw/j04/aopp/GuyDagan/barbedos_3deg_2/barbedos_3deg/barbedos_3deg_{}08_{}CDNC_2/*2D*.nc'
    
    # Start clock
    start_time = time.time()
    print("Processing hourly data for var={}, campaign={}, CDNC={}.".format(varname, campaign, cdnc))

    """Open the 2D files, get hrly means and take domain mean"""
    _ = xr.open_mfdataset(fpath.format(day, cdnc), combine='by_coords')
    da = _[varname].resample(time='1H').mean().mean('ncells')
    
    """Fix attrs and accumulation, if necessary"""
    # Only have cumulative precip, so need to use time-differencing
    if varname=='tot_prec' and accumulated==False:
        da=arr.diff('time')
        da.attrs['Manipulations'] = "Calculated hourly accumulated precip by taking discrete difference along time axis."

    t_elapsed=time.time()-start_time
    print("Done! Time elapsed = {} seconds.".format(np.round(t_elapsed, 1)))
        
    return da



def load_BC_dataarray(varname, campaign, day):
    """
    ==============================================================================================
    Load domain-mean data for the BCs --- MOSTLY JUST SST !!
    
    NB Because the BCs are actually daily data, we just loop 24 times
    ==============================================================================================
    """
    
    if (campaign=='NARVAL1' and day=='24') or (campaign=='NARVAL2' and day=='01'): # Missing day 24 in N1 and 01 in N2
        raise Exception("Shit, sorry! We are missing that data at the moment. Try a different day?")
    
    days = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', 
            '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31']
    
    if day not in days:
        raise Exception("Oops! That's not a legit day! Try another, eg. '03', '23'.")
    
    basepath='/gws/nopw/j04/aopp/GuyDagan/barbedos_3deg_2/barbedos_3deg'
    
    if campaign=='NARVAL1':
        fpath=basepath+'/NARVAL1/barbedos_3deg_12{}_20CDNC/initial_and_buondery/dei4_NARVALI_201312{}00_fg_DOM01_ML_0012.nc'
        
    if campaign=='NARVAL2':
        fpath=basepath+'/barbedos_3deg_{}08_20CDNC/initial_and_buondery/dei4_NARVALII_201608{}00_fg_DOM01_ML_0012.nc'

    _ = xr.open_dataset(fpath.format(day, day))
    da =  _[varname].mean('ncells')
                
    return da