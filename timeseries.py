"""
Andrew Williams -- April 2020
Useful functions for manipulating timeseries with xarray.
"""

def monthly_anomaly(da):
    """
    Function to calculate the anomaly timeseries 
    of a DataArray from the local monthly mean.
    
    https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes
    """ 
    da     = da.assign_coords(year_month=da.time.dt.strftime("%Y-%m"))
    result = da.groupby("year_month") - da.groupby("year_month").mean("time")
    
    return result


def weekly_anom(da):
    """
    Function to calculate the anomaly timeseries 
    of a DataArray from the local monthly mean.
    
    https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes
    """ 
    da = da.assign_coords(year_week=da.time.dt.strftime("%Y-%W"))
    result = da.groupby("year_week") - da.groupby("year_week").mean("time")
    return result


def season_selector(da, season):
    """
    Rough and ready way of selecting seasonal
    time indexes. It currently doesn't work 
    for DJF properly, because would take
    Jan+Feb+Dec from the same year! FIX THIS!
    """
    if season not in ['DJF', 'JJA', 'MAM', 'SON']: 
           raise Exception('Not a valid season! Please choose from [DJF, JJA, MAM, SON].')

    # Assertions to check that first month is Jan, last is Dec
    da_mnth = da.groupby('time.month').groups
    
    # Check that first month is Jan and last is Dec
    if da_mnth[1][0]!=0 or da_mnth[12][-1]!=(np.shape(da.time)[0]-1):
        raise Exception('First month not Jan and/or last month not Dec!')
    
    # Check that all other months are there too
    for idx in range(1,13):
        if idx not in da_mnth.keys():
            raise Exception(f'Month {idx} not in data!'+'We require all months.')
    
    season_dict = da.groupby('time.season').groups

    if season!='DJF':
        return da[season_dict[season]]
    
    elif season=='DJF':
        """
        DJF is 'year-skipping'
        so need to treat separately
        """
        # Select Jan, Feb, Dec from each year
        da_djf = da[season_dict['DJF']]
        
        # Remove Jan, Feb from first year
        # and Dec from last year
        da_djf = da_djf.assign_coords(year_month=da.time.dt.strftime("%Y-%m"))  
        
        first_yr = str(da_djf.time.dt.year.values[0])
        last_yr  = str(da_djf.time.dt.year.values[-1])

        condition = np.logical_and(
            da_djf.year_month!=first_yr+'-01', 
            np.logical_and(da_djf.year_month!=first_yr+'-02',
                           da_djf.year_month!=last_yr+'-12'))
        
        da_djf = da_djf.where(condition, drop=True)        
        
        return da_djf