import numpy as np

from scipy.stats import norm
from scipy.stats.mstats import theilslopes

import iris
import iris.coord_categorisation as cat

from misc import find_nearest

def add_band_membership(cube, coord, name, dlat):
    """
    Add a custom coord in latitude dimension which allows averaging over N-degree lat bands, eg. N=15 deg
    -------
    N.B. 
    1) Requires dlat is greater than or equal to the lat spacing of the native grid
    2) Requires that native grid is regularly spaced in latitude (I think??)
    
    TODO: Automatic re-gridding to 1deg lat-spacing??
    -------
    Args:
    * cube (:class:`iris.cube.Cube`):
        the cube containing 'from_coord'.  The new coord will be added into it.
    * name (string):
        name of the created coordinate
    * from_coord (:class:`iris.coords.Coord` or string):
        coordinate in 'cube', or the name of one
    * dlat:
        specified width of the aggregated latitude bands. 
        N.B. dlat should be >= the lat spacing of the native grid.
    """
    
    # Checks
    assert dlat == int(dlat), "dlat needs to be an integer!"
    assert dlat > np.mean(np.diff(cube.coord('latitude').points)), "dlat should be >= the lat spacing of the native grid!"
        
    lats = cube.coord('latitude')

    test_lats = np.linspace(-90,90,181)[::dlat]

    l = []
    for i in range(len(test_lats)):
        l.append(find_nearest(lats.points, test_lats[i]))    
    
    # For some reason the first value gets looped over twice, so need to start at -2
    # Also not really needed, this is more aesthetic than anything
    global ticker 
    ticker = -2
    
    """Category function"""
    def cat_func(coord, value):
        """
        Returning a category value for a coordinate point-value.
        ----------
        N.B. This function needs updating! 
        It's okay for rough latitude bands but sometimes has spurious end-points or round-off errors
        ----------
        Args:
        * coord (:class:`iris.coords.Coord` or string):
            coordinate in 'cube', or the name of one
        * value:
            coordinate point-value
        """
        
        global ticker 
        
        if value in l[:-1]: 
            ticker += 1
    
        return ticker
    
    cat.add_categorised_coord(cube, name, coord, cat_func)
