import numpy as np

def find_nearest(array, test_value):
    """
    Function to return the value of a 1d array which is closest to the specified input test value.
    -------
    Used in add_band_membership and lots of other functions.
    TODO: Upgrade function to be able to take N-dimensional inputs and return scalar output
    -------
    Args:
    * array (array_like):
        1d array
    * test_value (float):
        user-input test value
    
    Output:
    * array[idx] (float):
        the value in the 1D input array which is closest to the test value.
    """
    array = np.asarray(array)
    idx = (np.abs(array - test_value)).argmin()
    return array[idx]