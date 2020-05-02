import numpy as np

from scipy.stats import norm
from scipy.stats.mstats import theilslopes

def draw_bs_pairs_linreg(x, y, s, reg):
    """
    Perform pairs bootstrap for linear regression.
    1) Take [s] pairs of samples of length [len(x)] from (x,y)
    2) Use specified regression method to generate []
    -----
    Args:
    * (x,y) {1D arrays}:
        1D x,y arrays from which to draw the pairwise-samples
    * s: 
        number of bootstrap samples to draw
    * reg: 
        type of regression, 'Theil-Sen' or 'OLS'
    
    Outputs:
    * bs_slope_reps: 
        's' realizations of the slope
    * bs_intercept_reps: 
        's' realizations of the intercept (Not as useful)
    """
    
    # Set up array of indices to sample from: inds
    inds = np.arange(len(x))

    # Initialize replicates: bs_slope_reps, bs_intercept_reps
    bs_slope_reps     = np.empty(s)
    bs_intercept_reps = np.empty(shape=s)

    # Generate replicates
    for idx in range(s):
        bs_inds = np.random.choice(inds, size=len(inds)) # Re-sampling the indices (1d array requirement)
        bs_x, bs_y = x[bs_inds], y[bs_inds] # Choose the new (x,y) sample based of the random index sample above
        if reg == 'OLS':
            bs_slope_reps[idx], bs_intercept_reps[idx] = np.polyfit(bs_x, bs_y, 1)
        if reg == 'theil-sen':
            bs_slope_reps[idx]     = theilslopes(bs_y, bs_x)[0]
            bs_intercept_reps[idx] = theilslopes(bs_y, bs_x)[1]
            
    return bs_slope_reps, bs_intercept_reps

def bootstrap_resample(X, n=None):
    """ Bootstrap resample an array_like
    Parameters
    ----------
    X : array_like
      data to resample
    n : int, optional
      length of resampled array, equal to len(X) if n==None
    Results
    -------
    returns X_resamples
    """
    if n == None:
        n = len(X)

    resample_i = np.floor(np.random.rand(n)*len(X)).astype(int)
    X_resample = X[resample_i]
    return X_resample

def mk_test(x, alpha=0.05):
    """
    This function is derived from code originally posted by Sat Kumar Tomer
    (satkumartomer@gmail.com)
    See also: http://vsp.pnnl.gov/help/Vsample/Design_Trend_Mann_Kendall.htm
    The purpose of the Mann-Kendall (MK) test (Mann 1945, Kendall 1975, Gilbert
    1987) is to statistically assess if there is a monotonic upward or downward
    trend of the variable of interest over time. A monotonic upward (downward)
    trend means that the variable consistently increases (decreases) through
    time, but the trend may or may not be linear. The MK test can be used in
    place of a parametric linear regression analysis, which can be used to test
    if the slope of the estimated linear regression line is different from
    zero. The regression analysis requires that the residuals from the fitted
    regression line be normally distributed; an assumption not required by the
    MK test, that is, the MK test is a non-parametric (distribution-free) test.
    Hirsch, Slack and Smith (1982, page 107) indicate that the MK test is best
    viewed as an exploratory analysis and is most appropriately used to
    identify stations where changes are significant or of large magnitude and
    to quantify these findings.
    Input:
        x:   a vector of data
        alpha: significance level (0.05 default)
    Output:
        trend: tells the trend (increasing, decreasing or no trend)
        h: True (if trend is present) or False (if trend is absence)
        p: p value of the significance test
        z: normalized test statistics
    Examples
    --------
      >>> x = np.random.rand(100)
      >>> trend,h,p,z = mk_test(x,0.05)
    """
    n = len(x)

    # calculate S
    s = 0
    for k in range(n-1):
        for j in range(k+1, n):
            s += np.sign(x[j] - x[k])

    # calculate the unique data
    unique_x = np.unique(x)
    g = len(unique_x)

    # calculate the var(s)
    if n == g:  # there is no tie
        var_s = (n*(n-1)*(2*n+5))/18
    else:  # there are some ties in data
        tp = np.zeros(unique_x.shape)
        for i in range(len(unique_x)):
            tp[i] = sum(x == unique_x[i])
        var_s = (n*(n-1)*(2*n+5) - np.sum(tp*(tp-1)*(2*tp+5)))/18

    if s > 0:
        z = (s - 1)/np.sqrt(var_s)
    elif s < 0:
        z = (s + 1)/np.sqrt(var_s)
    else: # s == 0:
        z = 0

    # calculate the p_value
    p = 2*(1-norm.cdf(abs(z)))  # two tail test
    h = abs(z) > norm.ppf(1-alpha/2)

    if (z < 0) and h:
        trend = 'decreasing'
    elif (z > 0) and h:
        trend = 'increasing'
    else:
        trend = 'no trend'

    return trend, h, p, z

def autocorr(x):
    """
    Function to calculate autocorrelation on a timeseries.
    Uses lags over the whole time-series length.
    """
    lags = range(len(x))
    corr=[1. if l==0 else np.corrcoef(x[l:],x[:-l])[0][1] for l in lags]
    return np.array(corr)