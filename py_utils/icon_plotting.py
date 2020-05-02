"""
===================================================================
April 2020 --- Andrew Williams
===================================================================
Utility functions for plotting output and correlation matrices
for the ensemble of daily runs that Guy Dagan did in his ACP paper.
===================================================================
"""

import numpy as np
import xarray as xr
import pandas as pd
import seaborn as sns

from scipy.stats import ttest_ind, linregress

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches

"""
====================================================================================================================================
FUNCTIONS TO PLOT CORRELATION MATRICES
====================================================================================================================================
"""

def plot_corr_matrices(df, cloud_type, cdnc, lag=0, freq='hrly', save=False, plot_sig=True):
    """
    ==============================================================================================
    Plot three correlation matrices for N1+N2, N1 and N2. 
    
    Inputs:
    df: pd.DataFrame with columns for precip and all predictors, also require a column for 
        'NARVAL1' and 'NARVAL2', with `n_hours` rows
    
    cloud_type: 'shal', 'deep', 'high'
    cdnc: '20', '200'
    lag: int, +ve or -ve
    ==============================================================================================
    """
    
    # Load modules
    import numpy as np
    import xarray as xr
    import pandas as pd
    import seaborn as sns

    from scipy.stats import ttest_ind, linregress

    import matplotlib.pyplot as plt
    from matplotlib import cm
    import matplotlib.gridspec as gridspec
    import matplotlib.patches as patches
    
    """Subsample by cloud_type"""    
    
    if cloud_type != 'ALL':
        from loading_utils import load_cloudmask

        mask_N1 = load_cloudmask('NARVAL1', cdnc, cloud_type, freq)
        mask_N2 = load_cloudmask('NARVAL2', cdnc, cloud_type, freq)

        full_mask = np.concatenate([mask_N1, mask_N2])
        df=df[full_mask]

    
    #"""If required, lag the cloud mask and ensemble var dataframes"""
    
    """ Plot figure """
    fig = plt.figure(figsize=(10, 5), dpi=150)
    
    fig.suptitle(f'Correlation matrix, {freq} mean: {cloud_type} Clouds', 
                 y=.93, x=0.5, fontsize=15, fontweight='bold')

    """Gen axes"""
    gs = gridspec.GridSpec(nrows=1, ncols=4, 
                           width_ratios=[1, 1,1,0.1], 
                           height_ratios=[1]) # ncols, nrows

    """COLORBAR axis"""
    axcb = fig.add_subplot(gs[0, 3])
    
    """ Plot matrices for the different campaigns """
    for isubplot, campaign in enumerate(['NARVAL1+NARVAL2', 
                                         'NARVAL1', 'NARVAL2']):
        
        # Set up axis and make title
        ax = fig.add_subplot(gs[0, isubplot])
        ax.set_title(f'{campaign}', y=1.15, x=0.3, fontweight='bold')

        # Extract relevant campaigns from DataFrame
        if campaign=='NARVAL1+NARVAL2':
            df_plot=df.drop(columns='Campaign')
        else:
            df_plot=df[df['Campaign']==campaign].drop(columns='Campaign')
        
        # Perform cross-correlation
        # TODO: Introduce lag-correlation!!
        corr = df_plot.corr()
        mask = np.triu(np.ones_like(corr, dtype=np.bool))

        # Draw the heatmap with the mask and correct aspect ratio
        # *** If final axis, add a colorbar ***
        if isubplot==2:
            sns.heatmap(corr, mask=mask, cmap=cm.coolwarm, vmin=-1, vmax=1, 
                        center=0, square=True, linewidths=1, cbar_ax=axcb, 
                        cbar_kws={"shrink": 0.0}, ax=ax)
        else:
            sns.heatmap(corr, mask=mask, cmap=cm.coolwarm, cbar=False, vmin=-1,
                        vmax=1, center=0, square=True, linewidths=1, ax=ax)
        
        # Clean up fig params
        ax.tick_params(axis='x', rotation=60, labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        
        # How many valid data points?
        points = df_plot.shape[0]
        ax.text(4.3,4.5,f'{points} data points.', fontweight='bold')

        # Legend for points
        if isubplot==0:
            ax.text(3.3, 1.0, r'Dots: '    + r'$p<0.01$', bbox=dict(facecolor='lightblue', linewidth=1, edgecolor='k', alpha=0.5, pad=1))
            ax.text(3.3, 2.5, r'Crosses: ' + r'$p<0.05$', bbox=dict(facecolor='lightblue', linewidth=1, edgecolor='k', alpha=0.5, pad=1))
            #ax.text(3,2,r'$p<0.05$', bbox=dict(facecolor='lightblue', linewidth=1, edgecolor='k', alpha=0.5))
            #n = df_plot.shape[0]
            #rect = patches.Rectangle((2.8, 0.3),3.3+n/550,2,linewidth=1,
            #                         edgecolor='k',facecolor='lightblue', 
            #                         alpha=0.5)
            # Add the patches to the Axes
            #ax.add_patch(rect)

        # Plot a marker in box center if p<0.05 and plot_sig==True
        if plot_sig==True:
            for i in range(corr.shape[0]):
                for j in range(i+1, corr.shape[0]):
                    x=df_plot.iloc[:, i].values
                    y=df_plot.iloc[:, j].values
                    #print(x,y,campaign)
                    mask = ~np.isnan(x+y)
                    p = linregress(x[mask],y[mask])[3]
                    if p<0.01:
                        ax.scatter(np.array([i+0.5]), np.array([j+0.5]), marker='.', color='k', linewidth=0.5, alpha=0.7)
                    elif p<0.05:
                        ax.scatter(np.array([i+0.5]), np.array([j+0.5]), marker='x', color='k', linewidth=0.5, alpha=0.7)


    """# Adjust layout, shrink colorbar aspect ratio"""
    fig.tight_layout()
    # Change box layout for colorbar
    box = axcb.get_position()
    axcb.set_position([box.x0, box.y0+0.17, box.width * 1. , box.height * 0.6])
    
    if save==True:
        plt.savefig(f'Figs/corr_matrix_{freq}_{cloud_type}_{cdnc}CDNC_CampaignDecomp_p_values.pdf', dpi=300)
    
def scatter_subplots(x, y, df, freq, cloud_type, cdnc, save_fig=False):
    """
    Requires full df for N1+N2
    
    Plots 3 subplots of x vs y for NARVAL1+NARVAL2, NARVAL1 and NARVAL2
    
    Can also save figs...
    """
    
    """Subsample by cloud_type"""    
    
    if cloud_type != 'ALL':
        from loading_utils import load_cloudmask

        mask_N1 = load_cloudmask('NARVAL1', cdnc, cloud_type, freq)
        mask_N2 = load_cloudmask('NARVAL2', cdnc, cloud_type, freq)

        full_mask = np.concatenate([mask_N1, mask_N2])
        df_plt=df[full_mask]
        
    """Plot scatter graph"""
    units_dict = {'PRECIP': 'mm/hr', 'CAPE': 'J/kg', 'CIN': 'J/kg', 
                  'WS10': 'm/s', 'T2m': 'K', 'SST': 'K', 'LTS': 'K',
                  'RH850': '%', 'RH700': '%', 'w500': 'm/s', 'T850': 'K'}
    
    fig, ax = plt.subplots(nrows=1,ncols=3, dpi=200, figsize=(12, 4))
    
    ax[0].scatter(df_plt[x], df_plt[y], s=10)
    ax[1].scatter(df_plt[df_plt['Campaign']=='NARVAL1'][x], df_plt[df_plt['Campaign']=='NARVAL1'][y], s=10)
    ax[2].scatter(df_plt[df_plt['Campaign']=='NARVAL2'][x], df_plt[df_plt['Campaign']=='NARVAL2'][y], s=10)
    
    """Set labels and ax limits based off N1+N2 plot"""
    xlim, ylim = ax[0].get_xlim(), ax[0].get_ylim()
    ax[1].set_xlim(xlim)
    ax[1].set_ylim(ylim)
    ax[2].set_xlim(xlim)
    ax[2].set_ylim(ylim)
    
    ax[0].set_ylabel(y + f' [{units_dict[y]}]', color='red')
    ax[0].set_xlabel(x + f' [{units_dict[x]}]', color='red')
    ax[1].set_xlabel(x + f' [{units_dict[x]}]', color='red')
    ax[2].set_xlabel(x + f' [{units_dict[x]}]', color='red')
    
    ax[0].set_title('NARVAL1+NARVAL2')
    ax[1].set_title('NARVAL1')
    ax[2].set_title('NARVAL2')
    
    fig.tight_layout()
    
    """Add some explanatory text"""
    ax[0].text(0.65, 0.85, 
               f'{cdnc} CDNC' + '\n' + f'{freq} means', 
               transform=ax[0].transAxes, bbox=dict(facecolor='lightblue', linewidth=1, edgecolor='k', alpha=0.5))
    #plt.show()
    if save_fig == True:
        plt.savefig(f'Figs/April_2020/{freq}/scatterplot_{freq}_{y}_v_{x}_{cloud_type}clds_{cdnc}CDNC.pdf', dpi=300)
     