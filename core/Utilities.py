from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import matplotlib
import glob as glob
from cmocean import cm as cm
import numpy as np


def ThinGrid(lons_gcm,lats_gcm,gcm_data, gcm_mask, thinby=2, plot=False):
    """
    We don't want to predict at all locations.
    There are 41184 locations in the GCM grid -
    - too many to want to produce a full covariance matrix for.
    Of these, 27186 are ocean, the others are land, but that is still too many.

    As a simple fix, let's just subset by taking every nth value.
    We will also ignore the land and not predict there.

    gcm_mask should give the location of all the ocean grid cells.


    This approach reduces the number of grid points by approx 1-1/thinby**2

    """
    # create the GCM grid
    lon_grid_gcm, lat_grid_gcm = np.meshgrid(lons_gcm, lats_gcm)

    yplot = np.zeros(lon_grid_gcm.size)-10000.
    yplot[gcm_mask-1] = gcm_data # IS THIS RIGHT?
    gcm_grid = yplot.reshape(lats_gcm.size,lons_gcm.size)

    keep_lats= np.arange(0,lats_gcm.size,thinby)
    keep_lons= np.arange(0,lons_gcm.size,thinby)

    lon_grid_pred = lon_grid_gcm[keep_lats,:][:,keep_lons]
    lat_grid_pred = lat_grid_gcm[keep_lats,:][:,keep_lons]
    gcm_grid_pred = gcm_grid[keep_lats,:][:,keep_lons]

    # create an array of Falses, change the ocean values to true, then thin
    gcm_mask_TF = np.zeros(lon_grid_gcm.shape, dtype=bool)
    tmp = gcm_mask_TF.flatten()
    tmp[gcm_mask-1]=True
    gcm_mask_TF = tmp.reshape(lon_grid_gcm.shape)
    gcm_mask_TF_pred = gcm_mask_TF[keep_lats,:][:,keep_lons]

    if plot:
        mp2 = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
        llcrnrlon=-180,urcrnrlon=180,resolution='c')
        mp2.drawcoastlines()
        mp2.drawparallels(np.arange(-90.,91.,30.))
        mp2.drawmeridians(np.arange(-180.,181.,60.))
        mp2.drawmapboundary(fill_color='white')
        #mp2.scatter(lon_grid_pred.flatten(), lat_grid_pred.flatten())
        mp2.scatter(lon_grid_pred.flatten()[gcm_mask_TF_pred.flatten()], lat_grid_pred.flatten()[gcm_mask_TF_pred.flatten()])
        plt.show()

    # create the X locations for the prediction grid - thinned and with land removed
    X_pred =np.column_stack((lon_grid_pred.flatten()[gcm_mask_TF_pred.flatten()], lat_grid_pred.flatten()[gcm_mask_TF_pred.flatten()]))

    # return the thinned GCM output.
    gcm_grid_pred_S = gcm_grid_pred.flatten()[gcm_mask_TF_pred.flatten()]
    if gcm_grid_pred_S.min()<-100.:
        print('Error we have not remvoved all the land successfully')
    return X_pred, gcm_grid_pred_S[:,None]

def plot_map(lon_grid, lat_grid, levels, vals=None, mask=None, obs=None):

    plt.figure(figsize=(5.5,3.0),dpi=300)
    mp2 = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180)
    plt.xlabel('lon')
    plt.ylabel('lat')
    mp2.drawmapboundary(fill_color='white')

    if  vals is not None:
        ygrid = vals.reshape(lat_grid.shape)
  #      colours=plt.cm.get_cmap('RdYlBu_r')
        colours=plt.cm.get_cmap('RdYlBu_r')
        colours.set_bad('k', 1.0)
        mp2.contourf(lon_grid, lat_grid, ygrid, levels=levels, cmap = colours, extend='both')

    if obs is not None:
        # scatter does not easily support a discrete colormap, which would match the contour plots
	# https://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar
#        mp2.scatter(obs[:,0], obs[:,1], c=obs[:,2], cmap = 'RdYlBu_r', vmin=levels[0], vmax=levels[-1], s=40, edgecolors='k',zorder=20)
        mp2.scatter(obs[:,0], obs[:,1], c=obs[:,2], cmap = 'RdYlBu_r', vmin=levels[0], vmax=levels[-1], s=40, edgecolors='k',zorder=20)

    # this kind of works
    # need to blank out some of the level labels
    # when obs are included, there is only half the colorbar
    
    # colours=plt.cm.get_cmap('RdYlBu_r')
    # bounds = np.linspace(0,levels.size-1,levels.size)
    # norm = cols.BoundaryNorm(bounds, colours.N)
    # mp2.colorbar(cmap=colours, norm=norm, spacing='proportional', ticks=levels, boundaries=bounds)
    
    mp2.colorbar()

    if  mask is not None:
        mask_grid = mask.reshape(lat_grid.shape)
        mp2.contour(lon_grid, lat_grid, mask_grid, levels=0.5, colors='grey')

    return(mp2)
