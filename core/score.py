
from Cholesky import *
import numpy as np
from HaversineDist import Exponentialhaversine
def score(X_pred, gcm_out, X_obs, y_obs, model=None, method='gp_var'):
    """
    Code for calculating score. inputs
    - X_pred = a  n x 2 location array. One location of form (long, lat) per row. n is the number of gcm grid cells used.

    - gcm_out = a N x n array - one GCM run per row. N = number of GCMs used.

    - X_obs = a d x 2 array of locations (long, lat) of the observations. d = number of observations.

    - y_obs = a d x 1 array of observations.

    - model = a trained GP model from GPy

    - method =  one of 'gp_mean', 'gp_var', 'gp_full_cov', 'RMSE', (and I'll add RMSE_weighted later)

    gp_var computes the kriged mean and variance at all the GCM output locations. Then the loglikelihood is computed of the GCM output. Note this assumes independence between the outputs (no covariance).

    gp_full_cov is the same as gp_var, but now the full multivariate distribution is used, and the covariance estimates from the GP are used. Note that we have found this to be extremely sensitive to the fineness of the grid used, and so should probably not be used.

    gp_rmse computes the kriged mean at all the GCM output locations. Then the RMSE is computed between this and the GCM output.

    gp_rmse_cent computes the kriged mean at all the GCM output locations. Then the centred RMSE is computed between this and the GCM output.

    gp_corr computes the pattern correlation coefficient between the GP and each GCM

    gp_std computes the standard deviation coefficient of the GP and of each GCM (note first element = observation SD, for use in Taylor Diagram)

    RMSE computes the GCM grid cell closest to the observation, and then computes the RMSE between that output point and the data point it corresponds to. No scaling is done to account for the difference in measurement variance on the data points.

    RMSE_cent computes the GCM grid cell closest to the observation, and then computes the centred RMSE between that output point and the data point it corresponds to. No scaling is done to account for the difference in measurement variance on the data points.

    point_corr computes the pattern correlation coefficient between the point observations and each GCM

    point_std computes the standard deviation coefficient of the point observations and of each GCM (note first element = observation SD, for use in Taylor Diagram)

    """

    # Setup Nx3 confidence interval array for gp_rmse score (NaN for other cases)
    # Column 1 is lower limit
    # Column 2 is mean of lower and upper limit
    # Column 3 is upper limit
    # N row for GCMs
    
    n_GCM    = gcm_out.shape[0]


    #
    ### Gridded scores ###
    #
    if method == 'gp_full_cov':
        mu, Cov = model.predict_noiseless(X_pred, full_cov=True)
        Chol = np.linalg.cholesky(Cov)
        scores = dlogmvnorm(gcm_out.T, mu, Chol)
        orderings = np.argsort(-scores) #minus sign so that max is first
        relative = np.round(scores-np.max(scores),1)
    elif method == 'gp_var':
        mu, Var = model.predict_noiseless(X_pred, full_cov=False)
        Chol = np.linalg.cholesky(np.diagflat(Var))
        scores = dlogmvnorm(gcm_out.T, mu, Chol)
        orderings = np.argsort(-scores) #minus sign so that max is first
        relative = np.round(scores-np.max(scores),1)
    elif method == 'gp_RMSE':
        mu, Cov = model.predict_noiseless(X_pred, full_cov=True)
        scores = np.sqrt(np.mean((mu-gcm_out.T)**2, axis=0))
        orderings = np.argsort(scores) 
        relative = np.round(scores-np.min(scores),1)
    elif method == 'gp_RMSE_cent':
        mu, Var = model.predict_noiseless(X_pred, full_cov=False)
        scores = np.sqrt(np.mean((mu-gcm_out.T)**2, axis=0) - (np.mean(mu,axis=0) - np.mean(gcm_out.T,axis=0))**2)
        orderings = np.argsort(scores) 
        relative = np.round(scores-np.min(scores),1)
    #
    ### for gridded Taylor diagram ###
    #
    elif method == 'gp_corr':
        mu, Var = model.predict_noiseless(X_pred, full_cov=False)
	z=np.concatenate((mu.T,gcm_out))
        scores = np.corrcoef(z)[0,1:]
        orderings = np.argsort(-scores) #minus sign so that max is first
        relative = np.round(scores-np.max(scores),1)
    elif method == 'gp_std':
        mu, Var = model.predict_noiseless(X_pred, full_cov=False)
	z=np.concatenate((mu.T,gcm_out))
        scores = np.std(z,axis=1)
        orderings = np.argsort(scores) #not really relevant to SD, needed for return
        relative = np.round(scores/scores[0],1)
    #
    ### Point-based scores ###
    #
    elif method == 'RMSE':
        k = Exponentialhaversine(2)
        index = np.argmin(k._unscaled_dist(X_obs, X_pred), axis=1)
        scores = np.zeros(gcm_out.shape[0])
        for ii in range(gcm_out.shape[0]):
            y_gcm_grid = gcm_out.T[index,ii]
            scores[ii] = np.sqrt(np.mean((y_gcm_grid-y_obs)**2))
        orderings = np.argsort(scores) 
        relative = np.round(scores-np.min(scores),1)
    elif method == 'RMSE_cent':
        k = Exponentialhaversine(2)
        index = np.argmin(k._unscaled_dist(X_obs, X_pred), axis=1)
        scores = np.zeros(gcm_out.shape[0])
	for ii in range(gcm_out.shape[0]):
            y_gcm_grid = gcm_out.T[index,ii]
            scores[ii] = np.sqrt(np.mean((y_gcm_grid-y_obs)**2)-(np.mean(y_gcm_grid) - np.mean(y_obs))**2)
        orderings = np.argsort(scores) 
        relative = np.round(scores-np.min(scores),1)
    elif method == 'RMSE_weighted':
        pass
    #
    ### for point-based Taylor diagram ###
    #
    elif method == 'point_corr':
        k = Exponentialhaversine(2)
        index = np.argmin(k._unscaled_dist(X_obs, X_pred), axis=1)
        scores = np.zeros(gcm_out.shape[0])
        for ii in range(gcm_out.shape[0]):
            y_gcm_grid = gcm_out.T[index,ii]
            scores[ii] = np.corrcoef(y_obs.T,y_gcm_grid)[0,1]
        orderings = np.argsort(-scores) #minus sign so that max is first
        relative = np.round(scores-np.max(scores),1)
    elif method == 'point_std':
        k = Exponentialhaversine(2)
        index = np.argmin(k._unscaled_dist(X_obs, X_pred), axis=1)
        scores = np.zeros(gcm_out.shape[0]+1)
	scores[0]=np.std(y_obs)
        for ii in range(gcm_out.shape[0]):
            y_gcm_grid = gcm_out.T[index,ii]
            scores[ii+1] = np.std(y_gcm_grid)
        orderings = np.argsort(scores) #not really relevant to SD, needed for return
        relative = np.round(scores/scores[0],1)

    return(scores, orderings, relative)
