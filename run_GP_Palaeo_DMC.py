# v1.0 September 2019 F.J. Bragg Initial version
# v1.1 October   2019 F.J. Bragg Write GP out to netcdf
##                               Check existence of paths and files
##                               Add option to generate GP from observations and stop if no GCM files found
# v1.2 November  2019 S. Steinig Modified to estimate Eocene temperatures from DeepMIP compilation (see Inglis et al., 2020)

# Calling:
# python Palaeo_DMC.py <GCM_set> <Obs_set> <thinby>

# Directory structure within module_dir

# ./core - containing PEN DMC python code

# ./Model_Data/<GCM_set> - containing netcdf files of GMC data (*.nc) and mask data (mask.nc)
## dimensions 'lon', 'lat'; data variable 'var', mask variable 'mask'
## if only mask.nc supplied, code will generate GP from observations using mask.nc dimensions, then exit

# ./Observation_Data/<Obs_set> - containing ascii files (*.txt)
## 4 columns of observations with header: x, y, data, SD
## header row is ignored and provided solely for the convenience of the user

# Output directory will be created: ./Output/O-<Obs_set>_M-<GCM_set>_Th-<thinby>

#thinby=2 # this controls the degree of thinning which may be needed if memory problems occur.
## thinby=2 takes every second value
## thinby=1 uses the original GCM grid.

# Let's begin by loading the libraries we'll need.

#########################
### BASE PYTHON LIBRARIES
#########################
import sys
import os
import copy
from datetime import datetime
import re

import os
REF_FRAME = os.environ["REF_FRAME"]

startTime = datetime.now()

##################################################################
#### STANDARD PYTHON LIBRARIES THAT NEED TO BE INSTALLED FROM PyPI
##################################################################
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
# GPy is the Sheffield Gaussian Process library
import GPy

import glob as glob
#from cmocean import cm as cmo

##########################################
### LIBRARIES FROM PALAEO_DMC DISTRIBUTION
##########################################
# Add Palaeo_DMC code to path
module_dir = os.getcwd()
print("current directory is : " + module_dir)
sys.path.append(module_dir+"/core/")

from HaversineDist import Exponentialhaversine

####################
### INITIAL SETTINGS
####################

# missing value for output netcdf file
missing = -99999.0

########################
# SPECIFY DATA LOCATIONS
########################

GCM_set=str(sys.argv[1])
Obs_set=str(sys.argv[2])
thinby=int(sys.argv[3])

print('GCM_set = '+GCM_set)
GCM_dir = module_dir+'/Model_Data/'+GCM_set+'/'
if not os.path.exists(GCM_dir):
   sys.exit('STOPPING: GCM data directory '+GCM_dir+' not found')

print('Obs_set = '+Obs_set)
Obs_dir = module_dir+'/Observation_Data/'+Obs_set+'/'
if not os.path.exists(Obs_dir):
   sys.exit('STOPPING: Observations data directory '+Obs_dir+' not found')

print('Thinby = '+str(thinby))

Out_dir = module_dir+'/Output/O-'+Obs_set+'_M-'+GCM_set+'_Th-'+str(thinby)+'/'
print('Output will be sent to directory '+Out_dir)
if not os.path.exists(Out_dir):
   os.makedirs(Out_dir)

####################
# READ MASK
####################

# Need mask for grid dimensions and defintion of comparison region and continental outline
# If mask is supplied without any GCM data files, the script will stop after fitting the GP

# nc_m = GCM_dir+'mask_36x18.nc'  
# nc_m = GCM_dir+'mask_72x36.nc'  
nc_m = GCM_dir+'mask_144x72.nc'  
# nc_m = GCM_dir+'mask_144x72.nc'  
#nc_m = GCM_dir+'mask.nc'  

if not os.path.isfile(nc_m):
   sys.exit('STOPPING: Mask file '+nc_m+' not found')

nc_mid = Dataset(nc_m, 'r')  
lon_gcm = nc_mid.variables['lon'][:]
lat_gcm = nc_mid.variables['lat'][:]
mask_gcm = nc_mid.variables['mask'][:]

# Gridded lon/lat arrays for plotting, predicting GP.
lon_gcm_grid, lat_gcm_grid = np.meshgrid(lon_gcm, lat_gcm)

# only points in mask area are used for comparison
mask_gcm_ind=np.where(mask_gcm.flatten())[0]
mask_gcm_ind_nan=np.where(mask_gcm.flatten() == 0)[0]
n_grid = mask_gcm_ind.size
nc_mid.close()

########################
# LIST OBSERVATION FILES
########################

# We'll extract the observations into their coordinates (lon and lat) which I'll refer to as X, 
# the measurement (Y), and our estimate of the ratio of the variance of the measurement error at each point 

obs_files = glob.glob(Obs_dir+'DeepMIP_'+REF_FRAME+'*.txt')

if len(obs_files) == 0:
   sys.exit("STOPPING: No observation files: *.txt found in "+Obs_dir)
   
obs_files.sort()
obs_label = copy.deepcopy(obs_files)
n_obs=len(obs_files)

################
# LIST GCM FILES
################

gcm_files = glob.glob(GCM_dir+'*.nc')
gcm_files.remove(GCM_dir+'mask.nc')

# dmc is switch to indicate data-model comparison is to be performed
# if no GCM files are found, GP will be created from observations and save to netcdf
dmc = 1
if len(gcm_files) == 0:
   dmc = 0
   print("No GCM files *.nc found: Continuing with GP inference only")
   print("Data-model comparison will not be performed")
    
###################
# READ OBSERVATIONS
###################

for count_o in range(n_obs):
    obs_file = obs_files[count_o].split(Obs_dir)[-1]
    obs_label[count_o]=obs_file.split(".txt")[0]
    print('Reading observations from '+obs_file)
    observations = np.genfromtxt(Obs_dir+obs_file, skip_header=1)
    X_obs = observations[:,0:2]
    y_obs = observations[:,2].reshape(-1,1)
    # center anomalies around mean anomaly -> calculate residuals
    y_obs_simple_mean  = np.mean(y_obs)
    rad       = 4.0*np.arctan(1.0)/180.0
    clat      = np.cos(observations[:,1]*rad)
    y_obs_weighted_mean  = np.average(y_obs.flatten(),weights=clat)
    y_obs       = y_obs-y_obs_weighted_mean
    var_ratios = observations[:,3].reshape(-1,1)**2
    obs_SD=np.column_stack((observations[:,0:2],observations[:,3]))

    #####################################
    # ## FITTING A GAUSSIAN PROCESS MODEL
    #####################################
     
    k3 = Exponentialhaversine(2, lengthscale=6000)
    m3 = GPy.models.GPHeteroscedasticRegression(X=X_obs, Y=y_obs, kernel=k3)
    m3['.*het_Gauss.variance'] = var_ratios #Set the noise parameters to the error in Y
    m3.het_Gauss.variance.fix() #We can fix the noise term, since we already know it
    m3.Exponential.lengthscale.constrain_bounded(2000,10000) 

    print("start optimize")
    m3.optimize()

    print("start predict")
    X_plot=np.column_stack((lon_gcm_grid.flatten(), lat_gcm_grid.flatten())) # specifies the prediction locations
    mu3_v_mean,V3_v_mean = m3.predict_noiseless(X_plot)
    mu3_v,V3_v           = m3.predict_noiseless(X_plot, full_cov=True)
    # sample the posterior for GMST uncertainty
    print("start sampling")
    N                = 10000 # number used for publication results (slow)
    Ysample          = np.random.multivariate_normal(mu3_v.flatten(), V3_v, N)
    Ysample_mean     = np.mean(Ysample, axis=0)
    Ysample_std      = np.std(Ysample, axis=0)
    Ysample_gm_std   = np.std(np.mean(Ysample, axis=1), axis=0)

    print(m3) 


    print("start mask")
    # Mask data
    Ysample_masked             = Ysample    
    Ysample_masked = ma.masked_array(Ysample_masked)
    Ysample_masked[:,mask_gcm.flatten() == observations[0,4].astype(int)] = ma.masked
    Ysample_masked_reshape    = np.reshape(Ysample_masked,(N,len(lat_gcm),len(lon_gcm)))
    Ysample_masked_reshape_zm = np.nanmean(Ysample_masked_reshape, axis=2)
    
    print("start writing to disk")
    # Save the GP mean and SD to netcdf
    nc_gp = Out_dir+'GP_predict_'+obs_label[count_o]+'.nc'
    nc_gpid = Dataset(nc_gp, 'w', clobber=True)
    
    # Define dimensions
    dim_lon = nc_gpid.createDimension("lon", len(lon_gcm))
    dim_lat = nc_gpid.createDimension("lat", len(lat_gcm))
    dim_N   = nc_gpid.createDimension("N", N)

    # Add dimension variables and attributes
    var_lon = nc_gpid.createVariable("lon","f4",("lon",))
    var_lon[:] = lon_gcm
    var_lon.units = "degrees east"
    
    var_lat = nc_gpid.createVariable("lat","f4",("lat",))
    var_lat[:] = lat_gcm
    var_lat.units = "degrees north"
    
    # Add data variables and values
    # The basic GP is a global field
    var_GP_mean    = nc_gpid.createVariable("GP_mean","f4",("lat","lon"))
    var_GP_mean[:] = np.reshape(mu3_v_mean,(len(lat_gcm),len(lon_gcm)))

    var_GP_SD      = nc_gpid.createVariable("GP_SD","f4",("lat","lon"))
    var_GP_SD[:]   = np.reshape(np.sqrt(V3_v_mean),(len(lat_gcm),len(lon_gcm)))

    var_GP_Ysample    = nc_gpid.createVariable("Ysample_masked_reshape","f4",("N","lat","lon"))
    var_GP_Ysample[:] = Ysample_masked_reshape

    var_GP_Ysample_zm    = nc_gpid.createVariable("Ysample_masked_reshape_zm","f4",("N","lat"))
    var_GP_Ysample_zm[:] = Ysample_masked_reshape_zm

    var_y_obs_weighted_mean    = nc_gpid.createVariable("y_obs_weighted_mean","f4")
    var_y_obs_weighted_mean[:] = y_obs_weighted_mean

    nc_gpid.close()
    print('GP data written to netcdf '+nc_gp+'\n')


