"""
Code to test the implementation of the Ensemble Adjustment Kalman Filter (Anderson, 2003). 
This will assimimlate satellite observations of sea ice concentration in both the Arctic and Antarctic,
and then write out the increments (corrections) that are applied to each sub-grid sea ice concentration category

Author: Will Gregory
"""

import os
import glob
import numpy as np
import xarray as xr
import pickle
from global_land_mask import globe
from datetime import datetime,timedelta
from sklearn.metrics.pairwise import haversine_distances

print('Starting DA script:',datetime.now(),flush=True)

def preprocess(lon_mod,lat_mod,lon_obs,lat_obs,localization_radius=0.06,save='./grid_locs.pkl'):
    """
    pad the domain of the prior state variables and observations
    to appropriately handle the chosen localization distance.
    Then compute the localization matrix for this padded domain.

    Return the padded state variables, observations and
    localization matrix, and the indices of the original (unpadded) domain 
    """
    xT,yT = lon_obs.shape
    grid_info = {}
    
    mod_pos = np.deg2rad([lat_mod.ravel(),lon_mod.ravel()]).T
    m = 0
    for ix in range(xT):
        for jx in range(yT):
            if globe.is_ocean(lat_obs[ix,jx],lon_obs[ix,jx]):
                obs_pos = np.deg2rad([lat_obs[ix,jx],lon_obs[ix,jx]])[None]
                D = np.squeeze(haversine_distances(obs_pos,mod_pos))
                mod_loc = np.where(D==np.min(D))[0][0]
                halo = np.squeeze(np.where(D<=localization_radius))
                D = D[halo]
                W = np.copy(D)
                ratio = D/(localization_radius/2)
                cond1 = ratio<=1
                cond2 = ((ratio>1) & (ratio<=2))
                W[cond1] = -(ratio[cond1]**5) / 4 + ratio[cond1]**4 / 2 + 5 * ratio[cond1]**3 / 8 - 5 * ratio[cond1]**2 / 3 + 1
                W[cond2] = ratio[cond2]**5 / 12 - ratio[cond2]**4 / 2 + 5 * ratio[cond2]**3 / 8 + 5 * ratio[cond2]**2 / 3 - 5 * ratio[cond2] + 4 - 2 / 3 / ratio[cond2]
                W[D==0] = 1
                W[D>localization_radius] = 0

                grid_info['cell'+str(m)+'_locMatrix'] = W
                grid_info['cell'+str(m)+'_halo'] = halo
                grid_info['cell'+str(m)+'_trim'] = mod_loc
            
            m += 1
    
    with open(save, 'wb') as f:
        pickle.dump(grid_info, f)

def adjust_negatives(X):
    negatives = X<0
    dists = np.abs(X*negatives).sum(axis=0)
    positives = np.sum(X > 0, axis=0)
    positives[positives == 0] = 1
    subtract_values = dists/positives
    subtract_values_expanded = np.expand_dims(subtract_values, axis=0)
    X_adjusted = X - (X > 0) * subtract_values_expanded
    X_adjusted[negatives] = 0
    return X_adjusted
    
def postprocess(x):
    """
    post-processing to ensure updated sea ice concentration
    is bounded between 0 and 1.
    """
    x = x.transpose(1,0,2,3)
    while np.any(x<0):
        x = adjust_negatives(x)
    SIC = np.nansum(x,0)
    high = SIC>1
    ratio = 1/SIC[high]
    for CAT in range(x.shape[0]):
        x[CAT,high] = x[CAT,high]*ratio
    return x.transpose(1,0,2,3)

date = str(input('Please provide a date to assimilate observations in the format %Y%m%d:\n'))

### FILES ###
obs_file_NH = 'NSIDC0051_SEAICE_PS_N25km_'+date+'_v2.0.nc'
obs_file_SH = 'NSIDC0051_SEAICE_PS_S25km_'+date+'_v2.0.nc'
hemispheres = [obs_file_NH,obs_file_SH]
hem_labels = ['N','S']
grid = xr.open_dataset('ice.static.nc') #model grid
lon = grid.GEOLON.to_numpy()
lat = grid.GEOLAT.to_numpy()
ice_restarts = sorted(glob.glob('ice_model.res*')) #prior model states (RESTART files)
prior = xr.open_mfdataset(ice_restarts,concat_dim='ens',combine='nested',decode_times=False).part_size.to_numpy()[:,0,1:] #first index of the SIS2 model is open-water fraction
prior[np.isnan(prior)] = 0
nmembers,nCat,xT,yT = prior.shape #nCat is the number of sub-grid categories (each of which will be updated via DA)
prior_original = np.copy(prior)
prior = prior.reshape(nmembers,nCat,xT*yT)

### PARAMETERS ###
localization_radius = 0.06 #radians
obs_var = 0.01 #variance of observed SIC

results = []
for h,hem in enumerate(hemispheres):
    prior_temp = np.copy(prior)
    if os.path.exists(hem):
        loc_fp = 'gridinfo_PS_'+hem_labels[h]+'25km_for_OM4grid_locrad'+str(localization_radius)+'.pkl'
        obs = xr.open_dataset(hem)
        obs_lon = obs.lon.to_numpy()
        obs_lat = obs.lat.to_numpy()
        obs_sic = obs.sic.to_numpy().ravel()

        if os.path.exists(loc_fp):
            with open(loc_fp, 'rb') as f:
                grid_info = pickle.load(f)
        else:
            preprocess(lon,lat,obs_lon,obs_lat,localization_radius,loc_fp)
            with open(loc_fp, 'rb') as f:
                grid_info = pickle.load(f)

        obs_lonr = obs_lon.ravel()
        obs_latr = obs_lat.ravel()
        
        for m in range(len(obs_sic)):
            if (globe.is_ocean(obs_latr[m],obs_lonr[m])) & (~np.isnan(obs_sic[m])):
                priorH = np.sum(prior[:,:,grid_info['cell'+str(m)+'_trim']],1)
                priorH_mean = np.mean(priorH)
                priorH_var = np.var(priorH,ddof=1)
                if priorH_var != 0:
                    var_ratio = obs_var / (priorH_var + obs_var)
                    new_mean = var_ratio * (priorH_mean  + priorH_var*obs_sic[m] / obs_var)
                    a = np.sqrt(var_ratio)
                    obs_inc = a * (priorH - priorH_mean) + new_mean - priorH
                    for k in range(nCat):
                        if prior[:,k,grid_info['cell'+str(m)+'_halo']].shape[1] == 0:
                            prior_cov = np.cov(prior[:,k,grid_info['cell'+str(m)+'_trim']],priorH)[0,1]
                            ens_increment = (prior_cov/priorH_var) * obs_inc
                            prior[:,k,grid_info['cell'+str(m)+'_trim']] = prior[:,k,grid_info['cell'+str(m)+'_trim']] + ens_increment
                        else:
                            prior_cov = grid_info['cell'+str(m)+'_locMatrix'] * np.cov(prior[:,k,grid_info['cell'+str(m)+'_halo']].T,priorH)[:-1,-1]
                            ens_increment = (prior_cov/priorH_var)[:,np.newaxis] * obs_inc[np.newaxis,:]
                            prior[:,k,grid_info['cell'+str(m)+'_halo']] = prior[:,k,grid_info['cell'+str(m)+'_halo']] + ens_increment.T
    results.append(prior - prior_temp)
    
increments = np.nansum(results,0).reshape(nmembers,1,nCat,xT,yT)
posterior = np.zeros((nmembers,1,nCat+1,xT,yT))
posterior[:,0,1:] = postprocess(prior_original + increments[:,0])
posterior[:,0,0] = 1 - np.nansum(posterior[:,0,1:],1)
posterior[posterior<0] = 0
posterior[posterior>1] = 1

### SAVE INCREMENTS ###
ds = xr.Dataset(data_vars=dict(part_size=(['members','time', 'ct', 'yT', 'xT'], increments)), coords=dict(yT=grid['yT'], xT=grid['xT']))
ds.part_size.attrs['long_name'] = 'category_sea_ice_concentration_increments'
ds.part_size.attrs['units'] = 'area_fraction'
ds['time'] = [date]
ds.mean('members').to_netcdf(date+'.EAKF_increment.ens_mean.nc')
ds.to_netcdf(date+'.EAKF_increment.nc')

print('Finished DA:',datetime.now(),flush=True)
