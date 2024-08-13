import os
import glob
import numpy as np
import xarray as xr
import pickle
from global_land_mask import globe
from datetime import datetime,timedelta
from sklearn.metrics.pairwise import haversine_distances

print('Starting DA script:',datetime.now(),flush=True)

def enthalpy_ice(zTin,zSin):
    """
    compute enthalpy of ice based on liquidus temperature
    and salinity of mushy ice. Used to create a new 'sea ice
    profile' in the case the EnKF adds ice to grid cell which 
    was previously ice-free.
    """
    cp_wtr  = 4200
    cp_ice  = 2100
    Lfresh  = 3.34e5
    MIU = 0.054
    Tm = -MIU*zSin
    return cp_wtr*zTin + cp_ice*(zTin - Tm) + (cp_wtr - cp_ice)*Tm*np.log(zTin/Tm) + Lfresh*(Tm/zTin-1)

def liquidus_temperature_mush(Sbr):
    """
    compute liquidus temp of ice based on salinity of mushy ice. 
    Used to create a new 'sea ice profile' in the case the EnKF adds
    ice to grid cell which was previously ice-free.
    """
    Sb_liq =  123.66702800276086
    c1      = 1.0
    c1000   = 1000
    az1_liq = -18.48
    bz1_liq =   0.0
    az2_liq = -10.3085
    bz2_liq =  62.4
    az1p_liq = az1_liq / c1000
    bz1p_liq = bz1_liq / c1000
    az2p_liq = az2_liq / c1000
    bz2p_liq = bz2_liq / c1000
    M1_liq = az1_liq
    N1_liq = -az1p_liq
    O1_liq = -bz1_liq / az1_liq
    M2_liq = az2_liq
    N2_liq = -az2p_liq
    O2_liq = -bz2_liq / az2_liq

    if Sbr<=Sb_liq:
        t_high = 1.0
    else:
        t_high = 0.0

    return ((Sbr / (M1_liq + N1_liq * Sbr)) + O1_liq) * t_high + ((Sbr / (M2_liq + N2_liq * Sbr)) + O2_liq) * (1.0 - t_high)

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

def postprocess(x):
    """
    post-processing to ensure updated sea ice concentration
    is bounded between 0 and 1.
    """
    x = x.transpose(1,0,2,3)
    x[x<0] = 0    
    SIC = np.nansum(x,0)
    high = SIC>1
    ratio = 1/SIC[high]
    for CAT in range(x.shape[0]):
        x[CAT,high] = x[CAT,high]*ratio
    return x.transpose(1,0,2,3)

### PATHS ###
user = os.popen('whoami').read().split('\n')[0]
experiment = os.getcwd().split('/')[-2].split('.')[0]
savepath = '/gpfs/f5/gfdl_o/scratch/'+user+'/ENKF/increments/'+experiment+'/'
if os.path.exists(savepath)==False:
    os.makedirs(savepath)

### CURRENT MODEL TIME ###
y,m,d = np.genfromtxt('coupler.res',skip_header=1)[1,:3].astype(np.int32)
date = datetime(y,m,d).strftime('%Y%m%d')

### FILES ###
obs_file_NH = '/gpfs/f5/gfdl_o/scratch/William.Gregory/NTSIC/raw_data/NH/NSIDC0051_SEAICE_PS_N25km_'+date+'_v2.0.nc'
obs_file_SH = '/gpfs/f5/gfdl_o/scratch/William.Gregory/NTSIC/raw_data/SH/NSIDC0051_SEAICE_PS_S25km_'+date+'_v2.0.nc'
hemispheres = [obs_file_NH,obs_file_SH]
hem_labels = ['N','S']
grid = xr.open_dataset('/ncrc/home2/William.Gregory/dart_manhattan/ice.static.nc')
lon = grid.GEOLON.to_numpy()
lat = grid.GEOLAT.to_numpy()
ice_restarts = sorted(glob.glob('ice_model.res*')) #prior model states (RESTART files)
ocn_restarts = sorted(glob.glob('MOM.res.*')) #get ocean states for salinity-dependent freezing point
prior = xr.open_mfdataset(ice_restarts,concat_dim='ens',combine='nested',decode_times=False).part_size.to_numpy()[:,0,1:]
prior[np.isnan(prior)] = 0
nmembers,nCat,xT,yT = prior.shape
prior_original = np.copy(prior)
prior = prior.reshape(nmembers,nCat,xT*yT)
### PARAMETERS ###
localization_radius = 0.06 #radians
obs_var = 0.01 #variance of observed SIC

results = []
for h,hem in enumerate(hemispheres):
    prior_temp = np.copy(prior)
    if os.path.exists(hem):
        loc_fp = '/gpfs/f5/gfdl_o/scratch/William.Gregory/ENKF/tiles/gridinfo_PS_'+hem_labels[h]+'25km_for_OM4grid_locrad'+str(localization_radius)+'.pkl'
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
            if globe.is_ocean(obs_latr[m],obs_lonr[m]):
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
ds.mean('members').to_netcdf(savepath+date+'.EAKF_increment.ens_mean.nc')
ds.to_netcdf(savepath+date+'.EAKF_increment.nc')

### UPDATE SEA ICE RESTARTS & CREATE NEW SEA ICE PROFILES ###
rho_ice = 905.
rho_snow = 330.
phi_init = 0.75 #initial liquid fraction of frazil ice
Si_new = 5 #salinity of mushy ice
Ti = min(liquidus_temperature_mush(Si_new/phi_init),-0.1)
qi_new = enthalpy_ice(Ti, Si_new)
hlim = [1.0e-10, 0.1, 0.3, 0.7, 1.1, 1.5]
hmid = np.array([0.5*(hlim[n]+hlim[n+1]) for n in range(nCat)])
i_thick = np.tile((hmid*rho_ice)[None,:,None,None],(1,1,xT,yT))
for member,file in enumerate(ice_restarts):
    fr = xr.open_dataset(file)
    SSS = xr.open_dataset(ocn_restarts[member],decode_times=False).Salt.to_numpy()[:,0]
    prior_m = fr.part_size.to_numpy()
    post_m = posterior[member]

    cond1 = np.where((prior_m[:,1:]<=0) & (post_m[:,1:]>0)) #where original state was ice-free, but EnKF has added ice
    cond2 = np.where((prior_m[:,1:]>0) & (post_m[:,1:]<=0)) #where original state contained ice, but EnKF has made ice-free

    h_ice = fr.h_ice.to_numpy()
    h_ice[cond1] = i_thick[cond1]
    h_ice[cond2] = 0

    h_snow = fr.h_snow.to_numpy()
    h_snow[cond1] = 0
    h_snow[cond2] = 0

    enth_ice = fr.enth_ice.to_numpy()
    for layer in range(4):
        enth_ice[:,layer][cond1] = qi_new
        enth_ice[:,layer][cond2] = 0

    enth_snow = fr.enth_snow.to_numpy()
    enth_snow[0][cond1] = 0
    enth_snow[0][cond2] = 0

    T_skin = fr.T_skin.to_numpy()
    T_skin[cond1] = Ti
    T_skin[cond2] = np.maximum(-0.0539*(np.tile(SSS[:,None],(1,nCat,1,1))[cond2]),-2.)

    sal_ice = fr.sal_ice.to_numpy()
    for layer in range(4):
        sal_ice[:,layer][cond1] = Si_new
        sal_ice[:,layer][cond2] = 0

    h_pond = fr.h_pond.to_numpy()
    h_pond[cond1] = 0
    h_pond[cond2] = 0

    fr.part_size.loc[:] = post_m
    fr.h_ice.loc[:] = h_ice
    fr.h_snow.loc[:] = h_snow
    fr.h_pond.loc[:] = h_pond
    fr.enth_ice.loc[:] = enth_ice
    fr.enth_snow.loc[:] = enth_snow
    fr.T_skin.loc[:] = T_skin
    fr.sal_ice.loc[:] = sal_ice

    fr.to_netcdf(file,mode='a')
print('Finished DA:',datetime.now(),flush=True)
