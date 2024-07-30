import os
import pickle
import glob
import itertools
import numpy as np
import xarray as xr
from mpi4py import MPI
from datetime import datetime,timedelta
from sklearn.metrics.pairwise import haversine_distances

COMM = MPI.COMM_WORLD
if COMM.rank == 0:
    print('Number of cores used:',COMM.size)
    print('Starting DA script:',datetime.now(),flush=True)

def split(container, count):
    """
    function for dividing the number of tasks (container)
    across the number of available compute nodes (count)
    """
    return [container[_i::count] for _i in range(count)]

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

def preprocess(lon,lat,localization_radius,save):
    """
    pad the domain of the prior state variables and observations
    to appropriately handle the chosen localization distance.
    Then compute the localization matrix for this padded domain.

    Return the padded state variables, observations and
    localization matrix, and the indices of the original (unpadded) domain 
    """
    xT,yT = lon.shape
    xdiv = xT//20
    ydiv = yT//20

    X = np.deg2rad(np.array([lat.ravel(),lon.ravel()]).T)
    tiling = {}
    count = 0
    for ix in range(0,xT,xdiv):
        for jx in range(0,yT,ydiv):
            latr = lat[ix:ix+xdiv,jx:jx+ydiv].ravel()
            lonr = lon[ix:ix+xdiv,jx:jx+ydiv].ravel()
            Y = np.deg2rad(np.array([latr,lonr]).T)
            c = haversine_distances(X,Y)
            halo = np.unique(np.where(c<=localization_radius)[0])
            Z = np.deg2rad(np.array([lat.ravel()[halo],lon.ravel()[halo]]).T)
            trim = np.array([np.squeeze(np.where((Z == y).all(axis=1))) for y in Y])
            D = haversine_distances(Z,Z)
            W = np.copy(D)
            ratio = D/(localization_radius/2)
            cond1 = ratio<=1
            cond2 = ((ratio>1) & (ratio<=2))
            W[cond1] = -(ratio[cond1]**5) / 4 + ratio[cond1]**4 / 2 + 5 * ratio[cond1]**3 / 8 - 5 * ratio[cond1]**2 / 3 + 1
            W[cond2] = ratio[cond2]**5 / 12 - ratio[cond2]**4 / 2 + 5 * ratio[cond2]**3 / 8 + 5 * ratio[cond2]**2 / 3 - 5 * ratio[cond2] + 4 - 2 / 3 / ratio[cond2]
            W[D==0] = 1
            W[D>localization_radius] = 0
                        
            tiling['tile'+str(count)+'_locMatrix'] = W
            tiling['tile'+str(count)+'_halo'] = halo
            tiling['tile'+str(count)+'_trim'] = trim
            count += 1
    with open(save, 'wb') as f:
        pickle.dump(tiling, f)

def Kfilter(prior,obs,W,trim,reshape_dims,obs_error=0.01):
    """
    ensemble adjustment Kalman filter
    prior state variables of size: E, C, N
    E is the number of ensemble members
    C is the number of model states (categories)
    N is the number of grid points

    returns: posterior state variables and increments
    """

    E,C = prior.shape[0], prior.shape[1]
    dX,dY = reshape_dims

    if ((prior==0).all() & (obs==0).all()):
        tmp = np.zeros((E,C,dX,dY))
        return tmp,tmp
    elif np.isnan(obs).all():
        return prior[:,:,trim].reshape(E,C,dX,dY),np.zeros((E,C,dX,dY))
    else:
        valid_obs = np.atleast_1d(np.squeeze(np.where(~np.isnan(obs))))
        
        prior_mean = np.mean(prior,0)
        prior_anom = prior-prior_mean
        priorH = np.sum(prior,1)[:,valid_obs]
        priorH_mean = np.mean(priorH,0)
        priorH_anom = priorH - priorH_mean
	    
        innov = obs[valid_obs] - priorH_mean #bias of ensemble mean
        innovE = obs[valid_obs] - priorH #bias of each ensemble member
        N = priorH.shape[1]
        
        Bm = W[:,valid_obs] * (np.einsum('ijk,il->jkl', prior_anom, priorH_anom) / (E - 1))
        Bo = W[np.ix_(valid_obs,valid_obs)] * np.cov(priorH.T)
        Bo_i = np.linalg.inv(Bo + np.eye(N)*obs_error)
	
        increments = np.zeros((E,C,prior.shape[2]))
        for k in range(C):
            K = Bm[k] @ Bo_i #Kalman gain of ice category k
            posterior_mean = prior_mean[k] + np.dot(K, innov.T)
            increments[:,k] = np.dot(K, innovE.T).T + posterior_mean - prior[:,k]
	    
        posterior = (prior + increments)[:,:,trim]
        
        return postprocess(posterior.reshape(E,C,dX,dY)), increments[:,:,trim].reshape(E,C,dX,dY)

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
if COMM.rank == 0:
    experiment = os.getcwd().split('/')[-2].split('.')[0]
    savepath = '/gpfs/f5/gfdl_o/scratch/'+user+'/ENKF/increments/'+experiment+'/'
    if os.path.exists(savepath)==False:
        os.makedirs(savepath)

### CURRENT MODEL TIME ###
y,m,d = np.genfromtxt('coupler.res',skip_header=1)[1,:3].astype(np.int32)
date = datetime(y,m,d).strftime('%Y%m%d')

### FILES ###
obs_file = '/gpfs/f5/gfdl_o/scratch/William.Gregory/NTSIC/regrid_OM4/NSIDC0051_SEAICE_SPEAR1deg_'+date+'_v2.0.nc'
if os.path.exists(obs_file):
    grid = xr.open_dataset('/ncrc/home2/William.Gregory/dart_manhattan/ice.static.nc')
    lon = grid.GEOLON.to_numpy()
    lat = grid.GEOLAT.to_numpy()
    ice_restarts = sorted(glob.glob('ice_model.res*')) #prior model states (RESTART files)
    ocn_restarts = sorted(glob.glob('MOM.res.*')) #get ocean states for salinity-dependent freezing point
    prior = xr.open_mfdataset(ice_restarts,concat_dim='ens',combine='nested',decode_times=False).part_size.to_numpy()[:,0,1:]
    prior[np.isnan(prior)] = 0
    obs = xr.open_dataset(obs_file).sic.to_numpy()
    nmembers,nCat,xT,yT = prior.shape

    ### PARAMETERS ###
    localization_radius = 0.06 #radians
    tile_fp = '/gpfs/f5/gfdl_o/scratch/William.Gregory/ENKF/tiles/Tiling_OM4grid_locrad'+str(localization_radius)+'.pkl'
    obs_error = 0.01 #variance of observed SIC
    xdiv = xT//20
    ydiv = yT//20

    prior = prior.reshape(nmembers,nCat,xT*yT)
    obs = obs.ravel()
    if os.path.exists(tile_fp):
        with open(tile_fp, 'rb') as f:
            tiling = pickle.load(f)
    else:
        if COMM.rank == 0:
            preprocess(lon,lat,localization_radius,tile_fp)
        with open(tile_fp, 'rb') as f:
            tiling = pickle.load(f)

    Ntiles = np.arange(400)
    selected_variables = range(len(Ntiles))
    if COMM.rank == 0:
        splitted_jobs = split(selected_variables, COMM.size)
    else:
        splitted_jobs = None
    scattered_jobs = COMM.scatter(splitted_jobs, root=0)

    results = []
    for ix in scattered_jobs:
        tile = Ntiles[ix]
        prior_tile = prior[:,:,tiling['tile'+str(tile)+'_halo']]
        obs_tile = obs[tiling['tile'+str(tile)+'_halo']]
            
        outputs = Kfilter(prior_tile,obs_tile,tiling['tile'+str(tile)+'_locMatrix'],tiling['tile'+str(tile)+'_trim'],[xdiv,ydiv],obs_error=obs_error)
        results.append(outputs)

    results = COMM.gather(results, root=0)
    if COMM.rank == 0:
        posterior = np.zeros((nmembers,1,nCat+1,xT,yT))
        increments = np.zeros((nmembers,1,nCat,xT,yT))
        results = list(itertools.zip_longest(*results))
        results = [x for xs in results for x in xs]
        
        count = 0
        for ix in range(0,xT,xdiv):
            for jx in range(0,yT,ydiv):
                posterior[:,0,1:,ix:ix+xdiv,jx:jx+ydiv] = results[count][0]
                increments[:,0,:,ix:ix+xdiv,jx:jx+ydiv] = results[count][1]
                count += 1

        posterior[:,0,0] = 1 - np.nansum(posterior[:,0,1:],1)
        posterior[posterior<0] = 0
        posterior[posterior>1] = 1
        
        ### SAVE INCREMENTS ###
        ds = xr.Dataset(data_vars=dict(part_size=(['members','time', 'ct', 'yT', 'xT'], increments)), coords=dict(yT=grid['yT'], xT=grid['xT']))
        ds.part_size.attrs['long_name'] = 'category_sea_ice_concentration_increments'
        ds.part_size.attrs['units'] = 'area_fraction'
        ds['time'] = [date]
        ds.mean('members').to_netcdf(savepath+date+'.EnKF_increment.ens_mean.nc')
        ds.to_netcdf(savepath+date+'.EnKF_increment.nc')

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
else:
    if COMM.rank == 0:
        print('No observation file found for date',date,flush=True)
