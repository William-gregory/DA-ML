import os
import glob
import itertools
import numpy as np
import xarray as xr
from datetime import datetime,timedelta
from sklearn.metrics.pairwise import haversine_distances

def pp(x):
    """
    post-processing to ensure updated sea ice concentration
    is bounded between 0 and 1.
    """
    x = x.transpose(1,0,2,3)
    x[x<0] = 0    
    SIC = np.nansum(x,0)
    high = SIC>1
    ratio = 1/SIC[high]
    for CAT in range(5):
        x[CAT,high] = x[CAT,high]*ratio
    return x.transpose(1,0,2,3)

def enthalpy_ice(zTin,zSin):
    """
    compute enthalpy of ice based on liquidus temperature
    and salinity of mushy ice. Used to create a new 'sea ice
    profile' in the case the CNN adds ice to grid cell which 
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
    Used to create a new 'sea ice profile' in the case the CNN adds
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

def preprocess(prior,obs,lon1,lat1,lon2,lat2,localization_radius):
    """
    pad the domain of the prior state variables and observations
    to appropriately handle the chosen localization distance.
    Then compute the localization matrix for this padded domain.
    
    Return the padded state variables, observations and
    localization matrix, and the indices of the original (unpadded) domain 
    """
    E,C,dX,dY = prior.shape
    prior = prior.reshape(E,C,dX*dY)
    obs = obs.ravel()
    lon1 = lon1.ravel()
    lat1 = lat1.ravel()
    lon2 = lon2.ravel()
    lat2 = lat2.ravel()
    prior[np.isnan(prior)] = 0
    obs[np.isnan(obs)] = 0
    
    X = np.deg2rad(np.array([lat1,lon1]).T)
    Y = np.deg2rad(np.array([lat2,lon2]).T)
    c = haversine_distances(X,Y)                                                                                                                                                                                                                  
    pad_halo = np.unique(np.where(c<=localization_radius)[0])
    Z = np.deg2rad(np.array([lat1[pad_halo],lon1[pad_halo]]).T)
    trim_halo = np.array([np.squeeze(np.where((Z == y).all(axis=1))) for y in Y])
    W = haversine_distances(Z,Z)
    W[W<=localization_radius] = 1
    W[W!=1] = 0
    
    return prior[:,:,pad_halo], obs[pad_halo], W, trim_halo

def Kfilter(prior,obs,lon,lat,lon_sub,lat_sub,loc_rad=1,obs_error=0.1):
    """
    ensemble Kalman filter (still need to do 'adjustment' part)
    prior state variables of size: E, C, N
    E is the number of ensemble members
    C is the number of model states (categories)
    N is the number of grid points

    returns: posterior state variables and increments
    """
    prior, obs, W, trim_halo = preprocess(prior,obs,lon,lat,lon_sub,lat_sub,localization_radius=loc_rad)
    priorH = np.nansum(prior,1)
    E,C,N = prior.shape
    dX,dY = lon_sub.shape

    if ((priorH==0).all() & (obs==0).all()):
        tmp = np.zeros((E,C,dX,dY))
        return tmp,tmp
    else:
        Bm = np.einsum('ijk,il->jkl', prior-np.nanmean(prior,0), priorH-np.nanmean(priorH,0)) / (E - 1)
        Bo = np.cov(priorH.T) + np.eye(N)*obs_error
        K = W * (Bm @ np.linalg.inv(Bo))
        posterior = np.array([prior[x] + (K @ (obs - priorH[x])) for x in range(E)])[:,:,trim_halo]

        return pp(posterior.reshape(E,C,dX,dY)), (posterior-prior[:,:,trim_halo]).reshape(E,C,dX,dY)

### PATHS ###
experiment = os.getcwd().split('/')[-2].split('.')[0]
user = os.popen('whoami').read().split('\n')[0]
savepath = '/gpfs/f5/gfdl_o/scratch/'+user+'/EnKF_increments/'+experiment+'/'
if os.path.exists(savepath)==False:
    os.makedirs(savepath)

### CURRENT MODEL TIME ###
y,m,d = np.genfromtxt('coupler.res',skip_header=1)[1,:3].astype(np.int32)
date = datetime(y,m,d).strftime('%Y%m%d')

### FILES ###
grid = xr.open_dataset('../Data/ice.static.nc')
lon = grid.GEOLON.to_numpy()
lat = grid.GEOLAT.to_numpy()
ice_restarts = sorted(glob.glob('ice_model.res*')) #prior model states (RESTART files)
ocn_restart = sorted(glob.glob('MOM.res.*')) #get ocean states for salinity-dependent freezing point
fi = xr.open_mfdataset(ice_restarts,concat_dim='ens',combine='nested',decode_times=False).part_size.to_numpy()[:,0,1:]
fo = xr.open_mfdataset(ocn_restarts,concat_dim='ens',combine='nested',decode_times=False).s_surf.to_numpy()
obs = xr.open_dataset('/gpfs/f5/gfdl_o/proj-shared/obs/NTSIC/NSIDC0051_SEAICE_SPEAR1deg_'+date+'_v2.0.nc').sic.to_numpy()
nmembers,nCat,xT,yT = fi.shape

### PARAMETERS ###
rho_ice = 905.
rho_snow= 330.
phi_init = 0.75 #initial liquid fraction of frazil ice
Si_new = 5 #salinity of mushy ice
localization_radius = 0.03 #radians
observation_error = 0.1 #10% SIC observation error
Ti = min(liquidus_temperature_mush(Si_new/phi_init),-0.1)
qi_new = enthalpy_ice(Ti, Si_new)
hlim = [1.0e-10, 0.1, 0.3, 0.7, 1.1, 1.5]
hmid = np.array([0.5*(hlim[n]+hlim[n+1]) for n in range(nCat)])
i_thick = np.tile((hmid*rho_ice)[None,:,None,None],(1,1,xT,yT))
xdiv = xT//20
ydiv = yT//20

posterior = np.zeros((nmembers,1,nCat+1,xT,yT//2))
increments = np.zeros((nmembers,1,nCat,xT,yT//2))
hems = [np.arange(0,180),np.arange(180,360)]
for hem in range(2):
    lon_sub = lon[:,hems[hem]]
    lat_sub = lat[:,hems[hem]]
    for ix in range(0,xT,xdiv):
        for jx in range(0,yT//2,ydiv):
            outputs = Kfilter(prior[:,:,:,hems[hem]],obs[:,hems[hem]],\
                              lon_sub,lat_sub,lon_sub[ix:ix+xdiv,jx:jx+ydiv],\
                              lat_sub[ix:ix+xdiv,jx:jx+ydiv],loc_rad=localization_radius)
            posterior[hem,:,0,1:,ix:ix+xdiv,jx:jx+ydiv] = outputs[0]
            increments[hem,:,0,:,ix:ix+xdiv,jx:jx+ydiv] = outputs[1]
posterior = np.concatenate((posterior[0],posterior[1]),axis=4)
increments = np.concatenate((increments[1],increments[1]),axis=4)

ds = xr.Dataset(data_vars=dict(part_size=(['member','time', 'ct', 'yT', 'xT'], increments)), coords=dict(yT=f['yT'], xT=f['xT']))
ds.part_size.attrs['long_name'] = 'category_sea_ice_concentration_increments'
ds.part_size.attrs['units'] = 'area_fraction'
ds['time'] = [date]
ds.to_netcdf(savepath+date+'.EnKF_increment.nc')
posterior[:,0,0] = 1 - np.nansum(posterior[:,0,1:],1)
posterior[posterior<0] = 0
posterior[posterior>1] = 1

for member,file in enumerate(ice_restarts):
    fr = xr.open_dataset(file)
    prior = fr.part_size.to_numpy()
    post = posterior[member]

    cond1 = np.where((prior[:,1:]<=0) & (post[:,1:]>0)) #where original state was ice-free, but EnKF has added ice
    cond2 = np.where((prior[:,1:]>0) & (post[:,1:]<=0)) #where original state contained ice, but EnKF has made ice-free

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
    T_skin[cond2] = -0.054*fo[member,cond2]

    sal_ice = fr.sal_ice.to_numpy()
    for layer in range(4):
        sal_ice[:,layer][cond1] = Si_new
        sal_ice[:,layer][cond2] = 0

    h_pond = fr.h_pond.to_numpy()
    h_pond[cond1] = 0
    h_pond[cond2] = 0

    fr.part_size.loc[:] = post
    fr.h_ice.loc[:] = h_ice
    fr.h_snow.loc[:] = h_snow
    fr.h_pond.loc[:] = h_pond
    fr.enth_ice.loc[:] = enth_ice
    fr.enth_snow.loc[:] = enth_snow
    fr.T_skin.loc[:] = T_skin
    fr.sal_ice.loc[:] = sal_ice

    fr.to_netcdf(file,mode='a')