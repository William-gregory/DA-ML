import os
import glob
import itertools
import numpy as np
import xarray as xr
from mpi4py import MPI
from datetime import datetime,timedelta

COMM = MPI.COMM_WORLD

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

def haversine(coords,localization_radius=0.03):
    """
    compute the localization matrix based on Haversine distance.
    coords: Nx2 array of longitude and latitude points

    returns: NxN binary matrix with zeros outside of localization radius
    """
    lon = np.deg2rad(coords[:, 0])
    lat = np.deg2rad(coords[:, 1])
    lon1, lon2 = np.meshgrid(lon, lon)
    lat1, lat2 = np.meshgrid(lat, lat)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    a = np.sin(dlat / 2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0)**2
    c = 2 * np.arcsin(np.sqrt(a))
    c[c<=localization_radius] = 1
    c[c!=1] = 0
    
    return c

def split(container, count):
    """
    function for dividing the number of tasks (container)
    across the number of available compute nodes (count)
    """
    return [container[_i::count] for _i in range(count)]

def Kfilter(prior,obs,lon,lat,loc_rad,obs_error):
    """
    ensemble Kalman filter (still need to do 'adjustment' part)
    prior state variables of size: E, C, dX, dY
    E is the number of ensemble members
    C is the number of model states (categories)
    dX is the number of x grid points
    dY is the number of y grid points

    returns: posterior state variables and increments
    """
    E,C,dX,dY = prior.shape
    prior = prior.reshape(E,C,dX*dY)
    obs = obs.ravel()
    obs[np.isnan(obs)] = 0
    prior[np.isnan(prior)] = 0
    priorH = np.nansum(prior,1) #prior passed through obs operator
    if ((priorH==0).all() & (obs==0).all()):
        tmp = np.zeros((E,C,dX,dY))
        return tmp,tmp
    else:
        coords = np.array([lat.ravel(),lon.ravel()]).T
        W = haversine(coords,loc_rad)

        Bm = np.einsum('ijk,il->jkl', prior-np.nanmean(prior,0), priorH-np.nanmean(priorH,0)) / (E - 1)
        Bo = np.cov(priorH.T) + np.eye(len(W))*obs_error
        K = W * (Bm @ np.linalg.inv(Bo))
        posterior = np.array([prior[x] + (K @ (obs - priorH[x])) for x in range(E)])

        return pp(posterior.reshape(E,C,dX,dY)), (posterior-prior).reshape(E,C,dX,dY)

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
xdiv = xT//10
ydiv = yT//10
xindices,yindices = np.meshgrid(np.arange(0,xT,xdiv),np.arange(0,yT,ydiv))
xindices = xindices.ravel()
yindices = yindices.ravel()

selected_variables = range(len(yindices)) #divide globe into domains of size 32x36
if COMM.rank == 0:
    splitted_jobs = split(selected_variables, COMM.size)
else:
    splitted_jobs = None
scattered_jobs = COMM.scatter(splitted_jobs, root=0) #scatter the tasks to each of the computer nodes.

results = []
for ix in scattered_jobs: #each computer node will execute this loop and send its segement of data to the ENKF
    outputs = Kfilter(prior[:,:,xindices[ix]:xindices[ix]+xdiv,yindices[ix]:yindices[ix]+ydiv],\
                          obs[xindices[ix]:xindices[ix]+xdiv,yindices[ix]:yindices[ix]+ydiv],\
                          lon[xindices[ix]:xindices[ix]+xdiv,yindices[ix]:yindices[ix]+ydiv],\
                          lat[xindices[ix]:xindices[ix]+xdiv,yindices[ix]:yindices[ix]+ydiv],loc_rad=localization_radius,obs_error=observation_error)
    results.append(outputs)
results = COMM.gather(results, root=0)

if COMM.rank == 0: #tell the master node to compile the results into their own respective arrays and map back to the 2D domain
    posterior = np.zeros((nmembers,nCat,xT,yT))
    increments = np.zeros((nmembers,nCat,xT,yT))
    results = list(itertools.zip_longest(*results))
    ix = 0
    for r1 in results:
        for r2 in r1:
            if r2 is not None:
                posterior[:,:,xindices[ix]:xindices[ix]+xdiv,yindices[ix]:yindices[ix]+ydiv] = r2[0]
                increments[:,:,xindices[ix]:xindices[ix]+xdiv,yindices[ix]:yindices[ix]+ydiv] = r2[1]
            ix += 1

    ds = xr.Dataset(data_vars=dict(part_size=(['member','time', 'ct', 'yT', 'xT'], increments[:,np.newaxis])), coords=dict(yT=f['yT'], xT=f['xT']))
    ds.dSICN.attrs['long_name'] = 'category_sea_ice_concentration_increments'
    ds.dSICN.attrs['units'] = 'area_fraction'
    ds['time'] = [date]
    ds.to_netcdf(savepath+date+'.EnKF_increment.nc')

"""
for member,file in enumerate(ice_restarts): #compute increment and add to each ensemble member    
    
    ### ADD TO RESTART FILE ###                                                                                                                                                                                                    
    fr = xr.open_dataset(file)
    prior = fr.part_size.to_numpy()
    post = np.zeros((1,6,320,360))
    post[0,1:] = pp(prior[0,1:] + dSICN[member,0])
    post[0,0] = 1 - np.nansum(post[0,1:],0)
    post[post<0] = 0
    post[post>1] = 1

    cond1 = np.where((prior[:,1:]<=0) & (post[:,1:]>0)) #where original state was ice-free, but CNN has added ice
    cond2 = np.where((prior[:,1:]>0) & (post[:,1:]<=0)) #where original state contained ice, but CNN has made ice-free
    
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
    T_skin[cond2] = -1.8

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
"""
