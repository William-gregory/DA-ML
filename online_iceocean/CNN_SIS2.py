import os
import numpy as np
import xarray as xr
import glob
from preprocessing import pad
from NNetwork import *
from datetime import datetime,timedelta

def pp(x):
    """
    post-processing to ensure updated sea ice concentration
    is bounded between 0 and 1.
    """
    x[x<0] = 0
    SIC = np.nansum(x,0)
    high = SIC>1
    ratio = 1/SIC[high]
    for CAT in range(5):
        x[CAT,high] = x[CAT,high]*ratio
    return x

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

#dictionaries containing architecture/hyperparams for CNN
argsA = {
'kernel_size':3,
'zero_padding':0,
'h_channels':[32,64,128],
'n_classes':1,
'stride':1,
'bias':False,
'seed':711,
}

argsB = {
'kernel_size':1,
'zero_padding':0,
'h_channels':[32,64,128],
'n_classes':5,
'stride':1,
'bias':False,
'seed':711,
}

NetworkA_weights = '/ncrc/home2/William.Gregory/DA-ML/CNN_weights/NetworkA_weights_CNNopt.pt'
NetworkB_weights = '/ncrc/home2/William.Gregory/DA-ML/CNN_weights/NetworkB_weights_CNNopt.pt'

NetworkA_stats = np.load('/ncrc/home2/William.Gregory/DA-ML/data_files/NetworkA_statistics_1982-2017_allsamples.npz')
NetworkB_stats = np.load('/ncrc/home2/William.Gregory/DA-ML/data_files/NetworkB_statistics_1982-2017_allsamples.npz')

experiment = os.getcwd().split('/')[-2].split('.')[0]
user = os.popen('whoami').read().split('\n')[0]
savepath = '/lustre/f2/dev/'+user+'/CNN_increments/'+experiment+'/'
if os.path.exists(savepath)==False:
    os.mkdir(savepath)

files = sorted(glob.glob('../*ice_daily*')) #history files containing states to generate prediction                                                                                                                                
f = xr.open_mfdataset(files,combine='nested',concat_dim='ens')
y,m,d = np.genfromtxt('coupler.res',skip_header=1)[1,:3].astype(np.int32)
date_in = files[0].split('..')[1].split('.')[0] #get start date of history file                                                                                                                                                    
date_out = datetime(y,m,d).strftime('%Y%m%d') #get date at which correction is applied   

states = f.mean('time')
tend = f.diff('time').mean('time')
nmembers = len(f.ens)
yT = len(f.yT)
xT = len(f.xT)
pad_size = 4
scaling = len(f.time)/5 #applied to the increments at the end in case the correction is applied at different frequencies, e.g., 2-day vs 5-day etc.                                                                                
                        #CNN was originally trained on data from a 5-day DA cycle, so we just linearly scale. 

dSICN = np.zeros((nmembers,1,argsB['n_classes'],yT,xT)) #compute an increment for every ensemble member
inputs = ['siconc','SST','UI','VI','HI','SW','TS','SSS']

restarts = sorted(glob.glob('ice_model.res*')) #prior model states (raw RESTART files)                                                                                                                                             
rho_ice = 905.
rho_snow= 330.
phi_init = 0.75 #initial liquid fraction of frazil ice                                                                                                                                                                             
Si_new = 5 #salinity of mushy ice                                                                                                                                                                                                  
Ti = min(liquidus_temperature_mush(Si_new/phi_init),-0.1)
qi_new = enthalpy_ice(Ti, Si_new)
hlim = [1.0e-10, 0.1, 0.3, 0.7, 1.1, 1.5]
hmid = np.array([0.5*(hlim[n]+hlim[n+1]) for n in range(5)])
i_thick = np.tile((hmid*rho_ice)[None,:,None,None],(1,1,yT,xT))

for member,file in enumerate(restarts): #compute increment and add to each ensemble member    
    ### NETWORK A ### 
    X = []
    for label in inputs:
        X.append(pad(states[label].isel(ens=member).to_numpy()[None],label,pad_size))
        X.append(pad(tend[label].isel(ens=member).to_numpy()[None],label,pad_size))
    X = np.transpose(X,(1,0,2,3))

    land_mask = np.copy(X[:,0])
    land_mask[~np.isnan(land_mask)] = 1
    land_mask[np.isnan(land_mask)] = 0
    land_mask[:,:pad_size] = 0

    X = np.hstack((X,land_mask[:,None]))
    X[np.isnan(X)] = 0

    for N in range(X.shape[1]-1): 
        X[:,N] = (X[:,N]-NetworkA_stats['mu'][N])/NetworkA_stats['sigma'][N] #standardize inputs
        X[:,N][land_mask==0] = 0

    argsA['n_channels'] = X.shape[1]
    dSIC = Net(X,argsA,path=NetworkA_weights)[:,0] #generate aggregate SIC increment prediction
    dSIC[land_mask[:,pad_size:-pad_size,pad_size:-pad_size]==0] = 0
    
    ### NETWORK B ###                                                                                                                                                                  
    X = [dSIC]
    for CAT in range(5):
        X.append(states['CN'].isel(ct=CAT,ens=member).to_numpy()[None])
        X.append(tend['CN'].isel(ct=CAT,ens=member).to_numpy()[None])
    X = np.transpose(X,(1,0,2,3))

    X = np.hstack((X,land_mask[:,None,pad_size:-pad_size,pad_size:-pad_size]))
    X[np.isnan(X)] = 0

    for N in range(X.shape[1]-1):
        X[:,N] = (X[:,N]-NetworkB_stats['mu'][N])/NetworkB_stats['sigma'][N] #standardize inputs
        X[:,N][land_mask[:,pad_size:-pad_size,pad_size:-pad_size]==0] = 0

    argsB['n_channels'] = X.shape[1]
    dSICN_pred = Net(X,argsB,path=NetworkB_weights) #generate category SIC increment prediction
    for CAT in range(5):
        dSICN_pred[:,CAT][land_mask[:,pad_size:-pad_size,pad_size:-pad_size]==0] = 0
    dSICN[member] = scaling*dSICN_pred

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

ds = xr.Dataset(data_vars=dict(dSICN=(['time', 'ct', 'yT', 'xT'], np.nanmean(dSICN,0))), coords=dict(yT=f['yT'], xT=f['xT']))
ds.dSICN.attrs['long_name'] = 'category_sea_ice_concentration_increments'
ds.dSICN.attrs['units'] = 'area_fraction'
ds['time'] = [date_out]
ds.to_netcdf(savepath+date_in+'.CNN_increment.'+date_out+'.nc') #save ensemble mean increment
