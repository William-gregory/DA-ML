import numpy as np
import xarray as xr
import glob
from preprocessing import pad
from NNetwork import *

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

NetworkA_weights = '../CNN_weights/NetworkA_weights_1982-2017_allsamples.pt'
NetworkB_weights = '../CNN_weightsNetworkB_weights_1982-2017_allsamples.pt'

NetworkA_stats = np.load('../data_files/NetworkA_statistics_1982-2017_allsamples.npz')
NetworkB_stats = np.load('../data_files/NetworkB_statistics_1982-2017_allsamples.npz')

files = sorted(glob.glob('*ice_daily*'))
f = xr.open_mfdataset(files,combine='nested',concat_dim='ens')
states = f.mean('time')
tend = f.diff('time').mean('time')
nmembers = len(f.ens)
yT = len(f.yT)
xT = len(f.xT)
scaling = len(f.time)/5 #applied to the increments at the end in case the correction is applied at different frequencies, e.g., 2-day vs 5-day etc.
                        #CNN was originally trained on data from a 5-day DA cycle, so we just linearly scale.

dSICN = np.zeros((nmembers,1,argsB['n_classes'],yT,xT)) #compute an increment for every ensemble member
                                                                                                                                                                   
inputs = ['siconc','SST','UI','VI','HI','SW','TS','SSS']
for member in range(nmembers):
    
    ### NETWORK A ### 
    X = []
    for label in inputs:
        X.append(pad(states[label].isel(ens=member).to_numpy()[None],label,4))
        X.append(pad(tend[label].isel(ens=member).to_numpy()[None],label,4))
    X = np.transpose(X,(1,0,2,3))

    land_mask = np.copy(X[:,0])
    land_mask[~np.isnan(land_mask)] = 1
    land_mask[np.isnan(land_mask)] = 0
    land_mask[:,:4] = 0

    X = np.hstack((X,land_mask[:,None]))
    X[np.isnan(X)] = 0

    for N in range(X.shape[1]-1): 
        X[:,N] = (X[:,N]-NetworkA_stats['mu'][N])/NetworkA_stats['sigma'][N] #standardize inputs
        X[:,N][land_mask==0] = 0

    dSIC = Net(X,argsA,path=NetworkA_weights)[:,0] #generate aggregate SIC increment prediction
    dSIC[land_mask[:,4:-4,4:-4]==0] = 0
    
    ### NETWORK B ###                                                                                                                                                                  
    X = [dSIC]
    for CAT in range(5):
        X.append(states['CN'].isel(ct=CAT,ens=member).to_numpy()[None])
        X.append(tend['CN'].isel(ct=CAT,ens=member).to_numpy()[None])
    X = np.transpose(X,(1,0,2,3))

    X = np.hstack((X,land_mask[:,None,4:-4,4:-4]))
    X[np.isnan(X)] = 0

    for N in range(X.shape[1]-1):
        X[:,N] = (X[:,N]-NetworkB_stats['mu'][N])/NetworkB_stats['sigma'][N] #standardize inputs
        X[:,N][land_mask[:,4:-4,4:-4]==0] = 0

    dSICN_pred = Net(X,argsB,path=NetworkB_weights) #generate category SIC increment prediction
    for CAT in range(5):
        dSICN_pred[:,CAT][land_mask[:,4:-4,4:-4]==0] = 0
    dSICN[member] = dSICN_pred

np.save('dSICN_increment.npy',scaling*dSICN)
