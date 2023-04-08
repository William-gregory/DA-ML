import os
import numpy as np
import xarray as xr
from preprocessing import pad
from NNetwork import *

if os.path.exists('seaice_DA-ML_inputs_1982-2017.nc')==False:
    os.system('wget -nv ftp://sftp.gfdl.noaa.gov/perm/William.Gregory/seaice_DA-ML_inputs_1982-2017.nc')
if os.path.exists('seaice_DA-ML_outputs_1982-2017.nc')==False:
    os.system('wget -nv ftp://sftp.gfdl.noaa.gov/perm/William.Gregory/seaice_DA-ML_outputs_1982-2017.nc')

def LossA(outputs, targets):
    return torch.mean((outputs-targets)**2)
  
def LossB(outputs, targets):
    return torch.sum(torch.mean((outputs-targets)**2,(0,2,3)) + 5*torch.mean((torch.sum(outputs,1)-torch.sum(targets,1))**2)

argsA = {
'kernel_size':3,
'epochs':150,
'lr':0.001,
'batch_size':10,
'zero_padding':0,
'h_channels':[32,64,128],
'n_classes':1,
'stride':1,
'loss':LossA,
'wd':1e-7,
'bias':False,
'seed':711,
}

argsB = {
'kernel_size':1,
'epochs':125,
'lr':0.001,
'batch_size':10,
'zero_padding':0,
'h_channels':[32,64,128],
'n_classes':5,
'stride':1,
'loss':LossB,
'wd':1e-7,
'bias':False,
'seed':711,
}

NetworkA_weights = '../CNN_weights/NetworkA_weights_1982-2017_allsamples.pt'
NetworkB_weights = '../CNN_weightsNetworkB_weights_1982-2017_allsamples.pt'
                                                                       
inputs = ['SIC','SST','SIU','SIV','SIT','SW','TS','SSS']

forecasts = xr.open_dataset('seaice_DA-ML_inputs_1982-2017.nc')
increments = xr.open_dataset('seaice_DA-ML_outputs_1982-2017.nc')
dSIC = increments.dSIC.to_numpy()
dSICN = increments.dSICN.to_numpy()
lat = xr.open_dataset('../data_files/ice.static.nc').GEOLAT.to_numpy()

### NETWORK A ### 
X = []
for label in inputs:
    X.append(pad(forecasts[label].isel(n=0).to_numpy()[None],label,4))
    X.append(pad(forecasts[label].isel(n=1).to_numpy()[None],label,4))
X = np.transpose(X,(1,0,2,3))

land_mask = np.copy(X[:,0])
land_mask[~np.isnan(land_mask)] = 1
land_mask[np.isnan(land_mask)] = 0
land_mask[:,:4] = 0
                     
#compute statistics over ocean grid cells above 40 degrees latitude:
NH = np.where((lat>40) & (land_mask[0,4:-4,4:-4]==1))
SH = np.where((lat<-40) & (land_mask[0,4:-4,4:-4]==1))
normIDs = (np.concatenate((NH[0],SH[0])),np.concatenate((NH[1],SH[1]))) 

X = np.hstack((X,land_mask[:,None]))
X[np.isnan(X)] = 0

if os.path.exists('../data_files/NetworkA_statistics_1982-2017_allsamples.npz')
    stats = np.load('../data_files/NetworkA_statistics_1982-2017_allsamples.npz')
    for N in range(X.shape[1]-1): 
        X[:,N] = (X[:,N]-stats['mu'][N])/stats['sigma'][N] #standardize inputs
        X[:,N][land_mask==0] = 0
else:
    mu = np.zeros(X.shape[1]-1)
    sigma = np.zeros(X.shape[1]-1)
    for N in range(X.shape[1]-1):
        X_nopad = X[:,N,pad_size:-pad_size,pad_size:-pad_size]
        mu[N] = np.nanmean(X_nopad[:,normIDs[0],normIDs[1]])
        sigma[N] = np.nanstd(X_nopad[:,normIDs[0],normIDs[1]])
        X[:,N] = (X[:,N]-mu[N])/sigma[N]
        X[:,N][land_mask==0] = 0
    np.savez('../data_files/NetworkA_statistics_1982-2017_allsamples.npz',mu=mu,sigma=sigma)

argsA['n_channels'] = X.shape[1]
dSIC = Net(X,argsA,y_train=dSIC,x_valid=X,y_valid=dSIC,path=NetworkA_weights)[:,0] #generate aggregate SIC increment prediction
dSIC[land_mask[:,4:-4,4:-4]==0] = 0

### NETWORK B ###                                                                                                                                                                  
X = [dSIC]
for CAT in range(5):
    X.append(forecasts['SICN'].isel(n=0,ct=CAT).to_numpy()[None])
    X.append(forecasts['SICN'].isel(n=1,ct=CAT).to_numpy()[None])
X = np.transpose(X,(1,0,2,3))

X = np.hstack((X,land_mask[:,None,4:-4,4:-4]))
X[np.isnan(X)] = 0

if os.path.exists('../data_files/NetworkB_statistics_1982-2017_allsamples.npz')
    stats = np.load('../data_files/NetworkB_statistics_1982-2017_allsamples.npz')
    for N in range(X.shape[1]-1): 
        X[:,N] = (X[:,N]-stats['mu'][N])/stats['sigma'][N] #standardize inputs
        X[:,N][land_mask==0] = 0
else:   
    mu = np.zeros(X.shape[1]-1)
    sigma = np.zeros(X.shape[1]-1)
    for N in range(X.shape[1]-1):
        mu[N] = np.nanmean(X[:,N,normIDs[0],normIDs[1]])
        sigma[N] = np.nanstd(X[:,N,normIDs[0],normIDs[1]])
        X[:,N] = (X[:,N]-mu[N])/sigma[N]
        X[:,N][land_mask[:,4:-4,4:-4]==0] = 0
    np.savez('../data_files/NetworkB_statistics_1982-2017_allsamples.npz',mu=mu,sigma=sigma)

argsB['n_channels'] = X.shape[1]
dSICN_pred = Net(X,argsB,y_train=dSICN,x_valid=X,y_valid=dSICN,path=NetworkB_weights) #generate category SIC increment prediction
for CAT in range(5):
    dSICN_pred[:,CAT][land_mask[:,4:-4,4:-4]==0] = 0
dSICN[member] = dSICN_pred

np.save('dSICN_increment_1982-2017_allsamples.npy',dSICN)
