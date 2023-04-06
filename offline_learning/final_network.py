import numpy as np
import xarray as xr
from preprocessing import pad
from NNetwork import *

argsA = {
'kernel_size':3,
'epochs':150,
'lr':0.001,
'batch_size':10,
'zero_padding':0,
'h_channels':[32,64,128],
'n_classes':1,
'stride':1,
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
'bias':False,
'seed':711,
}

NetworkA_weights = '../CNN_weights/NetworkA_weights_1982-2017_allsamples.pt'
NetworkB_weights = '../CNN_weightsNetworkB_weights_1982-2017_allsamples.pt'

NetworkA_stats = np.load('../data_files/NetworkA_statistics_1982-2017_allsamples.npz')
NetworkB_stats = np.load('../data_files/NetworkB_statistics_1982-2017_allsamples.npz')
                                                                                                                                                                   
inputs = ['SIC','SST','SIU','SIV','SIT','SW','TS','SSS']

forecasts = xr.open_dataset('seaice_DAML_inputs_1982-2017.nc')
increments = xr.open_dataset('seaice_DAML_outputs_1982-2017.nc')
dSIC = increments.dSIC.to_numpy()
dSICN = increments.dSICN.to_numpy()

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

X = np.hstack((X,land_mask[:,None]))
X[np.isnan(X)] = 0

for N in range(X.shape[1]-1): 
    X[:,N] = (X[:,N]-NetworkA_stats['mu'][N])/NetworkA_stats['sigma'][N] #standardize inputs
    X[:,N][land_mask==0] = 0

dSIC = Net(X,argsA,y_train=dSIC,x_valid=X,y_valid=dSIC,path=NetworkA_weights)[:,0] #generate aggregate SIC increment prediction
dSIC[land_mask[:,4:-4,4:-4]==0] = 0

### NETWORK B ###                                                                                                                                                                  
X = [dSIC]
for CAT in range(5):
    X.append(forecasts['SICN'].isel(n=0,ct=CAT).to_numpy()[None])
    X.append(forecasts['SICN'].isel(n=1,ct=CAT).to_numpy()[None])
X = np.transpose(X,(1,0,2,3))

X = np.hstack((X,land_mask[:,None,4:-4,4:-4]))
X[np.isnan(X)]

for N in range(X.shape[1]-1):
    X[:,N] = (X[:,N]-NetworkB_stats['mu'][N])/NetworkB_stats['sigma'][N] #standardize inputs
    X[:,N][land_mask[:,4:-4,4:-4]==0] = 0

dSICN_pred = Net(X,argsB,y_train=dSICN,x_valid=X,y_valid=dSICN,path=NetworkB_weights) #generate category SIC increment prediction
for CAT in range(5):
    dSICN_pred[:,CAT][land_mask[:,4:-4,4:-4]==0] = 0
dSICN[member] = dSICN_pred

np.save('dSICN_increment_1982-2017_allsamples.npy',dSICN)
