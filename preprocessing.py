import numpy as np

def pad(data,label,size):
    """
    Pad a 3D array by 'size' along its spatial dimensions.
    If data are velocity files (UI or VI), first compute
    2-pt moving average so that fields are defined on the 
    same tracer grid as scalar fields
    """
    if label == 'SIU':
        data = np.nansum([data[:,:,1:],data[:,:,:-1]],0)/2
        sign = -1
    elif label == 'SIV':
        data = np.nansum([data[:,1:],data[:,:-1]],0)/2
        sign = -1
    else:
        sign = 1
        
    n,dX,dY = data.shape
    data = np.hstack((data,sign*np.flip(data[:,-size:,:],axis=(1,2))))
    return np.hstack((np.zeros((n,size,dY+2*size)),np.dstack((data[:,:,-size:],np.dstack((data,data[:,:,:size]))))))
