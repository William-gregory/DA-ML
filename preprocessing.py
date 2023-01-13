import numpy as np
import xarray as xr
from datetime import datetime, timedelta

def date_range(start, end):
    """
    for a given start and end date, return all dates in between
    """
    delta = end - start
    days = [(start + timedelta(days=i)).strftime('%Y%m%d') for i in range(delta.days + 1)]
    return np.array(days)

def pad(data,label,size):
    """
    Pad a 3D array by 'size' along all edges.
    If data are velocity files (UI or VI), first compute
    2-pt moving average so that fields are defined on the 
    same tracer grid as scalar fields
    """
    if label == 'UI':
        data = np.nansum([data[:,:,1:],data[:,:,:-1]],0)/2
        sign = -1
    elif label == 'VI':
        data = np.nansum([data[:,1:],data[:,:-1]],0)/2
        sign = -1
    else:
        sign = 1
        
    n,dX,dY = data.shape
    data = np.hstack((data,sign*np.flip(data[:,-size:,:],axis=(1,2))))
    return np.hstack((np.zeros((n,size,dY+2*size)),np.dstack((data[:,:,-size:],np.dstack((data,data[:,:,:size]))))))

def states(data,subset):
    """
    Compute 5-day means of DA foreast states
    """
    return xr.concat([data.isel(time=np.arange(x-5,x)).mean('time') for x in subset],dim='time')

def tendencies(data,subset):
   """
   Compute 5-day means of DA forecast tendencies
   """
   return xr.concat([data.isel(time=np.arange(x-5,x)).diff('time').mean('time') for x in subset],dim='time')
