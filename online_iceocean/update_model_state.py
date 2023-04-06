import xarray as xr
import glob
import numpy as np

def pp(x):
    """
    post-processing to ensure updated sea ice concentration
    is bounded between 0 and 1.
    """
    x[x<0] = 0
    SIC = np.nansum(x,0)
    high = SIC>1
    ratio = 1/SIC[high]
    for CN in range(5):
        x[CN,high] = x[CN,high]*ratio
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
  
  
dSICN = np.load('dSICN_increment.npy') #category increments computed during simulation
files = sorted(glob.glob('ice_model.res*')) #prior model states (raw RESTART files)
rho_ice = 905.
rho_snow= 330.
phi_init = 0.75 #initial liquid fraction of frazil ice                                                                                                                                 
Si_new = 5 #salinity of mushy ice                                                                                                                                                      
Ti = min(liquidus_temperature_mush(Si_new/phi_init),-0.1)
qi_new = enthalpy_ice(Ti, Si_new)
hlim = [1.0e-10, 0.1, 0.3, 0.7, 1.1, 1.5]
hmid = np.array([0.5*(hlim[n]+hlim[n+1]) for n in range(5)])
i_thick = np.tile((hmid*rho_ice)[None,:,None,None],(1,1,320,360))

for member,file in enumerate(files): #add increment to each ensemble member
    f = xr.open_dataset(file)
    prior = f.part_size.to_numpy()
    post = np.zeros((1,6,320,360))
    post[0,1:] = pp(prior[0,1:] + dSICN[member,0])
    post[0,0] = 1 - np.nansum(post[0,1:],0)
    post[post<0] = 0
    post[post>1] = 1

    cond1 = np.where((prior[:,1:]<=0) & (post[:,1:]>0)) #where original state was ice-free, but CNN has added ice
    cond2 = np.where((prior[:,1:]>0) & (post[:,1:]<=0)) #where original state contained ice, but CNN has made ice-free

    h_ice = f.h_ice.to_numpy()
    h_ice[cond1] = i_thick[cond1]
    h_ice[cond2] = 0

    h_snow = f.h_snow.to_numpy()
    h_snow[cond1] = 0
    h_snow[cond2] = 0

    enth_ice = f.enth_ice.to_numpy()
    for layer in range(4):
        enth_ice[:,layer][cond1] = qi_new
        enth_ice[:,layer][cond2] = 0
        
    enth_snow = f.enth_snow.to_numpy()
    enth_snow[0][cond1] = 0
    enth_snow[0][cond2] = 0

    T_skin = f.T_skin.to_numpy()
    T_skin[cond1] = Ti
    T_skin[cond2] = -1.8

    sal_ice = f.sal_ice.to_numpy()
    for layer in range(4):
        sal_ice[:,layer][cond1] = Si_new
        sal_ice[:,layer][cond2] = 0

    h_pond = f.h_pond.to_numpy()
    h_pond[cond1] = 0
    h_pond[cond2] = 0

    f.part_size.loc[:] = post
    f.h_ice.loc[:] = h_ice
    f.h_snow.loc[:] = h_snow
    f.h_pond.loc[:] = h_pond
    f.enth_ice.loc[:] = enth_ice
    f.enth_snow.loc[:] = enth_snow
    f.T_skin.loc[:] = T_skin
    f.sal_ice.loc[:] = sal_ice
    
    f.to_netcdf(file,mode='a')
