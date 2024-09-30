from netCDF4 import Dataset
import numpy as np

#Coordinates (lat,lon) of the CLAMPS deployment sites
coords_C1 = [33.363853, -91.26238]
coords_C2 = [32.9104, -90.3815]

#Dates (mmdd) and times (hhmm) of arrival of convection for the nine cases at LVA-C1
days_C1 = ['0307','0322','0330','0405','0411','0412','0413','0415','0416']
times_C1 = ['1112','1100','1944','0641','0707','2007','1938','2214','1419']

#Dates (mmdd) and times (hhmm) of arrival of convection for the nine cases at YCM-C2
days_C2 = ['0307','0322','0330','0405','0411','0413','0415','0416','0417']
times_C2 = ['1309','1642','2040','0825','1022','1950','2229','1624','2220']

#Function that returns the index of the value that is closest to a desired value in an array
def CVI(variable, desired_value):
    arr = np.asarray(variable)
    location = (np.abs(arr - desired_value)).argmin()
    return location

#Function that returns a netcdf dataset of processed WINDoe data for a given CLAMPS site and case number
def grab_WINDoe_nc(clamps,idx):
    """
    Parameters
    ----------
    clamps : int
        1 for CLAMPS1 or 2 for CLAMPS2
    idx : int
        case number (1-9) in time order

    Returns
    -------
    windoe : Dataset
        netcdf dataset with output variables from WINDoe

    """
        
    if clamps==1:
        days=days_C1
        times=times_C1
    else:
        days=days_C2
        times=times_C2
    try:
        file = f'/Users/matthew.ammon/Documents/Paper_Code/data/WINDoe/WINDoe_C{clamps}_5minres.2022{days[idx]}.nc'
        windoe = Dataset(file)
    except:
        pass
    return windoe

#Generates u, v, height, and the components of mean 0-6 km wind for a given clamps site and case number
def array_gen(clamps,time_idx):
    """
    Parameters
    ----------
    clamps : int
        1 for CLAMPS1 or 2 for CLAMPS2
    time_idx : int
        time step (0-13) of desired data, each timestep is 10 minutes starting from 2 hours prior to convection

    Returns
    -------
    all_u : array
        values of u component of wind for all cases at desired time step (time_idx)
    all_v : array
        values of v component of wind for all cases at desired time step (time_idx)
    all_hgt : array
        values of height levels for all cases at desired time step (time_idx)
    u_06mean : array
        values of u component of mean 0-6 km wind for all cases at desired time step (time_idx)
    v_06mean : array
        values of v component of mean 0-6 km wind for all cases at desired time step (time_idx)

    """
    u_06mean = []
    v_06mean = []
    i = 0
    while i < len(days_C1):
        nc = grab_WINDoe_nc(clamps,i)
        hgt = nc['height'][:]
        u = nc['u_wind'][time_idx,:]
        v = nc['v_wind'][time_idx,:]
        
        u_06 = np.nanmean(u[0:CVI(hgt,6)])
        u_06mean.append(u_06)
        v_06 = np.nanmean(v[0:CVI(hgt,6)])       
        v_06mean.append(v_06)

        if i==0:
            all_u = u
            all_v = v
            all_hgt = hgt
        else:
            all_u = np.vstack((all_u,u))
            all_v = np.vstack((all_v,v))
            all_hgt = np.vstack((all_hgt,hgt))
        i+=1
    return all_u,all_v,all_hgt,u_06mean,v_06mean