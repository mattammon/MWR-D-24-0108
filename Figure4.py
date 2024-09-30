####### Imports #######
from netCDF4 import Dataset
import numpy as np
import utils as CF
import metpy.calc as mpc
from metpy.units import units

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.axes as ax
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from datetime import datetime, timedelta


dt = datetime(2022, 3, 1)
c=0

"""

The following "while" statement loops through every available CLAMPS lidar observation
from both CLAMPS sites during PERiLS-2022 and distributes the data into arrays that contain
mean 0-1 km u, v, and w wind (also wind direction in degrees and total magnitude) 
components from each individual observation. These arrays are used to construct the "wind compass"
and "wind rose"

"""

while dt.month < 5:
    if dt.day < 10:
        date = f'20220{dt.month}0{dt.day}'
    else:
        date = f'20220{dt.month}{dt.day}'
    try:
        nc = Dataset(f'data/CLAMPS1/VAD_1_{date}.nc')
        hgt = nc['height'][:]
        top_idx = CF.CVI(hgt,1)
        
        nc2 = Dataset(f'data/CLAMPS2/VAD_2_{date}.nc')
        hgt2 = nc2['height'][:]
        top_idx2 = CF.CVI(hgt2,1)

        wspd = nc['wspd'][:,:top_idx+1]
        wdir = nc['wdir'][:,:top_idx+1]
        w = nc['w'][:,:top_idx+1]
        
        wspd2 = nc2['wspd'][:,:top_idx2+1]
        wdir2 = nc2['wdir'][:,:top_idx2+1]
        w2 = nc2['w'][:,:top_idx2+1]

        wspd[np.where(wspd>=100)]=np.nan
        wdir[np.where(wspd>=100)]=np.nan
        w[np.where(wspd>=100)]=np.nan
        w[np.where(abs(w)>=100)]=np.nan
        
        wspd2[np.where(wspd2>=100)]=np.nan
        wdir2[np.where(wspd2>=100)]=np.nan
        w2[np.where(wspd2>=100)]=np.nan
        w2[np.where(abs(w2)>=100)]=np.nan

        u,v = mpc.wind_components(units.Quantity(wspd,'m/s'),units.Quantity(wdir,'deg'))
        mean_u = np.nanmean(u,axis=1)
        mean_v = np.nanmean(v,axis=1)
        mean_w = np.nanmean(w,axis=1)
        
        u2,v2 = mpc.wind_components(units.Quantity(wspd2,'m/s'),units.Quantity(wdir2,'deg'))
        mean_u2 = np.nanmean(u2,axis=1)
        mean_v2 = np.nanmean(v2,axis=1)
        mean_w2 = np.nanmean(w2,axis=1)

        if c == 0:
            all_mean_u = mean_u
            all_mean_v = mean_v
            all_mean_w = mean_w
            
            all_mean_u2 = mean_u2
            all_mean_v2 = mean_v2
            all_mean_w2 = mean_w2

        else:
            all_mean_u = np.concatenate((all_mean_u,mean_u))
            all_mean_v = np.concatenate((all_mean_v,mean_v))
            all_mean_w = np.concatenate((all_mean_w,mean_w))
            
            all_mean_u2 = np.concatenate((all_mean_u2,mean_u2))
            all_mean_v2 = np.concatenate((all_mean_v2,mean_v2))
            all_mean_w2 = np.concatenate((all_mean_w2,mean_w2))

        all_mean_wdir = mpc.wind_direction(units.Quantity(all_mean_u,'m/s'),units.Quantity(all_mean_v,'m/s'))
        all_mean_wspd = mpc.wind_speed(units.Quantity(all_mean_u,'m/s'),units.Quantity(all_mean_v,'m/s'))
        
        all_mean_wdir2 = mpc.wind_direction(units.Quantity(all_mean_u2,'m/s'),units.Quantity(all_mean_v2,'m/s'))
        all_mean_wspd2 = mpc.wind_speed(units.Quantity(all_mean_u2,'m/s'),units.Quantity(all_mean_v2,'m/s'))

        print(all_mean_u.shape,all_mean_v.shape,all_mean_w.shape)
        print(all_mean_u2.shape,all_mean_v2.shape,all_mean_w2.shape)
        
        c+=1
    except:
        pass
    dt = dt + timedelta(days=1)

#####################################################
    
"""

The following section of code calculates the bin means, medians, and relative frequencies of vertical 
wind that are directly plotted on the wind rose and wind compass

"""
    
d = 0 
dd = 0

bin_wid = 4
freq_wid = 4

bin_meds = []
bin_means = []
bin_size = []

bin_meds2 = []
bin_means2 = []
bin_size2 = []

while d <= 360:
    bin_idxs = np.where((all_mean_wdir.m >=d) & (all_mean_wdir.m < d+bin_wid))
    freq_idxs = np.where((all_mean_wdir.m >=dd) & (all_mean_wdir.m < dd+freq_wid))
    bin_means.append(np.nanmean(all_mean_w[bin_idxs]))
    bin_meds.append(np.nanmedian(all_mean_w[bin_idxs]))
    bin_size.append(len(freq_idxs[0]))
    
    bin_idxs2 = np.where((all_mean_wdir2.m >=d) & (all_mean_wdir2.m < d+bin_wid))
    freq_idxs2 = np.where((all_mean_wdir2.m >=dd) & (all_mean_wdir2.m < dd+freq_wid))
    bin_means2.append(np.nanmean(all_mean_w2[bin_idxs2]))
    bin_meds2.append(np.nanmedian(all_mean_w2[bin_idxs2]))  
    bin_size2.append(len(freq_idxs2[0]))
    d+=bin_wid
    dd+=freq_wid

rel_freq = []
rel_freq2 = []

total = np.sum(bin_size)
total2 = np.sum(bin_size2)
i=0
while i < (360/freq_wid):
    rel_freq.append(bin_size[i]/total)
    rel_freq2.append(bin_size2[i]/total2)
    i+=1
i=0

#######################################################

"""

The following code handles the plotting of Figure 4

"""

fig = plt.figure(figsize = (10,5), facecolor = 'white')
bins = np.radians(np.arange(0,361,bin_wid))
bins_freq = np.radians(np.arange(0,360,freq_wid))
#bins_freq = [0,45,90,135,180,225,270,315]

ax = plt.subplot((121), projection = 'polar')
ax.plot(bins, bin_meds, color = 'blue',linewidth=1.5)
ax.plot(bins, bin_meds2, color = 'red',linewidth=1.5)
ax.fill_between(bins, bin_meds2, bin_meds, where = np.array(bin_meds2) > np.array(bin_meds), interpolate = True, color = 'red', alpha = 0.45)
ax.fill_between(bins, bin_meds2, bin_meds, where = np.array(bin_meds2) < np.array(bin_meds), interpolate = True, color = 'blue', alpha = 0.45)

ax.plot(np.radians(np.array([194, 14])), [0.25,0.25], color = 'k', alpha = 1, linewidth = 2)

ax.set_theta_zero_location('N')
ax.set_title(f'Median Vertical Motion (m/s) vs. Wind Direction', fontsize = 12)
ax.set_theta_direction(-1)
ax.set_ylim(-0.1,0.1)
ax.set_yticks(np.arange(-0.1,0.15,0.05))

ax.plot(bins, np.zeros_like(bins), color = 'black',linewidth=2,linestyle='--')

ax2 = plt.subplot((122), projection = 'polar')
ax2.bar(bins_freq, rel_freq, width=np.radians(freq_wid), bottom=0,align='edge',color='blue',edgecolor='blue',label='LVA-C1 (6673 Total Observations)',alpha=0.45)
ax2.bar(bins_freq, rel_freq, width=np.radians(freq_wid), bottom=0,align='edge',fill=False,edgecolor='blue',alpha=1)

ax2.bar(bins_freq, rel_freq2, width=np.radians(freq_wid), bottom=0,align='edge',color='red',edgecolor='red',label='YCM-C2 (6832 Total Observations)',alpha=0.45)
ax2.bar(bins_freq, rel_freq2, width=np.radians(freq_wid), bottom=0,align='edge',fill=False,edgecolor='red',alpha=1)

ax2.set_theta_zero_location('N')
ax2.set_title('Relative Frequency of Observance', fontsize = 12)
ax2.set_theta_direction(-1)
ax2.set_ylim(0,0.04)
ax2.set_yticks(np.arange(0,0.041,0.01))

legend_elements = [Patch(facecolor='red',edgecolor='red',label = 'YCM-C2 (6832 Total Observations)'),
                   Patch(facecolor='blue',edgecolor='blue',label = 'LVA-C1 (6673 Total Observations)'),
                  Line2D([0],[0],color='black',lw=2.5,label='Approximate Orientation of YBR'),
                  Line2D([0],[0],color='black',lw=2.5,ls='--',label='Median Vertical Motion of Zero')]

fig.legend(handles=legend_elements,loc='lower center')
plt.savefig('NEW_FIGS/MED_FREQ_BOTH.png')
