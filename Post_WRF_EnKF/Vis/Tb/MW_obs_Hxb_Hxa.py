#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from cartopy import crs as ccrs

#F_Obs = 'microwave_d03_201708221300_so.txt'
F_Hxb_mean = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/HARVEY/MW_THO/Obs_Hx/MW/201708221200/enkf/wrf_enkf_input_d03_mean.tb.ssmis_f17.crtm.conv.txt'
F_Hxa_mean = '/work2/06191/tg854905/stampede2/Pro2_PSU_MW/HARVEY/MW_THO/Obs_Hx/MW/201708221200/enkf/wrf_enkf_output_d03_mean.tb.ssmis_f17.crtm.conv.txt'
#def SSMIS(F_Obs, F_Hxb_mean, F_Hxa_mean):
tmp = np.fromfile(F_Hxb_mean, sep=' ')
ncol = 8

Lat_all = tmp[0::ncol] #latitude
Lon_all = tmp[1::ncol] #longitude
Chnum_all = tmp[2::ncol] #channel
Yo_all  = tmp[3::ncol] #observedTB
Yb_all  = tmp[4::ncol] #BackgroundTB    

tmp = np.fromfile(F_Hxa_mean, sep=' ') 
Ya_all = tmp[4::ncol]

#f, ax=plt.subplots(2, 3, dpi=300) #, linewidth=0.5, sharex='all', sharey='all',  dpi=300)
f, ax=plt.subplots(2, 3, subplot_kw={'projection': ccrs.PlateCarree()}, linewidth=0.5, sharex='all', sharey='all',  dpi=300)

# customize colormap
max_T = 300
min_T = 80
min_Jet = 150
jetLength = max_T - min_Jet + 1
notJetLength = (min_Jet - 1) - (min_T + 1) + 1 

jet = cm.get_cmap('jet', max_T-min_T)
Myjet = jet(np.linspace(0,1,jetLength))

jet_red = Myjet[:,0]
jet_green = Myjet[:,1]
jet_blue = Myjet[:,2]

jetAdd_red = np.linspace(1.0, jet_red[0], notJetLength)
jetAdd_green = np.linspace(1.0, jet_green[0], notJetLength)
jetAdd_blue  = np.linspace(1.0, jet_blue[0], notJetLength) 

cm_red = np.concatenate([np.insert(jetAdd_red,0,0),np.append(jet_red,0)])
cm_green = np.concatenate([np.insert(jetAdd_green,0,0),np.append(jet_green,0)])
cm_blue = np.concatenate([np.insert(jetAdd_blue,0,0),np.append(jet_blue,0)])

newJet_value = np.column_stack([cm_red,cm_green,cm_blue, np.ones(len(cm_red))])
newJet = ListedColormap(newJet_value)



ch_num = [13, 9]
scatter_size = [2.5, 2.5]

lat = np.array(Lat_all)
lon = np.array(Lon_all)
lat_min = np.amin(lat)
lat_max = np.amax(lat)
lon_min = np.amin(lon)
lon_max = np.amax(lon)

for i in range(2):
    idx = Chnum_all == ch_num[i]
    ax[i,0].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[i,0].coastlines(resolution='10m', color='black',linewidth=0.5)
    ax[i,0].scatter(Lon_all[idx],Lat_all[idx],scatter_size[i],c=Yo_all[idx],edgecolors='none', cmap=newJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())
    
    ax[i,1].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[i,1].coastlines(resolution='10m', color='black',linewidth=0.5)
    ax[i,1].scatter(Lon_all[idx],Lat_all[idx],scatter_size[i],c=Yb_all[idx],edgecolors='none', cmap=newJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree()) 

    ax[i,2].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax[i,2].coastlines(resolution='10m', color='black',linewidth=0.5)
    cs = ax[i,2].scatter(Lon_all[idx],Lat_all[idx],scatter_size[i],c=Ya_all[idx],edgecolors='none', cmap=newJet, vmin=min_T, vmax=max_T, transform=ccrs.PlateCarree())

#plt.xlim([-91,-84])
#plt.ylim([15,24])

# Colorbar
caxes = f.add_axes([0.2, 0.97, 0.4, 0.02])
cbar = f.colorbar(cs, orientation="horizontal", cax=caxes)
cbar.ax.tick_params(labelsize=5)
#plt.text( 0.8, 0.7, 'Brightness Temperature/K', fontsize=4)


#subplot title
font = {'size':6,}

ax[0,0].set_title('Yo', font, fontweight='bold')
ax[0,1].set_title('HXb', font, fontweight='bold')
ax[0,2].set_title('HXa', font, fontweight='bold')

plt.subplots_adjust(wspace=0, hspace=0)

plt.savefig('201708221200.png', dpi=300)
