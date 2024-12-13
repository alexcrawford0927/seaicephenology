#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 16:39:09 2023

Purpose: Create a figure of SLP and GPH during each of a set of events
"""
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import CycloneModule_13_2 as md
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

sicpath = '/Volumes/Miranda/SeaIce/G02202_V4/Daily'

eventpath = '/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses/Pauses_combined_final_bytime_276-400.csv'
eventpath2 = '/Users/telekineticturtle/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Pauses/Pauses_combined_final_bytime_276-400.csv'

sicoutpath = '/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses/Figures/Fig_SeaIceExtentExpansionPauses_V4.png'

sicvar = 'cdr_seaice_conc'
sicmin = 0.15

prj = ccrs.NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)
cmap = ListedColormap(["#C4E9FF", "magenta", "darkgreen", "white", "0.8"])

########

try:
    events = pd.read_csv(eventpath)
except: 
    events = pd.read_csv(eventpath2)
    
events = events.loc[(events['Use'] == 1) & (events['DOY'] < 369)]

sifiles = md.listdir(sicpath)

#### SEA ICE FIGURE ####
fig = plt.figure(figsize=(8,8))
axs = []

for row in range(len(events)):
    enddate = list(events.iloc[row][['year','month','day']].astype(int))
    startdate = md.timeAdd(enddate,[0,0,-6])

    # Load start and end date
    sicstt = xr.open_dataset(sicpath+"/"+str(startdate[0])+"/"+[f for f in md.listdir(sicpath+"/"+str(startdate[0])) if str(startdate[0])+md.dd[startdate[1]-1]+md.dd[startdate[2]-1] in f][0])
    sicend = xr.open_dataset(sicpath+"/"+str(enddate[0])+"/"+[f for f in md.listdir(sicpath+"/"+str(enddate[0])) if str(enddate[0])+md.dd[enddate[1]-1]+md.dd[enddate[2]-1] in f][0])

    sicsttarr = sicstt[sicvar][0,:].data.astype(float)
    sicendarr = sicend[sicvar][0,:].data.astype(float)
    
    # Identify extent on each date and shared
    sicsttarr[sicsttarr > 1] = np.nan
    sicendarr[sicendarr > 1] = np.nan
    
    sie = np.where( (sicsttarr >= sicmin) & (sicendarr >= sicmin), 3, np.where(sicendarr >= sicmin, 2, np.where(sicsttarr >= sicmin, 1, np.where(np.isfinite(sicsttarr*sicendarr), 0, 5))))
    xarr, yarr = np.meshgrid(sicstt['xgrid'][:].data, sicstt['ygrid'][:].data)

    # Plot result
    axs.append( fig.add_subplot(3,5,row+1, projection=prj) )
    cf = axs[row].pcolormesh(xarr, yarr, sie, cmap=cmap, transform=prj)    
    axs[-1].set_title(md.abc[row]+". "+str(startdate[0])+"-"+md.dd[startdate[1]-1]+"-"+md.dd[startdate[2]-1] + "\nto " + str(enddate[0])+"-"+md.dd[enddate[1]-1]+"-"+md.dd[enddate[2]-1], size=8, weight='bold')

plt.tight_layout(rect=[0,0,1,1])


cbar_ax = fig.add_axes([0.62, 0.14, 0.35, 0.015])
cbar1 = fig.colorbar(cf,cbar_ax,orientation='horizontal', ticks=[0.5,1.5,2.5,3.5,4.5])
cbar1.ax.set_xticklabels(['open water','ice loss','ice gain','ice cover','mask'],fontsize=7, ha='center', va='top')

plt.savefig(sicoutpath,dpi=300)

