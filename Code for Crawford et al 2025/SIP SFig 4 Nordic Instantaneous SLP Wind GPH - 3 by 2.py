#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 11:26:35 2024

@author: alex
"""

import pandas as pd
import xarray as xr
import xesmf as xe
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import CycloneModule_13_2 as md
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import os

'''*******************************************
Declare Variables
*******************************************'''
# Path Variables
path = '/media/alex/Datapool'
sicpath = path+'/SeaIce/G02202_V4/Daily'

outpath = '/home/alex/OneDrive/CMIP6 OWP/Sea Ice Extent Expansion Pauses/Figures'

dver, tver, subset, typ = '13_2E5', 'R', 'BBox10', 'Cyclone'
cycpath = path+"/CycloneTracking/tracking"+dver+tver+"/"+subset+"/"+typ+"Tracks"
cfpath = path+"/CycloneTracking/detection"+dver+"/CycloneFields"

e5path1 = path+'/ERA5/SLP'
e5path2  = path+'/ERA5/HorizontalWind_10m'
e5path3 = path+'/ERA5/Geopotential_MidLevel'
# e5path4 = path+'/ERA5/ArcticAmpMidLatLinkages/Blocking/ERA5_Blocking_AnomalyMethod_500hPa_1979-2020.nc'

eraprjpath = path+'/ERA5/DailyClimatology_1981-2010/ERA5_SLP_DailyClimatology_1981-2010.nc'
psnprjpath = path+'/Projections/psn_projection.nc'
# ease2prjpath = path+'/Projections/EASE2_N0_25km_Projection_uv.nc'

# Cyclone Variables
# tids = np.array([624,661,673,696,740]) # 1990 # [255,301,349,360,381] # 1983 # [390,399] # 2010 #
# cyclabels = np.arange(1,len(tids)+1)

# Date Variables
# times = [ [2015,12,26,0,0,0], [2015,12,28,0,0,0], [2015,12,30,0,0,0]]
# times = [ [2016,11,13,0,0,0], [2016,11,15,0,0,0], [2016,11,17,0,0,0]]
times = [ [1999,12,26,0,0,0], [1999,12,28,0,0,0], [1999,12,30,0,0,0]]

# times = [ [1984,11,25,0,0,0], [1984,11,27,0,0,0], [1984,11,29,0,0,0]]
# times = [ [1982,12,26,0,0,0], [1982,12,28,0,0,0], [1982,12,30,0,0,0]]
# times = [ [1980,12,16,0,0,0], [1980,12,18,0,0,0], [1980,12,20,0,0,0]]
# times = [ [1990,12,25,0,0,0], [1990,12,28,0,0,0], [1990,12,29,12,0,0] ]
dateref = [1900,1,1,0,0,0]

# Plotting Variables
prj = ccrs.NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)
extent = [-2800000,3500000,-4400000,700000] # Left, Right, Bottom, Top
cbbox = [55,90,-140,45]
# spres = 25000 # meters

prjsic = ccrs.NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)

colorbarlabel = ['Sea-Level Pressure (hPa)', 'Geopotential Height (m) at 500 hPa']
cmap = [ plt.cm.PuOr, plt.cm.PRGn_r]

slpbreaks = np.arange(95400,106000,200)
slpcbbreaks = np.arange(95400,106000,800)
slplabels = np.arange(954,1060,8)

gphbreaks = np.arange(4700,6100,60)*9.81
gphcbbreaks = np.arange(4700,6100,180)*9.81
gphlabels = np.arange(4700,6100,180)

'''*******************************************
Main Analysis
*******************************************'''

# Load projections
# ease2prj = xr.open_dataset(ease2prjpath)
era5prj = xr.open_dataset(eraprjpath)
psnprj = xr.open_dataset(psnprjpath)
regridder = xe.Regridder(era5prj, psnprj, 'bilinear')

# ext = (np.array(extent)/spres).astype(int)+360

# Start figure
fig = plt.figure(figsize=(7.5,5.5))

for i in range(3):
    # Define time
    time = times[i]
    time2 = md.timeAdd(time, [0,1,0,0,0,0])
    timeref = md.daysBetweenDates(dateref, time)
    timestr = str(time[0])+"-"+md.dd[time[1]-1]+"-"+md.dd[time[2]-1] + " " + md.hhmm[time[3]]

    # Define Axes
    ax1 = fig.add_subplot(2, 3, i+1, projection=prj)
    ax1.add_feature(cfeature.COASTLINE, zorder=6)
    ax1.set_extent(extent, prj)
    ax1.set_title(timestr + ' UTC', fontsize=8, weight='bold')
    ax1.annotate(md.abc[i],xy=(0.02,0.98), xycoords='axes fraction', ha='left', va='top', size=8, weight='bold')

    ax2 = fig.add_subplot(2, 3, i+4, projection=prj)
    ax2.add_feature(cfeature.COASTLINE, zorder=6)
    ax2.set_extent(extent, prj)
    ax2.annotate(md.abc[i+3],xy=(0.02,0.98), xycoords='axes fraction', ha='left',va='top',size=8, weight='bold')
    ax2.gridlines(xlocs=range(-180,180,45), ylocs=range(0,90,15),
                          draw_labels = False, linestyle = 'dashed', linewidth = 0.5, 
                          color = '0.3')
    # Plot SLP    
    xr1 = xr.open_dataset(e5path1+"/ERA5_SLP_Hourly_"+str(time[0])+md.dd[time[1]-1]+".nc")
    plotter = xr1.sel(time=timestr)['msl']
    cf1 = ax1.contourf(era5prj['longitude'], era5prj['latitude'], plotter, 
                       levels=slpbreaks, cmap=cmap[0], extend='both', transform=ccrs.PlateCarree())

    # Plot Wind
    xr2 = xr.open_dataset(e5path2+"/ERA5_HorizontalWind_10m_3h_"+str(time[0])+md.dd[time[1]-1]+".nc")
    uvpsn = regridder(xr2.sel(time=timestr))
    urot, vrot = md.rotateCoordsAroundOrigin(uvpsn['u10'], uvpsn['v10'], (psnprj['lon']-(-45))*np.pi/-180)
    vec = ax1.quiver(psnprj['x'], psnprj['y'], urot.data, vrot.data, transform=prjsic, regrid_shape=21, width=0.005, scale=300)

    # Plot GPH
    # xr3 = xr.open_dataset(e5path3+"/ERA5_Geopotential_MidLevel_3h_"+str(time[0])+md.dd[time[1]-1]+".nc")
    # plotter = xr3.sel(time=timestr).sel(level=500)['z']
    # cf3 = ax2.contourf(era5prj['longitude'], era5prj['latitude'], plotter, 
    #                    levels=gphbreaks, cmap=cmap[1], extend='both', transform=ccrs.PlateCarree())


    # Plot SIC
    sicnc = xr.open_dataset(sicpath+"/"+str(time[0])+"/"+[f for f in md.listdir(sicpath+"/"+str(time[0])) if str(time[0])+md.dd[time[1]-1]+md.dd[time[2]-1] in f][0])
    sicarr = sicnc['nsidc_nt_seaice_conc'][0,:,:].data
    sicarr[(sicarr > 1) | (psnprj['reg0780'] > 20) | (psnprj['reg0780'] < 0)] = np.nan
    
    ax1.contour(psnprj['x'], psnprj['y'], sicarr, levels=[0.15], transform=prjsic, colors='red', linewidths=1.5)
    # ax2.contour(psnprj['x'], psnprj['y'], sicarr, levels=[0.15], transform=prj, colors='gold', linewidths=1.5)


    # Plot Cyclone Tracks & Centers
    trs = pd.read_pickle(cycpath+"/"+str(time[0])+"/"+subset+typ.lower()+"tracks"+str(time[0])+md.dd[time[1]-1]+".pkl")
    trs2 = pd.read_pickle(cycpath+"/"+str(time2[0])+"/"+subset+typ.lower()+"tracks"+str(time2[0])+md.dd[time2[1]-1]+".pkl")
    # cfs = pd.read_pickle(cfpath+"/CF"+str(time[0])+md.dd[time[1]-1]+".pkl")
    # cf = np.array(cfs)[timeref == np.array([cf.time for cf in cfs])][0]
    
    # cflons, cflats, cftypes = np.array(cf.lons()), np.array(cf.lats()), np.array(cf.centerType())
   
    # ax1.plot(cflons,cflats,
    #          marker='o', markerfacecolor='white', markeredgecolor='k',
    #          transform=ccrs.Geodetic(), linewidth=0, zorder=10)
    # ax2.plot(cflons,cflats,
    #          marker='o', markerfacecolor='white', markeredgecolor='k',
    #          transform=ccrs.Geodetic(), linewidth=0, zorder=10)  

    for tr in trs+trs2:
        # Plot the storm track if it exists at this time
        if timeref in tr.data.time.values:
            lonx = tr.data.loc[tr.data['time'] == timeref,'lon'].values[0]
            laty = tr.data.loc[tr.data['time'] == timeref,'lat'].values[0]
            
            # ax1.plot(tr.data.lon, tr.data.lat, color='blue', alpha=0.5,
            #          linewidth=1.3,transform=ccrs.Geodetic())
            ax1.plot(lonx, laty, marker='o', markerfacecolor='yellow', 
                     markeredgecolor='k', transform=ccrs.Geodetic(), zorder=10)
            # ax1.plot(tr.data.iloc[0]['lon'],tr.data.iloc[0]['lat'],
            #          marker='x', color='blue',
            #          transform=ccrs.Geodetic(), zorder=10)
            
            if ((lonx > cbbox[2]) & (lonx < cbbox[3])) & (laty > cbbox[0]) & (laty < cbbox[1]):            
                ax2.plot(tr.data.iloc[0]['lon'],tr.data.iloc[0]['lat'],
                         marker='x', color='blue',
                         transform=ccrs.Geodetic(), zorder=10)
                ax2.plot(tr.data.lon, tr.data.lat, color='blue', alpha=0.5,
                         linewidth=1.3,transform=ccrs.Geodetic())
            
            ax2.plot(lonx, laty, marker='o', markerfacecolor='yellow', 
                     markeredgecolor='k', transform=ccrs.Geodetic(), zorder=10)
        
        # if tr.tid in tids:
        #     txt = ax2.annotate('  '+str(cyclabels[np.where(tids == tr.tid)[0]][0]), xy=(lonx,laty),
        #                  transform=ccrs.Geodetic(), size=7)
        #     txt.set_path_effects([pe.withStroke(linewidth=3, foreground='w')])

pb = 0.09 # Bottom of plottable area as figure fraction
plt.tight_layout(rect=[0.03,pb,1,1])

# Annotation and Legends
ax1.annotate('Sea-Level Pressure & 10-m Wind', xy=(0.01,1-((1-pb)/4)), xycoords='figure fraction',
            ha='left',va='center',size=8,weight='bold',rotation=90)
ax2.annotate('   500-hPa Geopotential Height', xy=(0.01, pb+((1-pb)/4) ), xycoords='figure fraction',
            ha='left',va='center',size=8,weight='bold',rotation=90)

ax1.plot([], color='red', linewidth=1.5, zorder=0, label='Sea Ice Edge')
ax1.plot([], marker='o', color='w', markerfacecolor='yellow', markeredgecolor='k',label='Cyclone Centers')
ax2.plot([], color='blue', marker='x', label = 'Cyclone Genesis')
ax2.plot([],color='blue',linewidth=1.3,label='Cyclone Tracks')
leg = ax1.legend(fontsize=7,loc='lower left',bbox_to_anchor=(-2,-1.38), framealpha=0, ncol=2)
leg2 = ax2.legend(fontsize=7,loc='lower left',bbox_to_anchor=(-0.4,-0.18), framealpha=0, ncol=2)

ax1.quiverkey(vec, 0.475, 0.14, 20, r'20 $\frac{m}{s}$', coordinates='figure', labelpos='E', fontproperties={'size':7})
cbar_ax = fig.add_axes([0.02, 0.10, 0.44, 0.025])
cb = fig.colorbar(cf1,cbar_ax,orientation='horizontal', ticks=slpcbbreaks)
cb.ax.set_xticklabels(slplabels, fontsize=7, rotation=45, ha='center')
cb.set_label(colorbarlabel[0], fontsize=7)

# cbar_ax2 = fig.add_axes([0.52, 0.10, 0.44, 0.025])
# cb2 = fig.colorbar(cf3,cbar_ax2,orientation='horizontal', ticks=gphcbbreaks)
# cb2.ax.set_xticklabels(gphlabels, fontsize=7, rotation=45, ha='center')
# cb2.set_label(colorbarlabel[1], fontsize=7)

plt.savefig(outpath+"/SFig6_InstantaneousSLP_Wind_GPH_"+str(times[0][0])+md.dd[times[0][1]-1]+".png", dpi=300)
