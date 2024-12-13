#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 27 May 2024
Modified: 03 Jul 2024

Author: Alex Crawford

Purpose: Create bar plots showing the trends of each atmospheric variable
for each 6-day period and each region
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import CycloneModule_13_2 as md
import matplotlib.patheffects as pe

path = '/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses/'
inpath = path+'/ERA5_and_SST_Avgs_forbboxes_seasonal_extremecounts.csv'
outpath = path+"/Figures/Fig9_AtmosphericExtremeCounts_Seasonal_1979-2023.png"

varlist2 = ['seb','seb','msl','u10','u10','msl','v10','v10']
perclist2 = ['90','10','90','90','10','10','90','10']
varlist3 = ['Extreme Downward Anomalies\nSurface Energy Balance ','Extreme Upward Anomalies\nSurface Energy Balance','Extreme Positive Anomalies\nSea-Level Pressure',
            'Extreme Westerly Anomalies\n10-m Zonal Wind','Extreme Easterly Anomalies\n10-m Zonal Wind','Exteme Negative Anomalies\nSea-Level Pressure',
            'Extreme Southerly Anomalies\n10-m Meridional Wind','Extreme Northerly Anomalies\n10-m Meridional Wind']

region = ("Barents & Kara Seas", "East Greenland Sea", "Hudson Bay", "Bering Sea", "Sea of Okhotsk")
regionlabels = ("Barents &\nKara Seas", "E. Greenland\nSea", "Hudson\nBay", "Bering\nSea", "Sea of\nOkhotsk")
# doys = np.array([312, 324, 336, 352, 354, 355, 357, 364, 366])
colorsmap = ['#0076ba','#7980f7','#ed8032','#4aa12f','#2abf77']
prj = ccrs.NorthPolarStereo(central_longitude=-45,true_scale_latitude=70)

basep = 42 # the percentage of days after the breakpoint -- i.e., the baseline percentage for determining if the number
## of extreme values after the breakpoint is above/below expectations

'''#####################
Main Analysis
#####################'''

tdf = pd.read_csv(inpath)
tdf = tdf.loc[np.in1d(tdf['bbox'],region)]
tdf2 = tdf.groupby(by=['bbox','var']).mean().reset_index()
tdf2 = tdf2.loc[np.in1d(tdf2['var'],varlist2)]

tdf2['n90diff'] = tdf2['n90aft'] / (tdf2['n90aft'] + tdf2['n90bef']) * 100
tdf2['n10diff'] = tdf2['n10aft'] / (tdf2['n10aft'] + tdf2['n10bef'])  * 100
tdf2['sig90'] = np.where((tdf2['n90p'] < 5) | (tdf2['n90p'] > 95), 1, 0)
tdf2['sig10'] = np.where((tdf2['n10p'] < 5) | (tdf2['n10p'] > 95), 1, 0)


### FIGURE ###
fig = plt.figure(figsize=(7.5,8.5))

for v, var in enumerate(varlist2):
    
    ax = fig.add_subplot(3,3,v+1)
    ax.grid(axis='x', linestyle='dashed', linewidth=0.5)
    ax.set_axisbelow(True)
    
    aftcount = tdf2.loc[(tdf2['var'] == var),'n'+perclist2[v]+'diff'].values
    sigs = tdf2.loc[(tdf2['var'] == var),'sig'+perclist2[v]].values
    
    x = np.arange(len(region))  # the label locations
    
    # difference (reported as percentage of counts in the later period #
    ax.barh(x, np.array(aftcount)-basep, color=colorsmap, left=0)
    ax.scatter(np.array(aftcount)-basep, x,  marker='X', s=sigs*15, color='k')

    ax.set_title(md.abc[v]+". "+varlist3[v],weight='bold',size=8)
    
    if v%3 == 0:
        ax.set_yticks(np.arange(0,len(regionlabels)))
        ax.set_yticklabels(regionlabels,size=7)
    else:
        ax.set_yticks([])
    
    if v > 4:
        ax.set_xlabel('Percentage of Extremes After 2004', size=7)
        
    ax.vlines(x=0,ymin=-0.6,ymax=len(aftcount)-0.4,color='k', linewidth=1)
    ax.set_ylim(ymin=-0.6,ymax=len(aftcount)-0.4)
    ax.invert_yaxis()
    
    ax.set_xticks([-(basep-20),-(basep-30),-(basep-40),(50-basep),(60-basep),(70-basep),(80-basep)])
    ax.set_xticklabels([20,30,40,50,60,70,80],size=7)
    ax.set_xlim(xmin=-10-basep/2,xmax=10+basep/2)
    ax.tick_params(axis='both', which='major', labelsize=7)
plt.tight_layout()

### Study Area ###

mapax = fig.add_axes([0.70,0.05,0.3,0.25], projection=prj)
mapax.add_feature(cfeature.COASTLINE)
mapax.add_feature(cfeature.LAND, fc='0.85')
mapax.set_extent([-3800000,3400000,-4500000,4500000], prj)

bboxname = ['Barents & Kara Seas','East Greenland Sea','Bering Sea','Sea of Okhotsk']
bbox = [[72,83,20,80],[70,80,-15,10],[55,65,-85,-75],[56,64,-180,-159],[50,65,130,165]]
for b in range(len(bbox)):
    mapax.plot((bbox[b][2],bbox[b][3],bbox[b][3],bbox[b][2],bbox[b][2]),
               (bbox[b][0],bbox[b][0],bbox[b][1],bbox[b][1],bbox[b][0]),
               color=colorsmap[b],linewidth=2,
               transform=ccrs.Geodetic())
mapax.set_title(md.abc[v+1] + '. Bounding Box Definitions', weight='bold', size=8)

# Add labels
txts = [ mapax.annotate('Okhotsk',xy=[0.3,0.9],xycoords='axes fraction',size=7,weight='bold'),
        mapax.annotate('Bering',xy=[0.05,0.58],xycoords='axes fraction',size=7,weight='bold'),
        mapax.annotate('Barents\n& Kara',xy=[0.76,0.5],xycoords='axes fraction',size=7,weight='bold'),
        mapax.annotate('E. Greenland',xy=[0.6,0.25],xycoords='axes fraction',size=7,weight='bold'),
        mapax.annotate('Hudson',xy=[0.05,0.08],xycoords='axes fraction',size=7,weight='bold')]

[txt.set_path_effects([pe.withStroke(linewidth=2.5, foreground='w')]) for txt in txts]

plt.savefig(outpath, dpi=300)
