#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 09:21:52 2023

@author: acrawfora
"""

import pandas as pd
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy.stats import linregress
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.patches import Rectangle
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import CycloneModule_13_2 as md

'''*******************************************
Declare Variables
*******************************************'''
print("Declaring Variables")

# Variables for all inputs
path = '/Volumes/Miranda/SeaIce'
# pbpath = '/Users/telekineticturtle/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Polar Bear Paper/PolarBearTelemetry_HB.csv'
# outpath = '/Users/telekineticturtle/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Hudson Bay Figures/ThicknessPhenologyFigures/Fig_SIEMethod_v_GridcellMethod_PMWAvg.png'

pbpath = '/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Polar Bear Paper/PolarBearTelemetry_HB.csv'
outpath ='/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Hudson Bay Figures/ThicknessPhenologyFigures/Fig_SIEMethod_v_GridcellMethod_PMWAvg_V2.png'

ymin, ymax = 1979, 2021
regvar = 'regHB3'
reg = 42

# Variables for SIC method
minc_ret = [10,10] # in %
minc_adv = minc_ret # in %

# Variables for SIE method
ct, et = 30, 30

ylabs = ['Period (days)','','Day of Year', 'Day of Year']

# Paths for other inputs
concmods = ['NASATeam','Bootstrap'] # ['Bootstrap','Bootstrap','Bootstrap','PIOMAS','PIOMAS','PIOMAS']
prjpath = '/Volumes/Miranda/Projections/psn_projection.nc'
colors = ['magenta','darkgreen','gold'] #['#D4382A','#A578F0','#00c9b2'] #['#A578F0','#297B3C','#00c9b2']#  ['#D4382A','#297B3C','#FF8D01']
linestyles = ['solid','solid','solid']
lwd = 1.2

prjplot = ccrs.LambertAzimuthalEqualArea(central_longitude=-80, central_latitude=90)
bbox = [-1050000,500000,-2600000,-4400000] # left, right, top, bottom
lines = [(63,-91,63,-88),(63,-88,56,-88)] # ,(60,-88,60,-77)

'''*******************************************
Pre-Processing
*******************************************'''
print("Processing Inputs")
years = np.arange(ymin,ymax+1)
YY = str(ymin)+"-"+str(ymax)

conclist = []
for i in range(len(concmods)):
    pdf = pd.read_csv(path+'/'+concmods[i]+'/AdvanceRetreat2/C'+str(minc_ret[i])+'/'+concmods[i]+'_Regionalized_'+regvar+'_Historical_C'+str(minc_ret[i])+'.csv')
    conclist.append( pdf.loc[pdf['Region'] == reg] )

molnarlist = [] 
for i in range(len(concmods)):
    pdf = pd.read_csv(path+"/"+concmods[i]+"/AdvanceRetreat-Molnar/siphenology-Molnar_C"+str(ct)+"-E"+str(et)+"_"+concmods[i]+"_"+regvar+"_"+YY+".csv")
    molnarlist.append( pdf.loc[pdf['Region'] == reg] )

# Average by method
obdf = conclist[0].loc[:,('Region','Year')]
mldf = molnarlist[0].loc[:,('Region','Year')]
madf = molnarlist[0].loc[:,('Region','Year')]

for col in ('OPC','LRD','FAD'):
    obdf[col] = np.array([df[col+'avg'] for df in conclist]).mean(axis=0)
    madf[col] = np.array([df[col+'adj'] for df in molnarlist]).mean(axis=0)
    mldf[col] = np.array([df[col] for df in molnarlist]).mean(axis=0)

pbdf = pd.read_csv(pbpath)

'''*******************************************
Plotting
*******************************************'''
print("Plotting")

fig = plt.figure(figsize=(179/25.4,129/25.4))
    
### OPC GRAPH ###
ax1 = fig.add_subplot(2,2,1)
ax1.plot(obdf['Year'], obdf['OPC'], c=colors[0], linestyle=linestyles[0])
ax1.plot(mldf['Year'], mldf['OPC'], c=colors[1], linestyle=linestyles[1])
ax1.plot(madf['Year'], madf['OPC'], c=colors[2], linestyle=linestyles[2])
ax1.plot(pbdf['Year'], pbdf['OPC'], c='k', linestyle='solid', marker='o', markersize=3)

### LRD GRAPH ###
ax3 = fig.add_subplot(2,2,3)
ax3.plot(obdf['Year'], obdf['LRD'], c=colors[0], linestyle=linestyles[0])
ax3.plot(mldf['Year'], mldf['LRD'], c=colors[1], linestyle=linestyles[1])
ax3.plot(madf['Year'], madf['LRD'], c=colors[2], linestyle=linestyles[2])
ax3.plot(pbdf['Year'], pbdf['LRD'], c='k', linestyle='solid', marker='o', markersize=3)

### FAD GRAPH ###
ax4 = fig.add_subplot(2,2,4)
ax4.plot(obdf['Year'], obdf['FAD'], c=colors[0], linestyle=linestyles[0])
ax4.plot(mldf['Year'], mldf['FAD'], c=colors[1], linestyle=linestyles[1])
ax4.plot(madf['Year'], madf['FAD'], c=colors[2], linestyle=linestyles[2])
ax4.plot(pbdf['Year'], pbdf['FAD'], c='k', linestyle='solid', marker='o', markersize=3)

for a, ax in enumerate([ax1,'',ax3,ax4]):
    if a != 1:
        txt = ax.annotate(md.abc[a], xy=(0.03,0.98),xycoords='axes fraction',weight='bold',size=10,ha='left',va='top')
        txt.set_path_effects([pe.withStroke(linewidth=5, foreground='w')])
        
        ax.set_ylabel(ylabs[a], size=8)
        ax.tick_params(axis='x', labelsize=8)        
        ax.tick_params(axis='y', labelsize=8)        

plt.tight_layout()

### Study Area Map ###
mapax = plt.axes([0.55,0.6,0.2,0.34], projection = prjplot)
mapax.set_extent(bbox, crs = prjplot)
mapax.add_feature(cfeature.COASTLINE)
mapax.add_feature(cfeature.OCEAN)
mapax.add_feature(cfeature.LAND)

# Add Study Area Lines
for l in range(len(lines)):
    if lines[l][0] == lines[l][2]:
        mapax.plot( np.arange(np.min([lines[l][1],lines[l][3]]),np.max([lines[l][1],lines[l][3]])+0.5,0.5), np.repeat(lines[l][0],np.ceil(1+(abs(lines[l][3]-lines[l][1]))/0.5)), color='k', transform=ccrs.Geodetic())
    else:
        mapax.plot(np.repeat(lines[l][1],np.ceil(1+(abs(lines[l][2]-lines[l][0]))/0.5)), np.arange(np.min([lines[l][2],lines[l][0]]),np.max([lines[l][2],lines[l][0]])+0.5,0.5), color='k', transform=ccrs.Geodetic())

txt = mapax.annotate('b', xy=(0.02,0.98),xycoords='axes fraction',weight='bold',size=10,ha='left',va='top')
txt.set_path_effects([pe.withStroke(linewidth=5, foreground='w')])

txt = mapax.annotate('WHB', xy=(0.02,0.68),xycoords='axes fraction',weight='bold',size=9,ha='left',va='top')
txt.set_path_effects([pe.withStroke(linewidth=5, foreground='w')])

# Add gridlines
grid = mapax.gridlines(xlocs=range(-180,180,10),ylocs=range(50,80,5),
    draw_labels = True, linestyle = 'dashed', linewidth = 0.5, color = '0.3')
grid.top_labels = False
grid.right_labels = False
grid.xformatter = LONGITUDE_FORMATTER
grid.yformatter = LATITUDE_FORMATTER
grid.xlabel_style = {'size': 8, 'color': 'k'}
grid.ylabel_style = {'size': 8, 'color': 'k'}


### Add Legend for Graphs ###
leg = plt.axes([0.77, 0.7, 0.19, 0.1])
leg.annotate('SIC < 10%',xy=(2,4),va='center',ha='left',size=7)
leg.annotate('SIE < 30%',xy=(2,3),va='center',ha='left',size=7)
leg.annotate('SIE < 30%, Adjusted',xy=(2,2),va='center',ha='left',size=7)
leg.annotate('Polar Bear Telemetry',xy=(2,1),va='center',ha='left',size=7)

leg.plot([1.0,1.7],[4,4], c=colors[0], linewidth=lwd, linestyle=linestyles[0])
leg.plot([1.0,1.7],[3,3], c=colors[1], linewidth=lwd, linestyle=linestyles[1])
leg.plot([1.0,1.7],[2,2], c=colors[2], linewidth=lwd, linestyle=linestyles[2])
leg.plot([1.0,1.7],[1,1], c='k', linewidth=lwd, linestyle='solid', marker='o', markersize=3)
leg.set_ylim(ymin=0.5,ymax=4.5)
leg.set_xlim(xmin=0.8,xmax=5)
leg.set_yticks([])
leg.set_xticks([])
for edge in ['top','bottom','left','right']:
    leg.spines[edge].set_visible(False)

plt.savefig(outpath,dpi=300)
