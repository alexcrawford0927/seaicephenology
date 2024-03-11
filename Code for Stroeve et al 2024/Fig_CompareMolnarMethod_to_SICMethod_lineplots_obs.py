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
outpath = '/Users/acrawfora/OneDrive - University of Manitoba/CMIP6 OWP/Hudson Bay Figures/ThicknessPhenologyFigures/Fig_SIEMethod_v_GridcellMethod.png'

ymin, ymax = 1979, 2021
regvar = 'regHB3'
reg = 42

# Variables for SIC method
minc_ret = [10,10] # in %
minc_adv = minc_ret # in %

# Variables for SIE method
ct, et = 30, 30

ylabs = ['Ice-Free (Fasting) Period','','Retreat (Onshore) Day', 'Advance (Offshore) Day']

# Paths for other inputs
concmods = ['NASATeam','Bootstrap'] # ['Bootstrap','Bootstrap','Bootstrap','PIOMAS','PIOMAS','PIOMAS']
prjpath = '/Volumes/Miranda/Projections/psn_projection.nc'
colors = ['#D4382A','#A578F0'] #['#A578F0','#297B3C','#00c9b2']#  ['#D4382A','#297B3C','#FF8D01']
linestyles = ['solid','dashed']
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

'''*******************************************
Plotting
*******************************************'''
print("Plotting")

fig = plt.figure(figsize=(179/25.4,129/25.4))

for i in range(len(concmods)):
    
    ### OPC GRAPH ###
    ax1 = fig.add_subplot(2,2,1)
    ax1.plot(conclist[i]['Year'], conclist[i]['OPCavg'], c=colors[i], linestyle='solid')
    ax1.plot(molnarlist[i]['Year'], molnarlist[i]['OPC'], c=colors[i], linestyle='dashed')
    
    ### LRD GRAPH ###
    ax3 = fig.add_subplot(2,2,3)
    ax3.plot(conclist[i]['Year'], conclist[i]['LRDavg'], c=colors[i], linestyle='solid')
    ax3.plot(molnarlist[i]['Year'], molnarlist[i]['LRD'], c=colors[i], linestyle='dashed')
       
    ### FAD GRAPH ###
    ax4 = fig.add_subplot(2,2,4)
    ax4.plot(conclist[i]['Year'], conclist[i]['FADavg'], c=colors[i], linestyle='solid')
    ax4.plot(molnarlist[i]['Year'], molnarlist[i]['FAD'], c=colors[i], linestyle='dashed')

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
leg = plt.axes([0.8, 0.6, 0.19, 0.1])
leg.annotate(str(minc_ret[0])+'% SIC',xy=(1.35,3.5),va='bottom',ha='left',size=7,rotation=45)
# leg.annotate('Molnár et al. (2020)',xy=(2.35,3.5),va='bottom',ha='left',size=7,rotation=45)
leg.annotate('SIE ('+str(et)+'%)',xy=(2.35,3.5),va='bottom',ha='left',size=7,rotation=45)
leg.annotate(concmods[0],xy=(3,3),va='center',ha='left',size=7)
leg.annotate(concmods[-1],xy=(3,1),va='center',ha='left',size=7)
leg.plot([1.0,1.7],[3,3], c=colors[0], linewidth=lwd, linestyle=linestyles[0])
leg.plot([1.0,1.7],[1,1], c=colors[-1], linewidth=lwd, linestyle=linestyles[0])
leg.plot([2.0,2.7],[3,3], c=colors[0], linewidth=lwd, linestyle=linestyles[-1])
leg.plot([2.0,2.7],[1,1], c=colors[-1], linewidth=lwd, linestyle=linestyles[-1])
leg.set_ylim(ymin=0.5,ymax=4)
leg.set_xlim(xmin=0.8,xmax=5)
leg.set_yticks([])
leg.set_xticks([])
for edge in ['top','bottom','left','right']:
    leg.spines[edge].set_visible(False)

# plt.savefig(outpath,dpi=300)
