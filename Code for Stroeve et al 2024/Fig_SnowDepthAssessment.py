#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 10:28:00 2023

@author: acrawfora
"""

'''*******************************************
Set up Modules
*******************************************'''
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import CycloneModule_13_2 as md
import matplotlib.patheffects as pe
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

def flatten(l):
    '''Courtesy Donnie Cameron'''
    if not isinstance(l, list):
        return [l]
    flat = []
    for sublist in l:
        flat.extend(flatten(sublist))
    return flat

'''*******************************************
Define Variables
*******************************************'''
exp = 'historical-ssp585'
YY = '1981-2021'
obsmod = 'SnowModel-LG'
var, obvar = 'sisnthick', 'snod'
multiplier = 100 # for m --> cm
varpath = var#+"-noJB"
tcoords = [2,4]
mon = 4

regnames =['regHB','regHB3']
regs = [[40],[41,42]] # [4,44,45,46] # [4,47,48,49] #  [4,41,42,43] # [4,40,41,42]
REGS = [['HB'],['SHB','WHB']] # ['HBC','HS','FB','Nar','NW','W','Cen','S','E','JB'] #

# Path Variables
obpath = '/Volumes/Miranda/SeaIce/SnowModel-LG/SM_snod_MERRA2_EASE2_MonthlyClimatology_1980-2021.nc'
c6path1 = '/Volumes/Cassandra/CMIP6'
csvpath = c6path1+"/RegionalStats/RawVariables/Avg"+YY
ease2path = "/Volumes/Miranda/Projections/EASE2_N0_25km_Projection.nc"
figpath = '/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Hudson Bay Figures/ThicknessPhenologyFigures/'

# Plotting Variables
prjdata = ccrs.LambertAzimuthalEqualArea(central_longitude=0, central_latitude=90)
prjplot = ccrs.LambertAzimuthalEqualArea(central_longitude=-80, central_latitude=90)
bbox = [-1000000,500000,-2700000,-4400000] # left, right, top, bottom

levs = [6,9,12,15,18,21,24]

lines = [(63,-91,63,-88),(63,-88,56,-88),(60,-88,60,-77)]

maplabel = ['SnowModel-LG', 'CMIP6 Multi-Model Mean', '', 'CMIP6 at 2°C', 'CMIP6 at 4°C']
ylabel = 'April snow depth (cm)'
styles = pd.read_csv('/Users/acrawfora/OneDrive - University of Manitoba/CMIP6 OWP/CMIP6_plot_styles.csv')

'''*******************************************
Plotting
*******************************************'''

# Load NetCDF Files for Maps
ease2 = xr.open_dataset(ease2path)
dso = xr.open_dataset(obpath)
ds1 = xr.open_dataset(c6path1+'/Climatology/ByMonth_'+YY+'/'+var+'/MultiModelMean_'+var+'_MonthlyClimatology_ssp585_'+YY+'.nc')
ds2 = xr.open_dataset(c6path1+"/data/bias-corrected/CMIP6_bias-corrected-by-t2mamly_"+exp+"_baselineyears_"+YY+"_"+var+"_multimodelmean.nc")

# Load Obs for scatter plot
obdf = pd.DataFrame()

for c in range(0,len(regs)):
    # All
    obdf1 = pd.read_csv(csvpath+'/'+varpath+'_ObsAvg-'+md.mmm[mon-1]+'-'+YY+'_'+regnames[c]+'.csv')
    obdf1 = obdf1.loc[np.in1d(obdf1['Region'],regs[c])]
    obdf1['REG'] = np.array( [ np.array(REGS[c])[np.where(np.array(regs[c]) == r)[0]][0] for r in list(obdf1['Region']) ] )
    obdf = obdf.append(obdf1, ignore_index=True)

obdf = obdf.sort_values(by=['Region'])

# Load CMIP6 Average
c6df = pd.read_csv(csvpath+'/'+varpath+'_CMIP6Avg-'+md.mmm[mon-1]+'-'+YY+'_'+regnames[0]+'_'+exp+'.csv')
c6df = c6df.loc[np.in1d(c6df['Region'],regs[0])]
c6df['REG'] = np.array( [ np.array(REGS[0])[np.where(np.array(regs[0]) == r)[0]][0] for r in list(c6df['Region']) ] )

for c in range(1,len(regs)):
    c6df1 = pd.read_csv(csvpath+'/'+varpath+'_CMIP6Avg-'+md.mmm[mon-1]+'-'+YY+'_'+regnames[c]+'_'+exp+'.csv')
    c6df1 = c6df1.loc[np.in1d(c6df1['Region'],regs[c])]
    c6df1['REG'] = np.array( [ np.array(REGS[c])[np.where(np.array(regs[c]) == r)[0]][0] for r in list(c6df1['Region']) ] )
    c6df = c6df.append(c6df1, ignore_index=True)

c6df = c6df[((c6df['Family'] != 'CESM2') & (c6df['Member'] == 1)) | ((c6df['Family'] == 'CESM2') & (c6df['Member'] == 4))]
c6df = c6df.sort_values(by=['Region','Family'])

c6avgdf = pd.read_csv(csvpath+'/'+varpath+'_CMIP6Avg-'+md.mmm[mon-1]+'-'+YY+'_MultiModelMean_'+regnames[0]+'_'+exp+'.csv')
c6avgdf = c6avgdf[np.in1d(c6avgdf['Region'], regs[0])]

for c in range(1,len(regs)):
    c6avgdf1 = pd.read_csv(csvpath+'/'+varpath+'_CMIP6Avg-'+md.mmm[mon-1]+'-'+YY+'_MultiModelMean_'+regnames[c]+'_'+exp+'.csv')
    c6avgdf1 = c6avgdf1[np.in1d(c6avgdf1['Region'], regs[c])]
    c6avgdf1['REG'] = np.array( [ np.array(REGS[c])[np.where(np.array(regs[c]) == r)[0]][0] for r in list(c6avgdf1['Region']) ] )
    c6avgdf = c6avgdf.append(c6avgdf1, ignore_index=True)

c6avgdf = c6avgdf.sort_values(by=['Region'])

# Load Max Internal Variability Average
c6sedf = pd.read_csv(csvpath+'/'+varpath+'_CMIP6Avg-'+md.mmm[mon-1]+'-'+YY+'_MultiModelMaxStdDev_'+regnames[0]+'_'+exp+'.csv')
c6sedf = c6sedf[np.in1d(c6sedf['Region'], regs[0])]

for c in range(1,len(regs)):
    c6sedf1 = pd.read_csv(csvpath+'/'+varpath+'_CMIP6Avg-'+md.mmm[mon-1]+'-'+YY+'_MultiModelMaxStdDev_'+regnames[c]+'_'+exp+'.csv')
    c6sedf1 = c6sedf1[np.in1d(c6sedf1['Region'], regs[c])]
    c6sedf1['REG'] = np.array( [ np.array(REGS[c])[np.where(np.array(regs[c]) == r)[0]][0] for r in list(c6sedf1['Region']) ] )
    c6sedf = c6sedf.append(c6sedf1, ignore_index=True)

c6sedf = c6sedf.sort_values(by=['Region'])

##### PLOTTING #####
fig = plt.figure(figsize=(7.5,5.75))
axs = []

### Observed Map ###
arr0 = np.flipud( dso.sel({'month':mon})[obvar].data )
arr0[ease2['z'] > 0] = np.nan

axs.append(fig.add_subplot(2,3,1, projection=prjplot) )
axs[0].add_feature(cfeature.COASTLINE)
axs[0].set_extent(bbox, crs=prjplot)
axs[0].contourf(ease2['x'], ease2['y'], np.where(np.isfinite(arr0), np.nan, 1), [0,1,2,3,4], cmap=plt.cm.gray_r,transform=prjdata)
axs[0].contourf(ease2['x'], ease2['y'], arr0*multiplier, levs, cmap=plt.cm.viridis, transform=prjdata, extend='both')

### CMIP6 Multi-Model Mean ###
arr1 = ds1.sel({'month':mon})[var].data
arr1[ease2['z'] > 0] = np.nan

axs.append( fig.add_subplot(2,3,2, projection=prjplot) )
axs[1].add_feature(cfeature.COASTLINE)
axs[1].set_extent(bbox, crs=prjplot)
axs[1].contourf(ease2['x'], ease2['y'], np.where(np.isfinite(arr1), np.nan, 1), [0,1,2,3,4], cmap=plt.cm.gray_r,transform=prjdata)
cf1 = axs[1].contourf(ease2['x'], ease2['y'], arr1*multiplier, levs, cmap=plt.cm.viridis, transform=prjdata, extend='both')

### Scatter Plot ###
axs.append( fig.add_subplot(2,3,3) )

# Internal Variability
axs[-1].errorbar(np.concatenate(REGS), obdf[var]*multiplier, yerr = 2*c6sedf[var+"_Avg"]*multiplier, fmt='o', ms=0.1, color='k', capsize=7, capthick=1.5, elinewidth=2, zorder=4)

# Obsservational Mean
axs[-1].scatter(np.concatenate(REGS),obdf[var]*multiplier,fc='w',ec='k', marker='o', s=60, zorder=5)

# Multi-Model Mean
axs[-1].scatter(np.arange(len(np.concatenate(REGS)))-0.05,c6avgdf[var+"_Avg"]*multiplier,fc='r',ec='k', marker='*', s=130, zorder=3)

# First Ensemble Members
mods = np.unique(c6df['Family'])
xvals = np.arange(-0.1,0.1,(0.1--0.1)/len(mods))

for m in range(len(mods)):
    axs[-1].scatter(x = xvals[m]+np.arange(len(np.concatenate(REGS))),
        y = c6df[c6df['Family'] == mods[m]][var+"_Avg"]*multiplier,zorder=2,
        color = list(styles.loc[styles['Model'] == mods[m]]['Color'])[0],
        marker = list(styles.loc[styles['Model'] == mods[m]]['Shape'])[0])

# Y Label
axs[-1].set_ylabel(ylabel, size=8)
axs[-1].tick_params(axis = 'y', labelsize = 8)
axs[-1].tick_params(axis = 'x', labelsize = 8)

### Temperatures ###
for ti in range(len(tcoords)):
    arr = ds2.sel({'tanom':tcoords[ti]})[var].data
    arr[ease2['z'] > 0] = np.nan
    
    axs.append( fig.add_subplot(2,3,4+ti, projection=prjplot) )
    axs[-1].add_feature(cfeature.COASTLINE)
    axs[-1].set_extent(bbox, crs=prjplot)
    axs[-1].contourf(ease2['x'], ease2['y'], np.where(np.isfinite(arr), np.nan, 1), [0,1,2,3,4], cmap=plt.cm.gray_r,transform=prjdata)
    axs[-1].contourf(ease2['x'], ease2['y'], arr*multiplier, levs, cmap=plt.cm.viridis, transform=prjdata, extend='both')

### Add Annotaton ###
for i in range(len(axs)):
    if i == 2:
        txt = axs[i].annotate(md.abc[i]+". "+maplabel[i],xy=(0.02,0.98),xycoords='axes fraction',ha='left',va='top',weight='bold',size=8)
    else:
        txt = axs[i].annotate(md.abc[i]+". "+maplabel[i],xy=(0.02,0.02),xycoords='axes fraction',ha='left',va='bottom',weight='bold',size=8)
    txt.set_path_effects([pe.withStroke(linewidth=3, foreground='w')])
    
    if i != 2:
        
        # Gridlines
        gdl = axs[i].gridlines(xlocs=range(-180,180,15),ylocs=range(0,90,5),linestyle='dashed',linewidth=0.5,draw_labels=True)
        gdl.top_labels = False
        gdl.right_labels = False
        # if i > 2:
        #     gdl.bottom_labels = False
        gdl.xformatter = LongitudeFormatter()
        gdl.yformatter = LatitudeFormatter()
        gdl.xlabel_style = {'size':8}
        gdl.ylabel_style = {'size':8}
        
        # Add Study Area Lines
        for l in range(len(lines)):
            if lines[l][0] == lines[l][2]:
                axs[i].plot( np.arange(np.min([lines[l][1],lines[l][3]]),np.max([lines[l][1],lines[l][3]])+0.5,0.5), np.repeat(lines[l][0],np.ceil(1+(abs(lines[l][3]-lines[l][1]))/0.5)), color='w', transform=ccrs.Geodetic())
            else:
                axs[i].plot(np.repeat(lines[l][1],np.ceil(1+(abs(lines[l][2]-lines[l][0]))/0.5)), np.arange(np.min([lines[l][2],lines[l][0]]),np.max([lines[l][2],lines[l][0]])+0.5,0.5), color='w', transform=ccrs.Geodetic())
        
        # Labels
        txt = [ axs[i].annotate('SHB', xy=(0.50,0.50), xycoords='axes fraction',va='center',ha='center',size=8),
               axs[i].annotate('WHB', xy=(0.25,0.65), xycoords='axes fraction',va='center',ha='center',size=8)]
        [ t.set_path_effects([pe.withStroke(linewidth=2, foreground='w')]) for t in txt ]
        
plt.tight_layout(rect=[0,0.07,1,1])

# Color bar
cbar_ax1 = fig.add_axes([0.02, 0.07, 0.40, 0.015]) # xloc, yloc, xsize, ysize
cbar1 = fig.colorbar(cf1,cbar_ax1,orientation='horizontal')
cbar1.ax.tick_params(labelsize=8)
cbar1.set_label(ylabel,size=8)

### First Legend ###
leg = plt.axes([0.67, 0.17, 0.19, 0.3])
leg.set_ylim(ymin=3.5,ymax=8)
leg.set_xlim(xmin=0.9, xmax=2.2)

leg.errorbar(x=1, y=4.5, yerr=0.5, fmt='o', ms=0.1, color='k', capsize=7,zorder=1)
leg.annotate('Internal\nVariability', xy=(1.3,4.5), ha='left', va='center', size=7)

leg.scatter(x=1, y=7, fc='w', ec='k', marker='o', s=60)
leg.annotate(obsmod, xy=(1.2,7), ha='left', va='center', size=7)

leg.scatter(x=1, y=6, fc='r', ec='k', marker='*', s=110)
leg.annotate('CMIP6 Multi-\nModel Mean', xy=(1.2,6), ha='left', va='center', size=7)

# Remove axes from legend
leg.spines['bottom'].set_color('w')
leg.spines['top'].set_color('w')
leg.spines['right'].set_color('w')
leg.spines['left'].set_color('w')
leg.set_xticks([]), leg.set_yticks([])

### Second Legend ###
leg2 = plt.axes([0.84,0.04,0.19,0.46]) # ll x, ll y, width, height
leg2.set_xticks( [] )
leg2.set_yticks( [] )
plt.setp(leg2.spines.values(), color='w')

xlocs = np.repeat(0,int(len(mods)))
ylocs = np.arange(int(len(mods)),0,-1)
leg2.set_xlim(xmin=-0.1,xmax=2)

for m, mod in enumerate(mods):
    leg2.scatter(x = xlocs[m],
                y = ylocs[m],
                color = list(styles.loc[styles['Model'] == mods[m]]['Color'])[0],
                marker = list(styles.loc[styles['Model'] == mods[m]]['Shape'])[0] )
    leg2.annotate(mod, xy=(xlocs[m]+0.14,ylocs[m]),xycoords="data",size=7,va='center')

plt.savefig(figpath+"/Fig_SnowDepth_ComparetoObs_and_Future-bias-corrected_V2.png",dpi=300)
