#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 26 Apr 2024
Modified: 03 Jul 2024

Author: Alex Crawford

Purpose: Basic count of SIE pause events and timing of their occurrence
"""

'''*******************************************
Load Modules
*******************************************'''
print(' ***Load Modules***')

import xarray as xr
import numpy as np
import pandas as pd
import CycloneModule_13_2 as md
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

'''*******************************************
Declare Variables
*******************************************'''
print(' ***Define Variables***')

inpath = '/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses/'
sicpath = '/Volumes/Miranda/SeaIce/G02202_V4/Daily'
sicpath2 = '/Volumes/Miranda/SeaIce/G10016_V2/Daily'
siepath = inpath+'/SIE_Daily_ArcticWide.csv'

outpath = inpath+"/Figures/Fig1_BasicCounts_and_Maps_bytime.png"

dmin, dmax = 276, 400
dmax2 = 369
ymin, ymax = 1979, 2022

ncvar = 'cdr_seaice_conc' # 'nsidc_nt_seaice_conc' # 
ncmin, ncmax = 0, 1
sicthresh = 0.15

colors=['sienna','brown','maroon','red','darkorange','gold','green','teal','blue','purple','violet','hotpink']
prj = ccrs.NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)
cmap = LinearSegmentedColormap.from_list('seaice',colors=[(0/256,0/256,0/256,1),(32/256,59/256,93/256,1),(87/256,123/256,192/256),(163/256,176/256,198/256,1),(256/256,256/256,256/256,1)],N=256)
doys = [1,32,60,91,121,152,182,213,244,274,305,335,366,397,425]
doyl = ["Jan 1","Feb 1","Mar 1","Apr 1","May 1","Jun 1","Jul 1","Aug 1","Sep 1","Oct 1","Nov 1","Dec 1","Jan 1","Feb 1","Mar 1"]

'''*******************************************
Process Variables
*******************************************'''
print(' ***Process Variables***')
siedfall = pd.read_csv(inpath+"/Pauses_combined_final_bytime_"+str(dmin)+"-"+str(dmax)+".csv")
siedfall['SMonth'] = np.where(siedfall['day'] > 3, siedfall['month'], siedfall['month']-1)
siedfall['SMonth'] = np.where(siedfall['SMonth'] < 1, 12, siedfall['SMonth'])
siedfmon = siedfall.loc[:,['SMonth','LowVal']].groupby(by=['SMonth']).count().reset_index()
siedfyear = siedfall.loc[:,['SYear','LowVal']].groupby(by=['SYear']).count().reset_index()

# xvals = [md.daysBetweenDates([1900,1,1], [y,10,1]) for y in np.arange(1980,2025,5)]
# xvallabs = ['' + str(y) for y in np.arange(1980,2025,5)]
xvals = np.arange(1980,2025,5)
xvallabs = xvals
yvals = [274,305,335,366,397]
ylab = ['Oct 1','Nov 1','Dec 1','Jan 1','Feb 1']

### Calculate the % of years that have sea ice present for each grid cell ###
siestartarrs, sieendarrs = [], []
for y in range(1979,2023):
    
    try:
        xrminfile = [ f for f in md.listdir(sicpath+'/'+str(y)) if '_'+str(y)+md.dd[10-1] in f][0]
        xrmin = xr.open_dataset( sicpath+'/'+str(y)+'/'+xrminfile)
    except:
        xrminfile = [ f for f in md.listdir(sicpath2+'/'+str(y)) if '_'+str(y)+md.dd[10-1] in f][0]
        xrmin = xr.open_dataset( sicpath2+'/'+str(y)+'/'+xrminfile)

    xrminarr = xrmin[ncvar][0,:,:].data
    xrminbool = (xrminarr >= sicthresh).astype(float)
    xrminbool[(xrminarr < ncmin) | (xrminarr > ncmax)] = np.nan
    
    siestartarrs.append( xrminbool )
    
    try:
        xrmaxfile = [ f for f in md.listdir(sicpath+'/'+str(y+1)) if '_'+str(y+1)+md.dd[1-1] in f][0]
        xrmax = xr.open_dataset( sicpath+'/'+str(y+1)+'/'+xrmaxfile)
    except:
        xrmaxfile = [ f for f in md.listdir(sicpath2+'/'+str(y+1)) if '_'+str(y+1)+md.dd[1-1] in f][0]
        xrmax = xr.open_dataset( sicpath2+'/'+str(y+1)+'/'+xrmaxfile)
        
    xrmaxarr = xrmax[ncvar][0,:,:].data    
    xrmaxbool = (xrmaxarr >= sicthresh).astype(float)
    xrmaxbool[(xrmaxarr < ncmin) | (xrmaxarr > ncmax)] = np.nan
    sieendarrs.append( xrmaxbool )  

probsie_start = np.mean(siestartarrs, axis=0)*100
probsie_end = np.mean(sieendarrs, axis=0)*100

# Percentiles
sie = pd.read_csv(siepath)
sie['SYear'] = np.where(sie['month'] > 8, sie['year'], sie['year']-1)
sie = sie.loc[sie['SYear'] > 1978]
sie['DOY'] = [md.daysBetweenDates([sie.loc[i,'SYear'],1,0],[sie.loc[i,'year'],sie.loc[i,'month'],sie.loc[i,'day']]) for i in sie.index]

### For each year, identify the day for which the SIE is closest to the siemin that is after the annual minimum extent and
# the year for which the SIE is closet to the siemax that is before the annual maximum extent ###
sie = sie.loc[((sie['month'] > 8) | (sie['month'] < 4)) &  (np.isfinite(sie['sie_nt']))]

sie_mean = sie.loc[:,['SYear','DOY','sie_cdr','cdr_delsie_5day']].groupby(by=['DOY']).mean()
sie_p50 = sie.loc[:,['SYear','DOY','sie_cdr','cdr_delsie_5day']].groupby(by=['DOY']).quantile(0.5)
sie_p05 = sie.loc[:,['SYear','DOY','sie_cdr','cdr_delsie_5day']].groupby(by=['DOY']).quantile(0.05)
sie_p95 = sie.loc[:,['SYear','DOY','sie_cdr','cdr_delsie_5day']].groupby(by=['DOY']).quantile(0.95)
sie_p25 = sie.loc[:,['SYear','DOY','sie_cdr','cdr_delsie_5day']].groupby(by=['DOY']).quantile(0.25)
sie_p75 = sie.loc[:,['SYear','DOY','sie_cdr','cdr_delsie_5day']].groupby(by=['DOY']).quantile(0.75)

sie_mean = sie_mean.reset_index()
sie_p50 = sie_p50.reset_index()
sie_p05 = sie_p05.reset_index()
sie_p95 = sie_p95.reset_index()
sie_p75 = sie_p75.reset_index()
sie_p25 = sie_p25.reset_index()

siedfall = siedfall.loc[siedfall['SYear'] > 1978]
sip = siedfall.loc[(siedfall['Use'] == 1) & (siedfall['DOY'] < dmax2)]
sipyears = np.unique(sip['SYear'])

# SIE at the point of pausing SIE growth
siepausearrs = []
for i in sip.index:
    row = sip.loc[i,:]
    
    try:
        xrfile = [ f for f in md.listdir(sicpath+'/'+str(int(row['year']))) if '_'+str(int(row['year']))+md.dd[int(row['month'])-1]+md.dd[int(row['day'])-1] in f][0]
        xrf = xr.open_dataset( sicpath+'/'+str(int(row['year']))+'/'+xrfile)
    except:
        xrfile = [ f for f in md.listdir(sicpath2+'/'+str(int(row['year']))) if '_'+str(int(row['year']))+md.dd[int(row['month'])-1]+md.dd[int(row['day'])-1] in f][0]
        xrf = xr.open_dataset( sicpath2+'/'+str(int(row['year']))+'/'+xrminfile)

    xrarr = xrf[ncvar][0,:,:].data
    xrbool = (xrarr >= sicthresh).astype(float)
    xrbool[(xrarr < ncmin) | (xrarr > ncmax)] = np.nan
    
    siepausearrs.append( xrbool )

probsie_pause = np.mean(siepausearrs, axis=0)*100

'''*******************************************
Plotting
*******************************************'''
print('***Plotting***')

fig = plt.figure(figsize=(7.5,6.00))
gs = GridSpec(8,4, figure=fig)

ax = fig.add_subplot(gs[1:4,0:4])
ax.grid(axis='x', linestyle='dashed', linewidth=0.5)
ax.set_axisbelow(True)
ax.set_xlim(xmin=siedfall.iloc[0]['SYear']-1, xmax=siedfall.iloc[-1]['SYear']+5)
ax.set_ylim(ymin=dmin-2, ymax=dmax+2)
# ax.vlines(x=xvals[-4]-0.5,ymin=yvals[0]-0.125,ymax=yvals[-1]+0.125, color='k', linewidth=0.6, linestyle='dashed',zorder=1)
ax.hlines(y=yvals,xmin=siedfall.iloc[0]['SYear']-1,xmax=siedfall.iloc[-1]['SYear'],color='0.8', linewidth=0.6, linestyle='dashed',zorder=1)

ax.set_xticks(xvals)
ax.set_xticklabels(xvallabs, size=7, rotation=0, ha='center', va='top')
ax.set_yticks(yvals)
ax.set_yticklabels(ylab, size=7)
ax.invert_yaxis()
   
ax.set_ylabel('Day of Year', size=8)
titl = ax.annotate('b. All SIE Expansion Pauses', xy=(0.02,0.98), xycoords='axes fraction', va='top', weight='bold', size=8)
titl.set_path_effects([pe.withStroke(linewidth=3, foreground='w')])

for r in range(len(siedfall)):
    x = siedfall.iloc[r].SYear
    y = siedfall.iloc[r].DOY
    
    ax.plot([x,x],[y,y-1], linewidth=2, color='b')
    
ax.vlines(x=[1987.5,2004.5],ymin=dmin-2, ymax=dmax+2, color='k', linestyle='dashed',linewidth=0.6)

# Monthly Totals
for i, m in enumerate([10,11,12,1]):
    print(m)
    if m in siedfmon['SMonth'].values:
        text = siedfmon.loc[siedfmon['SMonth']==m,'LowVal'].values[0]
    else:
        text = 0
        
    ax.annotate(text, xy=(siedfall.iloc[-1]['SYear']+2,yvals[i]+15),
                size=8,weight='bold', va='center')
ax.annotate('# Days per\nMonth', xy=(0.95,-0.01), xycoords='axes fraction', va='top', ha='center', size=7)

# Annual Totals
ax2 = fig.add_subplot(gs[0:1,0:4])
ax2.grid(axis='x', linestyle='dashed', linewidth=0.5)
ax2.set_axisbelow(True)
ax2.bar(siedfyear['SYear'],siedfyear['LowVal'], color='darkgreen')
ax2.set_ylim(ymin=0,ymax=max(siedfyear['LowVal'])+1)
ax2.set_xlim(xmin=siedfall.iloc[0]['SYear']-1, xmax=siedfall.iloc[-1]['SYear']+5)
ax2.set_xticks(xvals)
ax2.set_xticklabels([])
ax2.set_yticks([0,5])
ax2.set_yticklabels([0,5],size=7)
ax2.hlines(y=[5],xmin=siedfall.iloc[0]['SYear']-1,xmax=siedfall.iloc[-1]['SYear'],color='0.8', linewidth=0.6, linestyle='dashed',zorder=0)
ax2.vlines(x=[1987.5,2004.5],ymin=0, ymax=max(siedfyear['LowVal'])+1, color='k', linestyle='dashed',linewidth=0.6)
ax2.set_ylabel('# Days', size=8)
titl = ax2.annotate('a. Total Pause Days per Year', xy=(0.02,0.98), xycoords='axes fraction', va='top', weight='bold', size=8)
titl.set_path_effects([pe.withStroke(linewidth=3, foreground='w')])


# By Day of Year
ax5 = fig.add_subplot(gs[4:8,0:2])
ax5.plot(sie_p50['DOY'].values[3:-3],md.movingAverage(sie_p50['sie_cdr'].values,5)[3:-3],color='k')
ax5.fill_between(sie_p75['DOY'].values[3:-3], md.movingAverage(sie_p25['sie_cdr'].values,5)[3:-3], md.movingAverage(sie_p75['sie_cdr'].values,5)[3:-3], alpha=0.1, fc='k', ec=None)
ax5.fill_between(sie_p75['DOY'].values[3:-3], md.movingAverage(sie_p05['sie_cdr'].values,5)[3:-3], md.movingAverage(sie_p95['sie_cdr'].values,5)[3:-3], alpha=0.05, fc='k', ec=None)
ax5.vlines(x=366, ymin= 4, ymax=16, color='0.4', linewidth=1, linestyle='dashed')

# Indiviudal Years
for y in range(len(sipyears)):
    ax5.plot(sie.loc[sie['SYear'] == sipyears[y]]['DOY'],sie.loc[sie['SYear'] == sipyears[y]]['sie_cdr'], linewidth=1, alpha=0.5, color=colors[y], label=sipyears[y])
    
    sipsub = siedfall.loc[siedfall['SYear'] == sipyears[y]]
    
    for jd in sipsub['JulianDay']:
        siesub = sie.loc[ (sie['JulianDay'] >= jd-6) & (sie['JulianDay'] <= jd)]
        ax5.plot(siesub['DOY'],siesub['sie_cdr'], linewidth=1.4, alpha=1, color=colors[y])

ax5.tick_params(axis='both', which='major', labelsize=7)
ax5.set_ylabel('SIE (million km$^2$)', size=8, weight='bold')
ax5.annotate('c. Daily SIE during growth period', xy=(0.01,0.98), xycoords='axes fraction', ha='left', va='top', size=8, weight='bold')
ax5.set_xticks(doys)
ax5.set_xticklabels(doyl)
ax5.set_xlim(xmin=274,xmax=397)

plt.tight_layout()

# Typical Sea Ice Extent at end of 6-day pause #
ax3 = fig.add_axes([0.52,0.10,0.22,0.35], projection=prj)
ax3.set_facecolor('0.4')
ax3.add_feature(cfeature.LAND, facecolor='0.4',zorder=7)
ax3.contourf(xrmin['xgrid'], xrmin['ygrid'], probsie_pause, levels=np.arange(0,110,10), cmap=cmap, transform=prj)
ax3.annotate('d. After Pause in\nSIE Growth (Oct-Dec)', xy=(0.02,0.98), xycoords='axes fraction', 
              weight='bold', size=7, ha='left', va='top', color='w',zorder=12)


# Typical Sea Ice Extent at start and end of "growth period" #
ax4 = fig.add_axes([0.76,0.10,0.22,0.35], projection=prj)
ax4.set_facecolor('0.4')
ax4.add_feature(cfeature.LAND, facecolor='0.4',zorder=7)
pc4 = ax4.contourf(xrmax['xgrid'], xrmax['ygrid'], probsie_end, levels=np.arange(0,110,10), cmap=cmap, transform=prj)
ax4.contour(xrmin['xgrid'],xrmin['ygrid'], probsie_start, [50], colors='b', linewidths=1 , transform=prj, zorder=8)
ax4.annotate('e. On January 1', xy=(0.02,0.98), xycoords='axes fraction', 
              weight='bold', size=7, ha='left', va='top', color='w', zorder=12)
ax4.plot((1,2),(0,0), linewidth=1, color='b', zorder=1, label=' 50% Probability On Oct 1')

# Color Bar
cbar_ax = fig.add_axes([0.55, 0.06, 0.4, 0.025])
cbar1 = fig.colorbar(pc4,cbar_ax,orientation='horizontal')
cbar1.ax.tick_params(labelsize=7)
ax4.annotate('Probability of Sea Ice Presence (%)',xy=(0.75,0.48), xycoords='figure fraction', ha='center',va='top',weight='bold',fontsize=8)

# Legend
ax5.legend(loc='lower right',bbox_to_anchor=(0.97,0.01),fontsize=7, edgecolor='w', framealpha=0,ncols=3)
ax4.legend(loc='lower right',bbox_to_anchor=(1.04,-0.04),fontsize=7, edgecolor='w',facecolor='w', framealpha=1).set_zorder(16)





plt.savefig(outpath,dpi=300)
        



### STATS ###
# from scipy import stats
# import pymannkendall

# ### Count for Oct-Jan ###

# sy = pd.DataFrame({'SYear':np.arange(1979,2023)})
# sy['LowVal'] = 0
# for y in siedfyear['SYear']:
#     sy.loc[sy['SYear'] == y,'LowVal'] = siedfyear.loc[siedfyear['SYear'] == y,'LowVal'].values[0]
# sy['LowVal2'] = np.where(sy['SYear'] < 1988, sy['LowVal'] * 2, sy['LowVal'])

# stats.mannwhitneyu(sy.loc[sy['SYear']>2004,'LowVal'], sy.loc[(sy['SYear']>1987) & (sy['SYear']<2005),'LowVal'])


# md.theilsenslope(sy.loc[(sy['SYear']>1987)]['SYear'].values,sy.loc[(sy['SYear']>1987)]['LowVal'].values)
# md.theilsenslope(sy['SYear'].values,sy['LowVal'].values)
# md.theilsenslope(sy['SYear'].values,sy['LowVal2'].values)

# stats.linregress(sy['SYear'].values,sy['LowVal'].values)

# ## DOY ##

# siend = siedfall.loc[(siedfall['month'] > 9) & (siedfall['Use'] == 1)]
# siend.loc[:,['SYear','DOY','month','day']]
# stats.mannwhitneyu([357,355,336,364,365,364],[352,354,312,354,324])

# siend = siedfall.loc[(siedfall['month'] > 9)]
# stats.mannwhitneyu(siend.loc[siend['SYear']>2004,'DOY'],siend.loc[siend['SYear']<2005,'DOY'])
    
# stats.mannwhitneyu(siedfall.loc[(siedfall['SYear']>2004) & (siedfall['Use'] == 1),'DOY'], siedfall.loc[(siedfall['SYear']<2005) & (siedfall['Use'] == 1),'DOY'])
# stats.mannwhitneyu(siedfall.loc[siedfall['SYear']>2004,'DOY'],siedfall.loc[siedfall['SYear']<2005,'DOY'])

# ## Count for N-D ##
# synd = pd.DataFrame({'SYear':np.arange(1979,2023)})
# synd['LowVal'] = 0
# for y in siend['SYear']:
#     synd.loc[synd['SYear'] == y,'LowVal'] = len(siend.loc[siend['SYear'] == y])
# synd['LowVal2'] = np.where(synd['SYear'] < 1988, synd['LowVal'] * 2, synd['LowVal'])


# md.theilsenslope(synd['SYear'].values,synd['LowVal'].values)
# md.theilsenslope(synd['SYear'].values,synd['LowVal2'].values)
# stats.linregress(synd['SYear'].values,synd['LowVal'].values)

