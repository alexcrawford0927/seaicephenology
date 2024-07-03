#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 22 Apr 2024
Modified: 03 Jul 2024

Author: Alex Crawford

Purpose: Bar chart of regional contributions with reference map
"""

import matplotlib.pyplot as plt
import pandas as pd
import CycloneModule_13_2 as md
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.gridspec import GridSpec

inpath = '/home/alex/OneDrive/CMIP6 OWP/Sea Ice Extent Expansion Pauses/'
# inpath = '/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses/'
eventpath = inpath+'/Pauses_combined_final_bytime_276-400.csv'
outpath = inpath+"/Figures/Fig Regional Contributions and Frequnecy of SIE Expansion Pauses V4 bytime.png"

width = 0.5
labels1 = np.array(['Nordic', 'E. Canada', 'Pacific', 'Other'])
labels2 = np.array(['Barents', 'Kara', 'E. Greeland', 'Baffin', 'Labrador', 'St. Lawrence', 'Hudson', 'Chukchi', 'Bering', 'Okhotsk','Other'])

# Load Events
pdf = pd.read_csv(eventpath)
pdf = pdf.loc[(pdf['Use'] == 1) & (pdf['SYear'] > 1978) & (pdf['DOY'] < 369) ]
pdf = pdf.reset_index()

# Load Contexual SIE Data
siedf2 = pd.read_csv(inpath+"/SIE_Daily_Regional2_bytime.csv")
siedf3 = pd.read_csv(inpath+"/SIE_Daily_Regional3_bytime.csv")

fig = plt.figure(figsize=(7.5,5.5))
gs = GridSpec(2,4, figure=fig)

### Overall Regions ###
ax1 = fig.add_subplot(gs[0:,:2])
ax1.grid(axis='x', linestyle='dashed', linewidth=0.5)

# Extract SIE change for each region and event
ylabel, nor, can, pac, cao, ok, bgck, oth = [[] for i in range(8)]
for i in range(len(pdf)):
    year, month, day, jd = pdf.loc[i,['year','month','day','JulianDay']].astype(int)
    newdate = md.timeAdd([year,month,day],[0,0,-6])
    
    siesub = siedf3.loc[siedf3['JulianDay'] == jd]
    nor.append(siesub['BKG'].values[0]/1000000), can.append(siesub['Can'].values[0]/1000000) 
    pac.append(siesub['BgCkOk'].values[0]/1000000), cao.append(siesub['CAO'].values[0]/1000000)
    bgck.append( siesub['BgCk'].values[0]/1000000), ok.append( pac[-1] - bgck[-1] )
    oth.append( pdf.loc[i,'cdr_delsie_5day'] - (nor[-1] + pac[-1] + can[-1]) )
    ylabel.append( md.dd[newdate[2]-1] + " " + md.mmm[newdate[1]-1] + ' ' + str(newdate[0]) )

counts = np.array((nor,can,pac,oth))

# Calculate the "bottom" of each block to plot using cumulatie sums of negative and
# positive values separately

# Identify pos & neg
pos = np.where(counts > 0, counts, 0)
neg = np.where(counts < 0, counts, 0)

# Take cumulative sums
poscum = np.cumsum(pos,axis=0)
negcum = np.cumsum(neg,axis=0)

# Shift cumulatives so that they start at zero and don't inlcude the ending point
poscum1 = np.zeros(pos.shape)
negcum1 = np.zeros(neg.shape)
poscum1[1:] = poscum[:-1]
negcum1[1:] = negcum[:-1]

# Combine pos and neg bottoms
bottom = np.where(counts > 0, poscum1, negcum1)

for i in range(counts.shape[0]):
    ax1.barh(ylabel, counts[i,:], width, label=labels1[i], left=bottom[i,:])

ax1.vlines(x=0, ymin=-0.5,ymax=len(ylabel)-0.5, color='k')
ax1.set_ylim(ymin=-0.5, ymax=len(ylabel)-0.5)
ax1.tick_params(axis='both', which='major', labelsize=7)
ax1.set_xlabel('6-day SIE Change (millions km$^2$)', size=8)
ax1.set_title('a. Sectoral Contributions to Arctic SIE\nDuring 6-Day Pauses in Expansion', size=8, weight='bold')
ax1.set_axisbelow(True)
ax1.invert_yaxis()

### Frequency of Going Negative ###
neg2 = siedf2.loc[:,['region','Neg']].groupby(by=['region']).mean()*100
neg2 = neg2.reset_index()
neg2 = neg2.loc[np.in1d(neg2['region'],['Kara','Barents','E. Greenland','Baffin','Labrador','Hudson','Chukchi','Bering','Okhotsk'])]
neg2['Order'] = np.array([4,2,9,8,3,7,1,5,10])
neg2 = neg2.sort_values(by=['Order'])

neg3 = (siedf3 < 0).mean() *100
neg3['NorCan'] = np.mean((siedf3['Can'] < 0) & (siedf3['BKG'] < 0)) * 100
neg3['NorPac'] = np.mean((siedf3['BgCkOk'] < 0) & (siedf3['BKG'] < 0)) * 100
neg3['CanPac'] = np.mean((siedf3['Can'] < 0) & (siedf3['BgCkOk'] < 0)) * 100


ax3 = fig.add_subplot(gs[:1,2:])
ax3.grid(axis='y', linestyle='dashed', linewidth=0.5)

ax3.set_title('b. Distribution of 6-day Sectoral SIE Change', size=8, weight='bold')
bplot = ax3.boxplot(siedf3.loc[(siedf3['month'] >= 10) & np.isfinite(siedf3['BKG']+siedf3['Can']+siedf3['BgCkOk']),['BKG','Can','BgCkOk']].values/1000000,
                    vert=False, patch_artist=True, whis=(5,95), sym='')

colors = ['#0076ba','#ed8032','#4aa12f']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
for median in bplot['medians']:
    median.set_color('black')
ax3.set_axisbelow(True)
ax3.set_yticklabels(['Nordic','E. Canada','Pacific'], size=8)
ax3.set_xlabel('6-day Sectoral SIE Change (millions '+r'km$^2$'+')',size=8)
ax3.tick_params(axis='both', which='major', labelsize=7)
ax3.vlines(x=0, ymin=0.5,ymax=3.5, color='k')
ax3.set_ylim(ymin=0.5, ymax=3.5)
ax3.invert_yaxis()

ax4 = fig.add_subplot(gs[1:,2:3])
ax4.grid(axis='y', linestyle='dashed', linewidth=0.5)
ax4.set_title('.',color='w')

ax4.bar(('Nordic','E. Canada','Pacific'), 
        (neg3.loc[['BKG','Can','BgCkOk']].values), 
        color=['#0076ba','#ed8032','#4aa12f'])
ax4.tick_params(axis='both', which='major', labelsize=7)
ax4.set_ylabel('Frequency of Pauses\nin SIE Expansion (%)', size=8)
ax4.set_xticklabels(('Nordic','E. Canada','Pacific', 'Nor & Can'), rotation=45, ha='right')
ax4.set_yticks(np.arange(0,20,5))
ax4.set_axisbelow(True)

### Study Area Map ###
import xarray as xr
import cartopy.crs as ccrs

prjpath = '/media/alex/Datapool/Projections/psn_projection.nc'
# prjpath = '/Volumes/Miranda/Projections/psn_projection.nc'
prjnc = xr.open_dataset(prjpath)
prj = ccrs.NorthPolarStereo(central_longitude=-45,true_scale_latitude=70)

arr = prjnc['reg0780'][:].data
arr, lons, lats = prjnc['reg0780'][:].data, prjnc['lon'][:].data, prjnc['lat'][:].data
arr[(lons > -30) & (lons < 20) & (arr == 1)] = 8
arr[(lons > 20) & (lons < 65) & (arr == 1)] = 7
arr[(lons > 65) & (lons < 90) & (arr == 1)] = 6

arr[(arr == 0) | (arr == 35)] = 50
arr[(arr >= 6) & (arr <= 8)] = 51
arr[(arr == 9) | (arr == 10) | (arr == 11) | (arr == 19)] = 52
arr[(arr == 3) | (arr == 13) | (arr == 14)] = 53
arr[(arr > 0) & (arr < 30)] = 54
arr[(arr < 50)] = 55

# arr[(arr == 3) | (arr == 13) | (arr == 14) | (arr == 15) | (arr == 18)] = 3

cmap = ListedColormap(["lightskyblue","#0076ba", "#ed8032", "#4aa12f", "#e8001b", "0.7"])
ax5 = fig.add_subplot(gs[1:,3:], projection=prj)
ax5.pcolormesh(prjnc['x'], prjnc['y'], arr, cmap=cmap, transform=prj)

plt.tight_layout()

ax4.annotate('c. Frequency of 6-day Pauses\nin Expansion of Sectoral SIE', xy=(0.5,1.05), xycoords='axes fraction',size=8, weight='bold',ha='center')
ax5.set_title('d. Sector\nDefinitions', weight='bold', size=8)

ax1.legend(fontsize=7,loc='upper left',bbox_to_anchor=(0.75,0.55), framealpha=1)


plt.savefig(outpath,dpi=300)


# from scipy.stats import spearmanr
# siedf3 = siedf3.loc[np.isfinite(siedf3['Can'])]

# ## Atl + Pac ## --> 1.8%
# spearmanr(siedf3['Atl'],siedf3['Pac'])

# ## Atl + Can ## --> 0.2%
# spearmanr(siedf3['Atl'],siedf3['Can'])

# ## Pac + Can ## --> 1.1%
# spearmanr(siedf3['Pac'],siedf3['Can'])

# ## All ## --> 0.0%
# np.mean((siedf3['Can'] < 0) & (siedf3['Pac'] < 0) & (siedf3['Atl'] < 0)) * 100


### Sub Regions ###
# ax2 = fig.add_subplot(2,2,2)
# ax2.grid(axis='x', linestyle='dashed', linewidth=0.5, zorder=1)

# ver = 'cdr'
# # Extract SIE change for each region and event
# ylabel, bar, kar, gre, baf, lab, hud, gsl, chk, ber, okh, oth = [[] for i in range(12)]
# for i in range(len(pdf)):
#     year, month, day, jd = pdf.loc[i,['year','month','day','JulianDay']].astype(int)
    
#     siesub = siedf2.loc[siedf2['JulianDay'] == jd]
#     bar.append(siesub.loc[siesub['region'] == 'Barents',ver+'_delsie_5day'].values[0]/1000000)
#     kar.append(siesub.loc[siesub['region'] == 'Kara',ver+'_delsie_5day'].values[0]/1000000)
#     gre.append(siesub.loc[siesub['region'] == 'E. Greenland',ver+'_delsie_5day'].values[0]/1000000)
#     lab.append(siesub.loc[siesub['region'] == 'Labrador',ver+'_delsie_5day'].values[0]/1000000)
#     baf.append(siesub.loc[siesub['region'] == 'Baffin',ver+'_delsie_5day'].values[0]/1000000)
#     gsl.append(siesub.loc[siesub['region'] == 'St. Lawrence',ver+'_delsie_5day'].values[0]/1000000)
#     hud.append(siesub.loc[siesub['region'] == 'Hudson',ver+'_delsie_5day'].values[0]/1000000)
#     chk.append(siesub.loc[siesub['region'] == 'Chukchi',ver+'_delsie_5day'].values[0]/1000000)
#     ber.append(siesub.loc[siesub['region'] == 'Bering',ver+'_delsie_5day'].values[0]/1000000)
#     okh.append(siesub.loc[siesub['region'] == 'Okhotsk',ver+'_delsie_5day'].values[0]/1000000)

    
#     oth.append( pdf.loc[i,ver+'_delsie_5day'] - (bar[-1] + kar[-1] + gre[-1] + lab[-1] + baf[-1] + gsl[-1] + hud[-1] + chk[-1] + ber[-1] + okh[-1]) )
#     ylabel.append( md.dd[day-1] + " " + md.mmm[month-1] + ' ' + str(year) )

# counts = np.array((bar,kar,gre,baf,lab,gsl,hud,chk,ber,okh,oth))

# # Calculate the "bottom" of each block to plot using cumulative sums of negative and
# # positive values separately

# # Identify pos & neg
# pos = np.where(counts > 0, counts, 0)
# neg = np.where(counts < 0, counts, 0)

# # Take cumulative sums
# poscum = np.cumsum(pos,axis=0)
# negcum = np.cumsum(neg,axis=0)

# # Shift cumulatives so that they start at zero and don't inlcude the ending point
# poscum1 = np.zeros(pos.shape)
# negcum1 = np.zeros(neg.shape)
# poscum1[1:] = poscum[:-1]
# negcum1[1:] = negcum[:-1]

# # Combine pos and neg bottoms
# bottom = np.where(counts > 0, poscum1, negcum1)

# for i in range(counts.shape[0]):
#     ax2.barh(ylabel, counts[i,:], width, label=labels2[i], left=bottom[i,:])

# ax2.vlines(x=0, ymin=-0.5,ymax=len(ylabel)-0.5, color='k')
# ax2.set_ylim(ymin=-0.5, ymax=len(ylabel)-0.5)
# ax2.legend(fontsize=7,loc='lower left',bbox_to_anchor=(1.02,0.02), framealpha=1)
# ax2.tick_params(axis='both', which='major', labelsize=7)


# siann = siedf2.loc[:,['SYear','Neg','region']].groupby(by=['region','SYear']).sum().reset_index()
# siann['count'] = siedf2.loc[:,['SYear','region','Neg']].groupby(by=['region','SYear']).count()['Neg'].values
# siann['Negrel'] = siann['Neg'] / siann['count'] * 100

# siann.loc[siann['SYear'] < 1987,'Negrel'] = siann.loc[siann['SYear'] < 1987,'Negrel']*2

# from scipy.stats import mannwhitneyu

# for reg in ['Kara','Barents','E. Greenland','Baffin','Labrador','Hudson','Chukchi','Bering','Okhotsk']:
#     mannwhitneyu(siann.loc[(siann['SYear'] > 1987) & (siann['SYear'] < 2005) & (siann['region'] == reg),'Negrel'],
#                  siann.loc[(siann['SYear'] > 2004) & (siann['SYear'] < 2024) & (siann['region'] == reg),'Negrel'])
#     np.median(siann.loc[(siann['SYear'] > 1987) & (siann['SYear'] < 2005) & (siann['region'] == reg),'Negrel'])
#     np.median(siann.loc[(siann['SYear'] > 2004) & (siann['SYear'] < 2024) & (siann['region'] == reg),'Negrel'])
    
    
#     mannwhitneyu(siann.loc[(siann['SYear'] > 1978) & (siann['SYear'] < 2005) & (siann['region'] == reg),'Negrel'],
#                  siann.loc[(siann['SYear'] > 2004) & (siann['SYear'] < 2024) & (siann['region'] == reg),'Negrel'])
#     np.median(siann.loc[(siann['SYear'] > 1978) & (siann['SYear'] < 2005) & (siann['region'] == reg),'Negrel'])
#     np.median(siann.loc[(siann['SYear'] > 2004) & (siann['SYear'] < 2024) & (siann['region'] == reg),'Negrel'])

