#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 16:29:08 2024

@author: acrawfora
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import CycloneModule_13_2 as md

doymax = 369 # Maximum day of year for the end day of the pause

path = '/home/alex/OneDrive/CMIP6 OWP/Sea Ice Extent Expansion Pauses'
# path = '/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses'
atmpath = path+'/ERA5_and_SST_Avgs_forbboxes_forevents_V2.csv'
siipath = path+'/Pauses_combined_final_bytime_276-400.csv'
outpath = path+"/Figures/SFig_Atmosphere_Percentiles_Boxplots_doymax"+str(doymax)+" SLP.png"

atm = pd.read_csv(atmpath)
sii = pd.read_csv(siipath)
sii = sii.loc[(sii['Use'] == 1) & (sii['SYear'] > 1978) & (sii['DOY'] < doymax)]
sii = sii.reset_index()

bboxes = np.unique(atm['bbox'])
varlist = list(atm.columns[6:])
varlist2 = ['msl','t','seb','u10','v10']
varlist3 = ['Sea-Level Pressure','925-hPa Temperature','Surface Energy Balance','10-m Zonal Wind','10-m Meridional Wind']

region = ("Barents & Kara Seas", "East Greenland Sea", "Hudson Bay", "Bering Sea", "Sea of Okhotsk")
colors = ['#0076ba','#7980f7','#ed8032','#4aa12f','#2abf77']
prj = ccrs.NorthPolarStereo(central_longitude=-45,true_scale_latitude=70)


### RAW VALUES ###
enddatelist, ylist, mlist, dlist, blist = [[] for i in range(5)]
varlists = [[] for i in range(len(varlist))]

for i in range(len(sii)):
    
    year, month, day = sii.loc[i,['year','month','day']].astype(int)
    
    allrows = atm.loc[(atm['enddate'] == str(year)+'-'+md.dd[month-1]+'-'+md.dd[day-1])]

    row = allrows.loc[(allrows['year'] == year) & (allrows['shift'] == 0)]    

    for bbox in bboxes:
        
        enddatelist.append(str(year)+'-'+md.dd[month-1]+'-'+md.dd[day-1])
        ylist.append(year), mlist.append(month), dlist.append(day)
        blist.append(bbox)
        
        for vi, var in enumerate(varlist):
            val = row.loc[row['bbox'] == bbox,var].values[0]
            vals = allrows.loc[allrows['bbox'] == bbox,var].values
            vals = vals[np.isfinite(vals)]

            varlists[vi].append( np.mean(val > vals)*100 ) # percentile     
            
pdf = pd.DataFrame({'enddate':enddatelist, 'year':ylist, 'month':mlist,
                    'day':dlist, 'bbox': blist})
for vi in range(len(varlists)):
    pdf[varlist[vi]] = varlists[vi]

pdf['sst'] = np.where(pdf['year'] <= 1980, 50, pdf['sst'])
pdf = pdf.loc[np.in1d(pdf['bbox'],region)]

# Sort by preferred region order
pdf['Order'] = -1
for ri, reg in enumerate(region):
    pdf.loc[pdf['bbox'] == reg,'Order'] = ri
pdf = pdf.sort_values(by=['year','Order'])

# Identify all events
events = np.unique(pdf.enddate)
    
### FIGURE ###
fig = plt.figure(figsize=(7.5,5.5))

for v, var in enumerate(varlist2):
    ax = fig.add_subplot(2,3,v+1)
    ax.grid(axis='x', linestyle='dashed', linewidth=0.5)
    ax.set_axisbelow(True)

    x = np.arange(len(events))  # the label locations
    width = 1/(2+len(region))  # the width of the bars
    multiplier = 0

    a = [pdf.loc[pdf['Order'] == r,var].values for r in np.unique(pdf['Order'])]
    bplot = ax.boxplot(a, vert=False, patch_artist=True)
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
    for median in bplot['medians']:
        median.set_color('black')
        
    ax.set_xlim(xmin=-1,xmax=101)
    ax.set_xticks(np.arange(10,101,20))
    ax.set_xticklabels(np.arange(10,101,20), size=7)
    # ax.invert_yaxis()
    ax.set_title(md.abc[v]+". "+varlist3[v],weight='bold',size=8)
    
    if v%3 == 0:
        # ax.set_yticks(np.arange(0.33,len(events)+0.33))
        ax.set_yticklabels(region,size=7)
    else:
        ax.set_yticks([])
    
    ax.vlines(x=50, ymin=0.5, ymax=5.5, linewidth=1.2, color='0.3', linestyle='dashed')
    ax.set_ylim(ymin=0.5,ymax=5.5)
    
    if v >= 2:
        ax.set_xlabel('Percentile (wrt. 1979-2023)', size=8)
        
plt.tight_layout()

mapax = fig.add_axes([0.70,0.1,0.3,0.3], projection=prj)
mapax.add_feature(cfeature.COASTLINE)
mapax.add_feature(cfeature.LAND, fc='0.85')
mapax.set_extent([-3800000,3400000,-4500000,4500000], prj)

bboxname = ['Barents & Kara Seas','East Greenland Sea','Bering Sea','Sea of Okhotsk']
bbox = [[72,83,20,80],[70,80,-15,10],[55,65,-85,-75],[56,64,-180,-159],[50,65,130,165]]
for b in range(len(bbox)):
    mapax.plot((bbox[b][2],bbox[b][3],bbox[b][3],bbox[b][2],bbox[b][2]),
               (bbox[b][0],bbox[b][0],bbox[b][1],bbox[b][1],bbox[b][0]),
               color=colors[b],linewidth=2,
               transform=ccrs.Geodetic())
mapax.set_title(md.abc[v+1] + '. Bounding Box Definitions', weight='bold', size=8)

# ax.legend(fontsize=7,loc='lower left',bbox_to_anchor=(1.3,-0.05), framealpha=1)

plt.savefig(outpath, dpi=300)

### DETRENDED VALUES ###
