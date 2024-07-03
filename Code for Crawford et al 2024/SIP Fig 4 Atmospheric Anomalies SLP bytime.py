#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 22 Apr 2024
Modified: 03 Jul 2024

Author: Alex Crawford

Purpose: Compare atmospheric parameters for SIE pauses to distribution of
6-day periods for same time of year -- this verison includes SLP instead of SST
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import CycloneModule_13_2 as md

# path = '/home/alex/OneDrive/CMIP6 OWP/Sea Ice Extent Expansion Pauses'
path = '/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses'
atmpath = path+'/ERA5_and_SST_Avgs_forbboxes_forevents_V2.csv'
siipath = path+'/Pauses_combined_final_bytime_276-400.csv'
outpath = path+"/Figures/Fig_Atmosphere_Percentiles V3 SLP.png"

doymax = 369 # Maximum day of year for the end day of the pause

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
startdatelist, ylist, mlist, dlist, blist = [[] for i in range(5)]
varlists = [[] for i in range(len(varlist))]

for i in range(len(sii)):
    
    year, month, day = sii.loc[i,['year','month','day']].astype(int)
    startdate = md.timeAdd([year,month,day],[0,0,-6])
    
    allrows = atm.loc[(atm['enddate'] == str(year)+'-'+md.dd[month-1]+'-'+md.dd[day-1])]

    row = allrows.loc[(allrows['year'] == year) & (allrows['shift'] == 0)]    

    for bbox in bboxes:
        
        startdatelist.append(str(startdate[0])+'-'+md.dd[startdate[1]-1]+'-'+md.dd[startdate[2]-1])
        ylist.append(startdate[0]), mlist.append(startdate[1]), dlist.append(startdate[2])
        blist.append(bbox)
        
        for vi, var in enumerate(varlist):
            val = row.loc[row['bbox'] == bbox,var].values[0]
            vals = allrows.loc[allrows['bbox'] == bbox,var].values
            vals = vals[np.isfinite(vals)]

            varlists[vi].append( np.mean(val > vals)*100 ) # percentile     
            
pdf = pd.DataFrame({'startdate':startdatelist, 'year':ylist, 'month':mlist,
                    'day':dlist, 'bbox': blist})
for vi in range(len(varlists)):
    pdf[varlist[vi]] = varlists[vi]

pdf['sst'] = np.where(pdf['year'] <= 1980, 50, pdf['sst'])

events = np.unique(pdf.startdate)
    
### FIGURE ###
fig = plt.figure(figsize=(7.5,8.5))

for v, var in enumerate(varlist2):
    ax = fig.add_subplot(2,3,v+1)
    ax.grid(axis='x', linestyle='dashed', linewidth=0.5)
    ax.set_axisbelow(True)
    
    vals_means = {
        region[0]: list(pdf.loc[pdf['bbox'] == region[0],var].values), 
        region[1]: list(pdf.loc[pdf['bbox'] == region[1],var].values), 
        region[2]: list(pdf.loc[pdf['bbox'] == region[2],var].values), 
        region[3]: list(pdf.loc[pdf['bbox'] == region[3],var].values),
        region[4]: list(pdf.loc[pdf['bbox'] == region[4],var].values)
        }
    
    x = np.arange(len(events))  # the label locations
    width = 1/(2+len(region))  # the width of the bars
    multiplier = 0

    for attribute, measurement in vals_means.items():
        offset = width * multiplier
        ax.barh(x + offset, np.array(measurement)-50, width, color=colors[multiplier], label=attribute, left=0)
        multiplier += 1
        
    ax.set_xlim(xmin=-50,xmax=50)
    ax.set_xticks(np.arange(-40,51,20))
    ax.set_xticklabels(np.arange(10,101,20), size=7)
    ax.invert_yaxis()
    ax.set_title(md.abc[v]+". "+varlist3[v],weight='bold',size=8)
    
    if v%3 == 0:
        ax.set_yticks(np.arange(0.33,len(events)+0.33))
        ax.set_yticklabels(events,size=7)
    else:
        ax.set_yticks([])
    
    if v >= 2:
        ax.set_xlabel('Percentile (wrt. 1979-2023)', size=8)
        
plt.tight_layout()

mapax = fig.add_axes([0.70,0.15,0.3,0.3], projection=prj)
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
mapax.set_title(md.abc[v] + '. Bounding Box Definitions', weight='bold', size=8)

ax.legend(fontsize=7,loc='lower left',bbox_to_anchor=(1.3,-0.05), framealpha=1)

plt.savefig(outpath, dpi=300)

### DETRENDED VALUES ###
