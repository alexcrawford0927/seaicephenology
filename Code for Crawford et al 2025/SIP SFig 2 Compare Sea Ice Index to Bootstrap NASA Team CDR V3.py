#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 11:48:34 2024

@author: acrawfora

Make Figure that shows the sea ice declines based on four data sets
"""

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import pandas as pd
import numpy as np
import CycloneModule_13_2 as md

inpath = '/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses/'
outpath = '/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses/Figures/SFig2_Compare_SeaIceIndex_To_RegionalDatasetsV3.png'

siedfall = pd.read_csv(inpath+"/SIE_Daily_ArcticWide.csv")
siedfall = siedfall.loc[:,[col for col in siedfall.columns if col not in ['SeaIceIndex']]]
siidf = pd.read_csv(inpath+"/SIE_Daily_ArcticWide_sii.csv")
siedfall = pd.merge(siedfall,siidf.loc[:,['JulianDay','SeaIceIndex','sie_sii','sii_delsie_5day']],how='inner',on='JulianDay')

pauses = pd.read_csv(inpath+"/Pauses_combined_final_bytime_276-400.csv")
pauses0 = pauses.loc[(pauses['Use'] == 1) & (pauses['DOY'] < 369)]

sievar = 'sie_cdr'

# labels = ['NOAA/NSIDC CDR','CDR NASA Team', 'CDR Bootstrap','Goddard NASA Team', 'Goddard Bootstrap', 'Sea Ice Index']
# varlist = ['sie_cdr','sie_nnt','sie_nbt','sie_nt','sie_bt','SeaIceIndex']
labels = ['NOAA/NSIDC CDR','NASA Team', 'Bootstrap']
varlist = ['sie_cdr','sie_nt','sie_bt']

#### FIGURE ###

fig = plt.figure(figsize=(7.5,8.5))

for p in range(len(pauses0)):
    
    ax = fig.add_subplot(5,3,p+1)
    
    row = pauses0.iloc[p]
    allrows = pauses.loc[(pauses['JulianDay'] >= row['JulianDay'] -12) & (pauses['JulianDay'] <= row['JulianDay'] + 8),:]
    
    siesub = siedfall.loc[(siedfall['JulianDay'] >= row['JulianDay'] - 12) & (siedfall['JulianDay'] <= row['JulianDay'] + 8),:]
    
    for v in range(len(labels)):
        ax.plot(siesub.loc[np.isfinite(siesub[varlist[v]])]['JulianDay'],siesub.loc[np.isfinite(siesub[varlist[v]])][varlist[v]], label=labels[v])
    
    
    ax.add_patch(Rectangle( xy=(allrows['JulianDay'].min()-5,np.floor(row[sievar])-0.5), width=allrows['JulianDay'].max()-allrows['JulianDay'].min()+5,
                     height= (np.floor(row[sievar])+2.2)-(np.floor(row[sievar])), facecolor='black', alpha=0.1, zorder=1, label='Pause in SIE Expansion\n(for CDR)' ))
    # ax.vlines(x=allrows['JulianDay'],ymin=np.floor(row['Extent'])-0.5,ymax=np.floor(row['Extent'])+3, color='blue', linewidth=1, linestyle='dashed')
    
    # ax.scatter(allrows['JulianDay'],allrows['Extent'], label='End Date of 5-day\nPause in SIE Growth', s=20, color='grey')
    # ax.scatter(allrows['JulianDay'],allrows['Extent'], s=30)
    
    # for v in ['nt','bt','cdr']:
    #     ax.scatter(siesub.loc[siesub[v+'_delsie_5day'] < 0,'JulianDay'],siesub.loc[siesub[v+'_delsie_5day'] < 0,'sie_'+v], s=20)    
    
    ax.set_ylim(ymin=np.floor(row[sievar])-0.5,ymax=np.floor(row[sievar])+2.1)
    ax.set_yticks(np.arange(np.floor(row[sievar])-0.5,np.floor(row[sievar])+2.5,0.5))
    
    xticks = np.arange(siesub['JulianDay'].values[0],siesub['JulianDay'].values[-1]+1,2)
    ax.set_xticks(xticks)
    xticklabels = []
    for i in xticks:
        try:
            xticklabels.append(md.mmm[siesub.loc[siesub['JulianDay'] == i]['month'].values[0]-1] + ' ' + md.dd[siesub.loc[siesub['JulianDay'] == i]['day'].values[0]-1])
        except:
            xticklabels.append(md.mmm[siesub.loc[siesub['JulianDay'] == i-2]['month'].values[0]-1] + ' ' + md.dd[siesub.loc[siesub['JulianDay'] == i-2]['day'].values[0]+1])
    
    ax.set_xticklabels(xticklabels , rotation=90 )
    ax.tick_params(axis='both', which='major', labelsize=7)
    if p%3 == 0:
        ax.set_ylabel('Arctic Sea Ice Extent\n(million km$^2$)', size=8)
    
    newdate = md.timeAdd([int(row['year']),int(row['month']),int(row['day'])],[0,0,-6])

    ax.annotate(md.abc[p]+". "+ md.dd[newdate[2]-1] + " " +md.mmm[newdate[1] - 1] + " " + str(newdate[0]) , 
                xy=(0.02,0.98), xycoords='axes fraction',
                ha='left',va='top',size=8,weight='bold')

plt.tight_layout(rect=[0,0,1,1])

ax.legend(fontsize=7,loc='upper left',bbox_to_anchor=(1.3,0.75), edgecolor='w', framealpha=1,ncols=1)


### RANKS ###

# Where do each of these events rank in terms of slowest SIE growth for each method?
btall = siedfall.loc[(siedfall['bt_delsie_5day'] != 0) & ((siedfall['month'] > 9) | (siedfall['month'] < 3)) & 
             (siedfall['sie_bt'] >= 7.75) & (siedfall['sie_bt'] <= 13.25),'bt_delsie_5day'] 
ntall = siedfall.loc[(siedfall['nt_delsie_5day'] != 0) & ((siedfall['month'] > 9) | (siedfall['month'] < 3)) & 
             (siedfall['sie_nt'] >= 7.75) & (siedfall['sie_nt'] <= 13.25),'nt_delsie_5day'] 
cdrall = siedfall.loc[(siedfall['cdr_delsie_5day'] != 0) & ((siedfall['month'] > 9) | (siedfall['month'] < 3)) & 
             (siedfall['sie_cdr'] >= 7.75) & (siedfall['sie_cdr'] <= 13.25),'cdr_delsie_5day'] 
siiall = siedfall.loc[(siedfall['sii_delsie_5day'] != 0) & ((siedfall['month'] > 9) | (siedfall['month'] < 3)) & 
             (siedfall['SeaIceIndex'] >= 7.75) & (siedfall['SeaIceIndex'] <= 13.25),'sii_delsie_5day'] 
nbtall = siedfall.loc[(siedfall['nbt_delsie_5day'] != 0) & ((siedfall['month'] > 9) | (siedfall['month'] < 3)) & 
             (siedfall['sie_nbt'] >= 7.75) & (siedfall['sie_nbt'] <= 13.25),'nbt_delsie_5day'] 
nntall = siedfall.loc[(siedfall['nnt_delsie_5day'] != 0) & ((siedfall['month'] > 9) | (siedfall['month'] < 3)) & 
             (siedfall['sie_nnt'] >= 7.75) & (siedfall['sie_nnt'] <= 13.25),'nnt_delsie_5day'] 

ntrank, btrank, cdrrank, siirank, nntrank, nbtrank = [[] for i in range(6)]
for p in range(len(pauses0)):
    jd = pauses0.iloc[p]['JulianDay']
    
    btrank.append( np.sum(btall <= siedfall.loc[siedfall['JulianDay'] == jd,'bt_delsie_5day'].values[0]) )
    ntrank.append( np.sum(ntall <= siedfall.loc[siedfall['JulianDay'] == jd,'nt_delsie_5day'].values[0]) )
    cdrrank.append( np.sum(cdrall <= siedfall.loc[siedfall['JulianDay'] == jd,'cdr_delsie_5day'].values[0]) )
    siirank.append( np.sum(siiall <= siedfall.loc[siedfall['JulianDay'] == jd,'sii_delsie_5day'].values[0]) )
    nbtrank.append( np.sum(nbtall <= siedfall.loc[siedfall['JulianDay'] == jd,'nbt_delsie_5day'].values[0]) )
    nntrank.append( np.sum(nntall <= siedfall.loc[siedfall['JulianDay'] == jd,'nnt_delsie_5day'].values[0]) )
    
ranksize = [len(btall),len(ntall),len(cdrall), len(siiall), len(nbtall), len(nntall)]

# table = plt.table(cellText= [btrank, ntrank, cdrrank, siirank],
#                       bbox=[1.4,-0.3,0.85,0.5], fontsize=7.5, edges='open')

# ax.annotate('Bootstrap',xy=(1.38,0.13), xycoords='axes fraction', ha='right', fontsize=7)
# ax.annotate('NASA Team',xy=(1.38,0.0), xycoords='axes fraction', ha='right', fontsize=7)
# ax.annotate('NOAA/NSIDC CDR',xy=(1.38,-0.12), xycoords='axes fraction', ha='right', fontsize=7)
# ax.annotate('Sea Ice Index',xy=(1.38,-0.25), xycoords='axes fraction', ha='right', fontsize=7)

# txty = 0.2
# txtx = np.arange(1.36,4.0,0.105)
# for r in range(len(pauses0)):
#     ax.annotate( str(int(pauses0.iloc[r]['year'])) + " " +md.mmm[int(pauses0.iloc[r]['month'] - 1)] + " " + md.dd[int(pauses0.iloc[r]['day'] -1)],
#                 xy=(txtx[r],txty), xycoords='axes fraction',
#                 fontsize=7, rotation=90, ha='left')

# ax.legend(fontsize=7,loc='upper left',bbox_to_anchor=(-2,0), edgecolor='w', framealpha=1,ncols=3)
# ax.annotate(md.abc[p+1]+'. Rank (1 =\nMinimum\n5-day\nSIE Change)', xy=(1.22,0.3), xycoords='axes fraction',
#             ha='center',fontsize=7, weight='bold')


plt.savefig(outpath,dpi=300)



