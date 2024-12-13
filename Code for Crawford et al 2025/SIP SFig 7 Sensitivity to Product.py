#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date Created: 31 Jul 2024
Date Modified: 31 Jul 2024
Author: Alex Crawford
Purpose: Shows the senstivity of results to different choices of SIC algorithm or
their combination.
"""

'''*******************************************
Load Modules
*******************************************'''
print(' ***Load Modules***')

import numpy as np
import pandas as pd
import CycloneModule_13_2 as md
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

'''*******************************************
Declare Variables
*******************************************'''
print(' ***Define Variables***')

inpath = '/home/alex/OneDrive/CMIP6 OWP/Sea Ice Extent Expansion Pauses'
# inpath = '/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses/'
outpath = inpath+"/Figures/SFig_Sensitivity_to_SICDataset_and_SIEThreshold_V4_Combined_bytime.png"

siedfall = pd.read_csv(inpath+'/Pauses_combined_anyneg_bytime_276-400.csv')

dmin, dmax = 276, 369
ymin, ymax = 1979, 2022
plotlabs = ['CDR','NASA Team','Bootstrap']

'''*******************************************
Process Variables
*******************************************'''
print(' ***Process Variables***')
siedfall['time'] = siedfall['SYear'] + (siedfall['DOY']-dmin)/(365-dmin)


prelim = [siedfall.loc[(siedfall['cdr_delsie_5day'] < 0)  & (siedfall['LowVal'] < 3) & (siedfall['DOY'] > dmin) & ( (siedfall['DOY'] < dmax) )] ,
siedfall.loc[(siedfall['nt_delsie_5day'] < 0) & (siedfall['LowVal'] < 3)  & (siedfall['DOY'] > dmin) & ( (siedfall['DOY'] < dmax) )] ,
siedfall.loc[(siedfall['bt_delsie_5day'] < 0) & (siedfall['LowVal'] < 3)  & (siedfall['DOY'] > dmin) & ( (siedfall['DOY'] < dmax) )] ]

final = [siedfall.loc[(siedfall['cdr_delsie_5day'] < 0) & (siedfall['LowVal'] == 3) & (siedfall['DOY'] > dmin) & ( (siedfall['DOY'] < dmax) )] ,
siedfall.loc[(siedfall['nt_delsie_5day'] < 0) & (siedfall['LowVal'] == 3) & (siedfall['DOY'] > dmin) & ( (siedfall['DOY'] < dmax) )] ,
siedfall.loc[(siedfall['bt_delsie_5day'] < 0) & (siedfall['LowVal'] == 3) & (siedfall['DOY'] > dmin) & ( (siedfall['DOY'] < dmax) )] ]


# cdrall = siedfall.loc[(siedfall['cdr_delsie_5day'] < 0) & (siedfall['DOY'] > dmin) & ( (siedfall['DOY'] < dmax) )] 
# ntall = siedfall.loc[(siedfall['nt_delsie_5day'] < 0) & (siedfall['DOY'] > dmin) & ( (siedfall['DOY'] < dmax) )] 
# btall = siedfall.loc[(siedfall['bt_delsie_5day'] < 0) & (siedfall['DOY'] > dmin) & ( (siedfall['DOY'] < dmax) )] 

# ydf = pd.DataFrame({'SYear':np.arange(ymin,ymax+1)})

# cdrneg = cdrall.loc[:,['SYear','NegVal']].groupby(by=['SYear']).count().reset_index()
# cdrneg.columns = ['SYear','cdrneg']
# cdrfinal = cdrall.loc[cdrall['LowVal'] == 3,['SYear','NegVal']].groupby(by=['SYear']).count().reset_index()
# cdrfinal.columns = ['SYear','cdrfinal']

# btneg = btall.loc[:,['SYear','NegVal']].groupby(by=['SYear']).count().reset_index()
# btneg.columns = ['SYear','btneg']
# btfinal = btall.loc[btall['LowVal'] == 3,['SYear','NegVal']].groupby(by=['SYear']).count().reset_index()
# btfinal.columns = ['SYear','btfinal']

# ntneg = ntall.loc[:,['SYear','NegVal']].groupby(by=['SYear']).count().reset_index()
# ntneg.columns = ['SYear','ntneg']
# ntfinal = ntall.loc[ntall['LowVal'] == 3,['SYear','NegVal']].groupby(by=['SYear']).count().reset_index()
# ntfinal.columns = ['SYear','ntfinal']

# ydf = pd.merge(ydf,cdrneg,how='outer')
# ydf = pd.merge(ydf,cdrfinal,how='outer')
# ydf = pd.merge(ydf,btneg,how='outer')
# ydf = pd.merge(ydf,btfinal,how='outer')
# ydf = pd.merge(ydf,ntneg,how='outer')
# ydf = pd.merge(ydf,ntfinal,how='outer')

# for col in ydf.columns:
#     ydf[col] = np.where(np.isnan(ydf[col]), 0, ydf[col])

xvals = np.arange(1980,2025,5)
xvallabs = xvals
plab = [str(round(len(final[s]) / (len(final[s])+len(prelim[s])) * 100,1))+'%' for s in range(len(prelim))]

'''*******************************************
Plotting
*******************************************'''
print('***Plotting***')

fig = plt.figure(figsize=(7.5,2.75))

ax = fig.add_subplot(1,1,1)
ax.grid(axis='x', linestyle='dashed', linewidth=0.5)
# ax.vlines(x=xvals[-4]-0.5,ymin=yvals[0]-0.125,ymax=yvals[-1]+0.125, color='k', linewidth=0.6, linestyle='dashed',zorder=1)
# ax.hlines(y=[13.125,13.375],xmin=siedfall.iloc[0]['SYear'],xmax=siedfall.iloc[-1]['SYear'],color='k', linewidth=0.6, linestyle='dashed',zorder=1)

ax.set_xticks(xvals)
ax.set_xticklabels(xvallabs, size=7, rotation=90, ha='center', va='top')
ax.set_yticks([0,1,2])
ax.set_yticklabels([plotlabs[l] +'\n('+plab[l]+')' for l in range(len(plotlabs))], size=7)
ax.set_ylim(ymin=-0.5,ymax=2.5)
ax.set_xlim(xmin=1979,xmax=ymax+2)
# ax.set_title(md.abc[v]+". " + plotlabs[v], weight='bold', size=8)
ax.invert_yaxis()
   
for s, ss  in enumerate(['cdr','bt','nt']):
    
    # ax.scatter(ydf['SYear'], np.repeat(s,len(ydf)), marker='s', s=ydf[ss+'neg']*20, color='k', zorder=3)
    # ax.scatter(ydf['SYear'], np.repeat(s,len(ydf)), marker='s', s=ydf[ss+'final']*20, color='blue', zorder=4)

    ax.scatter(final[s]['time'], np.repeat(s,len(final[s])), marker='s', s=20, color='royalblue', zorder=3)
    ax.scatter(prelim[s]['time'], np.repeat(s,len(prelim[s])), marker='x', s=20, color='k', zorder=3)
        

scat1 = ax.scatter(prelim[s]['time'], np.repeat(s-10,len(prelim[s])), marker='x', s=20, color='k', label='SIE loss w/ this product, but\n∆SIE ≥ 0.1 million km$^2$ in at least one product', zorder=3)
scat2 = ax.scatter(final[s]['time'], np.repeat(s-10,len(final[s])), marker='s', s=20, color='royalblue', label='SIE loss w/ this product &\n∆SIE < 0.1 million km$^2$ w/ all products',zorder=4)


plt.tight_layout(rect=[0,0.1,1,1])

## Legend 1
handles, labels = plt.gca().get_legend_handles_labels()
legend1 = ax.legend([handles[i] for i in [0,1]], [labels[i] for i in [0,1]],
                    fontsize=7,loc='upper left',bbox_to_anchor=(0,-0.2), edgecolor='w', framealpha=1, ncols=2)

plt.savefig(outpath,dpi=300)
        
