#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 9 Aug 2023
Modified: 2 Feb 2024
Author: Alex Crawford
Purpose: Makes plot that show the retreat and advance 

"""

'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")
import netCDF4 as nc
import numpy as np
import xesmf as xe
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle

'''*******************************************
Declare Variables
*******************************************'''
print("Declaring Variables")

modsource = "CMIP6"
exp1, exp2 = 'historical', 'ssp585'
V = 'V9'
ct = '10cm'
regvar = 'regHB3'
refyears = '1979-2021'
weightver = 'V9.57_d0.49_s0.5'

doys = [1,32,60,91,121,152,182,213,244,274,305,335,366,397,425]
doyl = ["Jan 1","Feb 1","Mar 1","Apr 1","May 1","Jun 1","Jul 1","Aug 1","Sep 1","Oct 1","Nov 1","Dec 1","Jan 1","Feb 1","Mar 1"]

varlist = ['lrd','fad']
reglist = [41,42]
regname = ['Southern Hudson Bay','Western Hudson Bay']

verlist = ['Bias Adj. Average','Weighted Mean']
verlist2 = ['Bias Adj. Avg.','Weighted Avg.']
tanoms = np.arange(1.0,4.5,0.5)
color = ['blue', 'red']
yadj = [-0.1,0.1]
marker=['o','s']
ylabel = 'Global Mean Temperature Anomaly (°C)'
tint = 0.25
ymin, ymax = 0.5, 4.5
xmin, xmax = 91, 425

t_obs = 1.16 # 2012-2021 average temperature (from observations!)
t_obs_years = '2012-2021 Avg.'

# Path
snowpath = "/Volumes/Cassandra/CMIP6/RegionalStats/RawVariables/biascorrected1981-2021"

figpath = "/Users/acrawfora/OneDrive - University of Manitoba/CMIP6 OWP/Hudson Bay Figures/ThicknessPhenologyFigures"

#### Load Data ####
sisnthick = pd.read_csv(snowpath+"/sisnthick_AnnualBiasCorrected_by_GlobalTemperatureAnomaly_1981-2021_"+exp1+"-"+exp2+"_"+regvar+".csv")
siphen = pd.read_csv(figpath+"/Fig_compare_many_multi-model_means_Thickness-"+ct+"_"+regvar+"_"+exp1+"-"+exp2+"_Train"+refyears+"_"+weightver+'.csv')

#### Set Up Graphs ####
fig = plt.figure(figsize=(5,6))

ax1 = fig.add_subplot(2,1,1)
ax1.set_ylim(ymin=ymin,ymax=ymax)
ax1.set_xlim(xmin=xmin,xmax=xmax)
ax1.set_xticks( doys[4:])
ax1.set_xticklabels( doyl[4:], rotation=90, fontsize=8 )
ax1.grid(linestyle='dashed',axis='x')
ax1.set_ylabel(ylabel, size=8)
txt = ax1.annotate('a. ' + regname[0], xy=(95,ymax-0.04), xycoords='data', size=8, ha='left', va='top', weight='bold')
txt.set_path_effects([pe.withStroke(linewidth=2, foreground='w')])

ax2 = fig.add_subplot(2,1,2)
ax2.set_ylim(ymin=ymin,ymax=ymax)
ax2.set_xlim(xmin=xmin,xmax=xmax)
ax2.set_xticks( doys[4:])
ax2.set_xticklabels( doyl[4:], rotation=90, fontsize=8 )
ax2.grid(linestyle='dashed',axis='x')
ax2.set_ylabel(ylabel, size=8)
txt = ax2.annotate('b. ' + regname[1], xy=(95,ymax-0.04), xycoords='data', size=8, ha='left', va='top', weight='bold')
txt.set_path_effects([pe.withStroke(linewidth=2, foreground='w')])


#### Plot Ice-Free Period ####
for vi in range(len(verlist)):
    
    for t in tanoms:
        
        for ri in range(len(reglist)):

            lrdsub = siphen.loc[(siphen['Region'] == reglist[ri]) & (siphen['Version'] == verlist[vi]) & (siphen['taglobal'] == t)]['LRDavg']
            fadsub = siphen.loc[(siphen['Region'] == reglist[ri]) & (siphen['Version'] == verlist[vi]) & (siphen['taglobal'] == t)]['FADavg'] 
            
            if ri == 0:
                ax1.hlines(y=t, xmin=float(lrdsub), xmax=float(fadsub), colors=color[vi], alpha=0.5, linewidths=3)
                ax1.scatter((lrdsub,fadsub), (t,t),  c=color[vi], s=30, marker=marker[vi], zorder=4)
            else:
                ax2.hlines(y=t, xmin=float(lrdsub), xmax=float(fadsub), colors=color[vi], alpha=0.5, linewidths=3)
                if t == tanoms[0]:
                    ax2.scatter((lrdsub,fadsub), (t,t), c=color[vi], s=30, marker=marker[vi], zorder=4, label=verlist2[vi])
                else:
                    ax2.scatter((lrdsub,fadsub), (t,t), c=color[vi], s=30, marker=marker[vi], zorder=4)
    
# Add horizontal line for an observational temperature
ax1.plot((xmin, xmax),(t_obs,t_obs),linewidth=1, linestyle='dashed', color='0.6')
ax2.plot((xmin, xmax),(t_obs,t_obs),linewidth=1, linestyle='dashed', color='0.6')
txts = [ ax1.annotate(t_obs_years,xy=(xmax-32,t_obs-0.01), xytext=(xmax+10,t_obs-0.55),
                      xycoords='data',textcoords='data',size=7,ha='right',
                      arrowprops=dict(arrowstyle='->',facecolor='k')),
        ax2.annotate(t_obs_years,xy=(xmax-32,t_obs-0.01), xytext=(xmax+10,t_obs-0.55),
                              xycoords='data',textcoords='data',size=7,ha='right',
                              arrowprops=dict(arrowstyle='->',facecolor='k'))]
[txt.set_path_effects([pe.withStroke(linewidth=3, foreground='w')]) for txt in txts]
    

#### Plot April Snow Depth ####
sizelist = []
for ri in range(len(reglist)):
    for t in tanoms:
        meanbymod = sisnthick.loc[ (sisnthick['Region'] == reglist[ri]) & (sisnthick['taglobal'] >= t-tint) & (sisnthick['taglobal'] <= t+tint)].groupby(('Model')).mean()
        
        size = meanbymod.loc[(meanbymod['sisnthick'] > 0)]['sisnthick'].mean()

        if ri == 0:
            ax1.scatter(106, t, c='gray', s=2000*(size-0.05), marker='o')
        else:
            ax2.scatter(106, t, c='gray', s=2000*(size-0.05), marker='o')
        
        sizelist.append(size)


plt.tight_layout(rect=[0,0.05,1,1])

### Legends ###
leg = ax2.legend(bbox_to_anchor=(0.98,0.05), bbox_transform=fig.transFigure, ncol=2, fontsize=7)
leg.get_frame().set_edgecolor('w')
plt.annotate('Ice-Free Period', xy=(0.75,0.06), xycoords='figure fraction', size=8, va='top', ha='center')

# Legend for Snow
leg2 = plt.axes([0.05,0.01,0.45,0.045])

labels = [6,9,12,15,18]
sizes = [20*(label-5) for label in labels]

for i in range(len(sizes)):
    leg2.scatter(i*1.1, 1, c='gray', marker='o', s = sizes[i])
    leg2.annotate(labels[i], xy=(i+(sizes[i])/350+0.07, 0.99), size=7)

leg2.set_xlim(xmin=-0.25,xmax=len(sizes)+0.25)

# Remove axes from legend
leg2.spines['bottom'].set_color('w')
leg2.spines['top'].set_color('w')
leg2.spines['right'].set_color('w')
leg2.spines['left'].set_color('w')
leg2.set_xticks([]), leg2.set_yticks([])

# Add label
plt.annotate('April Snow\nDepth (cm)', xy=(0.065,0.085), xycoords='figure fraction', color='k', size=7, va='top', ha='center')
plt.annotate('', xy=(0.13,0.17), xytext=(0.06,0.09), xycoords='figure fraction',va='top',ha='center',
              arrowprops={'width':1,'headwidth':6,'facecolor':'k','shrink':0.05})

plt.savefig(figpath+"/AdvanceRetreat_LinePlot_by_GlobalTemperatureAnomaly_"+refyears+"_"+exp1+"-"+exp2+"_Thickness-"+ct+"_"+regvar+"_V2.png", dpi=300)


