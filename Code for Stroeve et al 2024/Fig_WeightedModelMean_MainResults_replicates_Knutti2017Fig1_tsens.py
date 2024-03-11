#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 10 Apr 2023
Modified: 10 Apr 2023

Author: Alex Crawford

Purpose: Produce a figure similar to Figure 1 in Knutti et al. (2017) -- three
time series:
    1) Spaghetti plot time series of regional T (color coded by model weight)
    2) Spaghatti plot time series of ice-free period (color coded by model weight)
    3) Line ± Shading plot time series of ice-free period in historical - ssp585
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import CycloneModule_13_2 as md

'''*******************************************
Define Inputs
*******************************************'''

V = 'V9'
V2 = V+'.57' # V+'.4' # V+'.47' #
typ = 'Thickness' # ''
thresh = '10cm'
regvar = 'regHB3'
sigma_d, sigma_s = 0.49, 0.50 # 0.51, 0.4 # 0.39, 0.5 #

ivar = 'OPCavg'
regs = [41,42]

exp1, exp2, exp3 = 'historical', 'ssp585', 'ssp245'
ymin, ymax = [1920, 2015], [2013, 2099]
YYT = "1979-2021" # Years used to train the weighting

path = '/Volumes/Cassandra/CMIP6/RegionalStats/'
opath = ['/Volumes/Miranda/SeaIce/Bootstrap/AdvanceRetreat2/C10/Bootstrap_Regionalized_regHB3_Historical_C10.csv','/Volumes/Miranda/SeaIce/NASATeam/AdvanceRetreat2/C10/NASATeam_Regionalized_regHB3_Historical_C10.csv']
tpath = '/Volumes/Theseus/SurfaceTemperature/BEST/BEST_Annual_ts_1850-2022.csv'
figpath = '/Users/acrawfora/OneDrive - University of Manitoba/CMIP6 OWP/Hudson Bay Figures/ThicknessPhenologyFigures'

mods1 = np.array(['ACCESS-CM2_1',  'AWI-CM-1-1-MR_1', 'BCC-CSM2-MR_1', 'CESM2-WACCM_1',
        'CESM2_4', 'CMCC-CM2-SR5_1', 'CMCC-ESM2_1','CNRM-CM6-1_1', 'CNRM-ESM2-1_1',
        'EC-Earth3-CC_1', 'IPSL-CM6A-LR_1', 'KIOST-ESM_1', 'MIROC-ES2L_1',
        'MIROC6_1','MPI-ESM1-2-HR_1','MPI-ESM1-2-LR_1','MRI-ESM2-0_1',
        'NESM3_1', 'NorESM2-LM_1', 'NorESM2-MM_1'])

wclass = np.array([0,0.01,0.03,0.10])
# wclass = np.array([0,0.005,0.01,0.05])
wclasscolor = ['0.8','yellow','orange','red']
xlabel = ['Year','Global Annual Temperature Anomaly (°C)','Global Annual Temperature Anomaly (°C)']
ylabel = ['SSP585 Global Annual\nTemperature Anomaly (°C)','Ice-Free Period (days)','Ice-Free Period (days)']
xlimits = [ [min(ymin),max(ymax)+1] , [0,5.5], [0,5.5]]
titles = ['Western Hudson Bay',' Southern Hudson Bay']

'''*******************************************
Load Universal Data
*******************************************'''
# Load Non-Averaged Files
pdf = pd.read_csv(path+"/"+V+"/Annual_"+exp1+"_"+str(ymin[0])+'_'+str(ymax[0])+"_"+typ+"-"+thresh+"_"+regvar+".csv")
pdf = pdf.append(pd.read_csv(path+"/"+V+"/Annual_"+exp2+"_"+str(ymin[1])+'_'+str(ymax[1])+"_"+typ+"-"+thresh+"_"+regvar+".csv"), ignore_index=True)
pdf['Model'] = pdf['Family']+"_"+pdf['Member'].astype(str)

# Load Averaged Files
wavgdf = pd.read_csv(path+"/"+V+"/ModelWeighting/AnnualWeightedMean_"+exp1+"_"+str(ymin[0])+'_'+str(ymax[0])+"_"+typ+"-"+thresh+"_"+regvar+"_"+exp1+'-'+exp2+"_Train"+YYT+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+".csv")
wavgdf = wavgdf.append( pd.read_csv(path+"/"+V+"/ModelWeighting/AnnualWeightedMean_"+exp2+"_"+str(ymin[1])+'_'+str(ymax[1])+"_"+typ+"-"+thresh+"_"+regvar+"_"+exp1+'-'+exp2+"_Train"+YYT+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+".csv"), ignore_index=True)

wsedf = pd.read_csv(path+"/"+V+"/ModelWeighting/AnnualWeightedSE_"+exp1+"_"+str(ymin[0])+'_'+str(ymax[0])+"_"+typ+"-"+thresh+"_"+regvar+"_"+exp1+'-'+exp2+"_Train"+YYT+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+".csv")
wsedf = wsedf.append(pd.read_csv(path+"/"+V+"/ModelWeighting/AnnualWeightedSE_"+exp2+"_"+str(ymin[1])+'_'+str(ymax[1])+"_"+typ+"-"+thresh+"_"+regvar+"_"+exp1+'-'+exp2+"_Train"+YYT+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+".csv"), ignore_index=True)

# Load weights
wdf = pd.read_csv(path+"/"+V+"/ModelWeighting/Weights_"+typ+"-"+thresh+"_"+regvar+"_"+exp1+'-'+exp2+"_Train"+YYT+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+".csv")
wdf['Model'] = wdf['Family']+"_"+wdf['Member'].astype(str)
wdf = wdf.sort_values(by=['Weight'])

# Load Observations - Sea Ice (taking multi-method mean)
odf = pd.DataFrame([])
for o in opath:
    odf = pd.concat((odf,pd.read_csv(o)))
odf = odf.groupby(by=['Year','Region']).mean()
odf = odf.reset_index()

tdf = pd.read_csv(tpath)

tmean = tdf.loc[(tdf['Year'] >= 1850) & (tdf['Year'] <= 1900)]['tglobal'].mean()

taglobal = np.array(tdf.loc[ (tdf['Year'] >= min(odf['Year'])) & (tdf['Year'] <= max(odf['Year']))].groupby(by=['Year']).mean()['tglobal'])
odf['taglobal'] = np.repeat(taglobal,len(np.unique(odf['Region'])))

tdf = tdf.loc[ (tdf['Year'] >= ymin[0])].groupby(by=['Year']).mean()
tdf = tdf.reset_index()

# Prep Figure
fig = plt.figure(figsize=(7,5.5))

'''*******************************************
Plot #2: Sea Ice v. Temperature - Spaghetti
*******************************************'''

# Sort by temperature
wavgdf = wavgdf.sort_values(by=['taglobal'])
wsedf = wsedf.sort_values(by=['taglobal'])
odf = odf.sort_values(by=['taglobal'])


for ri in range(len(regs)):

    ax3 = fig.add_subplot(2,2,1+ri)

    # Plot CMIP6 Models, shading by weight %
    for mod in np.array(wdf['Model']):
        # Subset to model
        wpdf = pdf.loc[(pdf['Model'] == mod) & (pdf['Region'] == regs[ri])]
        wpdf = wpdf.sort_values(by=['taglobal'])

        # Determine plotting color based on weight
        wc = np.max(np.where(float(wdf[wdf['Model'] == mod]['Weight']) > wclass))
        ax3.plot( wpdf[wpdf['Region'] == regs[ri]]['taglobal'],wpdf[wpdf['Region'] == regs[ri]][ivar], color=wclasscolor[wc], linewidth=0.6 )

    # Observations
    ax3.plot(odf[odf['Region'] == regs[ri]]['taglobal'] - tmean, odf.loc[odf['Region'] == regs[ri]][ivar], color='blue', linewidth=1.0)

    # Axis Settings
    ax3.set_ylabel(ylabel[1],fontsize=8)
    ax3.set_xlabel(xlabel[1],fontsize=8)
    plt.xticks(fontsize=8), plt.yticks(fontsize=8)
    ax3.set_xlim(xmin = xlimits[1][0], xmax = xlimits[1][1])
    ax3.annotate(md.abc[ri], xy=(0.02,0.98), xycoords='axes fraction', weight = 'bold', ha = 'left', va = 'top')
    ax3.set_yticks(np.arange(0,360,60))
    ax3.set_ylim(ymin=0, ymax=365)

    ax3.set_title(titles[ri], weight='bold', size=9)

    # Legend
    if ri == 0:
        legax = plt.axes([0.33,0.59,0.1,0.08])
        legax.plot([0,0.2],[0.985,0.985], color='blue', linewidth=2)
        legax.plot([0,0.2],[0.660,0.660], color='red', linewidth=2)
        legax.plot([0,0.2],[0.330,0.330], color='orange', linewidth=2)
        legax.plot([0,0.2],[0.015,0.015], color='yellow', linewidth=2)

        legax.set_xlim(xmin=0, xmax=1)
        legax.set_ylim(ymin=-0.02, ymax=1.02)
        legax.set_axis_off()
        legax.annotate(' PMW Avg.\n>10% weight\n>3% weight\n>1%weight', xy=(0.3,0.5), xycoords='data', size=8, va='center')


'''*******************************************
Plot #3: Sea Ice v. Time - Envelope
*******************************************'''


### Make Equal-Weighted Ensembles ###
# Mean
pdfavg1 = pdf.loc[np.in1d(pdf['Model'],mods1)].groupby(by=['Region','Year']).mean()
pdfavg1 = pdfavg1.reset_index()
pdfavg1 = pdfavg1.sort_values(by=['taglobal'])

# Standard Deviation
pdfsd1 = pdf.loc[np.in1d(pdf['Model'],mods1)].groupby(by=['Region','Year']).std(ddof=1)
pdfsd1 = pdfsd1.reset_index()
pdfsd1 = pdfsd1.sort_values(by=['taglobal'])

for ri in range(len(regs)):
    ax5 = fig.add_subplot(2,2,3+ri)

    # Unweighted Mean (Single Member Ensemble)
    minline = pdfavg1.loc[pdfavg1['Region'] == regs[ri]][ivar] - 1.645*pdfsd1.loc[pdfsd1['Region'] == regs[ri]][ivar]
    maxline = pdfavg1.loc[pdfavg1['Region'] == regs[ri]][ivar] + 1.645*pdfsd1.loc[pdfsd1['Region'] == regs[ri]][ivar]
    minline = np.where(minline < 0, 0, np.where(minline > 365, 365, minline))
    maxline = np.where(maxline < 0, 0, np.where(maxline > 365, 365, maxline))

    ax5.plot(pdfavg1.loc[pdfavg1['Region'] == regs[ri]]['taglobal'], pdfavg1.loc[pdfavg1['Region'] == regs[ri]][ivar], color='black', linewidth=0.8, label='CMIP6 Avg.')
    ax5.fill_between(pdfavg1.loc[pdfavg1['Region'] == regs[ri]]['taglobal'],
                     minline, maxline,
                     color='black', alpha=0.3, linewidth=0)

    # Weighted Mean
    minline = wavgdf.loc[wavgdf['Region'] == regs[ri]][ivar] - 1.645*wsedf.loc[wsedf['Region'] == regs[ri]][ivar]
    maxline = wavgdf.loc[wavgdf['Region'] == regs[ri]][ivar] + 1.645*wsedf.loc[wsedf['Region'] == regs[ri]][ivar]
    minline = np.where(minline < 0, 0, np.where(minline > 365, 365, minline))
    maxline = np.where(maxline < 0, 0, np.where(maxline > 365, 365, maxline))

    ax5.plot(wavgdf.loc[wavgdf['Region'] == regs[ri]]['taglobal'], wavgdf.loc[wavgdf['Region'] == regs[ri]][ivar], color='red', linewidth=0.8,label='Weighted Avg.')
    ax5.fill_between(wavgdf.loc[wavgdf['Region'] == regs[ri]]['taglobal'],
                     minline, maxline,
                     color='red', alpha=0.3, linewidth=0)

    # Observations
    ax5.plot(odf[odf['Region'] == regs[ri]]['taglobal'] - tmean, odf.loc[odf['Region'] == regs[ri]][ivar], color='blue', linewidth=0.8, label = 'PMW Avg.')

    # Axis Settings
    ax5.set_ylabel(ylabel[2], fontsize=8)
    ax5.set_xlabel(xlabel[2], fontsize=8)
    plt.xticks(fontsize=8), plt.yticks(fontsize=8)
    ax5.set_xlim(xmin = xlimits[2][0], xmax = xlimits[2][1])
    ax5.annotate(md.abc[2+ri], xy=(0.02,0.98), xycoords='axes fraction', weight = 'bold', ha = 'left', va = 'top')
    ax5.set_yticks(np.arange(0,360,60))
    ax5.set_ylim(ymin=0, ymax=365)

    if ri == 0:
        ax5.legend(loc='lower right',fontsize=8, frameon=False)

plt.tight_layout(rect=[0,0,1,1])

# plt.savefig(figpath+"/Fig_Replication_of_Knutti2017Fig1_tsens_"+typ+"-"+thresh+"_"+regvar+"_"+exp1+"-"+exp2+"_Train"+YYT+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+"_2by2.png",dpi=300)