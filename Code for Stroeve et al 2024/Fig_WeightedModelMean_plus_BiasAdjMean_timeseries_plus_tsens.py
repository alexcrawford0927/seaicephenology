#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 10 Apr 2023
Modified: 2 Feb 2024

Author: Alex Crawford

Purpose: Produce a figure similar to Figure 1 in Knutti et al. (2017) -- two
time series:
    1) Line ± Shading plot time series of ice-free period in historical - ssp585
    2) Line ± Shading plot v. global mean annual surface temperature anomaly
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import CycloneModule_13_2 as md
import matplotlib.patheffects as pe
from scipy import interpolate

'''*******************************************
Define Inputs
*******************************************'''

V = 'V9'
V2 = V+'.57' # V+'.1' # V+'.4'
typ = 'Thickness' # ''
thresh = '10cm'
regvar = 'regHB3'
sigma_d, sigma_s = 0.49, 0.5 # 0.23, 0.50 #  0.51, 0.4

ivar = 'OPCavg'
regs = [41,42]
REGS = ['Southern Hudson Bay', 'Western Hudson Bay']

exp1, exp2 = 'historical', 'ssp585'
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
mods2 = np.array(['ACCESS-CM2_1', 'BCC-CSM2-MR_1', 'CESM2-WACCM_1', 'CESM2_4',
         'CNRM-ESM2-1_1','IPSL-CM6A-LR_1','KIOST-ESM_1','MIROC6_1','MPI-ESM1-2-LR_1',
         'MRI-ESM2-0_1','NorESM2-LM_1'])
# modspecial = ['MRI-ESM2-0']

tamin = 0.25 # minimum temperature anomaly to plot

wclass = np.array([0,0.01,0.03,0.10])
# wclass = np.array([0,0.005,0.01,0.05])
# wclasscolor = ['0.8','yellow','orange','red']
obcol = 'gold' #'blue'
linecolors = ['k','k','black','b']
linestyles = ['solid','dashed','solid','dashed']
linelabels = ['Average - All Models','Average - Subset','Bias Adj. Average','Bias Adj. Avg. - Subset']
ylabel = ['Ice-Free Period (days)','Ice-Free Period (days)']
xlimits = [ [min(ymin)-1,max(ymax)+1] , [0,5.35]]
xlabel = ['Year','Global Annual Temperature Anomaly (°C wrt 1850-1900)']

'''*******************************************
Load Universal Data
*******************************************'''
# Load Non-Averaged Files
pdf = pd.read_csv(path+"/"+V+"/Annual_"+exp1+"_"+str(ymin[0])+'-'+str(ymax[0])+"_"+typ+"-"+thresh+"_"+regvar+".csv")
pdf = pdf.append(pd.read_csv(path+"/"+V+"/Annual_"+exp2+"_"+str(ymin[1])+'-'+str(ymax[1])+"_"+typ+"-"+thresh+"_"+regvar+".csv"), ignore_index=True)
pdf['Model'] = pdf['Family']+"_"+pdf['Member'].astype(str)
pdf = pdf.sort_values(by=['Family','Member','Region','Year'])

# Load Averaged Files
wavgdf = pd.read_csv(path+"/"+V+"/ModelWeighting/AnnualWeightedMean_"+exp1+"_"+str(ymin[0])+'_'+str(ymax[0])+"_"+typ+"-"+thresh+"_"+regvar+"_"+exp1+'-'+exp2+"_Train"+YYT+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+".csv")
wavgdf = wavgdf.append( pd.read_csv(path+"/"+V+"/ModelWeighting/AnnualWeightedMean_"+exp2+"_"+str(ymin[1])+'_'+str(ymax[1])+"_"+typ+"-"+thresh+"_"+regvar+"_"+exp1+'-'+exp2+"_Train"+YYT+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+".csv"), ignore_index=True)

wsedf = pd.read_csv(path+"/"+V+"/ModelWeighting/AnnualWeightedSE_"+exp1+"_"+str(ymin[0])+'_'+str(ymax[0])+"_"+typ+"-"+thresh+"_"+regvar+"_"+exp1+'-'+exp2+"_Train"+YYT+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+".csv")
wsedf = wsedf.append(pd.read_csv(path+"/"+V+"/ModelWeighting/AnnualWeightedSE_"+exp2+"_"+str(ymin[1])+'_'+str(ymax[1])+"_"+typ+"-"+thresh+"_"+regvar+"_"+exp1+'-'+exp2+"_Train"+YYT+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+".csv"), ignore_index=True)

# Load Bias-Corrected Files
pdfb = pd.read_csv(path+"/"+V+"/AnnualBiasCorrected_by_Time_"+YYT+"_"+exp1+"-"+exp2+'_'+typ+"-"+thresh+"_"+regvar+".csv")
pdfg = pd.read_csv(path+"/"+V+"/AnnualBiasCorrected_by_GlobalTemperatureAnomaly_"+YYT+"_"+exp1+"-"+exp2+'_'+typ+"-"+thresh+"_"+regvar+".csv")
pdft = pd.read_csv(path+"/"+V+"/AnnualBiasCorrected_by_RegionalTemperature_"+YYT+"_"+exp1+"-"+exp2+'_'+typ+"-"+thresh+"_"+regvar+".csv")

pdfb = pdfb.sort_values(by=['Family','Member','Region','Year'])
pdfg = pdfg.sort_values(by=['Family','Member','Region','Year'])
pdft = pdft.sort_values(by=['Family','Member','Region','Year'])

pdfb['taglobal'] = np.array(pdf['taglobal'])
pdfg['taglobal'] = np.array(pdf['taglobal'])
pdft['taglobal'] = np.array(pdf['taglobal'])

# Load Observations - Sea Ice (taking multi-method mean)
odf = pd.DataFrame([])
for o in opath:
    odf = pd.concat((odf,pd.read_csv(o)))
odf = odf.groupby(by=['Year','Region']).mean()
odf = odf.reset_index()

# Load Observatons - Temperature
tdf = pd.read_csv(tpath)

tmean = tdf.loc[(tdf['Year'] >= 1850) & (tdf['Year'] <= 1900)]['tglobal'].mean() # baseline period temperature

taglobal = np.array(tdf.loc[ (tdf['Year'] >= min(odf['Year'])) & (tdf['Year'] <= max(odf['Year']))].groupby(by=['Year']).mean()['tglobal']) - tmean
odf['taglobal'] = np.repeat(taglobal,len(np.unique(odf['Region'])))

'''*******************************************
Plot #1: Sea Ice v. Time - Envelope
*******************************************'''

# # Prep Figure
# fig = plt.figure(figsize=(179/25.4,8))

# ### Make Equal-Weighted & Bias-Corrected Ensembles ###
# pdfs = [pdf.loc[np.in1d(pdf['Model'],mods1)], pdf.loc[np.in1d(pdf['Model'],mods2)], pdfb.loc[np.in1d(pdfb['Model'],mods1)], pdfb.loc[np.in1d(pdfb['Model'],mods2)]]
# dfavg, dfse = [], []

# for df in pdfs:
#     # Mean
#     dfavg1 = df.groupby(by=['Region','Year']).mean()
#     dfavg1 = dfavg1.reset_index()
#     dfavg.append(dfavg1)

#     # Standard Deviation
#     dfsd1 = df.groupby(by=['Region','Year']).std(ddof=1)
#     dfsd1 = dfsd1.reset_index()
#     dfse.append(dfsd1)


# for ri in range(len(regs)):
#     ax1 = fig.add_subplot(2,1,1+ri)

#     for pi in [2]:
#         # Unweighted Mean (Single Member Ensemble)
#         minline = dfavg[pi].loc[dfavg[pi]['Region'] == regs[ri]][ivar] - 1.645*dfse[pi].loc[dfse[pi]['Region'] == regs[ri]][ivar]
#         maxline = dfavg[pi].loc[dfavg[pi]['Region'] == regs[ri]][ivar] + 1.645*dfse[pi].loc[dfse[pi]['Region'] == regs[ri]][ivar]
#         minline = np.where(minline < 0, 0, np.where(minline > 365, 365, minline))
#         maxline = np.where(maxline < 0, 0, np.where(maxline > 365, 365, maxline))

#         ax1.plot(dfavg[pi].loc[dfavg[pi]['Region'] == regs[ri]]['Year'], dfavg[pi].loc[dfavg[pi]['Region'] == regs[ri]][ivar], linestyle=linestyles[pi], color=linecolors[pi], linewidth=0.8, label=linelabels[pi])
#         ax1.fill_between(dfavg[pi].loc[dfavg[pi]['Region'] == regs[ri]]['Year'],
#                          minline, maxline,
#                          color=linecolors[pi], alpha=0.1, linewidth=0)

#     # Weighted Mean
#     minline = wavgdf.loc[wavgdf['Region'] == regs[ri]][ivar] - 1.645*wsedf.loc[wsedf['Region'] == regs[ri]][ivar]
#     maxline = wavgdf.loc[wavgdf['Region'] == regs[ri]][ivar] + 1.645*wsedf.loc[wsedf['Region'] == regs[ri]][ivar]
#     minline = np.where(minline < 0, 0, np.where(minline > 365, 365, minline))
#     maxline = np.where(maxline < 0, 0, np.where(maxline > 365, 365, maxline))

#     ax1.plot(wavgdf.loc[wavgdf['Region'] == regs[ri]]['Year'], wavgdf.loc[wavgdf['Region'] == regs[ri]][ivar], color='red', linewidth=0.8, label='Weighted Average')
#     ax1.fill_between(wavgdf.loc[wavgdf['Region'] == regs[ri]]['Year'],
#                      minline, maxline,
#                      color='red', alpha=0.3, linewidth=0)

#     # Single Model
#     # dfsp = pdf.loc[(pdf['Region'] == regs[ri]) &  (pdf['Family'] == modspecial[0])].groupby(by=['Year']).mean()
#     # dfsp = dfsp.reset_index()
#     # ax1.plot(dfsp['Year'], dfsp[ivar], color='purple', linestyle='dashed',linewidth=0.8, label=modspecial[0])

#     # Observations
#     ax1.plot(odf[odf['Region'] == regs[ri]]['Year'], odf.loc[odf['Region'] == regs[ri]][ivar], color=obcol, linewidth=0.8, label='PMW Average')

#     txt = ax1.annotate(md.abc[ri]+". "+REGS[ri], xy=(0.02,0.98), xycoords='axes fraction', weight = 'bold', ha = 'left', va = 'top')
#     txt.set_path_effects([pe.withStroke(linewidth=5, foreground='w')])

#     # Axis size/labels
#     ax1.set_ylabel(ylabel[0])
#     ax1.set_xlabel(xlabel[0])
#     ax1.set_xlim(xmin = xlimits[0][0], xmax = xlimits[0][1])
#     ax1.set_ylim(ymin=0, ymax=365)

#     # Legend
#     if ri == 0:
#         ax1.legend(loc='lower right', facecolor='0.98', fontsize=9)


#     # Grid
#     ax1.yaxis.grid(linestyle='dashed', linewidth=0.7, zorder=0)
#     ax1.set_yticks(np.arange(0,360,60))

# plt.tight_layout(rect=[0,0,1,1])

# plt.savefig(figpath+"/Fig_compare_many_multi-model_means_"+typ+"-"+thresh+"_"+regvar+"_"+exp1+"-"+exp2+"_Train"+YYT+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+"_timeseries_2by1_V2.png",dpi=300)

'''*******************************************
Plot #2: Sea Ice v. Temperature - Envelope
*******************************************'''

# Sort by temperature
wavgdf = wavgdf.sort_values(by=['taglobal'])
wsedf = wsedf.sort_values(by=['taglobal'])
odf = odf.sort_values(by=['taglobal'])

# Subset to tempeatures of interest
odf = odf.loc[(odf['taglobal'] >= tamin)]
wsedf = wsedf.loc[(wavgdf['taglobal'] >= tamin)]
wavgdf = wavgdf.loc[(wavgdf['taglobal'] >= tamin)]

### Make Equal-Weighted & Bias-Corrected Ensembles ###
# pdfs = [pdf.loc[np.in1d(pdf['Model'],mods1)], pdf.loc[np.in1d(pdf['Model'],mods2)], pdfg.loc[np.in1d(pdfg['Model'],mods1)], pdfg.loc[np.in1d(pdfg['Model'],mods2)], pdf.loc[np.in1d(pdf['Family'],modspecial)]]
pdfs = [pdf.loc[np.in1d(pdf['Model'],mods1)], pdf.loc[np.in1d(pdf['Model'],mods2)], pdfg.loc[np.in1d(pdfg['Model'],mods1)], pdfg.loc[np.in1d(pdfg['Model'],mods2)]]
dfavg, dfse = [], []

for df in pdfs:
    # Mean
    dfavg1 = df.groupby(by=['Region','Year']).mean()
    dfavg1 = dfavg1.reset_index()
    dfavg1 = dfavg1.sort_values(by=['taglobal'])
    dfavg.append(dfavg1.loc[(dfavg1['taglobal'] >= tamin)])

    # Standard Deviation
    dfsd1 = df.groupby(by=['Region','Year']).std(ddof=1)
    dfsd1 = dfsd1.reset_index()
    dfsd1 = dfsd1.sort_values(by=['taglobal'])
    dfse.append(dfsd1.loc[(dfavg1['taglobal'] >= tamin)])


# # Prep Figure
# fig = plt.figure(figsize=(179/25.4,8))


# for ri in range(len(regs)):
#     ax5 = fig.add_subplot(2,1,1+ri)

#     for pi in [2]:
#         minline = dfavg[pi].loc[dfavg[pi]['Region'] == regs[ri]][ivar] - 1.645*dfse[pi].loc[dfse[pi]['Region'] == regs[ri]][ivar]
#         maxline = dfavg[pi].loc[dfavg[pi]['Region'] == regs[ri]][ivar] + 1.645*dfse[pi].loc[dfse[pi]['Region'] == regs[ri]][ivar]
#         minline = np.where(minline < 0, 0, np.where(minline > 365, 365, minline))
#         maxline = np.where(maxline < 0, 0, np.where(maxline > 365, 365, maxline))

#         ax5.plot(dfavg[pi].loc[dfavg[pi]['Region'] == regs[ri]]['taglobal'], dfavg[pi].loc[dfavg[pi]['Region'] == regs[ri]][ivar], linestyle=linestyles[pi], color=linecolors[pi], linewidth=0.8, label=linelabels[pi],zorder=5)
#         ax5.fill_between(dfavg[pi].loc[dfavg[pi]['Region'] == regs[ri]]['taglobal'],
#                          minline, maxline,
#                          color=linecolors[pi], alpha=0.1, linewidth=0)

#     # Weighted Mean
#     minline = wavgdf.loc[wavgdf['Region'] == regs[ri]][ivar] - 1.645*wsedf.loc[wsedf['Region'] == regs[ri]][ivar]
#     maxline = wavgdf.loc[wavgdf['Region'] == regs[ri]][ivar] + 1.645*wsedf.loc[wsedf['Region'] == regs[ri]][ivar]
#     minline = np.where(minline < 0, 0, np.where(minline > 365, 365, minline))
#     maxline = np.where(maxline < 0, 0, np.where(maxline > 365, 365, maxline))

#     ax5.plot(wavgdf.loc[wavgdf['Region'] == regs[ri]]['taglobal'], wavgdf.loc[wavgdf['Region'] == regs[ri]][ivar], color='red', linewidth=0.8, label='Weighted Average',zorder=5)
#     ax5.fill_between(wavgdf.loc[wavgdf['Region'] == regs[ri]]['taglobal'],
#                      minline, maxline,
#                      color='red', alpha=0.3, linewidth=0)

#     # Single Model
#     # ax5.plot(dfavg[pi+1].loc[dfavg[pi+1]['Region'] == regs[ri]]['taglobal'], dfavg[pi+1].loc[dfavg[pi+1]['Region'] == regs[ri]][ivar], color='purple', linestyle='dashed',linewidth=0.8, label=modspecial[0])

#     # Observations
#     # interpolator = interpolate.interp1d(odf[odf['Region'] == regs[ri]]['taglobal'], odf.loc[odf['Region'] == regs[ri]][ivar], kind='linear')
#     # newx = np.arange(0.3,5.01,0.01)
#     # newy = md.movingAverage(interpolator(newx),n=5)
#     ax5.plot(odf[odf['Region'] == regs[ri]]['taglobal'], odf.loc[odf['Region'] == regs[ri]][ivar], color=obcol, linewidth=1.1, label='PMW Average')
#     # ax5.scatter(odf[odf['Region'] == regs[ri]]['taglobal'], odf.loc[odf['Region'] == regs[ri]][ivar], color=obcol, s=10, label='PMW Average',zorder=4)

#     # Note Times for certain temperatures
#     t1 = odf.loc[(odf['Year'] >= 1980) & (odf['Year'] <= 1989),'taglobal'].mean()
#     t2 = odf.loc[(odf['Year'] >= 2012) & (odf['Year'] <= 2021),'taglobal'].mean()
#     ax5.plot([t1,t1],[0,365], color='0.6', linestyle='dashed', linewidth=1)
#     ax5.plot([t2,t2],[0,365], color='0.6', linestyle='dashed', linewidth=1)
#     txt = [ax5.annotate('1980-1989 Avg.', xy=(t1-0.05,260), ha='right',va='center',size=8,rotation=90),
#            ax5.annotate('2012-2021 Avg.', xy=(t2-0.05,260), ha='right',va='center',size=8,rotation=90)]
#     [t.set_path_effects([pe.withStroke(linewidth=2, foreground='w')]) for t in txt]

#     # Annotation
#     ax5.set_ylabel(ylabel[1])
#     ax5.set_xlabel(xlabel[1])
#     ax5.set_xlim(xmin = xlimits[1][0], xmax = xlimits[1][1])
#     ax5.set_xticks(np.arange(0,xlimits[1][1], 0.5))

#     txt = ax5.annotate(md.abc[ri] + ". " + REGS[ri], xy=(0.02,0.98), xycoords='axes fraction', weight = 'bold', ha = 'left', va = 'top')
#     txt.set_path_effects([pe.withStroke(linewidth=5, foreground='w')])

#     ax5.set_ylim(ymin=0, ymax=365)
#     ax5.set_yticks(np.arange(0,360,60))

#     if ri == 0:
#         ax5.legend(loc='lower right', facecolor='0.98', fontsize=9)

#     # Grid
#     ax5.yaxis.grid(linestyle='dashed', linewidth=0.7, zorder=0)


# plt.tight_layout(rect=[0,0,1,1])

# plt.savefig(figpath+"/Fig_compare_many_multi-model_means_"+typ+"-"+thresh+"_"+regvar+"_"+exp1+"-"+exp2+"_Train"+YYT+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+"_tsens_2by1_V3.png",dpi=300)



########
ts = [1,1.5,2,2.5,3,3.5,4,4.5,5]
ivars = ['OPCavg','LRDavg','FADavg']

reglist, pilist, tlist, ilist = [], [], [], [[],[],[]]
for reg in regs:
    for t in ts:
        for pi in range(len(pdfs)):

            reglist.append(reg)
            pilist.append(linelabels[pi])
            tlist.append(t)
            
            for i, ivar in enumerate(ivars):
                ilist[i].append( np.mean( dfavg[pi].loc[(dfavg[pi]['Region'] == reg) & (dfavg[pi]['taglobal'] >= t-0.1) & (dfavg[pi]['taglobal'] <= t+0.1)][ivar] ) )

        reglist.append(reg)
        pilist.append('Weighted Mean')
        tlist.append(t)
        for i, ivar in enumerate(ivars):
            ilist[i].append( np.mean( wavgdf.loc[(wavgdf['Region'] == reg) & (wavgdf['taglobal'] >= t-0.1) & (wavgdf['taglobal'] <= t+0.1)][ivar] ) )

stats = pd.DataFrame({'Region':reglist,'Version':pilist,'taglobal':tlist})
for i, ivar in enumerate(ivars):
    stats[ivar] = ilist[i]

stats['LRDavg'] = np.where(stats['LRDavg'] > stats['FADavg'], stats['LRDavg']-365, stats['LRDavg'])
stats['OPCcalc'] = stats['FADavg'] - stats['LRDavg']

stats.to_csv(figpath+"/Fig_compare_many_multi-model_means_"+typ+"-"+thresh+"_"+regvar+"_"+exp1+"-"+exp2+"_Train"+YYT+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+".csv",index=False)
