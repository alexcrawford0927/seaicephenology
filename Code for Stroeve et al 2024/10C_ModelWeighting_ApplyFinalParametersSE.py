#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date Created: 4 Apr 2023
Date Modified: 31 May 2023

Applies the results of the model weighting script to regional sea ice phenology
data to calculate the uncertainty in the weighted multi-model means.
"""

'''*******************************************
Load Modules and Functions
*******************************************'''

import numpy as np
import pandas as pd
import CycloneModule_13_2 as md

'''*******************************************
Define Inputs
*******************************************'''

# All Data
V = 'V9'
V2 = V+'.57' # V+'.1' #
typ = 'Thickness' # ''
thresh = '10cm'
regvar = 'regHB3'

ivars = ['OPCavg','LRDavg','FADavg','CIPavg']
tvars = ['tas']

sigma_d, sigma_s = 0.49, 0.50 # 0.23, 0.5

experiment1, experiment2 = 'historical', 'ssp585'
ymin, ymax = [1920, 2015], [2013, 2099]
ytmin, ytmax = 1979, 2021
yamin, yamax = 1850, 1900

# Paths
path = '/Volumes/Cassandra/CMIP6/RegionalStats/'
inpath = path+"/ThicknessPhenology/"+regvar

'''*******************************************
Loading Data
*******************************************'''
YYT = str(ytmin)+'-'+str(ytmax)

### Load Weights ###
wdf = pd.read_csv(path+"/"+V+"/ModelWeighting/Weights_"+typ+"-"+thresh+"_"+regvar+"_"+experiment1+'-'+experiment2+"_Train"+YYT+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+".csv")
wdf['Model'] = wdf['Family'] + '_' + wdf['Member'].astype(str)
mods = np.array(wdf['Model'])

### Load Sea Ice Data ###
# Load CMIP6 Data And Subset
c6ilist = md.listdir(inpath)
c6ilist = [f for f in c6ilist if f.split('_')[3].split('Regional-')[-1] == regvar and f.split('_')[2] == thresh and f.split('_')[-2]+"_"+f.split('_r')[-1].split('i')[0] in mods]


### Load Temperature Data ###
# Identify the files to load
c6tlist = md.listdir(path+"/tas/")
c6glist = [f for f in c6tlist if f.split('Regional-')[1].split('_')[0] == 'global' and f.split('_')[4]+"_"+f.split('_r')[-1].split('i')[0] in mods]
c6tlist = [f for f in c6tlist if f.split('Regional-')[1].split('_')[0] == regvar and f.split('_')[4]+"_"+f.split('_r')[-1].split('i')[0] in mods]

### Idenitfy Regions ###
idf = pd.read_csv(inpath+"/"+c6ilist[0])
regs = np.unique(idf['Region'])
del idf

'''*******************************************
Apply Weighting Scheme -- BY YEAR
*******************************************'''

# Record the year, region, *weighted* global temperature & anomaly,
# regional temperature & anomaly, OPC, LRD, FAD, and CIP for the weighted multi-model mean

for e, exp in enumerate((experiment1,experiment2)):

    # Years
    years = np.arange(ymin[e],ymax[e]+1)

    # Load weighted means
    outdf = pd.read_csv(path+"/"+V+"/ModelWeighting/AnnualWeightedMean_"+exp+"_"+str(ymin[e])+'_'+str(ymax[e])+"_"+typ+"-"+thresh+"_"+regvar+"_"+experiment1+'-'+experiment2+"_Train"+YYT+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+".csv")

    # Initiate output data frame
    sddf = pd.DataFrame({'Region':outdf['Region'],'Year':outdf['Year'],'weight':outdf['weight']})

    # Add columns for weight and other variables with a default value of 0
    for var in ['tglobal','taglobal','treg','tareg']+ivars:
        sddf[var] = np.zeros_like(np.array(outdf['Year']))

    # Subset list of files to this experiment
    c6ilist2 = [f for f in c6ilist if f.split('_')[-3] == exp]
    c6glist2 = [f for f in c6glist if f.split('_')[-3] == exp]
    c6tlist2 = [f for f in c6tlist if f.split('_')[-3] == exp]

    c6gblist2 = [f for f in c6glist if f.split('_')[-3] == 'historical']
    c6tblist2 = [f for f in c6tlist if f.split('_')[-3] == 'historical']

    if len(c6ilist2) == len(c6tlist2):
        c6ilist2.sort()
        c6glist2.sort()
        c6tlist2.sort()
        c6gblist2.sort()
        c6tblist2.sort()

        # Process by Model
        for c in range(len(c6ilist2)):
            # Identify model weight
            w = float(wdf.loc[ (wdf['Family'] == c6ilist2[c].split('_')[-2]) & (wdf['Member'] == int(c6ilist2[c].split('_r')[-1].split('i')[0])),'Weight'])

            # Load CSVs
            idf = pd.read_csv(inpath+"/"+c6ilist2[c])
            gdf = pd.read_csv(path+"/tas/"+c6glist2[c])
            tdf = pd.read_csv(path+'/tas/'+c6tlist2[c])
            gbdf = pd.read_csv(path+"/tas/"+c6gblist2[c])
            tbdf = pd.read_csv(path+"/tas/"+c6tblist2[c])

            # Sort by Region, then Year
            idf = idf.sort_values(by=['Region','Year'])

            # Identify Baseline Temperature
            tdfavg = tbdf.loc[(tbdf['Year'] >= yamin) & (tbdf['Year'] <= yamax),['Region','tas']].groupby(by=['Region']).mean()
            gdfavg = gbdf.loc[(gbdf['Year'] >= yamin) & (gbdf['Year'] <= yamax),['Region','tas']].groupby(by=['Region']).mean()

            # Subset by years (and take annual averages for temperature)
            idf = idf.loc[(idf['Year'] >= ymin[e]) & (idf['Year'] <= ymax[e])]
            gdf2 = gdf.loc[(gdf['Year'] >= ymin[e]) & (gdf['Year'] <= ymax[e]),['Region','Year','tas']].groupby(by=['Region','Year']).mean()
            tdf2 = tdf.loc[(tdf['Year'] >= ymin[e]) & (tdf['Year'] <= ymax[e]),['Region','Year','tas']].groupby(by=['Region','Year']).mean()

            gdf2 = gdf2.reset_index()
            tdf2 = tdf2.reset_index()

            # Calculate annual anomalies for temperature data
            gdf2['tasa'] = gdf2['tas'] - np.float(gdfavg['tas'])
            tdf2['tasa'] = tdf2['tas'] - np.repeat( np.array(tdfavg['tas']), (1+ymax[e]-ymin[e]) )

            ### ADD TO OUTPUT PDF ###
            # Weighted Average of Global Temperature (which needs to be replicated for each region)

            sddf['tglobal'] += np.square(np.repeat( np.array(gdf2['tas']), len(regs) ).reshape(len(years), len(regs)).T.flatten() - np.array(outdf['tglobal'])) * w
            sddf['taglobal'] += np.square(np.repeat( np.array(gdf2['tasa']), len(regs) ).reshape(len(years), len(regs)).T.flatten() - np.array(outdf['taglobal'])) * w

            # Weighted Average of All Other Variables (no replications needed)
            sddf['treg'] += np.square( np.array(tdf2['tas']) - np.array(outdf['treg'])) * w
            sddf['tareg'] += np.square( np.array(tdf2['tasa']) - np.array(outdf['tareg'])) * w

            for ivar in ivars:
                sddf[ivar] += np.square( np.array(idf[ivar]) - np.array(outdf[ivar]) ) * w

    # Write to File

    for var in ['tglobal','taglobal','treg','tareg']+ivars:
        sddf[var] = np.sqrt( sddf[var] / ( sddf['weight'] * (len(mods)-1) / len(mods) ) )

    sddf.to_csv(path+"/"+V+"/ModelWeighting/AnnualWeightedSE_"+exp+"_"+str(ymin[e])+'_'+str(ymax[e])+"_"+typ+"-"+thresh+"_"+regvar+"_"+experiment1+'-'+experiment2+"_Train"+YYT+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+".csv",index=False)
