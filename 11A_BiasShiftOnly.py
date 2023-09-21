#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Alex Crawford
Date Created: 3 Aug 2023
Date Modified: 3 Aug 2023
    
Purpose: Applies a delta-shift bias correction on sea ice phenology data using
the time, regional temperature, or global temperature anomaly to define the period
of overlap.
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
V2 = ''
typ = 'Thickness' # ''
thresh = '10cm'
regvar = 'regHB3'

ivars = ['OPCavg','LRDavg','FADavg']
tvars = ['tas']


experiment1, experiment2 = 'historical', 'ssp585'
ymin, ymax = [1920, 2015], [2013, 2099]
yamin, yamax = 1850, 1900
ybmin, ybmax = 1979, 2021

# Paths
path = '/Volumes/Cassandra/CMIP6/RegionalStats/'
inpath = path+"/ThicknessPhenology/"+regvar

obpath = ['/Volumes/Miranda/SeaIce/Bootstrap/AdvanceRetreat2/C10/Bootstrap_Regionalized_regHB3_Historical_C10.csv','/Volumes/Miranda/SeaIce/NASATeam/AdvanceRetreat2/C10/NASATeam_Regionalized_regHB3_Historical_C10.csv']
otpath = '/Volumes/Theseus/SurfaceTemperature/BEST/RegionalStats/tas/BEST_tasRegionalized_'+regvar+'_1979-2021.csv'
ogpath = '/Volumes/Theseus/SurfaceTemperature/BEST/BEST_Annual_ts_1850-2022.csv'

mods = np.array(['ACCESS-CM2_1', 'ACCESS-CM2_2', 'ACCESS-CM2_3', 'AWI-CM-1-1-MR_1',
        'BCC-CSM2-MR_1', 'CESM2-WACCM_1', 'CESM2-WACCM_2', 'CESM2-WACCM_3',
        'CESM2_10', 'CESM2_11', 'CESM2_4', 'CMCC-CM2-SR5_1', 'CMCC-ESM2_1',
        'CNRM-CM6-1_1', 'CNRM-ESM2-1_1', 'EC-Earth3-CC_1',
        'IPSL-CM6A-LR_1', 'IPSL-CM6A-LR_14', 'IPSL-CM6A-LR_2',
        'IPSL-CM6A-LR_3', 'IPSL-CM6A-LR_33', 'IPSL-CM6A-LR_4',
        'IPSL-CM6A-LR_6', 'KIOST-ESM_1', 'MIROC-ES2L_1', 'MIROC-ES2L_10',
        'MIROC-ES2L_2', 'MIROC-ES2L_3', 'MIROC-ES2L_4', 'MIROC-ES2L_5',
        'MIROC-ES2L_6', 'MIROC-ES2L_7', 'MIROC-ES2L_8', 'MIROC-ES2L_9',
        'MIROC6_1', 'MIROC6_2', 'MIROC6_3', 'MPI-ESM1-2-HR_1',
        'MPI-ESM1-2-HR_2', 'MPI-ESM1-2-LR_1', 'MPI-ESM1-2-LR_10',
        'MPI-ESM1-2-LR_11', 'MPI-ESM1-2-LR_12', 'MPI-ESM1-2-LR_13',
        'MPI-ESM1-2-LR_14', 'MPI-ESM1-2-LR_15', 'MPI-ESM1-2-LR_16',
        'MPI-ESM1-2-LR_17', 'MPI-ESM1-2-LR_18', 'MPI-ESM1-2-LR_19',
        'MPI-ESM1-2-LR_2', 'MPI-ESM1-2-LR_20', 'MPI-ESM1-2-LR_21',
        'MPI-ESM1-2-LR_22', 'MPI-ESM1-2-LR_23', 'MPI-ESM1-2-LR_24',
        'MPI-ESM1-2-LR_25', 'MPI-ESM1-2-LR_26', 'MPI-ESM1-2-LR_27',
        'MPI-ESM1-2-LR_28', 'MPI-ESM1-2-LR_29', 'MPI-ESM1-2-LR_3',
        'MPI-ESM1-2-LR_30', 'MPI-ESM1-2-LR_4', 'MPI-ESM1-2-LR_5',
        'MPI-ESM1-2-LR_6', 'MPI-ESM1-2-LR_7', 'MPI-ESM1-2-LR_8',
        'MPI-ESM1-2-LR_9', 'MRI-ESM2-0_1', 'MRI-ESM2-0_2', 'MRI-ESM2-0_3',
        'MRI-ESM2-0_4', 'MRI-ESM2-0_5', 'NESM3_1', 'NESM3_2',
        'NorESM2-LM_1', 'NorESM2-MM_1'])

# mods = np.array(['ACCESS-CM2_1',  'AWI-CM-1-1-MR_1', 'BCC-CSM2-MR_1', 'CESM2-WACCM_1',
#        'CESM2_4', 'CMCC-CM2-SR5_1', 'CMCC-ESM2_1','CNRM-CM6-1_1', 'CNRM-ESM2-1_1',
#        'EC-Earth3-CC_1', 'IPSL-CM6A-LR_1', 'KIOST-ESM_1', 'MIROC-ES2L_1',
#        'MIROC6_1','MPI-ESM1-2-HR_1','MPI-ESM1-2-LR_1','MRI-ESM2-0_1',
#        'NESM3_1', 'NorESM2-LM_1', 'NorESM2-MM_1'])

'''*******************************************
Loading Data
*******************************************'''
### Load Sea Ice Data ###
# Load Observations
# Load Observations - Sea Ice (taking multi-method mean)
odf = pd.DataFrame([])
for o in obpath:
    odf = pd.concat((odf,pd.read_csv(o)))
odfavg = odf.groupby(by=['Region']).mean()
odfavg = odfavg.reset_index()

# Identify the files to load
c6ilist = md.listdir(inpath)
c6ilist = [f for f in c6ilist if f.split('_')[3].split('Regional-')[-1] == regvar and f.split('_')[2] == thresh and f.split('_')[-2]+"_"+f.split('_r')[-1].split('i')[0] in mods]

### Load Temperature Data ###
# Load Observations
otdf = pd.read_csv(otpath) # Regional Temperature
otdf = otdf.loc[(otdf['Year'] >= ybmin) & (otdf['Year'] <= ybmax)].groupby(by=['Region','Year']).mean()
otdf = otdf.reset_index()

ogdf = pd.read_csv(ogpath) # Global Temperature Anomaly
ogdf['taglobal'] = ogdf['tglobal'] - np.mean( ogdf[(ogdf['Year'] >= yamin) & (ogdf['Year'] <= yamax)]['tglobal'] )

# Identify the files to load
c6tlist = md.listdir(path+"/tas/")
c6glist = [f for f in c6tlist if f.split('Regional-')[1].split('_')[0] == 'global' and f.split('_')[4]+"_"+f.split('_r')[-1].split('i')[0] in mods]
c6tlist = [f for f in c6tlist if f.split('Regional-')[1].split('_')[0] == regvar and f.split('_')[4]+"_"+f.split('_r')[-1].split('i')[0] in mods]

### Identify Regions ###
regs = np.unique(odf['Region'])

'''*******************************************
Load Data By Year
*******************************************'''

# Record the year, region, *weighted* global temperature & anomaly,
# regional temperature & anomaly, OPC, LRD, FAD, and CIP for the weighted multi-model mean

c6df = pd.DataFrame()

for e, exp in enumerate((experiment1,experiment2)):

    # prep arrays for years and regions
    years = np.arange(ymin[e],ymax[e]+1)
    regout = np.repeat(regs,len(years))
    yearout = np.repeat(years, len(regs)).reshape(len(years), len(regs)).T.flatten()

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

        outarrs = [np.empty(shape=(0)) for i in range(8+len(ivars))]

        # Process by Model
        for c in range(len(c6ilist2)):
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

            ### Concatenate everything to output arrays ###
            # Model Family and Member
            outarrs[0] = np.concatenate( ( outarrs[0], np.repeat(mods[c].split('_')[0],len(regout)) ), axis = 0 )
            outarrs[1] = np.concatenate( ( outarrs[1], np.repeat(mods[c].split('_')[1],len(regout)) ), axis = 0 )

            # Regions and Years
            outarrs[2], outarrs[3] = np.concatenate( (outarrs[2],regout), axis = 0), np.concatenate( (outarrs[3],yearout), axis = 0 )

            # Global Temperature
            outarrs[4] = np.concatenate( (outarrs[4], np.repeat(np.array(gdf2['tas']), len(regs)).reshape(len(years), len(regs)).T.flatten()), axis = 0 )
            outarrs[5] = np.concatenate( (outarrs[5], np.repeat(np.array(gdf2['tasa']), len(regs)).reshape(len(years), len(regs)).T.flatten()), axis = 0 )

            # Regional Temperature
            outarrs[6] = np.concatenate( (outarrs[6], np.array(tdf2['tas'])), axis = 0 )
            outarrs[7] = np.concatenate( (outarrs[7], np.array(tdf2['tasa'])), axis = 0 )

            for i, ivar in enumerate(ivars):
                outarrs[8+i] = np.concatenate( (outarrs[8+i], np.array(idf[ivar]) ), axis = 0 )

        # Combine
        pdf = pd.DataFrame(outarrs).T
        pdf.columns = ['Family','Member','Region','Year','tglobal','taglobal','treg','tareg']+ivars

        c6df = c6df.append( pdf, ignore_index=True )

# Combine Experiments
c6df = c6df.sort_values(by=['Family','Member','Region','Year'])
for col in c6df.columns[1:]:
    c6df[col] = c6df[col].astype(float)

c6df['LRDavg'] = np.where(c6df['LRDavg'] > 450, c6df['LRDavg'] - 365, c6df['LRDavg'])

'''*******************************************
Apply Bias Correction
*******************************************'''

### IDENTIFY AVERAGES ###

# Identify averages for the time period of overlap
c6dfb = c6df.loc[(c6df['Year'] >= ybmin) & (c6df['Year'] <= ybmax)].groupby(by=['Family','Member','Region']).mean()
c6dfb = c6dfb.reset_index()

# Identify averages for the global temperature anomaly of overlap
c6dfg = c6df.loc[ ( c6df['taglobal'] >=  np.min(ogdf[(ogdf['Year'] >= ybmin) & (ogdf['Year'] <= ybmax)]['taglobal']) ) & ( c6df['taglobal'] <= np.max(ogdf[(ogdf['Year'] >= ybmin) & (ogdf['Year'] <= ybmax)]['taglobal']) )].groupby(by=['Family','Member','Region']).mean()
c6dfg = c6dfg.reset_index()

# Identify averages for the regional temperature overlap
c6dft = pd.DataFrame()
for reg in regs:
    c6dftr = c6df.loc[ ( c6df['Region'] == reg ) & ( c6df['treg'] >=  np.min(otdf[(otdf['Region'] == reg) & (otdf['Year'] >= ybmin) & (otdf['Year'] <= ybmax)]['tas']) ) & ( c6df['treg'] <= np.max(otdf[(otdf['Region'] == reg) & (otdf['Year'] >= ybmin) & (otdf['Year'] <= ybmax)]['tas']) )].groupby(by=['Family','Member','Region']).mean()
    c6dftr = c6dftr.reset_index()
    c6dft = c6dft.append(c6dftr, ignore_index=True)

for col in c6dft.columns[3:]:
    c6dft[col] = c6dft[col].astype(float)

### CALCULATE BIAS ###

# By Time
biasb = c6dfb.iloc[:,:3]
for ivar in ivars:
    biasb[ivar] = np.nan
    for reg in regs:
        biasb.loc[biasb['Region'] == reg,ivar] = np.array(c6dfb[c6dfb['Region'] == reg][ivar]) - float( odfavg[odfavg['Region'] == reg][ivar] )


# By Global Temperature Anomaly
biasg = c6dfb.iloc[:,:3]
for ivar in ivars:
    biasg[ivar] = np.nan
    for reg in regs:
        biasg.loc[biasg['Region'] == reg,ivar] = np.array(c6dfg[c6dfg['Region'] == reg][ivar]) - float( odfavg[odfavg['Region'] == reg][ivar] )

# By Regional Temperature
biast = c6dfb.iloc[:,:3]
for ivar in ivars:
    biast[ivar] = np.nan
    for reg in regs:
        biast.loc[biast['Region'] == reg,ivar] = np.array(c6dft[c6dft['Region'] == reg][ivar]) - float( odfavg[odfavg['Region'] == reg][ivar] )

### ADJUST FOR BIAS ###

c6df['Model'] = c6df['Family'] + "_" + c6df['Member'].astype(int).astype(str)

# By Time
c6adjb = c6df.iloc[:,:4]
c6adjb['Model'] = c6df['Model']
for ivar in ivars:
    c6adjb[ivar] = np.nan
    for mod in np.unique(c6df['Model']):
        for reg in regs:
            bias = float(biasb.loc[(biasb['Family'] == mod.split('_')[0]) & (biasb['Member'] == float(mod.split('_')[1])) & (biasb['Region'] == reg)][ivar])
            c6adjb.loc[(c6adjb['Region'] == reg) & (c6adjb['Model'] == mod),ivar] = np.array(c6df.loc[(c6df['Model'] == mod) & (c6df['Region'] == reg),ivar]) - bias

# By Global Temperature Anomaly
c6adjg = c6df.iloc[:,:4]
c6adjg['Model'] = c6df['Model']
for ivar in ivars:
    c6adjg[ivar] = np.nan
    for mod in np.unique(c6df['Model']):
        for reg in regs:
            bias = float(biasg.loc[(biasg['Family'] == mod.split('_')[0]) & (biasg['Member'] == float(mod.split('_')[1])) & (biasg['Region'] == reg)][ivar])
            c6adjg.loc[(c6adjg['Region'] == reg) & (c6adjg['Model'] == mod),ivar] = np.array(c6df.loc[(c6df['Model'] == mod) & (c6df['Region'] == reg),ivar]) - bias

# By Regional Temperature
c6adjt = c6df.iloc[:,:4]
c6adjt['Model'] = c6df['Model']
for ivar in ivars:
    c6adjt[ivar] = np.nan
    for mod in np.unique(c6df['Model']):
        for reg in regs:
            bias = float(biast.loc[(biast['Family'] == mod.split('_')[0]) & (biast['Member'] == float(mod.split('_')[1])) & (biast['Region'] == reg)][ivar])
            c6adjt.loc[(c6adjt['Region'] == reg) & (c6adjt['Model'] == mod),ivar] = np.array(c6df.loc[(c6df['Model'] == mod) & (c6df['Region'] == reg),ivar]) - bias


### Write to File ###
c6adjb.to_csv(path+"/"+V+"/AnnualBiasCorrected_by_Time_"+str(ybmin)+'-'+str(ybmax)+"_"+experiment1+"-"+experiment2+'_'+typ+"-"+thresh+"_"+regvar+".csv",index=False)
c6adjg.to_csv(path+"/"+V+"/AnnualBiasCorrected_by_GlobalTemperatureAnomaly_"+str(ybmin)+'-'+str(ybmax)+"_"+experiment1+"-"+experiment2+'_'+typ+"-"+thresh+"_"+regvar+".csv",index=False)
c6adjt.to_csv(path+"/"+V+"/AnnualBiasCorrected_by_RegionalTemperature_"+str(ybmin)+'-'+str(ybmax)+"_"+experiment1+"-"+experiment2+'_'+typ+"-"+thresh+"_"+regvar+".csv",index=False)

