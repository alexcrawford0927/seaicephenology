#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date Created: 4 April 2023
Date Modified: 31 May 2023

Creates a weighting scheme following the general method of Knutti et al. (2017)
that uses the average and trends in sea ice phenlogy variables and regional 
temperature as inputs. Model simulations receive greater weight if a) they better
match observations and b) they are independent from other models. This script
tests the sensitivity to two input parameters (sigma_d and sigma_s). Next, it
saves one chosen set of weights.
"""

'''*******************************************
Load Modules and Functions
*******************************************'''


import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy.stats import linregress
import CycloneModule_13_2 as md


def modelweights_knutti(models, obs, sigma_d=0.4, sigma_s=0.5):
    '''
    Calculates the weighting for models (intended for a multi-model mean
    projection), giving more weight to models that are close to observations
    and independent of other models for a series of relevant variables. The
    variables can be any physical property (temperature, sea ice concentration,
    sea-level pressure, etc.), can be regional averages or by grid cell, and
    can be a temporal average or trend. It can handle providing multiple
    simulations from the same model since it will just weight each one lower
    for being relatively similar.

    The three key issues to consider are a) which variables to use (including
    possible spatial averaging), b) the reletive importance of models matching
    observation, and c) the relative importance of model independence.

    If you make sigma_d too low, you risk over-fitting to observations. If you
    make sigma_s too high, you risk underweighting models that are similar
    because they both match observations (even though they are independent).
    On the other hand, if you make sigma_d too high or sigma_s too low, you
    risk making the weighting exercise meaningless -- converging toward a
    model democracy that cannot identify model dependence. (In other words,
    you want sigma_d to be as low as possible without over-fitting and you
    want both results to be insenstive to small changes in sigma_d and sigma_s
    values.)

    Method based on Knutti et al. (2017) - https://doi.org/10.1002/2016gl072012

    Parameters
    ----------
    models : numpy array with the shape [M,k]
        M = number of models; k = number of variables used for distance
        calculations -- order of k must match order in obs!
    obs : numpy array with the shape [k] or [1,k]
        k = number of variables used for distance calculations -- order of
        k must match order in models (axis 1)
    sigma_d : observational distance sensitivity, optional
        A value between 0 and 1 is typical -- with larger values converging to
        model democracy and smaller values weighting a few models that better
        match the observations. The default is 0.4.
    sigma_s : model independence sensitivity, optional
        A value between 0 and 1 is typical -- with larger values emphasizing a
        few independent models and smaller values converging toward treating
        all models as independent. The default is 0.5.

    Returns
    -------
    weight : numpy array of size [M]
        Final weight for each model, normalized so that sum of weights is 1.

    '''
    ### Calculate Distance from Observations (D) ###

    # Create a RMSD (root mean square difference) distance metric
    sqdiff =  np.square(models - obs) # squared differences
    rmsd = np.sqrt( np.mean( sqdiff / np.median(sqdiff, axis = 0), axis = 1 ) ) # normalize each difference by its median

    numerator = np.exp( -1 * np.square(rmsd) / np.square(sigma_d) ) # calculate numerator

    ### Calculate Distance from Other Models (S) ###
    denominator = np.zeros_like(models[:,0])
    for i in range(models.shape[0]): # For each model...
        # Create a RMSD (root mean square difference) distance metric
        models_single = np.delete(np.square(models - models[i]), i, axis = 0) # squared differences (then remove model in question)
        rmsd_single = np.sqrt( np.mean( models_single / np.median(models_single, axis = 0), axis = 1 ) ) # normalize each difference by its median

        denominator[i] = 1 + np.sum( np.exp( -1 * np.square(rmsd_single) / np.square(sigma_s) ) ) # calculate denominator (for single model)

    ### Calculate Final Weights (w) ###
    weight = numerator / denominator
    weight /= np.sum(weight)

    return weight

def sd_weighted(values, weights, ddof=1):
    '''
    Takes a weighted standard deviation

    Parameters
    ----------
    values : list of numpy array
        Values for which to the standard deviation.
    weights : list or numpy array (must have same shape as values)
        Weights for values for which to take the standard deviation.
    ddof : integer (optional)
        The difference between the number of values in a sample and the
        degrees of freedom (dof = n - ddof); 1 by default, meaning dof = n - 1

    Returns
    -------
    A single value for the weighted standard deviation

    '''
    mean_weighted = np.sum(values*weights) / np.sum(weights) # the weighted mean
    n_nonzero = np.sum(weights != 0) # The number of obs with non-zero weight

    return np.sqrt( np.sum( np.square(values - mean_weighted) * weights ) / (np.sum(weights)*(n_nonzero-ddof)/n_nonzero ) )

'''*******************************************
Define Inputs
*******************************************'''

# All Data
V = 'V9'
V2 = V+'.5'
typ, thresh = 'Thickness', '10cm' # '', 'C10' #
otyp, othresh, omod = '', 'C10', ['Bootstrap','NASATeam']

siregvar = 'regHB3'
tregvarA = 'regHB3'
tregvarB = 'global'

sigma_d, sigma_s = 0.51, 0.4
z_score = 1.645 #1.96 #

# Training Data
experiment = 'historical-ssp585'
ymin, ymax = 1979, 2021

# Validation Data
experiment2 = 'ssp585'
y2min, y2max = 2070, 2099

# Paths
path = '/Volumes/Cassandra/CMIP6/RegionalStats/'
path2 = '/Volumes/Theseus/SurfaceTemperature/BEST/RegionalStats/tas'
outpath = path+"/"+V+"/ModelWeighting"

modstoskip = np.array([])

modstoskip = np.array([ 'MPI-ESM1-2-HR_8', 'MPI-ESM1-2-HR_9', 'MPI-ESM1-2-HR_10',
  'MPI-ESM1-2-LR_8', 'MPI-ESM1-2-LR_9', 'MPI-ESM1-2-LR_10',
  'MPI-ESM1-2-LR_11', 'MPI-ESM1-2-LR_12', 'MPI-ESM1-2-LR_13', 'MPI-ESM1-2-LR_14',
  'MPI-ESM1-2-LR_15', 'MPI-ESM1-2-LR_16', 'MPI-ESM1-2-LR_17', 'MPI-ESM1-2-LR_18',
  'MPI-ESM1-2-LR_19', 'MPI-ESM1-2-LR_20', 'MPI-ESM1-2-LR_21', 'MPI-ESM1-2-LR_22',
  'MPI-ESM1-2-LR_23', 'MPI-ESM1-2-LR_24', 'MPI-ESM1-2-LR_25', 'MPI-ESM1-2-LR_26',
  'MPI-ESM1-2-LR_27', 'MPI-ESM1-2-LR_28', 'MPI-ESM1-2-LR_29', 'MPI-ESM1-2-LR_30',
  'MIROC-ES2L_8','MIROC-ES2L_9','MIROC-ES2L_10'])

'''*******************************************
Loading Data
*******************************************'''
YY2 = str(y2min)+"-"+str(y2max)
YY = str(ymin)+'-'+str(ymax)
years = np.arange(ymin, ymax+1)

### Load Sea Ice Data ###
# Training Data
sic6df = pd.read_csv(path+'/'+V+'/'+typ+'AvgTrend'+YY+'/CMIP6AvgTrend_'+siregvar+'_'+experiment+'_'+thresh+'_'+V+'.csv')
siobdf = pd.read_csv(path+'/'+V+'/'+otyp+'AvgTrend'+YY+'/ObsAvgTrend_'+siregvar+'_'+experiment+'_'+othresh+'_'+V+'.csv')
siobdf = siobdf.loc[np.in1d(siobdf['Model'], omod)].groupby(by=['Region']).mean()
siobdf = siobdf.reset_index()

sic6df['Model'] = sic6df['Family'] + '_' + sic6df['Member'].astype(str)
sic6df = sic6df.loc[np.in1d(np.array(sic6df['Model']),modstoskip) == 0]
sic6df = sic6df.sort_values(by=['Family','Member'])

# Identify the CMIP6 models
mods = [ f.split('_') for f in np.unique(sic6df['Family'] + '_' + sic6df['Member'].astype(str)) ]

# Validation Data
sivaldf = pd.read_csv(path+'/'+V+'/'+typ+'AvgTrend'+YY2+'/CMIP6AvgTrend_'+siregvar+'_'+experiment2+'_'+thresh+'_'+V+'.csv')
sivaldf['Model'] = sivaldf['Family'] + '_' + sivaldf['Member'].astype(str)
sivaldf = sivaldf.loc[np.in1d(sivaldf['Model'], sic6df['Model'])]
sivaldf = sivaldf.sort_values(by=['Family','Member'])

### Load Temperature Data ###
# Identify the files to load
c6tlist = md.listdir(path+"/tas/")
c6tlist = [f for f in c6tlist if [f.split('_')[4],f.split('_r')[-1].split('i')[0]] in mods]

c6tlistA1 = [f for f in c6tlist if f.split('_')[3] == experiment.split('-')[0] and f.split('Regional-')[1].split('_')[0] == tregvarA]
c6tlistA2 = [f for f in c6tlist if f.split('_')[3] == experiment.split('-')[1] and f.split('Regional-')[1].split('_')[0] == tregvarA]
c6tlistB1 = [f for f in c6tlist if f.split('_')[3] == experiment.split('-')[0] and f.split('Regional-')[1].split('_')[0] == tregvarB]
c6tlistB2 = [f for f in c6tlist if f.split('_')[3] == experiment.split('-')[1] and f.split('Regional-')[1].split('_')[0] == tregvarB]

# Load & Process Files
familylist, memberlist, regionlist, meanlist, trendlist, yintlist, plist, selist, r2list = [[] for i in range(9)]
for model in mods:

    # Load Data
    tdfA1 = pd.read_csv(path+"/tas/"+[f for f in c6tlistA1 if [f.split('_')[4],f.split('_r')[-1].split('i')[0]] == model][0])
    tdfA2 = pd.read_csv(path+"/tas/"+[f for f in c6tlistA2 if [f.split('_')[4],f.split('_r')[-1].split('i')[0]] == model][0])
    tdfB1 = pd.read_csv(path+"/tas/"+[f for f in c6tlistB1 if [f.split('_')[4],f.split('_r')[-1].split('i')[0]] == model][0])
    tdfB2 = pd.read_csv(path+"/tas/"+[f for f in c6tlistB2 if [f.split('_')[4],f.split('_r')[-1].split('i')[0]] == model][0])

    tdf = tdfA1.append(tdfA2,ignore_index=True).append(tdfB1,ignore_index=True).append(tdfB2,ignore_index=True)
    del tdfA1, tdfA2, tdfB1, tdfB2
    tdf = tdf.loc[(tdf['Year'] >= ymin) & (tdf['Year'] <= ymax)]

    # Aggregate by region and year
    tdf = tdf.groupby(by=['Region','Year']).mean()
    tdf = tdf.reset_index()

    # Calculate mean and trend
    for reg in np.unique(tdf['Region']):

        yvar = np.array(tdf[tdf['Region'] == reg]['tas'])

        # Average
        meanlist.append( np.mean(tdf[tdf['Region'] == reg]['tas']) )

        # Trend
        finite = np.isfinite(yvar)
        lm = linregress(years[finite], yvar[finite])

        trendlist.append( lm[0] )
        yintlist.append( lm[1] )
        r2list.append( lm[2]*lm[2] )
        plist.append( lm[3] )
        selist.append( lm[4] )

        # Acessories
        familylist.append(model[0]), memberlist.append(int(model[1]))
        regionlist.append(reg)

tc6df = pd.DataFrame({'Family':familylist,'Member':memberlist,
                        'Region':regionlist,'T_Avg':meanlist,'T_Trend':trendlist,
                        'T_Trend_yint':yintlist,'T_Trend_r2':r2list,
                        'T_Trend_p':plist,'T_Trend_se':selist})
tc6df = tc6df.sort_values(by=['Family','Member'])

# Load Observational Temperature Data
tobdfA = pd.read_csv(path2+'/'+[f for f in md.listdir(path2) if f.split('_')[-2] == tregvarA][0])
tobdfB = pd.read_csv(path2+'/'+[f for f in md.listdir(path2) if f.split('_')[-2] == tregvarB][0])
tobdf0 = tobdfA.append(tobdfB,ignore_index=True)

# Subset/aggregate to annual data
tobdf0 = tobdf0.loc[(tobdf0['Year'] >= ymin) & (tobdf0['Year'] <= ymax)]
tobdf0 = tobdf0.groupby(by=['Region','Year']).mean()
tobdf0 = tobdf0.reset_index()

# Calculate Mean/Trend
familylist, memberlist, regionlist, meanlist, trendlist, yintlist, plist, selist, r2list = [[] for i in range(9)]

for reg in np.unique(tdf['Region']):

    yvar = np.array(tobdf0[tobdf0['Region'] == reg]['tas'])

    # Average
    meanlist.append( np.mean(tobdf0[tobdf0['Region'] == reg]['tas']) )

    # Trend
    finite = np.isfinite(yvar)
    lm = linregress(years[finite], yvar[finite])

    trendlist.append( lm[0] )
    yintlist.append( lm[1] )
    r2list.append( lm[2]*lm[2] )
    plist.append( lm[3] )
    selist.append( lm[4] )

    # Acessories
    familylist.append('BEST'), memberlist.append(1)
    regionlist.append(reg)

tobdf = pd.DataFrame({'Family':familylist,'Member':memberlist,
                        'Region':regionlist,'T_Avg':meanlist,'T_Trend':trendlist,
                        'T_Trend_yint':yintlist,'T_Trend_r2':r2list,
                        'T_Trend_p':plist,'T_Trend_se':selist})


'''*******************************************
Design Weighting Scheme
*******************************************'''

### Set up inputs to weighting function ###

# The observational reference must be an array of shape [1,k], where k is the
# number of variables (dimensions/diagnostics) being used
if V2.split('.')[-1] in ['1','2']:
    obs = np.array([siobdf[siobdf['Region'] == 41]['OPCavg_Avg'], siobdf[siobdf['Region'] == 42]['OPCavg_Avg'],
                    siobdf[siobdf['Region'] == 41]['OPCavg_TSens'], siobdf[siobdf['Region'] == 42]['OPCavg_TSens'],
                    tobdf[tobdf['Region'] == 40]['T_Avg'], tobdf[tobdf['Region'] == -1]['T_Avg'],
                    tobdf[tobdf['Region'] == 40]['T_Trend'], tobdf[tobdf['Region'] == -1]['T_Trend'] ]).T

else:
    obs = np.array([siobdf[siobdf['Region'] == 41]['OPCavg_Avg'], siobdf[siobdf['Region'] == 42]['OPCavg_Avg'],
                    siobdf[siobdf['Region'] == 41]['OPCavg_TSens'], siobdf[siobdf['Region'] == 42]['OPCavg_TSens'],
                    siobdf[siobdf['Region'] == 41]['LRDavg_Avg'], siobdf[siobdf['Region'] == 42]['LRDavg_Avg'],
                    siobdf[siobdf['Region'] == 41]['LRDavg_TSens'], siobdf[siobdf['Region'] == 42]['LRDavg_TSens'],
                    siobdf[siobdf['Region'] == 41]['FADavg_Avg'], siobdf[siobdf['Region'] == 42]['FADavg_Avg'],
                    siobdf[siobdf['Region'] == 41]['FADavg_TSens'], siobdf[siobdf['Region'] == 42]['FADavg_TSens'],
                    tobdf[tobdf['Region'] == 41]['T_Avg'], tobdf[tobdf['Region'] == 42]['T_Avg'],tobdf[tobdf['Region'] == -1]['T_Avg'],
                    tobdf[tobdf['Region'] == 41]['T_Trend'], tobdf[tobdf['Region'] == 42]['T_Trend'], tobdf[tobdf['Region'] == -1]['T_Trend'] ]).T

# The CMIP6 data must be an array of M models by k variables (shape is [M,k])
if V2.split('.')[-1] in ['1','2']:
    cmip6 = np.array([ sic6df[sic6df['Region'] == 41]['OPCavg_Avg'], sic6df[sic6df['Region'] == 42]['OPCavg_Avg'],
                        sic6df[sic6df['Region'] == 41]['OPCavg_TSens'], sic6df[sic6df['Region'] == 42]['OPCavg_TSens'],
                        tc6df[tc6df['Region'] == 40]['T_Avg'], tc6df[tc6df['Region'] == -1]['T_Avg'],
                        tc6df[tc6df['Region'] == 40]['T_Trend'], tc6df[tc6df['Region'] == -1]['T_Trend'] ]).T
else:
    cmip6 = np.array([ sic6df[sic6df['Region'] == 41]['OPCavg_Avg'], sic6df[sic6df['Region'] == 42]['OPCavg_Avg'],
                        sic6df[sic6df['Region'] == 41]['OPCavg_TSens'], sic6df[sic6df['Region'] == 42]['OPCavg_TSens'],
                        sic6df[sic6df['Region'] == 41]['LRDavg_Avg'], sic6df[sic6df['Region'] == 42]['LRDavg_Avg'],
                        sic6df[sic6df['Region'] == 41]['LRDavg_TSens'], sic6df[sic6df['Region'] == 42]['LRDavg_TSens'],
                        sic6df[sic6df['Region'] == 41]['FADavg_Avg'], sic6df[sic6df['Region'] == 42]['FADavg_Avg'],
                        sic6df[sic6df['Region'] == 41]['FADavg_TSens'], sic6df[sic6df['Region'] == 42]['FADavg_TSens'],
                        tc6df[tc6df['Region'] == 41]['T_Avg'], tc6df[tc6df['Region'] == 42]['T_Avg'],tc6df[tc6df['Region'] == -1]['T_Avg'],
                        tc6df[tc6df['Region'] == 41]['T_Trend'], tc6df[tc6df['Region'] == 42]['T_Trend'], tc6df[tc6df['Region'] == -1]['T_Trend'] ]).T

# The future CMIP6 data used for validation should be in the same format as the training data, but use np.nan wherever you
## don't want to test that particular variable during validation (i.e., not all training variables are validation variables)
val6 = np.zeros_like(cmip6)*np.nan
val6[:,0] = sivaldf[sivaldf['Region'] == 41]['OPCavg_Avg']
val6[:,1] = sivaldf[sivaldf['Region'] == 42]['OPCavg_Avg']

### Sensitivity Test (using Perfect Model Setup) ###
sigma_ds, sigma_ss = np.round(np.arange(0.1,1.01,0.01),2), np.round(np.arange(0.1,1.01,0.1),2)

correlation, median_bias, median_se, within_range = [[] for i in range(4)]
sigmadlist, sigmaslist = [[] for i in range(2)]
for d in sigma_ds:
    print(d)
    for s in sigma_ss:

        predictions, predictions_se = [], []
        for m in range(cmip6.shape[0]):

            # Calculate Weights
            weights = modelweights_knutti(np.delete(cmip6, m, axis=0), cmip6[m,:], d, s)

            # Use weights to predict the "true" model for the validation period
            predictions.append( np.sum( np.delete(val6, m, axis=0).T * weights.T, axis=1) )
            predictions_se.append( np.apply_along_axis( sd_weighted, 0, np.delete(val6, m, axis=0), weights, ddof=1) )

        # Correlate the predictions and the "true" values
        corr_row = []
        for v in range(val6.shape[1]):
            try:
                corr_row.append( pearsonr(val6[:,v], np.array(predictions)[:,v] )[0] )
            except:
                corr_row.append( np.nan )
        correlation.append (np.array(corr_row) )

        # Calculate the median bias magnitude between predictions and "true" model
        median_bias.append( np.median(np.abs(np.array(predictions) - val6), axis = 0))

        # Calculate the median weighted standard deviation for predictions
        median_se.append( np.median(np.array(predictions_se), axis = 0) )

        # Note if the true value is within error range of prediction
        within_range.append( np.mean(  (val6 <= np.array(predictions) + z_score*np.array(predictions_se)) & (val6 >= np.array(predictions) - z_score*np.array(predictions_se)), axis = 0 ) )

        # Append accessories
        sigmadlist.append(d)
        sigmaslist.append(s)

# Compile data from sensitivity test
weights_sens_test = pd.DataFrame({'sigma_d':sigmadlist, 'sigma_s':sigmaslist })
for v in range(val6.shape[1]):
    if np.isfinite(np.sum(val6[:,v])):
        weights_sens_test['correlation_'+str(v)] = np.array(correlation)[:,v]
        weights_sens_test['median_bias_'+str(v)] = np.array(median_bias)[:,v]
        weights_sens_test['median_se_'+str(v)] = np.array(median_se)[:,v]
        weights_sens_test['within_range_'+str(v)] = np.array(within_range)[:,v]

weights_sens_test.to_csv(outpath+"/WeightsSensitivity_"+typ+"-"+thresh+"_"+siregvar+"_"+experiment+"_Train"+YY+"_Validate"+YY2+"_"+V2+".csv",index=False)

'''*******************************************
Save Weighting Scheme
*******************************************'''

# Run Perfect Model Experiment One More Time With Chosen Parameters
predictions, predictions_se = [], []
for m in range(cmip6.shape[0]):

    # Calculate Weights
    weights = modelweights_knutti(np.delete(cmip6, m, axis=0), cmip6[m,:], sigma_d, sigma_s)

    # Use weights to predict the "true" model for the validation period
    predictions.append( np.sum( np.delete(val6, m, axis=0).T * weights.T, axis=1) )
    predictions_se.append( np.apply_along_axis( sd_weighted, 0, np.delete(val6, m, axis=0), weights, ddof=1) )

predictions = np.array(predictions)

# Make PDF for output
output = pd.DataFrame({'Family': np.array(mods)[:,0], 'Member': np.array(mods)[:,1].astype(int)})
for v in range(val6.shape[1]):
    if np.isfinite(np.sum(val6[:,v])):
        output['truth_'+str(v)] = val6[:,v]
        output['predict_'+str(v)] = predictions[:,v]

# Write to file
output.to_csv(outpath+"/PerfectModelExperiment_"+typ+"-"+thresh+"_"+siregvar+"_"+experiment+"_Train"+YY+"_Validate"+YY2+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+".csv",index=False)


### Sensitivity Test (using Observations) ###
# sigma_ds, sigma_ss = np.arange(0.1,1.01,0.05), np.arange(0.1,1.01,0.05)
# weights_sensobs_test = pd.DataFrame({'Family':np.array([f[0] for f in mods]), 'Member':np.array([f[1] for f in mods])})
# for d in sigma_ds:
#     for s in sigma_ss:
#         weights_sensobs_test[str(np.round(d,2))+'-'+str(np.round(s,2))] = modelweights_knutti(cmip6, obs, d, s)

# Final Weights
weights = modelweights_knutti(cmip6, obs, sigma_d, sigma_s)
wdf = pd.DataFrame({'Family': np.array(mods)[:,0],
                    'Member': np.array(mods)[:,1].astype(int),
                    'Weight': weights})
wdf.to_csv(outpath+"/Weights_"+typ+"-"+thresh+"_"+siregvar+"_"+experiment+"_Train"+YY+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+".csv",index=False)
