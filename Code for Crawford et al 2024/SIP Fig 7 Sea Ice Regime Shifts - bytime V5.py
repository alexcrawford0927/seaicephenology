#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 25 Mar 2024
Modified: 03 Jul 2024

Author: Alex Crawford

Purpose: Trends / Regime shifts in sea ice metrics
"""

import pandas as pd
import xarray as xr
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import CycloneModule_13_2 as md
from matplotlib.colors import LinearSegmentedColormap, ListedColormap


siemin, siemax = 8, 13.5
dmin, dmax = 276, 400
dmax2 = 400
ncmin, ncmax = 0, 1
mmin, mmax = 10, 2
colors=['sienna','brown','maroon','red','darkorange','gold','green','teal','blue','purple','violet','hotpink']

doys = [1,32,60,91,121,152,182,213,244,274,305,335,366,397,425]
doyl = ["Jan 1","Feb 1","Mar 1","Apr 1","May 1","Jun 1","Jul 1","Aug 1","Sep 1","Oct 1","Nov 1","Dec 1","Jan 1","Feb 1","Mar 1"]

path = '/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses'
siepath = path+'/SIE_Daily_ArcticWide.csv'
eventpath = path+"/Pauses_combined_final_bytime_"+str(dmin)+"-"+str(dmax)+".csv"

'''*******************************************
Processing
*******************************************'''
sie = pd.read_csv(siepath)
sie['SYear'] = np.where(sie['month'] > 8, sie['year'], sie['year']-1)
sie = sie.loc[sie['SYear'] > 1978]
sie['DOY'] = [md.daysBetweenDates([sie.loc[i,'SYear'],1,0],[sie.loc[i,'year'],sie.loc[i,'month'],sie.loc[i,'day']]) for i in sie.index]

### For each year, identify the day for which the SIE is closest to the siemin that is after the annual minimum extent and
# the year for which the SIE is closet to the siemax that is before the annual maximum extent ###
sie = sie.loc[((sie['month'] > 8) | (sie['month'] < 4)) &  (np.isfinite(sie['sie_nt']))]

annual = pd.DataFrame({'SYear':np.unique(sie['SYear'])[:]})

# Identify the annual minimum and maximum SIE dates
annual['siemin'] = sie.loc[:,('SYear','sie_cdr')].groupby(by='SYear').min()['sie_cdr'].values
annual['siemax'] = sie.loc[:,('SYear','sie_cdr')].groupby(by='SYear').max()['sie_cdr'].values
annual['sieminJD'] = [sie.loc[(sie['SYear'] == annual.iloc[i]['SYear']) & (sie['month'] < 10) & (sie['month'] > 8) & (sie['sie_cdr'] == annual.iloc[i]['siemin']), 'JulianDay'].values[0] for i in range(len(annual))]
annual['siemaxJD'] = [sie.loc[(sie['SYear'] == annual.iloc[i]['SYear']) & (sie['month'] < 5) & (sie['month'] > 1) & (sie['sie_cdr'] == annual.iloc[i]['siemax']), 'JulianDay'].values[0] for i in range(len(annual))]
annual['siemaxmonth'] = [sie.loc[(sie['SYear'] == annual.iloc[i]['SYear']) & (sie['month'] < 5) & (sie['month'] > 1) & (sie['sie_cdr'] == annual.iloc[i]['siemax']), 'month'].values[0] for i in range(len(annual))]
annual['siemaxday'] = [sie.loc[(sie['SYear'] == annual.iloc[i]['SYear']) & (sie['month'] < 5) & (sie['month'] > 1) & (sie['sie_cdr'] == annual.iloc[i]['siemax']), 'day'].values[0] for i in range(len(annual))]

# Identify the first SIE value to exceed the analysis min/max SIE
annual['siestart'] = [sie.loc[(sie['SYear'] == i) & (sie['DOY'] >= 274),'sie_cdr'].values[0] for i in annual['SYear'].values]
annual['sieend'] = [sie.loc[(sie['SYear'] == i) & (sie['DOY'] >= dmax2-3),'sie_cdr'].values[0] for i in annual['SYear'].values]
   
# Identify the start and end dates for the analysis period -- as Julian Day and then as Y, M, D
annual['siestartJD'] = [sie.loc[(sie['SYear'] == annual.iloc[i]['SYear']) & (sie['JulianDay'] >= annual.iloc[i]['sieminJD']) & (sie['sie_cdr'] == annual.iloc[i]['siestart']), 'JulianDay'].values[0] for i in range(len(annual))]
annual['sieendJD'] = [sie.loc[(sie['SYear'] == annual.iloc[i]['SYear']) & (sie['JulianDay'] <= annual.iloc[i]['siemaxJD']) & (sie['sie_cdr'] == annual.iloc[i]['sieend']), 'JulianDay'].values[0] for i in range(len(annual))]
annual['siegrowthperiod'] = annual['sieendJD'] - annual['siestartJD']
annual['siegrowthrateSD'] = [np.std( sie.loc[(sie['SYear'] == annual.iloc[i]['SYear']) & (sie['JulianDay'] >= annual.iloc[i]['sieminJD']) & (sie['JulianDay'] <= annual.iloc[i]['siemaxJD']),'cdr_delsie_5day'] ) for i in range(len(annual))]
annual['siegrowthrateAvg'] = [np.mean( sie.loc[(sie['SYear'] == annual.iloc[i]['SYear']) & (sie['JulianDay'] >= annual.iloc[i]['sieminJD']) & (sie['JulianDay'] <= annual.iloc[i]['siemaxJD']),'cdr_delsie_5day'] ) for i in range(len(annual))]
annual['siegrowthtotal'] =  annual['sieend'] - annual['siestart']

annual['siestartyear'] =  [ sie.loc[sie['JulianDay'] == annual['siestartJD'][i],'year'].values[0] for i in range(len(annual))]
annual['siestartmonth'] =  [ sie.loc[sie['JulianDay'] == annual['siestartJD'][i],'month'].values[0] for i in range(len(annual))]
annual['siestartday'] =  [ sie.loc[sie['JulianDay'] == annual['siestartJD'][i],'day'].values[0] for i in range(len(annual))]

annual['sieendyear'] =  [ sie.loc[sie['JulianDay'] == annual['sieendJD'][i],'year'].values[0] for i in range(len(annual))]
annual['sieendmonth'] =  [ sie.loc[sie['JulianDay'] == annual['sieendJD'][i],'month'].values[0] for i in range(len(annual))]
annual['sieendday'] =  [ sie.loc[sie['JulianDay'] == annual['sieendJD'][i],'day'].values[0] for i in range(len(annual))]

annual['siedmin'] = sie.loc[(sie['month'] == mmin)].groupby(by='SYear').first()['sie_cdr'].values
annual['siedmax'] = sie.loc[(sie['month'] == mmax)].groupby(by='SYear').first()['sie_cdr'].values

annual['jitter'] = 1+np.random.randn(len(annual))/10

### ID First Date at siemin million  ###
# sie = sie.loc[(sie['month'] >= 9) | (sie['month'] < 3) ]

# Percentiles
sie_mean = sie.loc[:,['SYear','DOY','sie_cdr','cdr_delsie_5day']].groupby(by=['DOY']).mean()
sie_p50 = sie.loc[:,['SYear','DOY','sie_cdr','cdr_delsie_5day']].groupby(by=['DOY']).quantile(0.5)
sie_p05 = sie.loc[:,['SYear','DOY','sie_cdr','cdr_delsie_5day']].groupby(by=['DOY']).quantile(0.05)
sie_p95 = sie.loc[:,['SYear','DOY','sie_cdr','cdr_delsie_5day']].groupby(by=['DOY']).quantile(0.95)
sie_p25 = sie.loc[:,['SYear','DOY','sie_cdr','cdr_delsie_5day']].groupby(by=['DOY']).quantile(0.25)
sie_p75 = sie.loc[:,['SYear','DOY','sie_cdr','cdr_delsie_5day']].groupby(by=['DOY']).quantile(0.75)

sie_mean = sie_mean.reset_index()
sie_p50 = sie_p50.reset_index()
sie_p05 = sie_p05.reset_index()
sie_p95 = sie_p95.reset_index()
sie_p75 = sie_p75.reset_index()
sie_p25 = sie_p25.reset_index()

sip0 = pd.read_csv(eventpath)
sip0['SYear'] = np.where(sip0['month'] >= 9, sip0['year'], sip0['year']-1)
sip0 = sip0.loc[sip0['SYear'] > 1978]
sip = sip0.loc[(sip0['Use'] == 1) & (sip0['DOY'] < dmax2)]
sipyears = np.unique(sip['SYear'])

'''*******************************************
Plotting
*******************************************'''
fig = plt.figure(figsize=(7,3))
gs = GridSpec(4,4, figure=fig)

# Plot the SIE on the start and stop points, noting years of interest #
ax7 = fig.add_subplot(gs[:2,2:])
ax7var = 'siedmin' # 'siegrowthperiod' # 'siegrowthrateSD' # 
changeyears = [1999,2007]
ax7.plot(annual['SYear'],annual[ax7var], color='0.6')
ax7.hlines(y=annual.loc[(annual['SYear'] < changeyears[0]),ax7var].mean(), 
           xmin=1979, xmax=changeyears[0]-1, linestyle='dashed', linewidth=1, color='k')
ax7.hlines(y=annual.loc[(annual['SYear'] >= changeyears[1]),ax7var].mean(), 
           xmin=changeyears[1], xmax=2023, linestyle='dashed', linewidth=1, color='k')

# for y in range(len(sipyears)):
#     ax7.scatter(annual.loc[annual['SYear'] == sipyears[y],'SYear'], annual.loc[annual['SYear'] == sipyears[y],ax7var],
#                 color=colors[y], s=15, zorder=4, label=sipyears[y])
ax7.tick_params(axis='both', which='major', labelsize=7)
ax7.set_ylabel('SIE (million km$^2$)', size=8, weight='bold')
ax7.annotate('b. Oct 1 SIE', xy=(0.01,0.02), xycoords='axes fraction', ha='left', va='bottom', size=8, weight='bold')
ax7.set_ylim(ymin=4, ymax=8.7)

ax8 = fig.add_subplot(gs[2:,2:])
changeyears = [1995,2004]
ax8var = 'siedmax' # 'siegrowthperiod' # 'siegrowthrateSD' # 
ax8.plot(annual['SYear'],annual[ax8var], color='0.6')
ax8.hlines(y=annual.loc[(annual['SYear'] < changeyears[0]),ax8var].mean(), 
           xmin=1979, xmax=changeyears[0]-1, linestyle='dashed', linewidth=1, color='k')
ax8.hlines(y=annual.loc[(annual['SYear'] >= changeyears[1]),ax8var].mean(), 
           xmin=changeyears[1], xmax=2023, linestyle='dashed', linewidth=1, color='k')

# for y in range(len(sipyears)):
#     ax8.scatter(annual.loc[annual['SYear'] == sipyears[y],'SYear'], annual.loc[annual['SYear'] == sipyears[y],ax8var],
#                 color=colors[y], s=15, zorder=4, label=sipyears[y])
ax8.tick_params(axis='both', which='major', labelsize=7)
ax8.set_ylabel('SIE (million km$^2$)', size=8, weight='bold')
ax8.annotate('d. Feb 1 SIE', xy=(0.01,0.02), xycoords='axes fraction', ha='left', va='bottom', size=8, weight='bold')
ax8.set_ylim(ymin=12.5,ymax=17.2)

ax5 = fig.add_subplot(gs[:2,:2])
changeyear = 2005
ax5var = 'siegrowthrateSD'
ax5.plot(annual['SYear'],annual[ax5var], color='0.6')

ax5.hlines(y=annual.loc[(annual['SYear'] < changeyear),ax5var].mean(), 
           xmin=1979, xmax=changeyear-1, linestyle='dashed', linewidth=1, color='k')
ax5.hlines(y=annual.loc[(annual['SYear'] >= changeyear),ax5var].mean(), 
           xmin=changeyear, xmax=2023, linestyle='dashed', linewidth=1, color='k')

# for y in range(len(sipyears)):
#     ax5.scatter(annual.loc[annual['SYear'] == sipyears[y],'SYear'], annual.loc[annual['SYear'] == sipyears[y],ax5var],
#                 color=colors[y], s=15, zorder=4, label=sipyears[y])
ax5.tick_params(axis='both', which='major', labelsize=7)
ax5.set_ylabel('Std. Dev. of ∆SIE\n(million km$^2$ [6 days]$^{-1}$)', weight='bold', size=8)
ax5.annotate('a. Variability of SIE Growth Rate', xy=(0.01,0.98), xycoords='axes fraction', ha='left', va='top', size=8, weight='bold')

# Plot Growth Period Length #

ax6 = fig.add_subplot(gs[2:,:2])
changeyear = 2005
ax6var = 'siegrowthtotal' # 'siegrowthperiod' # 'siegrowthrateSD' # 
# ax6.boxplot(annual[ax6var], vert=False, medianprops={'color':'k'}, zorder=5, whis=[0,100])
ax6.plot(annual['SYear'],annual[ax6var], color='0.6')
# ax6.bar(annual['SYear'],annual['siegrowthperiod'], color='0.6')
ax6.hlines(y=annual.loc[(annual['SYear'] < changeyear),ax6var].mean(), 
           xmin=1979, xmax=changeyear-1, linestyle='dashed', linewidth=1, color='k')
ax6.hlines(y=annual.loc[(annual['SYear'] >= changeyear),ax6var].mean(), 
           xmin=changeyear, xmax=2023, linestyle='dashed', linewidth=1, color='k')


# for y in range(len(sipyears)):
#     ax6.scatter(annual.loc[annual['SYear'] == sipyears[y],'SYear'], annual.loc[annual['SYear'] == sipyears[y],ax6var],
#                 color=colors[y], s=15, zorder=4, label=sipyears[y])
    # ax6.scatter(x=annual.loc[annual['SYear'] == sipyears[y],ax6var],y=annual.loc[annual['SYear'] == sipyears[y],'jitter'],
    #             color=colors[y], zorder=4, label=sipyears[y])
ax6.tick_params(axis='both', which='major', labelsize=7)
ax6.set_ylabel('Total SIE Growth\n(million km$^2$)', size=8, weight='bold')
ax6.annotate('c. Total Growth', xy=(0.01,0.97), xycoords='axes fraction', ha='left', va='top', size=8, weight='bold')
# ax6.set_xticks([])

plt.tight_layout()

plt.savefig('/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses/Figures/Fig7_SeaIceTrends_bytime_V4C.png', dpi=300)


# ####### STATS #####
from scipy import stats
def rodionov_regimes(data, l, p):
    """
    Original Author: Trevor J Amestoy
    Date Created: July 2023
    Original Source: https://github.com/TrevorJA/Rodionov_regime_shifts/blob/main/rodionov.py
    Date Modified by Alex Crawford: 5 Apr 2024
    
    Implements Rodionov's (2004) regime shift detection algorithm:

    Rodionov, S. N. (2004). A sequential algorithm for testing climate regime shifts. 
    Geophysical Research Letters, 31(9).

    Args:
        data (array): Timeseries array of std values
        l (int): The assumed minimum regime length
        p (float): The singificance probability to use when assessing shifts

    Returns:
        list, list: Two lists: The regime-shift indices, the RSI values 
    """
    # Step 1: Set the cut-off length l of the regimes
    # l: Number of years for each regime
    # p: Probability level for significance
    n = len(data)
    regime_shift_indices = []
    rsi = np.zeros(n)

    # Step 2: Determine the difference diff for statistically significant mean values
    t_stat = np.abs(stats.t.ppf(p, (2*l-2)))
    avg_var = np.mean([np.var(data[i:(i+l)]) for i in range(n-l)])
    diff = t_stat * np.sqrt(2 * avg_var / l)

    # Step 3: Calculate initial mean for regime R1
    r1 = np.mean(data[:l])
    r1_lower = r1 - diff
    r1_upper = r1 + diff

    i = l + 1
    while i < n:
        
        # Step 4: Check if the value exceeds the range of R1 ± diff
        if data[i] < r1_lower or data[i] > r1_upper:
            j = i

            # Estimate the mean of regime 2 as the upper bound of r1 distribution
            test_r2 = r1_lower if data[i] < r1_lower else r1_upper
            
            # Step 5: Calculate regime shift index (RSI) across next l-window
            for k in range(j + 1, min(j + l,n)):
                if data[j] > r1:
                    rsi[j] += (data[k] - test_r2)/(avg_var*l)
                elif data[j] < r1:
                    rsi[j] += (test_r2 - data[k])/(avg_var*l)

                # Step 6: Test for a regime shift at year j
                if rsi[j] < 0:
                    rsi[j] = 0
                    break

            # Step 7: Confirm significant regime shift at year j
            if rsi[j] > 0:
                regime_shift_indices.append(j)
                r2 = np.mean(data[j:min(j+l,n)])
                r1 = r2
                r1_lower = r1 - diff
                r1_upper = r1 + diff

            i = j + 1
        else:
            # Recalculate average for R1 to include the new value
            r1 = ((l - 1) * r1 + data[i]) / l

        i += 1
    return regime_shift_indices, rsi


# rodionov_regimes(annual[ax7var], 10, 0.01)
# annual['SYear'][20]
# annual['SYear'][28]

# rodionov_regimes(annual[ax8var], 10, 0.01)
# annual['SYear'][16]
# annual['SYear'][25]

rodionov_regimes(annual[ax6var], 10, 0.01)
annual['SYear'][26]

# rodionov_regimes(annual[ax5var], 10, 0.01)
# annual['SYear'][26]

# annual['siegrowthamly'] = 0
# annual.loc[annual['SYear'] < 2005,'siegrowthamly'] = annual['siegrowthtotal'] - np.mean(annual.loc[annual['SYear'] < 2005,'siegrowthtotal'])
# annual.loc[annual['SYear'] > 2004,'siegrowthamly'] = annual['siegrowthtotal'] - np.mean(annual.loc[annual['SYear'] > 2004,'siegrowthtotal'])

# stats.mannwhitneyu(annual.loc[np.in1d(annual['SYear'],sip0['SYear']) == 0,'siegrowthamly'], annual.loc[np.in1d(annual['SYear'],sip0['SYear']),'siegrowthamly'])
# stats.ttest_ind(annual.loc[np.in1d(annual['SYear'],sip0['SYear']) == 0,'siegrowthamly'], annual.loc[np.in1d(annual['SYear'],sip0['SYear']),'siegrowthamly'])

