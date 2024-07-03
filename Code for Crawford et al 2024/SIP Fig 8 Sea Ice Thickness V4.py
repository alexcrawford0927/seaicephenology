#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 01 May 2024
Modified: 03 Jul 2024

Author: Alex Crawford

Purpose: Sea ice thickness time series
"""

import pandas as pd
import xarray as xr
import xesmf as xe
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import CycloneModule_13_2 as md
from matplotlib.colors import LinearSegmentedColormap


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

def pettitt_test(data, times=np.nan, nmin=5):
    '''
    Performs a Pettitt Test of homogeneity (for change-point detection in a 
        time series). Returns the index and averages for the two periods 
        of interest. Also returns the change year if years are provided. Note 
        that the index/year represent the first instance of the second period.
        Also note that it is possible for there to be multiple change points,
        in which case multiple means, indices, and years are returned. However,
        because this test is designed to identify a the *maximum* change point,
        this only happens if to potential change points hae equal intensity. 
        In such a case, the p value will be identical, hence only returning a
        single p-value.

    Parameters
    ----------
    data : numpy array or list
        1-D array of values of interest.
    times : numpy array or list, optional
        1-D array of times -- must be the same shape as data.
    nmin : integer, optional
        minimum size of periods. The default is 5

    Returns
    -------
    A tuple of: (period 1 mean, period 2 mean, p-value, index of change point, time of change point)
        If no time array provided, the index is repeated for the last member of the tuple
    '''
    data = np.array(data) # Ensure data is a numpy array
    
    taulist = []
    for tau in range(nmin,len(data)-nmin): # Loop through each possible change point
        
        isum = 0
        for i in range(tau): # Loop through each value in the first period
        
            # Add the sign of the difference
            isum += np.sign( data[tau:] - data[i] ).sum()
                    
        # Append to full list of potential break points
        taulist.append( isum )
        
    tauarr = np.abs(taulist)
    K = np.max(tauarr)
    Ktau = np.where(tauarr == K)[0]
    
    p = 2*np.exp((-6*K**2)/(len(data)**2 + len(data)**3))
    
    index = Ktau + nmin
    
    try:
        return [np.mean(data[:i]) for i in index], [np.mean(data[i:]) for i in index], p, index, np.array(times)[index]
    except:
        return [np.mean(data[:i]) for i in index], [np.mean(data[i:]) for i in index], p, index, index

def pettitt2_test(data, times=np.nan, nmin=5):
        '''
        Performs a Pettitt Test of homogeneity (for change-point detection in a 
            time series). Returns the index and averages for the two periods 
            of interest. Also returns the change year if years are provided. Note 
            that the index/year represent the first instance of the second period.
            Also note that it is possible for there to be multiple change points,
            in which case multiple means, indices, and years are returned. However,
            because this test is designed to identify a the *maximum* change point,
            this only happens if to potential change points hae equal intensity. 
            In such a case, the p value will be identical, hence only returning a
            single p-value.

        Parameters
        ----------
        data : numpy array or list
            1-D array of values of interest.
        times : numpy array or list, optional
            1-D array of times -- must be the same shape as data.
        nmin : integer, optional
            minimum size of periods. The default is 5

        Returns
        -------
        A tuple of: (period 1 mean, period 2 mean, p-value, index of change point, time of change point)
            If no time array provided, the index is repeated for the last member of the tuple
        '''
        data = np.array(data) # Ensure data is a numpy array
        
        taulistsum = []
        taulistmean = []
        for tau in range(nmin,len(data)-nmin): # Loop through each possible change point
            
            isum, imean = 0, 0
            for i in range(tau): # Loop through each value in the first period
            
                # Add the sign of the difference
                isum += np.sign( data[tau:] - data[i] ).sum()
                imean += np.sign( data[tau:] - data[i] ).mean()
                        
            # Append to full list of potential break points
            taulistmean.append( imean / (tau) )
            taulistsum.append( isum )
            
        tauarrsum, tauarrmean = np.abs(taulistsum), np.abs(taulistmean)
        Ksum, Kmean = np.max(tauarrsum), np.max(tauarrmean)
        Ktau = np.where(tauarrmean == Kmean)[0]
        
        p = 2*np.exp((-6*Ksum**2)/(len(data)**2 + len(data)**3))
        
        index = Ktau + nmin
                
        try:
            return [np.mean(data[:i]) for i in index], [np.mean(data[i:]) for i in index], p, index, np.array(times)[index]
        except:
            return [np.mean(data[:i]) for i in index], [np.mean(data[i:]) for i in index], p, index, index
        
'''*******************************************
Declare Variables
*******************************************'''

siemin, siemax = 7.75, 13.0 # 8, 13 #7.8, 13.4 # 7.75, 13.25
ncvar = 'nsidc_nt_seaice_conc' # 'cdr_seaice_conc'
sicthresh = 0.15
ncmin, ncmax = 0, 1
sitmin = 0
colors=['red','darkorange','gold','green','blue','purple','hotpink']

doys = [1,32,60,91,121,152,182,213,244,274,305,335,366,397,425]
doyl = ["Jan 1","Feb 1","Mar 1","Apr 1","May 1","Jun 1","Jul 1","Aug 1","Sep 1","Oct 1","Nov 1","Dec 1","Jan 1","Feb 1","Mar 1"]

path = '/Volumes/Miranda/SeaIce'
siepath = '/Volumes/Miranda/SeaIce/SeaIceIndex_NH_19781026-20240325_v3.0.csv'
sicpath = '/Volumes/Miranda/SeaIce/G02202_V4/Daily'
sicpath2 = '/Volumes/Miranda/SeaIce/G10016_V2/Daily'
siapath = '/Volumes/Miranda/SeaIce/Age'
sitpaths = [path+'/ISSITGR4-IceSat-Seasonal']
frampath = '/Volumes/Miranda/SeaIce/Ice_thickness_distribution_Fram_Strait/FramStrait_SIT_Summary.csv'
soriotpath = '/Volumes/Miranda/SeaIce/Soirot/Soriot_PMW_SIT_Avg_Summary.csv'

sitprjs = ['PSN']
sats = ['ICESat']
satvers = [['']]
satvs = [['']]
sitvars = [['ice_thickness_int']]
sitlatvars = ['latitude']
sitlonvars = ['latitude']
sitpreyear = ['ISSITGR4_01_ON']
sityearlen = [2,4]
sitmons = [['_ON','_ON']]
mons = [10,11]
bbox = [80,85,-30,100] # minlat, maxlat, minlon, maxlon

ease2path = '/Volumes/Miranda/Projections/Cryosat2-UBristol_Projection.nc'
psnpath = '/Volumes/Miranda/Projections/psn_projection.nc'
siaprjpath = siapath+"/iceage_nh_12.5km_20110101_20111231_v4.1.nc"

outpath = '/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses/Figures/Fig_SeaIceThickness_Oct-Dec_V4.png'

prj = ccrs.NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)
cmap = LinearSegmentedColormap.from_list('seaice',colors=[(0/256,0/256,0/256,1),(32/256,59/256,93/256,1),(87/256,123/256,192/256),(163/256,176/256,198/256,1),(256/256,256/256,256/256,1)],N=256)

'''*******************************************
Processing
*******************************************'''

### Calculate the % of years that have sea ice present for each grid cell ###
siestartarrs, sieendarrs = [], []
for y in range(1979,2023):
    
    try:
        xrminfile = [ f for f in md.listdir(sicpath+'/'+str(y)) if '_'+str(y)+md.dd[10-1] in f][0]
        xrmin = xr.open_dataset( sicpath+'/'+str(y)+'/'+xrminfile)
    except:
        xrminfile = [ f for f in md.listdir(sicpath2+'/'+str(y)) if '_'+str(y)+md.dd[10-1] in f][0]
        xrmin = xr.open_dataset( sicpath2+'/'+str(y)+'/'+xrminfile)

    xrminarr = xrmin[ncvar][0,:,:].data
    xrminbool = (xrminarr >= sicthresh).astype(float)
    xrminbool[(xrminarr < ncmin) | (xrminarr > ncmax)] = np.nan
    
    siestartarrs.append( xrminbool )
    
    try:
        xrmaxfile = [ f for f in md.listdir(sicpath+'/'+str(y+1)) if '_'+str(y+1)+md.dd[1-1] in f][0]
        xrmax = xr.open_dataset( sicpath+'/'+str(y+1)+'/'+xrmaxfile)
    except:
        xrmaxfile = [ f for f in md.listdir(sicpath2+'/'+str(y+1)) if '_'+str(y+1)+md.dd[1-1] in f][0]
        xrmax = xr.open_dataset( sicpath2+'/'+str(y+1)+'/'+xrmaxfile)
        
    xrmaxarr = xrmax[ncvar][0,:,:].data    
    xrmaxbool = (xrmaxarr >= sicthresh).astype(float)
    xrmaxbool[(xrmaxarr < ncmin) | (xrmaxarr > ncmax)] = np.nan
    sieendarrs.append( xrmaxbool )  

probsie_start = np.mean(siestartarrs, axis=0)*100
probsie_end = np.mean(sieendarrs, axis=0)*100

# Create Mask for PSN Grid #
maskpsn = np.where(probsie_start >= 0.5, 1, np.nan)


# Create Mask for EASE2 Grid #
psn = xr.open_dataset(psnpath)
ease2 = xr.open_dataset(ease2path)
psn_to_ease2 = xe.Regridder(psn,ease2,method='bilinear')
maskease2 = np.where(psn_to_ease2(probsie_start) >= 0.5, 1, np.nan)

# Write out Masks
psn = psn.assign(mask=(['y','x'],maskpsn))
ease2 = ease2.assign(mask=(['y','x'],maskease2))

psnout = psn[['lat','lon','mask','area']]
ease2out = ease2[['lat','lon','mask','area']]

psnout.to_netcdf('/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses/masks/psn_25km_projection.nc')
ease2out.to_netcdf('/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses/masks/ease2_25km_projection.nc')

# Define a mask for each Thickness Dataset Depending on Projection #
sitmasks = [maskpsn if p == 'PSN' else maskease2 for p in sitprjs]
areas = [psn['area'][:].data if p == 'PSN' else ease2['area'][:].data  for p in sitprjs]


ylist, mlist, satlist, vlist, verlist, avglist1, p75list1, avglist2, p75list2 = [[] for i in range(9)]
for s in range(len(sats)):
    print(sats[s])
    
    files = md.listdir(sitpaths[s])
    
    for m in range(len(mons)):
        
        mfiles = [f for f in files if sitmons[s][m] in f]
        
        for file in mfiles:
            xrd = xr.open_dataset(sitpaths[s]+"/"+file)
            
            # Define Bounding Box Mask
            lats, lons = xrd[sitlatvars[s]], xrd[sitlonvars[s]]
            lons = np.where(lons > 180, lons-360, lons)
            maskbb = np.where((lats >= bbox[0]) & (lats <= bbox[1]) & (lons >= bbox[2]) & (lons <= bbox[3]), 1, np.nan)
            
            for v in range(len(satvers[s])):
                arr = xrd[sitvars[s][v]].data
                arr[arr < sitmin] = np.nan
                
                # First Subset -- by bounding box
                avgsit1 = np.nansum( arr * maskbb * areas[s] ) / np.sum( np.isfinite(arr * maskbb) * areas[s])
                p75sit1 = np.nanpercentile( arr * maskbb, 75)
                
                # Second Subset -- by mask
                avgsit2 = np.nansum( arr * sitmasks[s] * areas[s] ) / np.sum( np.isfinite(arr * sitmasks[s]) * areas[s])
                p75sit2 = np.nanpercentile( arr * sitmasks[s], 75)
                
                # Append to lists
                mlist.append( mons[m] ), ylist.append(int(file.split(sitpreyear[s])[1][0:sityearlen[s]]) )
                satlist.append(sats[s]), verlist.append(satvers[s][v]), vlist.append(satvs[s][v])
                avglist1.append(avgsit1), avglist2.append(avgsit2)
                p75list1.append(p75sit1), p75list2.append(p75sit2)
            
sitdf = pd.DataFrame({'Year':ylist,'Month':mlist,'Satellite':satlist,
                      'Version':vlist,'LongName':verlist,'AvgSIT-BBox':avglist1,
                      'perc75SIT-BBox':p75list1,'AvgSIT-SIEMask':avglist2,
                      'perc75SIT-SIEMask':p75list2})
sitdf['Year'] = np.where(sitdf['Year'] < 1000, sitdf['Year']+2000, sitdf['Year'])
                
sitdf2 = sitdf.loc[:,['Year','Satellite','Version','LongName','AvgSIT-BBox','perc75SIT-BBox', 'AvgSIT-SIEMask', 'perc75SIT-SIEMask']].groupby(by=['Year','Satellite','Version']).mean()
sitdf2 = sitdf2.reset_index()

# Open Up Fram Data #
framdf = pd.read_csv(frampath)
framdf = framdf.groupby(by='time_start').mean()
framdf = framdf.reset_index()
framdf['year'] = [int(date[0:4]) for date in framdf['time_start']]
framdf['month'] = [int(date[5:7]) for date in framdf['time_start']]
framdf['syear'] = np.where(framdf['month'] > 8, framdf['year'], framdf['year']-1)
framdf = framdf[np.in1d(framdf['month'],[10,11,12])]
framdf = framdf.groupby(by=['syear']).mean()
framdf = framdf.reset_index()

# Open Soriot Data #
soriot = pd.read_csv(soriotpath)
soriot['SYear'] = np.where(soriot['Month'] > 8, soriot['Year'], soriot['Year']-1)
soriot1 = soriot.loc[np.in1d(soriot['Month'],[10,11,12])].groupby(by=['SYear']).mean()
soriot1 = soriot1.reset_index()
soriot2 = soriot.loc[np.in1d(soriot['Month'],[10,11])].groupby(by=['SYear']).mean()
soriot2 = soriot2.reset_index()

# Open Sea Ice Age Data #
xsia = xr.open_dataset(siaprjpath)
siatopsn = xe.Regridder(xsia, psn, 'nearest_s2d')

gt1, gt2, gt3, gt4 = [], [], [], []
for file in md.listdir(siapath,end='.nc'):
    
    # Open File
    xsia = siatopsn( xr.open_dataset(siapath+"/"+file) )
    
    # Use the Week of Oct 7 -- about Week 40
    age = xsia['age_of_sea_ice'][40,:,:].data
    
    # Subset to mask & Take percentage over a given age
    gt1.append( np.nanmean( (age > 1) * (age < 20) * psn['mask'] ) )
    gt2.append( np.nanmean( (age > 2) * (age < 20) * psn['mask'] ) )
    gt3.append( np.nanmean( (age > 3) * (age < 20) * psn['mask'] ) )
    gt4.append( np.nanmean( (age > 4) * (age < 20) * psn['mask'] ) )

agedf = pd.DataFrame({'Year':np.arange(1984,2023), 'gt1':gt1, 'gt2':gt2, 'gt3':gt3, 'gt4':gt4})
    

'''####################
Plotting 
####################'''

fig = plt.figure(figsize=(6,6))

ax1 = fig.add_subplot(2,2,1)

# Plot Soriot Data first #
ax1.vlines(x=2004.5, ymin=0, ymax=4.5, color='0.8', linestyle='dashed', linewidth=2, zorder=1)
ax1.vlines(x=2004.5, ymin=0, ymax=4.5, color='k', alpha=0.5, linestyle='dashed', linewidth=1.2, zorder=1)
# ax1.vlines(x=2003.5, ymin=0, ymax=4.5, color='red', alpha=0.5, linestyle='dashed', linewidth=1.2, zorder=1)

ax1.plot(soriot1['Year'],soriot1['AvgSIT-SIEMask'], color='k', label='PMW (Oct-Dec)')
# ax1.plot(soriot2['Year'],soriot2['AvgSIT-SIEMask'], color='red', label='PMW (Oct-Nov)')

df = sitdf2.loc[(sitdf2['Satellite'] == 'ICESat')]
ax1.plot(df['Year'],df['AvgSIT-SIEMask'], color='blue', label='IceSat (Oct-Nov)')

ax1.set_xlim(xmin=1989,xmax=2024)
ax1.set_ylim(ymin=0.3,ymax=4.5)
ax1.tick_params(axis='both', which='major', labelsize=7)
ax1.set_title('a. Sea Ice Thickness\nfrom Satellites',size=8, weight='bold')
ax1.set_ylabel('Sea Ice Thickness (m)', size=8)

# Plot Fram Strait Data
ax2 = fig.add_subplot(2,2,3)
ax2.vlines(x=2004.5, ymin=0, ymax=4.5, color='0.8', linestyle='dashed', linewidth=2, zorder=1)
ax2.vlines(x=2007.5, ymin=0, ymax=4.5, color='k', alpha=0.5, linestyle='dashed', linewidth=1.2, zorder=1)
ax2.vlines(x=2009.5, ymin=0, ymax=4.5, color='red', alpha=0.5, linestyle='dashed', linewidth=1.2, zorder=1)

ax2.plot(framdf['year'],framdf['modefitted_sit'], color='k', label='Second Mode')
ax2.plot(framdf['year'],framdf['p75_sit'], color='red', linestyle='solid', label='75th Percentile')
# ax2.plot(framdf['year'],framdf['mode2_sit'], color='gray')
# ax2.plot(framdf['year'],framdf['avg_sit'], color='blue', label='Average')
# ax2.plot(framdf['year'],framdf['median_sit'], color='green', label='Median')

ax2.set_xlim(xmin=1989,xmax=2024)
ax2.set_ylim(ymin=0.3,ymax=4.5)
ax2.tick_params(axis='both', which='major', labelsize=7)
ax2.set_title('b. Oct-Dec Sea Ice Thickness\nfrom Fram Strait Moorings', size=8, weight='bold')
ax2.set_ylabel('Sea Ice Thickness (m)', size=8)

### Sea Ice Age ###
ax4 = fig.add_subplot(2,2,2)
ax4.vlines(x=2004.5, ymin=0, ymax=4.5, color='0.8', linestyle='dashed', zorder=1)
ax4.vlines(x=2001.5, ymin=0, ymax=4.5, color='red', alpha=0.5, linestyle='dashed', linewidth=1.2, zorder=1)
ax4.vlines(x=2004.5, ymin=0, ymax=4.5, color='blue', alpha=0.5, linestyle='dashed', linewidth=1.2, zorder=1)
ax4.vlines(x=2004.6, ymin=0, ymax=4.5, color='magenta', alpha=0.5, linestyle='dashed', linewidth=1.2, zorder=1)
ax4.vlines(x=2005.5, ymin=0, ymax=4.5, color='green', alpha=0.5, linestyle='dashed', linewidth=1.2, zorder=1)


ax4.plot(agedf['Year'], agedf['gt1'], color='red', label='Age > 1 yr')
ax4.plot(agedf['Year'], agedf['gt2'], color='blue', label='Age > 2 yrs')
ax4.plot(agedf['Year'], agedf['gt3'], color='magenta', label='Age > 3 yrs')
ax4.plot(agedf['Year'], agedf['gt4'], color='green', label='Age > 4 yrs')


ax4.set_xticks(np.arange(1980,2030,5))
ax4.set_xlim(xmin=1983,xmax=2024)
ax4.set_ylim(ymin=0,ymax=1)
ax4.tick_params(axis='both', which='major', labelsize=7)
ax4.set_title('c. Week 40 Sea Ice Age Fraction', size=8, weight='bold')
ax4.set_ylabel('Fraction Exceeding Given Age', size=8)

plt.tight_layout()

### Location Map ###

# Typical Sea Ice Extent at siemin million km2 #
ax3 = fig.add_axes([0.58,0.04,0.25,0.43], projection=prj)
ax3.set_facecolor('0.4')
ax3.add_feature(cfeature.LAND, facecolor='0.4',zorder=7)
cf3 = ax3.contourf(xrmax['xgrid'], xrmax['ygrid'], probsie_end, levels=np.arange(0,110,10), cmap=cmap, transform=prj)
ax3.contour(xrmin['xgrid'],xrmin['ygrid'], probsie_start, [50], colors='b', linewidths=1 , transform=prj, zorder=8)
ax3.annotate('d. Sea Ice Probablity', xy=(0.02,0.98), xycoords='axes fraction', 
              weight='bold', size=7, ha='left', va='top', color='w', zorder=12)
ax3.plot((1,2),(0,0), linewidth=1, color='b', zorder=1, label=' On Oct 1 ')
ax3.plot( [-4.5], [78.8],  transform=ccrs.Geodetic(), color='red', markersize=15, marker='.', linewidth=0, zorder=10, label = 'Fram Strait Moorings')


### Legends ###
ax1.legend(fontsize=7, loc='lower right',bbox_to_anchor=(1,0), edgecolor='w', framealpha=1)
ax2.legend(fontsize=7, loc='upper left', bbox_to_anchor=(0.6,1), framealpha=1)
ax3.legend(fontsize=7, loc='upper left', bbox_to_anchor=(0,-0.02), edgecolor='w')
ax4.legend(fontsize=7, loc='upper right',bbox_to_anchor=(1,1), edgecolor='w', framealpha=1)


# Color Bar
cbar_ax = fig.add_axes([0.87, 0.08, 0.02, 0.37])
cbar1 = fig.colorbar(cf3,cbar_ax,orientation='vertical')
cbar1.ax.tick_params(labelsize=7)
cbar1.set_label('Probability of Sea Ice Presence (%) on Jan 1', size=7)

plt.savefig(outpath,dpi=300)


### RODIONOV ###
rodavgpmw, rodp75pmw, rodavgpmw2, rodmodfram, rodp75fram = [], [], [], [], []
ages = [[] for i in range(4)]
for l in range(5,15):
    
    breaks, ris = rodionov_regimes(soriot1['AvgSIT-SIEMask'], l, 0.01) # Identify breaks
    i = np.where(np.max(ris) == ris)[0][0] # Find largest break
    rodavgpmw.append( soriot1['SYear'].values[i] ) # Append

    breaks, ris = rodionov_regimes(soriot1['perc75SIT-SIEMask'], l, 0.01) # Identify breaks
    i = np.where(np.max(ris) == ris)[0][0] # Find largest break
    rodp75pmw.append( soriot1['SYear'].values[i] ) # Append

    breaks, ris = rodionov_regimes(soriot2['AvgSIT-SIEMask'], l, 0.01) # Identify breaks
    i = np.where(np.max(ris) == ris)[0][0] # Find largest break
    rodavgpmw2.append( soriot2['SYear'].values[i] ) # Append

    breaks, ris = rodionov_regimes(framdf['modefitted_sit'], l, 0.01) # Identify breaks
    i = np.where(np.max(ris) == ris)[0][0] # Find largest break
    rodmodfram.append( framdf['syear'].values[i] ) # Append    
        
    breaks, ris = rodionov_regimes(framdf['p75_sit'], l, 0.01) # Identify breaks
    i = np.where(np.max(ris) == ris)[0][0] # Find largest break
    rodp75fram.append( framdf['syear'].values[i] ) # Append   

    for a in range(1,5):
        breaks, ris =  rodionov_regimes(agedf['gt'+str(a)], l, 0.01) # Identify breaks
        i = np.where(np.max(ris) == ris)[0][0] # Find largest break
        ages[a-1].append( agedf['Year'].values[i] ) # Append   

# ### PETTITT TEST ###
# l = 0
# var, mean1, mean2, p, year = [[] for i in range(5)]
# for a in range(1,5):
#     pet = pettitt_test(agedf['gt'+str(a)], agedf['Year'].values, nmin=l)

#     var.append( 'Age_gt'+str(a) ) # Append  
#     mean1.append(pet[0]), mean2.append(pet[1])
#     p.append(pet[2]), year.append(pet[-1])


# pettitt_test(soriot1['AvgSIT-SIEMask'], soriot1['SYear'], nmin=l)
# pettitt_test(soriot1['perc75SIT-SIEMask'], soriot1['SYear'], nmin=l)
# pettitt_test(framdf['modefitted_sit'], framdf['syear'], nmin=l)
# pettitt_test(framdf['p75_sit'], framdf['syear'].values, nmin=l)
    
# ### PETTITT 2 TEST ###
# l = 5
# var, diff, mean1, mean2, p, year = [[] for i in range(6)]
# for a in range(1,5):
#     pet = pettitt2_test(agedf['gt'+str(a)], agedf['Year'].values, nmin=l)

#     var.append( 'Age_gt'+str(a) ) # Append  
#     mean1.append(pet[0]), mean2.append(pet[1])
#     p.append(pet[2]), year.append(pet[-1])
#     diff.append( np.array(pet[1]) -  np.array(pet[0]))

# [np.array(year[i])[ np.where( np.abs(diff[i]) == np.max(np.abs(diff[i])) ) ] for i in range(len(diff))]


# pettitt2_test(soriot1['AvgSIT-SIEMask'], soriot1['SYear'], nmin=l)
# pettitt2_test(soriot1['perc75SIT-SIEMask'], soriot1['SYear'], nmin=l)
# pettitt2_test(framdf['modefitted_sit'], framdf['syear'], nmin=l)
# pettitt2_test(framdf['p75_sit'], framdf['syear'].values, nmin=l)

        
            
            
    
        
        
        




