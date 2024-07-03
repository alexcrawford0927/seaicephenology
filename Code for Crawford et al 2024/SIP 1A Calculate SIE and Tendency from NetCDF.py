#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 16 Apr 2024
Modified: 03 Jul 2024

Author: Alex Crawford

Purpose: 
    1) Create a single CSV of regional daily sea ice extent for each NSIDC
    region and for each major sea ice concentration algorithm. 
    2) Create a single CSV file that has the 6-day tendency of sea ice extent 
    for each NSIDC region for each major sea ice concentration algorithm.
    3) Calculate % frequency of negative SIE tendency for each region/sector
"""

'''*******************************************
Load Modules
*******************************************'''
print(' ***Load Modules***')

import xarray as xr
import numpy as np
import pandas as pd
import CycloneModule_13_2 as md

'''*******************************************
Declare Variables
*******************************************'''
print(' ***Define Variables***')

path = '/Volumes/Miranda/'
sipath = '/Volumes/Miranda/SeaIce/G02202_V4/Daily'
sinrtpath = '/Volumes/Miranda/SeaIce/G10016_V2/Daily'

btpath = '/Volumes/Miranda/SeaIce/Bootstrap/Daily'
ntpath = '/Volumes/Miranda/SeaIce/NASATeam/Daily'
ntnrtpath = '/Volumes/Miranda/SeaIce/NASATeam/NearRealTime'

prjpath = path+"/Projections/psn_projection.nc"
areapath = path+'/Projections/NSIDC0771_CellArea_PS_N25km_v1.0.nc'

outpath = '/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses/'

tstep = 6
ymin, ymax = 1979, 2023
sicmin = 0.15
doymin, doymax = 276, 369

regvar = 'reg0780'
regs = np.arange(1,20)

varlist = ['cdr_seaice_conc','nsidc_bt_seaice_conc','nsidc_nt_seaice_conc']
vlist = ['sie_cdr','sie_nnt','sie_nbt','sie_bt','sie_nt']

ncvar = 'ICECON'

'''*******************************************
Step 1: Daily Regional SIE
*******************************************'''
print(' ***Start Phase 1***')


# Load projection file
prj = xr.open_dataset(prjpath)
regarr, lons, lats = prj['reg0780'][:].data, prj['lon'][:].data, prj['lat'][:].data
regarr[(lons > -30) & (lons < 20) & (regarr == 1)] = 8
regarr[(lons > 20) & (lons < 65) & (regarr == 1)] = 7
regarr[(lons > 65) & (lons < 90) & (regarr == 1)] = 6

area = xr.open_dataset(areapath)

# masks = [(prj[regvar].values == reg) * prj['area'].values for reg in regs]
masks = [(regarr == reg) * area['cell_area'].values / 1000 / 1000 for reg in regs] 

ylist, mlist, dlist, jdlist, reglist, ntlist, btlist, nntlist, nbtlist, cdrlist = [[] for i in range(10)]
for y in range(ymin,ymax+1):
    print('    ' + str(y))
    
    # Load File Names
    files = md.listdir(sipath+"/"+str(y))
    btfiles = md.listdir(btpath+"/"+str(y))
    ntfiles = md.listdir(ntpath+"/"+str(y)) 
    
    # Loop through CDR files
    for file in files:
        ## CDR File ##
        try:
            xrd = xr.open_dataset(sipath+"/"+str(y)+"/"+file)
        except:
            xrd = xr.open_dataset(sinrtpath+"/"+str(y)+"/"+file)
            
        # Calculate Sea Ice Cover Binary for each cell
        sics = [ np.where((xrd[var] >= sicmin) & (xrd[var] <= 1), 1, 0) for var in varlist]
        
        # Loop through each region & sea ice cocentration dataset
        for r in range(len(regs)):
            # store the regional SIE
            cdrlist.append( np.sum(sics[0] * masks[r]) )
            nbtlist.append( np.sum(sics[1] * masks[r]) )
            nntlist.append( np.sum(sics[2] * masks[r]) )
            reglist.append( prj['reg0780'].attrs['flag_names'][regs[r]] )
            ylist.append(y)
            mlist.append(int(file.split('_nh_')[1][4:6]))
            dlist.append(int(file.split('_nh_')[1][6:8]))
            jdlist.append(md.daysBetweenDates([1900,1,1],[y,mlist[-1],dlist[-1]]))
            # doylist.append(md.daysBetweenDates([y,1,0],[y,mlist[-1],dlist[-1]]))
        
        ## Bootstrap ##
        btfilesd = [btf for btf in btfiles if btf.split('_')[-2][4:8] == file.split('_nh_')[1][4:8]]
        if len(btfilesd) == 0:
            for r in range(len(regs)):
                btlist.append( np.nan )
        
        else:
            xrd = xr.open_dataset(btpath+"/"+str(y)+"/"+btfilesd[0])
            sic = xrd[[key for key in list(xrd.variables) if ncvar in key][0]][0,:].data
            
            # Fill the pole hole
            sic[(sic < 0) | (sic > 1)] = np.nan
            
            # Fill in Pole Hole
            holeR, holeC = np.where((lats >= 83) & (regarr < 20) & (np.isnan(sic)))
            fillR, fillC = np.where((lats >= np.min(lats[holeR,holeC])-1) & (np.isfinite(sic)))
            sic[holeR,holeC] = np.nanmean(sic[fillR,fillC])
            sic = np.where((sic >= sicmin), 1, 0)
        
            for r in range(len(regs)):
                btlist.append( np.sum( sic * masks[r] ) )
            
        ## NASA Team ##
        ntfilesd = [ntf for ntf in ntfiles if ntf.split('_')[-2][4:8] == file.split('_nh_')[1][4:8]]
        if len(ntfilesd) == 0:
            for r in range(len(regs)):
                ntlist.append( np.nan )
        
        else:
            try:
                xrd = xr.open_dataset(ntpath+"/"+str(y)+"/"+ntfilesd[0])
            except:
                xrd = xr.open_dataset(ntnrtpath+"/"+ntfilesd[0])

            sic = xrd[[key for key in list(xrd.variables) if ncvar in key][0]][0,:].data
            
            # Fill the pole hole
            sic[(sic < 0) | (sic > 1)] = np.nan
            
            # Fill in Pole Hole
            holeR, holeC = np.where((lats >= 83) & (regarr < 20) & (np.isnan(sic)))
            fillR, fillC = np.where((lats >= np.min(lats[holeR,holeC])-1) & (np.isfinite(sic)))
            sic[holeR,holeC] = np.nanmean(sic[fillR,fillC])
            sic = np.where((sic >= sicmin), 1, 0)
        
            for r in range(len(regs)):
                ntlist.append( np.sum( sic * masks[r] ) )            

siedf = pd.DataFrame({'year':ylist, 'month':mlist,'day':dlist,'JulianDay':jdlist,
                      'region':reglist, 'sie_cdr':cdrlist,'sie_nbt':nbtlist, 'sie_nnt':nntlist,
                      'sie_bt':btlist, 'sie_nt':ntlist})
siedf.to_csv(outpath+"/SIE_Daily_Regional.csv",index=False)

'''*******************************************
Step 2: Daily Regional SIE Tendency
*******************************************'''
print(' ***Start Phase 2***')
siedf = pd.read_csv(outpath+"/SIE_Daily_Regional.csv")

### Diff for any 5-day period (t=0 days minus t=-6 days) ###
if 'JulianDay' not in siedf.columns:
    jds = [md.daysBetweenDates([1900,1,1],[siedf['year'].values[i],siedf['month'].values[i],siedf['day'].values[i]]) for i in range(len(siedf))]
    siedf['JulianDay'] = jds
    siedf.to_csv(outpath+"/SIE_Daily_Regional.csv",index=False)
   
for var in vlist:
    siedf[var.split('_')[1]+'_delsie_5day'] = np.nan

for r in regs:
    reg = prj['reg0780'].attrs['flag_names'][r]
    print('    ' + reg)

    sierf = siedf.loc[(siedf['region'] == reg)]
    sierf = sierf.reset_index()
    
    for var in vlist:
        delsie5day = []
        for i in range(len(sierf)):
            try:
                delsie5day.append( sierf.loc[i,var] - sierf.loc[sierf['JulianDay'] == sierf.loc[i,'JulianDay']-tstep,var].values[0] )
            except:
                delsie5day.append( np.nan )
    
        sierf[var.split('_')[1]+'_delsie_5day'] = delsie5day
    
        # Reappend to main dataframe
        siedf.loc[siedf['region'] == reg,var.split('_')[1]+'_delsie_5day'] = sierf[var.split('_')[1]+'_delsie_5day'].values

siedf.to_csv(outpath+"/SIE_Daily_Regional.csv", index=False)

'''*******************************************
Step 3: Daily Arctic-Wide SIE & SIE Tendency
*******************************************'''
print(' ***Start Phase 3***')

# Load daily pan-Arctic SIE
sii = pd.read_csv('/Volumes/Miranda/SeaIce/SeaIceIndex_NH_19781026-20240325_v3.0.csv')

# Add a pan-Arctic SIE row to the siedf
siedf['SeaIceIndex'] = np.nan 
for row in range(len(sii)):
    rowdf = sii.loc[row]
    
    siedf.loc[(siedf['year'] == rowdf['Year']) & (siedf['month'] == rowdf['Month']) & (siedf['day'] == rowdf['Day']),'SeaIceIndex'] = rowdf['Extent']

siedf.to_csv(outpath+"/SIE_Daily_Regional.csv", index=False)

####

siedf = pd.read_csv(outpath+"/SIE_Daily_Regional.csv")

# Calculate the Arctic-wide SIE based on the regions
siedfall = siedf.loc[:,['JulianDay', 'sie_cdr', 'sie_bt','sie_nt', 'sie_nbt','sie_nnt','cdr_delsie_5day', 'bt_delsie_5day', 'nt_delsie_5day','nbt_delsie_5day', 'nnt_delsie_5day']].groupby(by=['JulianDay']).sum()/1000000
siedfall = siedfall.reset_index()
for var in ['year', 'month', 'day','SeaIceIndex']:
    siedfall[var] = siedf.loc[siedf['region'] == 'CAO',var].values
    
siedfall['sie_bt'] = np.where(siedfall['sie_bt'] == 0, np.nan, siedfall['sie_bt'])
siedfall['sie_nt'] = np.where(siedfall['sie_nt'] == 0, np.nan, siedfall['sie_nt'])

siedfall.to_csv(outpath+"/SIE_Daily_ArcticWide.csv",index=False)

'''*******************************************
Step 4: Regional Frequency of ∆SIE < 0 km2
*******************************************'''
### Frequency of Negative 5-day Tendency in any region ###
siedf2 = siedf.loc[( (siedf['DOY'] >= doymin) | (siedf['DOY'] <= doymax) )]
siedf2['Neg'] = (siedf['cdr_delsie_5day'] < 0)
siedf2.to_csv(outpath+"/SIE_Daily_Regional2_bytime.csv", index=False)

neg1 = siedf2.loc[:,['month','region','Neg']].groupby(by=['month','region']).mean()
neg1 = neg1.reset_index()

neg2 = siedf2.loc[:,['region','Neg']].groupby(by=['region']).mean()
neg2 = neg2.reset_index()

### Frequency of joint-negative 5-day Tendency in region pairs ###
# (Barents | Kara) --> 36.1%
100* np.mean(siedf2.loc[siedf2['region'] == "Kara",'Neg'].values | siedf2.loc[siedf2['region'] == "Barents",'Neg'].values)

# (Barents | Kara) & (Bering | Chukchi) --> 10.5%
100*np.mean( (siedf2.loc[siedf2['region'] == "Kara",'Neg'].values | siedf2.loc[siedf['region'] == "Barents",'Neg'].values) & (siedf2.loc[siedf2['region'] == "Chukchi",'Neg'].values | siedf2.loc[siedf2['region'] == "Bering",'Neg'].values) )

# (Barents | Kara | Greenland) & (Bering | Chukchi) --> 15.3%
100*np.mean( (siedf2.loc[siedf2['region'] == "Kara",'Neg'].values | siedf2.loc[siedf2['region'] == "Barents",'Neg'].values | siedf2.loc[siedf2['region'] == "E. Greenland",'Neg'].values) & (siedf2.loc[siedf2['region'] == "Chukchi",'Neg'].values | siedf2.loc[siedf2['region'] == "Bering",'Neg'].values) )

# (Greenland) & (Hudson) --> 1.9%
100*np.mean( siedf2.loc[siedf2['region'] == 'E. Greenland','Neg'].values & siedf2.loc[siedf2['region'] == 'Hudson','Neg'].values )

# (Greenland) & (Hudson) & (Gulf of St Law) --> 0.1%
100*np.mean( siedf2.loc[siedf2['region'] == 'E. Greenland','Neg'].values & siedf2.loc[siedf2['region'] == 'Hudson','Neg'].values & siedf2.loc[siedf2['region'] == 'St. Lawrence','Neg'].values)

# (Barents | Kara) & (G of St Law) --> 10.5%
100*np.mean( (siedf2.loc[siedf2['region'] == "Kara",'Neg'].values | siedf2.loc[siedf['region'] == "Barents",'Neg'].values) & (siedf2.loc[siedf2['region'] == "St. Lawrence",'Neg'].values) )

# (Barents & Kara) --> 8.4%
100*np.mean(siedf2.loc[siedf2['region'] == "Kara",'Neg'].values & siedf2.loc[siedf2['region'] == "Barents",'Neg'].values)

### Frequency of larger regions ###

siedf3 = siedf2.loc[siedf2['region'] == 'Kara',['year','month','day','JulianDay']]

## Combine Barents and Kara --> 20.0%
siedf3['BK'] = siedf2.loc[siedf2['region'] == 'Barents','cdr_delsie_5day'].values + siedf2.loc[siedf2['region'] == 'Kara','cdr_delsie_5day'].values
100*np.mean((siedf3['BK'] < 0))

## Combine All Atlantic --> 16.5%
siedf3['BKG'] = siedf3['BK'] + siedf2.loc[siedf2['region'] == 'E. Greenland','cdr_delsie_5day'].values
100*np.mean(siedf3['BKG'] < 0)

## Combine All Atlantic --> 16.5%
siedf3['Atl'] = siedf3['BK'] + siedf2.loc[siedf2['region'] == 'E. Greenland','cdr_delsie_5day'].values + siedf2.loc[siedf2['region'] == 'Baltic','cdr_delsie_5day'].values
100*np.mean(siedf3['Atl'] < 0)

## Combine Bering and Chukchi --> 14.8%
siedf3['BgCk'] = siedf2.loc[siedf2['region'] == 'Bering','cdr_delsie_5day'].values + siedf2.loc[siedf2['region'] == 'Chukchi','cdr_delsie_5day'].values
100*np.mean(siedf3['BgCk'] < 0)

# Combine Bering, Chukchi, Okhotsk --> 13.8%
siedf3['BgCkOk'] = siedf3['BgCk'] + siedf2.loc[siedf2['region'] == 'Okhotsk','cdr_delsie_5day'].values
np.mean(siedf3['BgCkOk'] < 0)

## Combine All Pacific --> 13.2%
siedf3['Pac'] = siedf3['BgCk'] + siedf2.loc[siedf2['region'] == 'Alaska','cdr_delsie_5day'].values + siedf2.loc[siedf2['region'] == 'Okhotsk','cdr_delsie_5day'].values + siedf2.loc[siedf2['region'] == 'Japan','cdr_delsie_5day'].values
100*np.mean(siedf3['Pac'] < 0)

## Combine Baffin & Labrador --> 8.5%
siedf3['BfLb'] = siedf2.loc[siedf2['region'] == 'Labrador','cdr_delsie_5day'].values + siedf2.loc[siedf2['region'] == 'Baffin','cdr_delsie_5day'].values
100*np.nanmean(siedf3['BfLb'] < 0)

## Combine Hudson, Baffin, Labrador, and GoSL # 5.0%
siedf3['Can'] = siedf2.loc[siedf2['region'] == 'Hudson','cdr_delsie_5day'].values + siedf2.loc[siedf2['region'] == 'Labrador','cdr_delsie_5day'].values + siedf2.loc[siedf2['region'] == 'Baffin','cdr_delsie_5day'].values + siedf2.loc[siedf2['region'] == 'St. Lawrence','cdr_delsie_5day'].values
100*np.nanmean(siedf3['Can'] < 0)

## Combine Lap & E Sib
siedf3['LES'] = siedf2.loc[siedf2['region'] == 'Laptev','cdr_delsie_5day'].values + siedf2.loc[siedf2['region'] == 'E. Siberian','cdr_delsie_5day'].values
100*np.mean(siedf3['LES'] < 0)

## Combine CAO, Beaufort, E Sib, Laptev, CAA # 26.7%
siedf3['CAO'] = siedf2.loc[siedf2['region'] == 'CAO','cdr_delsie_5day'].values + siedf2.loc[siedf2['region'] == 'CAA','cdr_delsie_5day'].values + siedf2.loc[siedf2['region'] == 'Beaufort','cdr_delsie_5day'].values + siedf2.loc[siedf2['region'] == 'Laptev','cdr_delsie_5day'].values + siedf2.loc[siedf2['region'] == 'E. Siberian','cdr_delsie_5day'].values
100*np.nanmean(siedf3['CAO'] < 0)

siedf3.to_csv(outpath+"/SIE_Daily_Regional3.csv", index=False)

from scipy.stats import spearmanr
siedf3 = siedf3.loc[np.isfinite(siedf3['Can'])]

## Atl + Pac ## --> 1.9%
np.mean((siedf3['Pac'] < 0) & (siedf3['Atl'] < 0)) * 100
spearmanr(siedf3['Atl'],siedf3['Pac'])

## Atl + Can ## --> 0.2%
np.mean((siedf3['Can'] < 0) & (siedf3['Atl'] < 0)) * 100
spearmanr(siedf3['Atl'],siedf3['Can'])

## Atl + CAO ##
np.mean((siedf3['CAO'] < 0) & (siedf3['Atl'] < 0)) * 100
spearmanr(siedf3['Atl'],siedf3['CAO'])

## Pac + Can ## --> 1.1%
np.mean((siedf3['Can'] < 0) & (siedf3['Pac'] < 0)) * 100
spearmanr(siedf3['Pac'],siedf3['Can'])

## Pac + CAO ##
np.mean((siedf3['CAO'] < 0) & (siedf3['Pac'] < 0)) * 100
spearmanr(siedf3['CAO'],siedf3['Pac'])

## Can + CAO ##
np.mean((siedf3['CAO'] < 0) & (siedf3['Can'] < 0)) * 100
spearmanr(siedf3['CAO'],siedf3['Can'])

## All ## --> 0.0%
np.mean((siedf3['Can'] < 0) & (siedf3['Pac'] < 0) & (siedf3['Atl'] < 0)) * 100


