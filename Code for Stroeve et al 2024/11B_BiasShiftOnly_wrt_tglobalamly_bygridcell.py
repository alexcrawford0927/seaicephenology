#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Alex Crawford
Date Created: 3 Aug 2023
Date Modified: 17 Aug 2023

Purpose: Performs a delta shift bias correction on sea ice phenology variables
by gridcell, defining the period of overlap between model and observations 
by the global temperature anomaly (rather than time).
"""

'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")
import xarray as xr
import numpy as np
import xesmf as xe
import pandas as pd
import CycloneModule_13_2 as md

'''*******************************************
Declare Variables
*******************************************'''
print("Declaring Variables")

modsource = "CMIP6"
exp1, exp2 = 'historical', 'ssp585'
obssource = ['Bootstrap','NASATeam']

obsfileyears = '1979-2021'
modyears2 = '2014-2099'
baseyears = [1979,2021]
refyears = [1850,1900]

othres = 'C15'
c6thres = '10cm'
varlist = ['cip','lrd','fad'] # ['lrd','fad'] # 
long_name = ['ice-free period (continuous open-water period)','last retreat day','first advance day'] # ['last retreat day','first advance day'] # 
nanval = 0

# Path Variables
obspath = ["/Volumes/Miranda/SeaIce/"+o+"/AdvanceRetreat2" for o in obssource]
modpath = "/Volumes/Cassandra/"+modsource+"/SeaIce/ThicknessPhenology/"
psnpath = "/Volumes/Miranda/Projections/psn_projection.nc"
ease2path = "/Volumes/Miranda/Projections/EASE2_N0_25km_Projection.nc"
ogpath = '/Volumes/Theseus/SurfaceTemperature/BEST/BEST_Annual_ts_1850-2022.csv'
c6tpath = '/Volumes/Cassandra/CMIP6/RegionalStats/V9/'

'''*******************************************
Main Analysis
*******************************************'''
print('Main Analysis')
# Load projections
psn = xr.open_dataset(psnpath)
ease2 = xr.open_dataset(ease2path)

regrid_psn_to_ease2 = xe.Regridder(psn, ease2, 'bilinear')

#### Read in sea ice observations -- average baseline period ####
oblist = [[] for v in varlist]
for oi in range(len(obssource)):
    # Load
    obsnc = xr.open_dataset(obspath[oi]+"/"+othres+"/siphenology"+othres+"_"+obssource[oi]+"_"+obsfileyears+".nc", decode_times = False)
    
    # Subset to baseline
    obsnc.sel(time=slice(baseyears[0],baseyears[1]))
    
    # Average each varible, keeping track of NA values
    for v in range(len(varlist)):
        if varlist[v] == 'cip':
            n = np.sum(obsnc['opc'] >= nanval, axis=0)
            obsnc['opc'].data[obsnc['opc'].data < nanval] = 0
            oblist[v].append( obsnc['opc'].sum(axis=0) / n ) 
        else:
            n = np.sum(obsnc[varlist[v]] >= nanval, axis=0)
            obsnc[varlist[v]].data[obsnc[varlist[v]].data < nanval] = 0
            oblist[v].append( obsnc[varlist[v]].sum(axis=0) / n )   
        
##### Take multi-method mean -- including reprojection ####
for v in range(len(varlist)):
    oblist[v] = regrid_psn_to_ease2( np.array(oblist[v]).mean(axis=0) )
    oblist[v][ease2['lat'].data < 0] = np.nan # set NaNs for SH

    if varlist[v] == 'cip':
        oblist[v] = 365 - oblist[v]

#### Read in observed temperature -- average anomaly for baseline period ####
ogdf = pd.read_csv(ogpath) # Global Temperature Anomaly
ogdf['taglobal'] = ogdf['tglobal'] - np.mean( ogdf[(ogdf['Year'] >= refyears[0]) & (ogdf['Year'] <= refyears[1])]['tglobal'] )

# Identify the max and min global temperature anoamly during the baseline period
tbmax = ogdf.loc[(ogdf['Year'] >= baseyears[0]) & (ogdf['Year'] <= baseyears[1])]['taglobal'].max()
tbmin = ogdf.loc[(ogdf['Year'] >= baseyears[0]) & (ogdf['Year'] <= baseyears[1])]['taglobal'].min()

#### Read in CMIP6 Temperature ####
tdf = pd.read_csv('/Volumes/Cassandra/CMIP6/CMIP6_Annual_tas_'+exp1+'.csv')
tdf = tdf.append(pd.read_csv('/Volumes/Cassandra/CMIP6/CMIP6_Annual_tas_'+exp2+'.csv'), ignore_index=False)
tdf['Model2'] = [m.split('_')[0] + '_' + m.split('_r')[1].split('i')[0] for m in tdf['Model'].values]
for model in np.unique(tdf.Model2):
    tmean = tdf[ (tdf['Model2'] == model) & (tdf['Year'] >= refyears[0]) & (tdf['Year'] <= refyears[1])]['tglobal'].mean()
    tdf.loc[(tdf['Model2'] == model),'taglobal'] = tdf.loc[(tdf['Model2'] == model),'tglobal'].values - tmean

#### Read in CMIP6 data -- reproject -- and calculate bias ####
files1 = md.listdir(modpath+"/"+exp1+"/"+c6thres)
files2 = md.listdir(modpath+"/"+exp2+"/"+c6thres)
files1 = [f for f in files1 if (int(f.split('_vl_r')[-1].split('i')[0]) == 1) and (f.split('_')[1] != 'CESM2') or ( (f.split('_')[1] == 'CESM2') and (int(f.split('_vl_r')[-1].split('i')[0]) == 4) )]
files2 = [f for f in files2 if (int(f.split('_vl_r')[-1].split('i')[0]) == 1) and (f.split('_')[1] != 'CESM2') or ( (f.split('_')[1] == 'CESM2') and (int(f.split('_vl_r')[-1].split('i')[0]) == 4) )]
mods1 = [f.split('_')[1] for f in files1]
mods2 = [f.split('_')[1] for f in files2]
mods = np.intersect1d(mods1,mods2)
mods = [mod for mod in mods if 'AWI' not in mod]

files1 = [f for f in files1 if f.split('_')[1] in mods]
files2 = [f for f in files2 if f.split('_')[1] in mods]

variants = [f.split('_')[5] for f in files1]

biascorrlist = [[] for v in varlist]
tcorrlist = []
for mi in range(len(files1)):
    print(files1[mi])
    
    # Load Models
    modfile1 = md.listdir(modpath+"/"+exp1+"/"+c6thres+"/"+files1[mi])[0]
    modfile2 = md.listdir(modpath+"/"+exp2+"/"+c6thres+"/"+files2[mi])[0]

    modnc1 = xr.open_dataset(modpath+"/"+exp1+"/"+c6thres+"/"+files1[mi]+"/"+modfile1, decode_times=False)
    modnc2 = xr.open_dataset(modpath+"/"+exp2+"/"+c6thres+"/"+files2[mi]+"/"+modfile2, decode_times=False)

    # Rename lat and lon dimensions to be 'lat' and 'lon'
    latname = [key for key in list(modnc1.variables) if key.startswith('lat') or key.startswith('nav_lat')][0]
    lonname = [key for key in list(modnc1.variables) if key.startswith('lon') or key.startswith('nav_lon')][0]
    modnc1 = modnc1.rename({latname:'lat',lonname:'lon'})
    modnc2 = modnc2.rename({latname:'lat',lonname:'lon'})
    
    # Combine
    modnc = xr.concat([modnc1,modnc2], dim='time')

    # Reprojection
    regridder = xe.Regridder(modnc1, ease2, 'bilinear',ignore_degenerate=True)
    modnc = regridder(modnc)

    # Identify years that have a temperature anomaly within the range of baseline period
    years = tdf.loc[(tdf['Family'] == modfile1.split('_')[1]) & (tdf['Member'] == int(modfile1.split('_vl_r')[1].split('i')[0])) & (tdf['taglobal'] >= tbmin) & (tdf['taglobal'] <= tbmax),'Year'].values
    
    # Append relevant taglobal data to output list
    tcorrlist.append( tdf.loc[(tdf['Family'] == modfile1.split('_')[1]) & (tdf['Member'] == int(modfile1.split('_vl_r')[1].split('i')[0])),'taglobal'] )
    
    # Find bias
    for v in range(len(varlist)):
        # Find average for the temperature years
        meanarr = np.nanmean(modnc[varlist[v]].data[np.in1d(modnc['time'],years),:,:], axis=0)
        
        # Find bias
        bias = meanarr - oblist[v]
        
        # Adjust by bias
        biascorrlist[v].append( modnc[varlist[v]] - bias )
        
#### Construct an xarray dataset out of the results ####
tarr = np.array([np.array(arr[:250]) for arr in tcorrlist])
yarr = np.arange(1850,2100)

# Add Data Variables
for v in range(len(varlist)):
    print("Writing " + varlist[v])
    xro = xr.Dataset(coords=dict(model=mods, variant=(['model'], variants),
                    time=yarr, y=ease2['y'], x=ease2['x'], 
                    lat=(['y','x'], ease2['lat'].data), lon=(['y','x'], ease2['lon'].data)),
               attrs=dict(description='''Bias-corrected ''' + long_name[v] + ''' values for
                          ''' + exp1 + ' and ' + exp2 +''' experiments relative to 
                          global annual mean surface temperature anomaly'''))
    xro['z'] = (('y','x'), ease2['z'].data)
    xro['tanom'] = (('model','time'), tarr)
    

    xro[varlist[v]] = (('model','time','y','x'), np.array(biascorrlist[v]).astype(np.float32))
    
    if varlist[v] == 'cip':
        xro[varlist[v]].data = np.where( xro[varlist[v]].data > 365, 365, xro[varlist[v]].data )
        xro[varlist[v]].data = np.where( xro[varlist[v]].data < 0, 0, xro[varlist[v]].data )    
        xro[varlist[v]].data = np.where( (xro['lat'] < 0), np.nan, xro[varlist[v]].data )    

        # xro[varlist[v]].data = np.where( np.isnan( xro[varlist[v]] ) & (xro['lat'] < 45) & (xro['z'] <= 0), 0, xro[varlist[v]].data )    
        # xro[varlist[v]].data = np.where( np.isnan( xro[varlist[v]] ) & (xro['lat'] < 60) & ((xro['lon'] < -100) | (xro['lon'] > -45) ) & (xro['z'] <= 0), 0, xro[varlist[v]].data )    
        # xro[varlist[v]].data = np.where( np.isnan( xro[varlist[v]] ) & (xro['lat'] < 75) & (xro['lon'] < 30) & (xro['lon'] > 0) & (xro['z'] <= 0), 0, xro[varlist[v]].data )    
        # xro[varlist[v]].data = np.where( np.isnan( xro[varlist[v]] ) & (xro['lat'] < 70) & (xro['lon'] < 0) & (xro['lon'] > -10) & (xro['z'] <= 0), 0, xro[varlist[v]].data )    
        # xro[varlist[v]].data = np.where( np.isnan( xro[varlist[v]] ) & (xro['lat'] < 65) & (xro['lon'] < -10) & (xro['lon'] > -15) & (xro['z'] <= 0), 0, xro[varlist[v]].data )    
        # xro[varlist[v]].data = np.where( np.isnan( xro[varlist[v]] ) & (xro['lat'] < 63) & (xro['lon'] < -15) & (xro['lon'] > -25) & (xro['z'] <= 0), 0, xro[varlist[v]].data )    
        # xro[varlist[v]].data = np.where( (xro[varlist[v]] == 0 ) & (xro['lat'] < 60) &  (xro['z'] <= 0), 0, xro[varlist[v]].data )    

    xro[varlist[v]].attrs = modnc[varlist[v]].attrs
    
    xro.to_netcdf(modpath+"/CMIP6_bias-corrected-by-t2mamly_"+exp1+"-"+exp2+"_baselineyears_"+str(baseyears[0])+"-"+str(baseyears[1])+"_"+varlist[v]+".nc")
