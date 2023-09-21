#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date Created: 4 Apr 2023
Date Modified: 8 Aug 2023

Applies the results of the model weighting script to grids of sea ice phenology
data to acquire weighted multi-model means for each grid cell.
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

obsfileyears = '1979-2021'
modyears2 = '2014-2099'
baseyears = [1979,2021]
refyears = [1850,1900]

othres = 'C15'
c6thres = '10cm'
varlist = ['cip','lrd','fad'] # ['lrd','fad'] # 
long_name = ['ice-free period (continuous open-water period)','last retreat day','first advance day'] # ['last retreat day','first advance day'] # 
nanval = 0

tcoords = np.arange(0,5.01,0.5)
tint  = 0.25

V2 = 'V9.57' # V+'.1'
sigma_d, sigma_s = 0.49, 0.50 # 0.23, 0.50
typ ='Thickness'
regvar = 'regHB3'

# Path Variables
modpath = "/Volumes/Cassandra/"+modsource+"/SeaIce/ThicknessPhenology/"
psnpath = "/Volumes/Miranda/Projections/psn_projection.nc"
ease2path = "/Volumes/Miranda/Projections/EASE2_N0_25km_Projection.nc"
ogpath = '/Volumes/Theseus/SurfaceTemperature/BEST/BEST_Annual_ts_1850-2022.csv'
c6tpath = '/Volumes/Cassandra/CMIP6/RegionalStats/V9/'

modstoskip = ['AWI-CM-1-1-MR']

'''*******************************************
Main Analysis
*******************************************'''
print('Main Analysis')
# Load projections
ease2 = xr.open_dataset(ease2path)

### Load Weights ###
wdf = pd.read_csv(c6tpath+"/ModelWeighting/Weights_"+typ+"-"+c6thres+"_"+regvar+"_"+exp1+'-'+exp2+"_Train"+str(baseyears[0])+"-"+str(baseyears[1])+"_"+V2+"_d"+str(sigma_d)+"_s"+str(sigma_s)+".csv")

# Limit models and re-calibrate weights to sum to 1
wdf = wdf.loc[np.in1d(wdf['Family'],modstoskip) == 0]
wdf['Wight'] = wdf['Weight'] / wdf['Weight'].sum()
wdf['Model'] = wdf['Family'] + '_' + wdf['Member'].astype(str)
mods = np.array(wdf['Model'])

#### Read in CMIP6 Temperature ####
tdf = pd.read_csv('/Volumes/Cassandra/CMIP6/CMIP6_Annual_tas_'+exp1+'.csv')
tdf = tdf.append(pd.read_csv('/Volumes/Cassandra/CMIP6/CMIP6_Annual_tas_'+exp2+'.csv'), ignore_index=False)
tdf['Model2'] = [m.split('_')[0] + '_' + m.split('_r')[1].split('i')[0] for m in tdf['Model'].values]
tdf = tdf.loc[np.in1d(tdf['Model2'],mods)]
for model in np.unique(tdf.Model2):
    tmean = tdf[ (tdf['Model2'] == model) & (tdf['Year'] >= refyears[0]) & (tdf['Year'] <= refyears[1])]['tglobal'].mean()
    tdf.loc[(tdf['Model2'] == model),'taglobal'] = tdf.loc[(tdf['Model2'] == model),'tglobal'].values - tmean

#### Read in CMIP6 data -- reproject -- bin by temperature ####
files1 = md.listdir(modpath+"/"+exp1+"/"+c6thres)
files2 = md.listdir(modpath+"/"+exp2+"/"+c6thres)
files1 = [f for f in files1 if (f.split('_')[1]+'_'+f.split('_vl_r')[-1].split('i')[0] in mods) ]
files2 = [f for f in files2 if (f.split('_')[1]+'_'+f.split('_vl_r')[-1].split('i')[0] in mods) ]

files1.sort()
files2.sort()

families = [f.split('_')[1] for f in files1]
variants = [f.split('_')[5] for f in files1]
models = [families[i]+"_"+variants[i].split('i')[0][1:] for i in range(len(families))]

weightedlist = [[] for mi in range(len(files1))]
weights = [[] for mi in range(len(files1))]
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
    
    # Combine and subset modnc vars to save memory
    modnc = xr.concat([modnc1[varlist],modnc2[varlist]], dim='time')
    
    # Define OPC
    # modnc['opc'] = 365 - modnc['cip']
    # modnc['opc'].data = np.where(modnc['opc'] > 365, 365, modnc['opc'])
    # modnc['opc'].data = np.where(modnc['opc'] < 0, 0, modnc['opc'])
    # modnc['opc'].attrs = dict(long_name='Continuous open-water period (aka ice-free period)', units='days')
    
    # Reprojection
    regridder = xe.Regridder(modnc1, ease2, 'bilinear', ignore_degenerate=True)
    modnc = regridder(modnc)
    
    for t in tcoords:
        years = tdf.loc[(tdf['Model2'] == models[mi]) & (tdf['taglobal'] > t-tint) & (tdf['taglobal'] < t+tint), 'Year' ].values
        arr = modnc.isel(time=np.in1d(modnc['time'].data, years)).mean(dim='time', skipna=False)
        for dvar in list(arr.data_vars):
            arr[dvar].data = np.mean( modnc.isel(time=np.in1d(modnc['time'].data, years))[dvar].data, axis=0 )
        
        weightedlist[mi].append( arr )

        if years.shape[0] == 0: # Assign 0 weight if this model does not experience the current temperature anomaly
            weights[mi].append( arr*0 )
        else: # Otherwise, use normal weight
            weights[mi].append(np.isfinite(arr)*wdf.loc[(wdf['Model'] == models[mi]),'Weight'].values)
            for dvar in list(arr.data_vars):
                weights[mi][-1][dvar].data = np.isfinite( arr[dvar] )
        
        weights[mi][-1] = weights[mi][-1].assign_coords({'model':np.array(mods[mi]),'tanom':np.array(t)})
            
#### Take Multi-model weighted mean ####
print('Multi Model Mean')
outlist = [[] for t in tcoords]
for v in range(len(varlist)):   
    
    for ti in range(len(tcoords)):
        
        varr = np.zeros_like( weightedlist[0][ti][varlist[v]] )
        narr = np.zeros_like( weightedlist[0][ti][varlist[v]] )
        
        for mi in range(len(weightedlist)):
            newarr = weightedlist[mi][ti][varlist[v]].data
            newarr[np.isnan(newarr)] = 0
            
            narr += np.isfinite( weightedlist[mi][ti][varlist[v]].data ) * weights[mi][ti][varlist[v]].data
            varr += newarr * weights[mi][ti][varlist[v]].data
        
        outlist[ti].append(varr / narr)
outarr = np.array(outlist)

weights = xr.concat([xr.concat(weights[mi], dim='tanom') for m in range(len(files1))], dim='model')

print('Writing File')
# Add Data Variables
xro = xr.Dataset(coords=dict(y=ease2['y'], x=ease2['x'], tanom=tcoords, model=models,
                lat=(['y','x'], ease2['lat'].data), lon=(['y','x'], ease2['lon'].data)),
            attrs=dict(description='''Weighted mean of sea ice phenology values for
                      ''' + exp1 + ' and ' + exp2 +''' experiments relative to 
                      global annual mean surface temperature anomaly'''))
                      
xro['z'] = (('y','x'), ease2['z'].data)

for v in range(len(varlist)):
    xro[varlist[v]] = (('tanom','y','x'), outarr[:,v,:,:].astype(np.float32))
    xro[varlist[v]].data = np.where(ease2['z'] > 0, np.nan, xro[varlist[v]].data)
    xro[varlist[v]].attrs = modnc[varlist[v]].attrs
    
    xro[varlist[v]+'_weight'] = (('model','tanom','y','x'), weights[varlist[v]].data)


v = 0
xro[varlist[v]].data = np.where( xro[varlist[v]].data > 365, 365, xro[varlist[v]].data )
xro[varlist[v]].data = np.where( xro[varlist[v]].data < 0, 0, xro[varlist[v]].data )    
# xro[varlist[v]].data = np.where( np.isnan( xro[varlist[v]] ) & (xro['lat'] < 60) & (xro['z'] <= 0), 365, xro[varlist[v]].data )    
# xro[varlist[v]].data = np.where( np.isnan( xro[varlist[v]] ) & (xro['lat'] < 75) & (xro['lon'] < 30) & (xro['lon'] > 0) & (xro['z'] <= 0), 0, xro[varlist[v]].data )    
# xro[varlist[v]].data = np.where( np.isnan( xro[varlist[v]] ) & (xro['lat'] < 70) & (xro['lon'] < 0) & (xro['lon'] > -10) & (xro['z'] <= 0), 0, xro[varlist[v]].data )    
# xro[varlist[v]].data = np.where( np.isnan( xro[varlist[v]] ) & (xro['lat'] < 65) & (xro['lon'] < -10) & (xro['lon'] > -15) & (xro['z'] <= 0), 0, xro[varlist[v]].data )    
# xro[varlist[v]].data = np.where( np.isnan( xro[varlist[v]] ) & (xro['lat'] < 63) & (xro['lon'] < -15) & (xro['lon'] > -25) & (xro['z'] <= 0), 0, xro[varlist[v]].data )    
# xro[varlist[v]].data = np.where( (xro[varlist[v]] == 0 ) & (xro['lat'] < 60) &  (xro['z'] <= 0), 0, xro[varlist[v]].data )    

xro.to_netcdf(modpath+"/CMIP6_weightedaverage_tanombins_"+exp1+"-"+exp2+"_baselineyears_"+str(baseyears[0])+"-"+str(baseyears[1])+".nc")


