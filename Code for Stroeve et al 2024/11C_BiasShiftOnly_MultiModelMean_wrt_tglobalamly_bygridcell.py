#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Alex Crawford
Date Created: 3 Aug 2023
Date Modified: 17 Aug 2023

Purpose: Calculates a multi-model mean (by grid-cell) for preivously bias-
corrected sea ice phenology fields.
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

modsource = 'CMIP6'
exp1, exp2 = 'historical', 'ssp585'
modpath = "/Volumes/Cassandra/"+modsource+"/SeaIce/ThicknessPhenology/"
ease2path = "/Volumes/Miranda/Projections/EASE2_N0_25km_Projection.nc"

tcoords = np.arange(0.5,4.01,0.5)
tint  = 0.25

baseyears = [1979,2021]

varlist = ['cip','lrd','fad'] # ['lrd','fad'] # 
long_name = ['ice-free period (continuous open-water period)','last retreat day','first advance day'] # ['last retreat day','first advance day'] # 


ease2 = xr.open_dataset(ease2path)


### Initiate Output ###
xro = xr.Dataset(coords=dict(y=ease2['y'], x=ease2['x'], tanom=tcoords,
                lat=(['y','x'], ease2['lat'].data), lon=(['y','x'], ease2['lon'].data)),
           attrs=dict(description='''Multi-model mean of bias-corrected sea ice 
                      phenology values for
                      ''' + exp1 + ' and ' + exp2 +''' experiments relative to 
                      global annual mean surface temperature anomaly'''))
xro['z'] = (('y','x'), ease2['z'].data)


for v in range(len(varlist)):
    print(varlist[v])
    
    # Load data for this variable
    ds = xr.open_dataset(modpath+"/CMIP6_bias-corrected-by-t2mamly_"+exp1+"-"+exp2+"_baselineyears_"+str(baseyears[0])+"-"+str(baseyears[1])+"_"+varlist[v]+".nc")

    tlist = []
    for t in tcoords:
        print(t)
    
        mlist = []
        for m in ds['model']:
            # Subset to model and to years of temperature bin
            dssub = ds.sel(model=m).isel(time=np.where((ds.sel(model=m)['tanom'] > t-tint) & (ds.sel(model=m)['tanom'] < t+tint))[0])
            
            # Take temporal average, store
            mlist.append( dssub[varlist[v]].mean(dim='time').data )
        
        # Take multi-model mean, store
        tlist.append ( np.nanmean( np.array(mlist), axis = 0) )

    # Store in ouput dataset
    xro[varlist[v]] = (('tanom','y','x'), np.array(tlist).astype(np.float32))
    xro[varlist[v]]['long_name'] = long_name[v]

xro.to_netcdf(modpath+"/"+modsource+"_bias-corrected-by-t2mamly_"+exp1+"-"+exp2+"_baselineyears_"+str(baseyears[0])+"-"+str(baseyears[1])+"_multimodelmean.nc")

        