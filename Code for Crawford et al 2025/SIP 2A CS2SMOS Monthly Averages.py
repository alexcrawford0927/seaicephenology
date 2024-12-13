#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 16:24:10 2024

@author: acrawfora
"""

import xarray as xr
import numpy as np
import CycloneModule_13_2 as md

inpath = '/Volumes/Miranda/SeaIce/CS2SMOS'
outpath = '/Volumes/Miranda/SeaIce/CS2SMOS-Monthly'

ymin, ymax = 2010, 2010

for y in np.arange(ymin,ymax+1):
    
    files = md.listdir(inpath+"/"+str(y))
    
    mons = np.unique([f.split('_')[5][4:6] for f in files])
    
    for mon in mons:
        mfiles = [f for f in files if f.split('_')[5][4:6] == mon]
        
        xrlist = []
        for file in mfiles:
            xrlist.append( xr.open_dataset(inpath+"/"+str(y)+"/"+file) )
        
        xr2 = xr.concat(xrlist,dim='time').mean(dim='time')
        
        xr2.to_netcdf(outpath+"/CS2SMOS_"+str(y)+mon+".nc")
            
            
