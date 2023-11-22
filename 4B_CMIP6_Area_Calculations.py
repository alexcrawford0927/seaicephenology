"""
Author: Alex Crawford
Date Created: 31 Jul 2020
Date Modified: 31 Jul 2020
Purpose: Creates a numpy array with area values (in sq. km) for each grid cell
 in CMIP6 models
"""

'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import numpy as np
import pandas as pd
import netCDF4 as nc
import os
import CycloneModule_13_2 as md

'''#####################
 Declare Variables
#####################'''


lonmin, lonmax = -180, 360
latmin, latmax = -90,90
minlat = -90 # -90 for atmosphere, 45 for NH sea ice

var = 'tas'
frequency = 'mon'
tableid = 'Amon' # may be identical to frequency, but not always
experiment = 'historical'

### Path Variables ###
path = "/Volumes/Troilus"
inpath = path+"/CMIP6/"+var+"_"+frequency+"_"+experiment
outpath = "/Volumes/Troilus/CMIP6/GridAreas/Atmosphere"

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

################
# Find Models! #
################
pretext = var+"_"+tableid+"_"
posttext = "_"+experiment

# Load model names
files = os.listdir(inpath)
files = [f for f in files if ( (f.startswith('.') == 0) and (pretext in f) and ("AWI" not in f))]

models = np.unique([f.split(pretext)[1].split(posttext)[0] for f in files])

cmplmods = os.listdir(outpath)
cmplmods = [f.split(posttext)[0][6:-4]  for f in cmplmods if f.startswith('.') == 0]

newmods = [f for f in models if f not in cmplmods]
newmods.sort()

################
# Calc Areas!! #
################

for mod in newmods:
    print(mod)
    
    # 1) Load lat & lon
    modfiles = [f for f in files if tableid+"_"+mod+"_"+experiment in f]
    ncf = nc.Dataset(inpath+"/"+modfiles[0])
    
    try:
        lons = ncf['lon'][:].data
        lats = ncf['lat'][:].data
    except:
        try:
            lats = ncf.variables['latitude'][:].data
            lons = ncf.variables['longitude'][:].data
        except:
            lats = ncf.variables['nav_lat'][:].data
            lons = ncf.variables['nav_lon'][:].data
        
    lons[(lons < lonmin) | (lons > lonmax)] = np.nan
    lons = np.where(lons > 180, lons-360, lons)
    lats[(lats < latmin) | (lats > latmax)] = np.nan
    
    # Create a grid mesh
    if len(lats.shape) == 1:
        lons, lats = np.meshgrid(lons,lats)
        
    # Trim based on minimum latitude
    irows = np.unique(np.where((np.isfinite(lats) == 1) & (lats >= minlat))[0])
    lats = lats[irows,:]
    lons = lons[irows,:]
    
    # 2) Make arrays that estimate grid cell corners & center
    nrow, ncol = lats.shape
    
    # Center
    latsC = np.hstack(( np.zeros((nrow+2, 1))*np.nan , \
            np.vstack(( np.zeros((1,ncol))*np.nan, lats ,np.zeros((1,ncol))*np.nan )), \
            np.zeros((nrow+2, 1))*np.nan ))
    lonsC = np.hstack(( np.zeros((nrow+2, 1))*np.nan , \
            np.vstack(( np.zeros((1,ncol))*np.nan, lons ,np.zeros((1,ncol))*np.nan )), \
            np.zeros((nrow+2, 1))*np.nan ))
        
    # Upper-left
    latsUL = np.hstack(( np.zeros((nrow+2, 0))*np.nan , \
            np.vstack(( np.zeros((0,ncol))*np.nan, lats ,np.zeros((2,ncol))*np.nan )), \
            np.zeros((nrow+2, 2))*np.nan ))
    lonsUL = np.hstack(( np.zeros((nrow+2, 0))*np.nan , \
            np.vstack(( np.zeros((0,ncol))*np.nan, lons ,np.zeros((2,ncol))*np.nan )), \
            np.zeros((nrow+2, 2))*np.nan ))
           
    # Upper-right
    latsUR = np.hstack(( np.zeros((nrow+2, 0))*np.nan , \
            np.vstack(( np.zeros((2,ncol))*np.nan, lats ,np.zeros((0,ncol))*np.nan )), \
            np.zeros((nrow+2, 2))*np.nan ))
    lonsUR = np.hstack(( np.zeros((nrow+2, 0))*np.nan , \
            np.vstack(( np.zeros((2,ncol))*np.nan, lons ,np.zeros((0,ncol))*np.nan )), \
            np.zeros((nrow+2, 2))*np.nan ))
    
    # Lower-right
    latsLR = np.hstack(( np.zeros((nrow+2, 2))*np.nan , \
            np.vstack(( np.zeros((2,ncol))*np.nan, lats ,np.zeros((0,ncol))*np.nan )), \
            np.zeros((nrow+2, 0))*np.nan ))
    lonsLR = np.hstack(( np.zeros((nrow+2, 2))*np.nan , \
            np.vstack(( np.zeros((2,ncol))*np.nan, lons ,np.zeros((0,ncol))*np.nan )), \
            np.zeros((nrow+2, 0))*np.nan ))
    
    # Lower-left
    latsLL = np.hstack(( np.zeros((nrow+2, 2))*np.nan , \
            np.vstack(( np.zeros((0,ncol))*np.nan, lats ,np.zeros((2,ncol))*np.nan )), \
            np.zeros((nrow+2, 0))*np.nan ))
    lonsLL = np.hstack(( np.zeros((nrow+2, 2))*np.nan , \
            np.vstack(( np.zeros((0,ncol))*np.nan, lons ,np.zeros((2,ncol))*np.nan )), \
            np.zeros((nrow+2, 0))*np.nan ))  

    # Calculate cell corner locations from shifted arrays
    lats1 = (latsC + latsUL)/2
    lats2 = (latsC + latsUR)/2
    lats3 = (latsC + latsLR)/2
    lats4 = (latsC + latsLL)/2
    
    lons1 = (lonsC + lonsUL)/2
    lons2 = (lonsC + lonsUR)/2
    lons3 = (lonsC + lonsLR)/2
    lons4 = (lonsC + lonsLL)/2
    
    # Band-aid where the longitude crosses +/- 180
    lons1 = np.where( (lonsC*lonsUL < 1) & (np.abs(lonsC) > 90), (360+lonsC+lonsUL)/2, lons1)
    lons2 = np.where( (lonsC*lonsUR < 1) & (np.abs(lonsC) > 90), (360+lonsC+lonsUR)/2, lons2)
    lons3 = np.where( (lonsC*lonsLR < 1) & (np.abs(lonsC) > 90), (360+lonsC+lonsLR)/2, lons3)
    lons4 = np.where( (lonsC*lonsLL < 1) & (np.abs(lonsC) > 90), (360+lonsC+lonsLL)/2, lons4)

    # Calculate areas
    area = np.zeros_like(lats)*np.nan
    
    for ri in range(1,lats1.shape[0]-1):
        for ci in range(1,lats1.shape[1]-1):
            lati = np.array([lats1[ri,ci],lats2[ri,ci],lats3[ri,ci],lats4[ri,ci]])
            loni = np.array([lons1[ri,ci],lons2[ri,ci],lons3[ri,ci],lons4[ri,ci]])
            
            # Adjust if there are both negative and positive longitudes over magnitude 90
            if np.min(loni) < -90 and np.max(loni) > 90:
                loni[loni < 0] = 360+loni[loni < 0]
            
            area[ri-1,ci-1] = md.area_polygon_latlong(lati,loni)
    
    # Estimate edges
    area[0,:] = area[1,:] + (area[1,:] - area[2,:])
    area[-1,:] = area[-2,:] + (area[-2,:] - area[-3,:])
    area[:,0] = area[:,1] + (area[:,1] - area[:,2])
    area[:,-1] = area[:,-2] + (area[:,-2] - area[:,-3])
    
    pd.to_pickle(area,outpath+"/Areas_"+mod+".pkl")

    
    