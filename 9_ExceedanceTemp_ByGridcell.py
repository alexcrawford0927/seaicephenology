"""
Author: Alex Crawford
Date Created: 16 Oct 2020
Date Modified: 16 Oct 2020
Purpose: Calculate the temperature at which a variable exceeds some threshold 
value in CMIP6 models. Models are also re-projected to a common grid in order
to make a grid-cell-by-grid-cell multi-model-mean when plotting. Only first
ensemble members are used in this script.

Inputs:
var = variable of interest (e.g., 'opc' for continuous open-water period)
varName = Human-readable version of variable name
tvar = variable for temperature ('tglobal' or 'tregion')

ct = concentration threshold for defining sea ice variables
experiment1, experiment2 = experiments to use for linking variable to temperature anomalies
exclude1, exclude2 = any models to exclude from this analysis

threslist = list of threshold values for the given variable
miny, maxy = inclusive range of years for which all models have data for both temperature 
    and sea ice variable -- note, some models only have sea ice data back to 1950
miny2, maxy2 = inclusive range of years used for the temperature baseline period
    (note that the IPCC uses 1850-1900)

Outputs: Netcdf file for each threshold with one dimension being the list of
CMIP6 models.

"""

'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import numpy as np
import pandas as pd
import os
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata

np.seterr(all='ignore')

'''#####################
 Declare Variables
#####################'''

ct = 15
experiment1 = 'historical'
experiment2 = 'ssp585'
exclude1 = ["CMCC-CM2-SR5_"+experiment1+"_r1i1p1f1_gn",'CanESM5_'+experiment1+'_r1i1p2f1_gn','MRI-ESM2-0_'+experiment1+'_r1i2p1f1_gn','MRI-ESM2-0_'+experiment1+'_r1i1000p1f1_gn']
exclude2 = ["CMCC-CM2-SR5_"+experiment2+"_r1i1p1f1_gn",'CanESM5_'+experiment2+'_r1i1p2f1_gn','MRI-ESM2-0_'+experiment2+'_r1i2p1f1_gn','MRI-ESM2-0_'+experiment2+'_r1i1000p1f1_gn']
var = 'opc'
tvar = 'tglobal'
varName = "Open Water Period"

threslist = [90,180,270,360] # The threshold that the variable must cross
miny, maxy = 1950, 2099 # Used to limit to common period for all models
miny2, maxy2 = 1850, 1900 # Used for temperature baseline

inpath1 = "/Volumes/Troilus/CMIP6/AdvanceRetreat/"+experiment1+"/C"+str(ct)
inpath2 = "/Volumes/Troilus/CMIP6/AdvanceRetreat/"+experiment2+"/C"+str(ct)
outpath = "/Volumes/Troilus/CMIP6/RegionalStats/ExceedanceTemp/C"+str(ct)

tpath1 = "/Volumes/Troilus/CMIP6/CMIP6_Annual_tas_"+experiment1+".csv"
tpath2 = "/Volumes/Troilus/CMIP6/CMIP6_Annual_tas_"+experiment2+".csv"

# Inputs for reprojection
minlat, maxlat = 45, 90
bb = [-89.99,-45,-89.99,135] # in degrees [ll lat, ll lon, ur lat, ur lon]
xsize, ysize = 50000, -50000 # in meters
nx, ny = int(180*(100000/xsize)), int(180*(100000/xsize)) # number of grid cells
lon_0 = 0 # Central Meridian (which longitude is at the 6 o'clock position)
lat_0 = 90 # Reference Latitude (center of projection)

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

# Load Temperature data and split up the model variable
tdf1 = pd.read_csv(tpath1)
tdf1 = tdf1.loc[~tdf1['Model'].isin(exclude1)] # The tilde is a "NOT" operator
tdf1['Family'] = tdf1['Model'].str.split('_'+experiment1,expand=True)[0]
tdf1['Member'] = tdf1['Model'].str.split('_'+experiment1,expand=True)[1].str.split('i1p',expand=True)[0].str.split('i2p',expand=True)[0].str.split('_r',expand=True)[1].astype(int)

tdf2 = pd.read_csv(tpath2)
tdf2 = tdf2.loc[~tdf2['Model'].isin(exclude2)] # The tilde is a "NOT" operator
tdf2['Family'] = tdf2['Model'].str.split('_'+experiment2,expand=True)[0]
tdf2['Member'] = tdf2['Model'].str.split('_'+experiment2,expand=True)[1].str.split('i1p',expand=True)[0].str.split('i2p',expand=True)[0].str.split('_r',expand=True)[1].astype(int)

# Load output projection data
projnc = nc.Dataset("/Volumes/Troilus/Projections/EASE2_N0_"+str(int(xsize/1000))+"km_Projection.nc")
outlons = projnc['lon'][:]
outlats = projnc['lat'][:]

# Create map
mp = Basemap(projection='laea',lat_0=lat_0,lon_0=lon_0,\
    llcrnrlat=bb[0], llcrnrlon=bb[1],urcrnrlat=bb[2], urcrnrlon=bb[3],resolution='c')

# Identify the aligned models (exist in both experiments)
models2 = os.listdir(inpath2)
models2.sort()
models2 = [model for model in models2 if model.startswith('.') == 0 and '_r1i' in model and model not in exclude2]
models22 = [model.split('_'+experiment2)[0] for model in models2]

models1 = os.listdir(inpath1)
models1.sort()
models1 = [model for model in models1 if model.startswith('.') == 0 and '_r1i' in model and model not in exclude1]
models12 = [model.split('_'+experiment1)[0] for model in models1]

models = np.intersect1d(models12,models22)

models1 = [models1[i] for i in range(len(models1)) if models12[i] in models]
models2 = [models2[i] for i in range(len(models2)) if models22[i] in models]

# For each model...
for thres in threslist:
    print(thres)
    arrlist = []
    for i in range(len(models)):
        ### Prepare Temperatures ###
        # Load Temperature Data (eliminate 2014 b/c we couldn't detect sea ice phenology for it)        
        t1 = tdf1.loc[(tdf1['Family'] == models[i]) & (tdf1['Member'] == 1) & (tdf1['Year'] >= miny) & (tdf1['Year'] != 2014) & (tdf1['Year'] <= maxy), tvar]
        t2 = tdf2.loc[(tdf2['Family'] == models[i]) & (tdf2['Member'] == 1) & (tdf2['Year'] >= miny) & (tdf2['Year'] != 2014) & (tdf2['Year'] <= maxy), tvar]
        t = t1.append(t2,ignore_index=1).sort_values()
        indices = np.array(t.index)
        
        # Convert to a temperature anomaly based on model's historical run so that baseline is 1850-1900
        tavg2 = tdf1.loc[(tdf1['Family'] == models[i]) & (tdf1['Member'] == 1) & (tdf1['Year'] >= miny2) & (tdf1['Year'] <= maxy2),tvar].mean() # Anomaly wrt 1979-1998
        tavg = np.array(t) - tavg2
        
        # Identify the min and max allowed values
        vmin = np.min(tavg) - np.std(tavg)/2
        vmax = np.max(tavg) + np.std(tavg)/2
        
        ### Load each model's netcdf file (1st ensemble members) ###
        ncpath1 = os.listdir(inpath1+"/"+models1[i])
        ncpath2 = os.listdir(inpath2+"/"+models2[i])
        
        file1 = [f for f in ncpath1 if f.startswith("siphenologyC"+str(ct)+"_"+models1[i])][0]
        file2 = [f for f in ncpath2 if f.startswith("siphenologyC"+str(ct)+"_"+models2[i])][0]
      
        ncf1 = nc.Dataset(inpath1+"/"+models1[i]+"/"+file1,'r')
        ncf2 = nc.Dataset(inpath2+"/"+models2[i]+"/"+file2,'r')
    
        ### Concatenate model arrays ###
        arr = np.concatenate((ncf1[var][:].data,ncf2[var][:].data),axis=0)
        times = np.concatenate((ncf1['time'][:].data,ncf2['time'][:].data),axis=0)
        arr2 = arr[times >= miny][indices,:,:]
        
        ### Identify the first temperature at which the variable is above the given threshold ###
        varexceed = np.zeros((arr.shape[1],arr.shape[2]))*np.nan
        
        jrows, jcols = np.where(np.isfinite(arr2[0,:,:])) 
        for j in range(len(jrows)): # Apply separately to each ocean grid cell
            if np.min(arr2[:,jrows[j],jcols[j]]) >= thres: # If the variable always exceeds the threshold
                varexceed[jrows[j],jcols[j]] = vmin
            elif arr2[-1,jrows[j],jcols[j]] < thres: # If the variable doesn't exceed the treshold at the warmest temperature
                varexceed[jrows[j],jcols[j]] = vmax
            else: # If sometimes the variable is above and sometimes below threshold, find the coldest temperature at or above which the variable is always above the threshold
                varexceed[jrows[j],jcols[j]] = tavg[ np.where(arr2[:,jrows[j],jcols[j]] < thres)[0][-1]+1 ]
        
        ### Reproject Output to a Common Grid - Step 1 - Convert to Regular Lat/Lon Grid ###    
        # Load Lats and Lons
        inlons = ncf1['lon'][:].data
        inlats = ncf1['lat'][:].data
        
        # Convert Lon to -180 to 180
        if np.nanmax(inlons) > 180:
            inlons = np.where(inlons > 180, inlons-360, inlons)
        if np.nanmin(inlons) < -180:
            inlons = np.where(inlons < -180, inlons+360, inlons)
        
        latinterval = (maxlat-minlat)/inlats.shape[0]
        lats2 = np.arange(minlat+latinterval/2,maxlat,latinterval)
    
        loninterval = 360/inlons.shape[1]
        lons2 = np.arange(-180+loninterval/2,180,loninterval)
        
        lons22, lats22 = np.meshgrid(lons2,lats2)
        
        varexceed_reg = griddata((inlons.flatten(),inlats.flatten()),varexceed.flatten(),(lons22,lats22),method='nearest')
        
        ### Reproject Output to a Common Grid - Step 2 - Convert to Regular Lat/Lon Grid ###    
        outArr = mp.transform_scalar(varexceed_reg,lons2,lats2,nx,ny,returnxy=False)
        outArr[outlats < minlat] = np.nan
        
        arrlist.append(outArr)
    
    ### Write to File ###
    ncf = nc.Dataset(outpath+"/ExceedanceTemp_C"+str(ct)+"_"+var+str(thres)+"_"+experiment2+".nc", 'w')
    ncf.description = 'Degree at or above which ' + varName + ' (concentration threshold = ' + str(ct)+'%) always exceeds a value of ' + str(thres)
    ncf.source = 'netCDF4 python module'
    
    ncf.createDimension('model', len(arrlist))
    ncf.createDimension('x', nx)
    ncf.createDimension('y', ny)
    ncft = ncf.createVariable('model', np.int, ('model',))
    ncfx = ncf.createVariable('x', np.float64, ('x',))
    ncfy = ncf.createVariable('y', np.float64, ('y',))
    ncfArr = ncf.createVariable(var, np.float64, ('model','y','x'))
    
    ncft.units = 'A CMIP6 model run'+str(models)
    ncfx.units = 'm'
    ncfy.units = 'm'
    ncfArr.units = 'Degrees C'
    
    # For x and y, note that the upper left point is the edge of the grid cell, but
    ## for this we really want the center of the grid cell, hence dividing by 2.
    ncft[:] = np.arange(len(models))
    ncfx[:] = np.arange(-xsize*(nx-1)/2, xsize*(nx-1)/2+xsize, xsize)
    ncfy[:] = np.arange(-ysize*(ny-1)/2, ysize*(ny-1)/2+ysize, ysize)
    ncfArr[:] = np.array(arrlist)
    
    ncf.close()