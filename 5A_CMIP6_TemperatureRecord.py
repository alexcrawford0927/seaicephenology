'''*********************************************
Authors: Alex Crawford
Date Created: 6 Aug 2020
Date Modified: 6 Aug 2020
Purpose: Makes a CSV file of the regional and global annual temperature 
    for CMIP6 models from monthly data

*********************************************'''
# Import clock:
from time import perf_counter as clock
# Start script stopwatch. The clock starts running when time is imported
start = clock()

'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import os
import netCDF4 as nc
import numpy as np
import pandas as pd

np.seterr(all='ignore')

'''*******************************************
Declare Variables
*******************************************'''
print("Declaring Variables")

var = 'tas' # CMIP6 variable name 
frequency = 'mon'
tableid = 'Amon' # May be identical to frequency
experiment = 'ssp245'

bbox = [60,-180,90,360] # Min Lat, Min Lon, Max Lat, Max Lon

path =  "/Volumes/Troilus/CMIP6"
areapath = path+"/GridAreas/Atmosphere"
inpath = path+"/"+var+"_"+frequency+"_"+experiment
outpath = path

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")
# Check if output already started; initiate if not
outname = "CMIP6_Annual_"+var+"_"+experiment+".csv"
try:
    pdf = pd.read_csv(outpath+"/"+outname)
    completed = np.unique(pdf['Model'])
except:
    pdf = pd.DataFrame()
    completed = []

# Identify valid models and their area files
files = os.listdir(inpath)
files = [f for f in files if (f.startswith(var+"_"+tableid) and (f.split(var+"_"+tableid+"_")[1][:-17] not in completed))]
files.sort()

areas = os.listdir(areapath)

#### Loop through each file ####
for f in files:
    print(f)
    # Identify model & member
    FAM = f.split(var+"_"+tableid+"_")[1].split("_"+experiment)[0]
    MOD = f.split(var+"_"+tableid+"_")[1][:-17]
    MEM = int(f.split(experiment+"_r")[1].split("i")[0])
    
    # Load file
    ncf = nc.Dataset(inpath+"/"+f,'r')
    
    # Load grid cell area
    area = pd.read_pickle(areapath+"/"+[a for a in areas if a[6:-4] == FAM][0])
    area = np.array([area,area,area, area,area,area, area,area,area, area,area,area])
    
    # Load the time dimension
    time = ncf['time'][:].data
    
    # Load the spatial dimensions
    try:
        lons = ncf['lon'][:].data
        lats = ncf['lat'][:].data
    except:
        try:
            lats = ncf.variables['latitude'][:]
            lons = ncf.variables['longitude'][:]
        except:
            lats = ncf.variables['nav_lat'][:]
            lons = ncf.variables['nav_lon'][:]
    
    # Create a grid mesh for lat and lon
    if len(lats.shape) == 1:
        lons, lats = np.meshgrid(lons,lats)
    
    # Create a mask for the regional calculation
    mask = np.where((lons >= bbox[1]) & (lons <= bbox[3]) & (lats >= bbox[0]) & (lats <= bbox[2]))
    
    if len(ncf[var].shape) == 2:
        print("Error:" + f + " has 2 dimensions, expecting 3")
    
    # For each year in the model...
    else:
        for i in [t for t in list(range(ncf['time'].shape[0])) if t%12 == 0]:
            # Take the annual average
            arr = ncf[var][i:(i+12),:,:]
            
            # Take spatial averages (with area-weighting)
            mean_global = np.sum(arr*area)/np.sum(area)
            mean_regional = np.sum(arr[:,mask[0],mask[1]]*area[:,mask[0],mask[1]])/np.sum(area[:,mask[0],mask[1]])
            
            # Append to data frame
            pdf = pdf.append(pd.DataFrame([{'Model':MOD,'Family':FAM,'Member':MEM,\
                'Year':int(f[-16:-12]) + int(i/12), 'tglobal':mean_global, 'tregion':mean_regional}]),\
                ignore_index=True)

    pdf.to_csv(outpath+"/"+outname,index=False)
