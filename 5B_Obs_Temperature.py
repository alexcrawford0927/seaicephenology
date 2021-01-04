'''*********************************************
Authors: Alex Crawford
Date Created: 6 Aug 2020
Date Modified: 6 Aug 2020
Purpose: Makes a CSV file of the regional and global annual temperature 
    for each observational dataset

*********************************************'''

'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import netCDF4 as nc
import numpy as np
import pandas as pd
import MERRA_Module as md

'''*******************************************
Declare Variables
*******************************************'''
print("Declaring Variables")

path = "/Volumes/Troilus/SurfaceTemperature/"
bbox = [60,-180,90,360] # Min Lat, Min Lon, Max Lat, Max Lon
arrmin, arrmax = -99, 99 # Any value equal to or higher magnitude than these limits is set to NA
ymax = 2019 # Maximum allowed year 

'''*******************************************
Main Analysis
*******************************************'''
###############
### GISTEMP ###
###############
print("GISTEMP")

ncf = nc.Dataset(path+"/GISTEMP/gistemp1200_GHCNv4_ERSSTv5.nc")
ncf.variables.keys()

# Explore Time
ncf['time'].units
ncf['time'].shape
startyear = md.timeAdd([1800,1,1,0,0,0],[0,0,ncf['time'][0],0,0,0],lys=1)[0]

# Load lats and lons and regional mask
lons, lats = np.meshgrid(ncf['lon'][:],ncf['lat'][:])
mask = np.where((lons >= bbox[1]) & (lons <= bbox[3]) & (lats >= bbox[0]) & (lats <= bbox[2]))

# Establish data frame
pdf = pd.DataFrame()

for i in [j for j in list(range(ncf['time'].shape[0])) if ((j%12 == 0) & (startyear+np.floor(j/12) <= ymax))]: # Loop via each January
    # Take the annual average
    marr= ncf['tempanomaly'][i:(i+12),:,:].data
    marr[(marr <= arrmin) | (marr >= arrmax)] = np.nan
    arr = np.apply_along_axis(np.mean,0,marr)

    # Take spatial averages (with area-weighting)
    mean_global = np.nansum(arr*np.cos(lats*np.pi/180))/np.sum(np.isfinite(arr)*np.cos(lats*np.pi/180))
    mean_regional = np.nansum(arr[mask]*np.cos(lats[mask]*np.pi/180))/np.sum(np.isfinite(arr[mask])*np.cos(lats[mask]*np.pi/180))
    
    # Append to data frame
    pdf = pdf.append(pd.DataFrame([{'Year':int(startyear+np.floor(i/12)), 'tglobal':mean_global, \
                                    'tregion':mean_regional}]),ignore_index=True)

pdf.to_csv(path+"/GISTEMP/GISTEMP_Annual_ts_"+str(startyear)+"-"+str(int(startyear+np.floor(i/12)))+".csv",index=False)

del pdf, arr, marr, mean_global, mean_regional, mask, lons, lats, startyear
ncf.close()

########################
### NOAA Global Temp ###
########################
print("NOAAGlobalTemp")
ncf = nc.Dataset(path+'/NOAAGlobalTemp/NOAAGlobalTemp_v5.0.0_gridded_s188001_e202006_c20200707T133318.nc')
ncf.variables.keys()

# Explore Time
ncf['time'].units
ncf['time'].shape
startyear = md.timeAdd([1800,1,1,0,0,0],[0,0,ncf['time'][0],0,0,0],lys=1)[0]

# Load lats and lons and regional mask
lons, lats = np.meshgrid(ncf['lon'][:],ncf['lat'][:])
mask = np.where((lons >= bbox[1]) & (lons <= bbox[3]) & (lats >= bbox[0]) & (lats <= bbox[2]))

# Establish data frame
pdf = pd.DataFrame()

for i in [j for j in list(range(ncf['time'].shape[0])) if ((j%12 == 0) & (startyear+np.floor(j/12) <= ymax))]: # Loop via each January
    # Take the annual average
    marr= ncf['anom'][i:(i+12),0,:,:].data
    marr[(marr <= arrmin) | (marr >= arrmax)] = np.nan
    arr = np.apply_along_axis(np.mean,0,marr)

    # Take spatial averages (with area-weighting)
    mean_global = np.nansum(arr*np.cos(lats*np.pi/180))/np.sum(np.isfinite(arr)*np.cos(lats*np.pi/180))
    mean_regional = np.nansum(arr[mask]*np.cos(lats[mask]*np.pi/180))/np.sum(np.isfinite(arr[mask])*np.cos(lats[mask]*np.pi/180))
    
    # Append to data frame
    pdf = pdf.append(pd.DataFrame([{'Year':int(startyear+np.floor(i/12)), 'tglobal':mean_global, \
                                    'tregion':mean_regional}]),ignore_index=True)

pdf.to_csv(path+"/NOAAGlobalTemp/NOAAGlobalTemp_Annual_ts_"+str(startyear)+"-"+str(int(startyear+np.floor(i/12)))+".csv",index=False)

del pdf, arr, marr, mean_global, mean_regional, mask, lons, lats, startyear
ncf.close()

###############
### HadCRUT ###
###############
print("HadCRUT")
ncf = nc.Dataset(path+'/HadCRUT/HadCRUT.4.6.0.0.median.nc')
ncf.variables.keys()

# Explore Time
ncf['time'].units
ncf['time'].shape
startyear = md.timeAdd([1850,1,1,0,0,0],[0,0,ncf['time'][0],0,0,0],lys=1)[0]

# Load lats and lons and regional mask
lons, lats = np.meshgrid(ncf['longitude'][:],ncf['latitude'][:])
mask = np.where((lons >= bbox[1]) & (lons <= bbox[3]) & (lats >= bbox[0]) & (lats <= bbox[2]))

# Establish data frame
pdf = pd.DataFrame()

for i in [j for j in list(range(ncf['time'].shape[0])) if ((j%12 == 0) & (startyear+np.floor(j/12) <= ymax))]: # Loop via each January
    # Take the annual average
    marr= ncf['temperature_anomaly'][i:(i+12),:,:].data
    marr[(marr <= arrmin) | (marr >= arrmax)] = np.nan
    arr = np.apply_along_axis(np.mean,0,marr)

    # Take spatial averages (with area-weighting)
    mean_global = np.nansum(arr*np.cos(lats*np.pi/180))/np.sum(np.isfinite(arr)*np.cos(lats*np.pi/180))
    mean_regional = np.nansum(arr[mask]*np.cos(lats[mask]*np.pi/180))/np.sum(np.isfinite(arr[mask])*np.cos(lats[mask]*np.pi/180))
    
    # Append to data frame
    pdf = pdf.append(pd.DataFrame([{'Year':int(startyear+np.floor(i/12)), 'tglobal':mean_global, \
                                    'tregion':mean_regional}]),ignore_index=True)

pdf.to_csv(path+"/HadCRUT/HadCRUT_Annual_ts_"+str(startyear)+"-"+str(int(startyear+np.floor(i/12)))+".csv",index=False)

del pdf, arr, marr, mean_global, mean_regional, mask, lons, lats, startyear
ncf.close()

############
### BEST ###
############
print("Berkeley (BEST)")
# ncf = nc.Dataset(path+'/BEST/BEST_Complete_TAVG_LatLong1.nc')
ncf = nc.Dataset(path+"/BEST/BEST_Land_and_Ocean_LatLong1.nc")
ncf.variables.keys()

ncf['time'].units
ncf['time'].shape
times = np.floor(ncf['time'][:].data)
startyear = int(times[0])
endyear = int(times[-1])
if endyear > ymax:
    endyear = ymax

# Load lats and lons and regional mask
lons, lats = np.meshgrid(ncf['longitude'][:],ncf['latitude'][:])
mask = np.where((lons >= bbox[1]) & (lons <= bbox[3]) & (lats >= bbox[0]) & (lats <= bbox[2]))

# Establish data frame
pdf = pd.DataFrame()

for y in range(startyear,endyear+1):
    # Take the annual average
    marr= ncf['temperature'][np.where(times == y)[0],:,:].data
    arr = np.apply_along_axis(np.mean,0,marr)

    # Take spatial averages (with area-weighting)
    mean_global = np.nansum(arr*np.cos(lats*np.pi/180))/np.sum(np.isfinite(arr)*np.cos(lats*np.pi/180))
    mean_regional = np.nansum(arr[mask]*np.cos(lats[mask]*np.pi/180))/np.sum(np.isfinite(arr[mask])*np.cos(lats[mask]*np.pi/180))
    
    # Append to data frame
    pdf = pdf.append(pd.DataFrame([{'Year':y, 'tglobal':mean_global, \
                                    'tregion':mean_regional}]),ignore_index=True)

pdf.to_csv(path+"/BEST/BEST_Annual_ts_"+str(startyear)+"-"+str(endyear)+".csv",index=False)

del pdf, arr, marr, mean_global, mean_regional, mask, lons, lats, startyear, endyear, times
ncf.close()
