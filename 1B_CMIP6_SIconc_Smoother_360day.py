'''*********************************************
Authors: Alex Crawford
Date Created: 18 May 2020
Date Modified: 18 May 2020

Purpose: To calculate the day of the year on which sea ice retreats and
advances in each grid cell each year for models with 360-day calendars 
(only relevant for UKESM1)

Inputs: 
    n = the moving average size
    minlat = minimum latitude (e.g., 45 means 45°N)
    ncvar = the variable name for sea ice concentration ('siconc' or 'siconca')
    experiment = the experiment name
    exclude = model simulations that should be avoided because of grid 
    incompatibility
    
Outputs: A series of decadal netcdf files with a smoothed time series.
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
import MERRA_Module as md

'''*******************************************
Declare Variables
*******************************************'''
print("Declaring Variables")

### Input Variables ###
n = 5 # Moving Average Size (# observations total, including the central
# value -- should always be an odd number)
minlat = 45 # in degrees
ncvar = 'siconca'
experiment = 'ssp245'

### Path Variables ###
path = "/Volumes/Prospero/CMIP6/" # "E:/CMIP6" 
inpath = path+"/"+ncvar+"_day_"+experiment
outpath = path+"/SmoothedMA"+str(n)+"/"+experiment

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")
# Time Set Up
mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
days = ["01","02","03","04","05","06","07","08","09","10","11","12","13",\
    "14","15","16","17","18","19","20","21","22","23","24","25","26","27",\
    "28","29","30","31"]

# Read in list of ALL file names
files = os.listdir(inpath)
files = [f for f in files if ((f.endswith('.nc') == 1) & (f.startswith("siconca_SIday_UKESM1")))]
files.sort()

for f in files[:]:
    print(f)
    ncf = nc.Dataset(inpath+"/"+f)
    
    try:
        os.chdir(outpath+"/"+f[14:-21])
    except:
        os.mkdir(outpath+"/"+f[14:-21])
        os.chdir(outpath+"/"+f[14:-21])
    
    # Determine Start Year & Leap Years
    timenc = ncf.variables['time']
    startyear = int(timenc.units[11:15])
    endtime = float(timenc[-1])
        
    # Generate Lat & Lon
    lats = ncf.variables['lat'][:]
    lons = ncf.variables['lon'][:]
    
    if lats.ndim == 1:
        lons, lats = np.meshgrid(lons,lats)

    irows = np.unique(np.where((np.isfinite(lats) == 1) & (lats >= minlat))[0])
    
    # Identify starting years for each decade
    knownyears = [int(f[-20:-16])] + [int(f[-11:-7])]
    start, end = np.min(knownyears), np.max(knownyears)
    
    starts = np.arange(int(np.floor(start/10)*10),end+1,10) # Create starting years for each decade
    starts[0] = start # Ensure that first starting year is the first overall year
    ends = np.arange(int(np.ceil((start+1)/10)*10)-1,int(np.ceil(end/10)*10),10) # Create ending years in the 9 for each decade
    if len(ends) < len(starts):
        ends = np.concatenate((ends,np.array([end])))
    ends[-1] = end # Ensure that the final ending year is the final overall year
        
    # Loop through the netcdf file by decade
    for d in range(len(starts)):  
        istart, iend = int(starts[d]-starts[0])*360, int(ends[d]-starts[0]+1)*360
        arr0 = ncf[ncvar][istart:iend,irows,:].data
        arr0[(arr0 > 100.1) | (arr0 < -0.1)] = np.nan # Set NaNs
        validrows, validcols = np.where( (np.isfinite(arr0[0,:,:]) == 1) )
        
        # Create a faux-time that matches a 365-day calendar
        fauxtime = []
        for ii in range( 0, int(ends[d]-starts[d])+1 ):
            # Create the times relative to startyear if it were a 365-day calendar
            fstart = md.daysBetweenDates([startyear,1,1], [int(starts[d])+ii,1,1], 0)
            fauxyear = list(range(fstart,fstart+365,1))
            
            # Remove DOYs 2, 4, 183, 362, 364 so that the length is back to 360
            fauxtime = fauxtime+fauxyear[0:1]+fauxyear[2:3]+fauxyear[4:182]+fauxyear[183:361]+fauxyear[362:363]+fauxyear[364:]

        # Smoothing
        arrMA = np.zeros(arr0.shape)*np.nan
        print(" -- Smoothing at " + str(d) + " of " + str(len(starts)))
        for i in range(len(validrows)):
            arrMA[:,validrows[i],validcols[i]] = md.movingAverage2(arr0[:,validrows[i],validcols[i]],n)
        
        # Write new netcdf file
        ts = str(starts[d])
        te = str(ends[d]+1)
        fsmooth = "siconcMA"+str(n)+f[7:-20]+ts+"0101-"+te+"0101.nc"
        
        ncf1 = nc.Dataset(outpath+"/"+f[14:-21]+"/"+fsmooth, 'w', format='NETCDF4')
        ncf1.createDimension('y', arrMA.shape[1])
        ncf1.createDimension('x', arrMA.shape[2])
        ncf1.createDimension('time', arrMA.shape[0])
        
        yNC = ncf1.createVariable('y', np.float32, ('y',))
        xNC = ncf1.createVariable('x', np.float32, ('x',))
        tNC = ncf1.createVariable('time', np.float32, ('time',))
        
        sicNC = ncf1.createVariable('siconc', np.float32, ('time','y','x',))
        latNC = ncf1.createVariable('lat', np.float32, ('y','x',))
        lonNC = ncf1.createVariable('lon', np.float32, ('y','x',))
        
        ncf1.description = 'Sea Ice Concentration with 5-Day Moving Average'
        ncf1.source = 'netCDF4 python module'
        tNC.units = ncf.variables['time'].units #e.g., 'days since 1979-01-01 00:00:00.0'
        latNC.units = 'degrees north'
        lonNC.units = 'degrees east'      
        sicNC.units = 'Percentage'

        tNC[:] = fauxtime
        latNC[:] = lats[irows,:]
        lonNC[:] = lons[irows,:]
        sicNC[:] = arrMA

        ncf1.close()
        
        del arrMA, arr0, validrows, validcols