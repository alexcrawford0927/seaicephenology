'''*********************************************
Authors: Alex Crawford
Date Created: 18 May 2020
Date Modified: 25 Aug 2020

Purpose: To calculate the day of the year on which sea ice retreats and
advances in each grid cell each year in CMIP6 simulations

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
import os
import netCDF4 as nc
import numpy as np
import MERRA_Module as md

'''*******************************************
Declare Variables
*******************************************'''
### Input Variables ###
n = 5 # Moving Average Size (# observations total, including the central
# value -- should always be an odd number)
minlat = 45 # in degrees
ncvar = 'siconc'
experiment = 'ssp'

### Path Variables ###
path = "/Volumes/Prospero/CMIP6/" # "E:/CMIP6" 
inpath = path+"/"+ncvar+"_day_"+experiment
outpath = path+"/SmoothedMA"+str(n)+"/"+experiment

exclude = ["AWI-CM-1-1-MR","AWI-ESM-1-1-LR","GISS-E2-1-G","r1i1p1f2_gr1_","MIROC-ES2L_"+experiment+"_r2i1p1f2_gr1","MIROC-ES2L_"+experiment+"_r3i1p1f2_gr1"]

'''*******************************************
Main Analysis
*******************************************'''
# Read in list of ALL file names
files = os.listdir(inpath)
files = [f for f in files if ((f.endswith('.nc') == 1) & (f.startswith(".") == 0) & ("i1p1f" in f))]
files.sort()

for ex in exclude: # Exclude abnormal models
    files = [f for f in files if (ex not in f)]

# Identify unique models
models = np.unique([f.split(ncvar+'_SIday_')[1][:-21] for f in files])

for model in models[2:10]: #
    print(model)
    # Identify the number of model files
    modfiles = [f for f in files if (model in f)]
    
    # Use the first file from that model to establish...
    # ... the outpath
    ncf = nc.Dataset(inpath+"/"+modfiles[0])
    
    try:
        os.chdir(outpath+"/"+model)
    except:
        os.mkdir(outpath+"/"+model)
        os.chdir(outpath+"/"+model)
    
    complete = os.listdir(outpath+"/"+model)
        
    # ... the lats/lons
    try:
        lats = ncf.variables['latitude'][:].data
        lons = ncf.variables['longitude'][:].data
    except:
        try:
            lats = ncf.variables['lat'][:].data
            lons = ncf.variables['lon'][:].data
        except:
            lats = ncf.variables['nav_lat'][:].data
            lons = ncf.variables['nav_lon'][:].data
    
    if lats.ndim == 1:
        lons, lats = np.meshgrid(lons,lats)
    
    lats[(lats > 90) | (lats < -90)] = np.nan
    lons[(lons > 360) | (lons < -360)] = np.nan
    irows = np.unique(np.where((np.isfinite(lats) == 1) & (lats >= minlat))[0])
    
    # Establish Time Units                
    timeslen = ncf['time'].shape[0]
    ncf.close()
        
    # Establish whether there are leap years
    for f in modfiles[1:]:
        ncf = nc.Dataset(inpath+"/"+f)
        timeslen += ncf['time'].shape[0]
        ncf.close()

    if timeslen%365 == 0 or "CESM" in model:
        lyb = 0
    else:
        lyb = 1

    # Identify starting years for each decade
    knownyears = [int(f[-20:-16]) for f in modfiles] + [int(f[-11:-7]) for f in modfiles]
    if int(modfiles[-1][-7:-3]) >= 1230: # If the last date is Dec 30 or 31, include all years
        start, end = np.min(knownyears), np.max(knownyears)
    else: # otherwise, truncate the final year
        start, end = np.min(knownyears), np.max(knownyears)-1
    
    starts = np.arange(int(np.floor(start/10)*10),end+1,10) # Create starting years for each decade
    starts[0] = start # Ensure that first starting year is the first overall year
    ends = np.arange(int(np.ceil((start+1)/10)*10)-1,int(np.ceil(end/10)*10),10) # Create ending years in the 9 for each decade
    if len(ends) < len(starts):
        ends = np.concatenate((ends,np.array([end])))
    ends[-1] = end # Ensure that the final ending year is the final overall year
    
    for i in list(range(len(starts)))[:]: # For each decade... 
        s, e = starts[i], ends[i]
        
        #If this decade has already been smoothed, skip it
        fsmooth = "siconcMA"+str(n)+modfiles[0].split(ncvar)[1][:-20]+str(s)+"0101-"+str(e+1)+"0101.nc"
        if fsmooth in complete:
            continue
        
        else:
            # Identify needed files (starts before the end, ends after the start)
            decfiles = [f for f in modfiles if ( (int(f[-11:-7]) >= s) & (int(f[-20:-16]) < e+1) )]
            
            arrs = [] # container for arrays
            for ff in decfiles: # For each file... load the array
                ncf = nc.Dataset(inpath+"/"+ff)
                
                # Identify valid times in the netcdf file
                ncstart = int(ncf['time'].units.strip('days since ')[0:4])
                times = ncf['time'][:].data
                # times = np.arange(60225,91615) # Used for KIOST-ESM 2015-2100 only
                ti = np.where( (times >= md.daysBetweenDates([ncstart,1,1],[s,1,1],lyb)) & (times < md.daysBetweenDates([ncstart,1,1],[e,12,31],lyb)+1) )[0]
                
                # Load Data
                if len(ncf[ncvar].shape) == 3:
                    arrs.append( ncf[ncvar][ti,irows,:].data )
                else:
                    arrs.append( ncf[ncvar][ti,irows].data )
            
            # Concatenate arrays for all files and standardize structure as 0-100 with NaNs
            arr0 = np.concatenate(tuple(arrs),0) # Converts to 2 or 3-D array
            arr0[(arr0 > 100.1) | (arr0 < -0.1)] = np.nan # Set NaNs
            if np.nanmax(arr0) <= 1.1: # Convert from decimal to percentage
                arr0 = arr0*100
            
            # Identify valid sea ice locations
            if len(arr0.shape) == 3:
                validrows, validcols = np.where( (np.isfinite(arr0[0,:,:]) == 1) )
            else:
                validrows = np.where( (np.isfinite(arr0[0,:]) == 1) )
    
            # Smoothing
            arrMA = np.zeros(arr0.shape)*np.nan
            print(" -- Smoothing at " + str(s) + ", ends at " + str(ends[-1]))
            if len(arr0.shape) == 3:
                for i in range(len(validrows)):
                    arrMA[:,validrows[i],validcols[i]] = md.movingAverage2(arr0[:,validrows[i],validcols[i]],n)
            else:
                 for i in range(len(validrows)):
                    arrMA[:,validrows[i]] = md.movingAverage2(arr0[:,validrows[i]],n)
    
            # Write new netcdf file
            ncf1 = nc.Dataset(outpath+"/"+model+"/"+fsmooth, 'w', format='NETCDF4')
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
            tNC.units =  'days since 1850-01-01 00:00:00.0'
            latNC.units = 'degrees north'
            lonNC.units = 'degrees east'      
            sicNC.units = 'Percentage'
            
            startday = int(md.daysBetweenDates([1850,1,1,0,0,0],[s,1,1,0,0,0],lyb))
            tNC[:] = np.arange(startday,startday+arr0.shape[0],1)
            latNC[:] = lats[irows,:]
            lonNC[:] = lons[irows,:]
            sicNC[:] = arrMA
    
            ncf1.close()
            
            del arrMA, arr0, validrows, validcols