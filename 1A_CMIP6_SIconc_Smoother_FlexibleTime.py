'''*********************************************
Authors: Alex Crawford
Date Created: 6/9/15
Date Modified: 1/28/19; 8/20/19 edited for Python 3
18 May 2020 --> parameterized the minimum latitude
1 Aug 2020 --> Finished conversion to a more flexible format that can
    handle more input structures (e.g., annual v. decadal)
3 Aug 2020 --> Added a NaN statement for NESM model, reorganized time handling
    for MRI model (starts in 1919 insted of a multiple of 10), and added a
    conversion from decimal to percent if needed (for NESM, as well, but only some members)
7 Aug 2020 --> simplified the leapyear test (should be faster now) and made more flexible to
    allow for "siconca" as a a netcdf variable as well as "siconc"
21 Aug 2020 --> Changed the order for which lat/lon variables are assessed so that the BCC models
    have the right grid and changed the limits for SIC from 0 to 100 to -0.1 to 100.1 to account
    for the BCC models.
24 Aug 2020 --> Changed the "decfiles" definition to be based on both "s" and "e"
    variables instead of just "s". This allows for starts that are not multiples of 10.
25 Aug 2020 --> Changed the output time units to always be relative to 1850 (the original method
    was only a problem for models that store dates relative to the start date and had
    starts that were not multiples of 10... e.g., BCC-CSM2-MR for the ssp585 experiment).

Purpose: To calculate the day of the year on which sea ice retreats and
advances in each grid cell of a particular sector each year.

Inputs:
    n = the moving average size -- b/c of bi-daily data in the 1980s in the
    observational record, 5 or 7 is recommended.
    minlat = minimum latitude (e.g., 45 means 45°N)

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
import CycloneModule_13_2 as md

'''*******************************************
Declare Variables
*******************************************'''
### Input Variables ###
n = 5 # Moving Average Size (# observations total, including the central
# value -- should always be an odd number)
minlat = 45 # in degrees
ncvar = 'siconc'
experiment = 'historical' # 'ssp585' #'ssp126' # 'ssp370' # 'ssp245' #'ssp585' # 

### Path Variables ###
path =  '/project/6061839/crawfora/CMIP6/' # "/Volumes/Cassandra/CMIP6/" # "E:/CMIP6"
inpath = path+"/data/e_"+experiment+"/v_"+ncvar
outpath = path+"/SeaIce/SmoothedMA"+str(n)+"/"+experiment

exclude = ["AWI-CM-1-1-MR","AWI-ESM-1-1-LR","GISS-E2-1-G","r1i1p1f2_gr1_","MIROC-ES2L_e_"+experiment+"_vl_r2i1p1f2_gr1","MIROC-ES2L_e_"+experiment+"_vl_r3i1p1f2_gr1","ICON-"]
# include = ['CNRM-CM6-1-HR']

'''*******************************************
Main Analysis
*******************************************'''
# Read in list of ALL file names
models = md.listdir(inpath,contains='i1p1f')
models.sort()

for ex in exclude: # Exclude abnormal models
    models = [f for f in models if (ex not in f)]

# models = [f for f in models if f in include]

for model in models: #
    print(model)
    
    # Identify the number of model files
    modfiles = md.listdir(inpath+"/"+model)

    # Use the first file from that model to establish...
    # ... the outpath
    ncf = nc.Dataset(inpath+"/"+model+"/"+modfiles[0])
    
    outmod = modfiles[0].split('_')[2]+"_"+experiment+"_"+modfiles[0].split('_')[4]+'_'+modfiles[0].split('_')[5]

    try:
        os.chdir(outpath+"/"+outmod)
    except:
        os.mkdir(outpath+"/"+outmod)
        os.chdir(outpath+"/"+outmod)

    complete = os.listdir(outpath+"/"+outmod)

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
        ncf = nc.Dataset(inpath+"/"+model+"/"+f)
        timeslen += ncf['time'].shape[0]
        ncf.close()

    if timeslen%365 == 0 or "CESM" in outmod:
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
                ncf = nc.Dataset(inpath+"/"+model+"/"+ff)

                # Identify valid times in the netcdf file
                ncstart = int(ncf['time'].units.strip('days since ')[0:4])
                times = ncf['time'][:].data
                # times = np.arange(60225,91615) # Used for KIOST-ESM 2015-2100 in ssp126 only
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
                    arrMA[:,validrows[i],validcols[i]] = md.movingAverage(arr0[:,validrows[i],validcols[i]],n)
            else:
                 for i in range(len(validrows)):
                    arrMA[:,validrows[i]] = md.movingAverage(arr0[:,validrows[i]],n)

            # Write new netcdf file
            ncf1 = nc.Dataset(outpath+"/"+outmod+"/"+fsmooth, 'w', format='NETCDF4')
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