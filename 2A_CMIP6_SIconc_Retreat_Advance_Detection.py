'''*********************************************
Authors: Alex Crawford
Date Created: 18 May 2020
Date Modified: 25 Aug 2020

Purpose: To calculate the day of the year on which sea ice retreats and
advances in each grid cell each year. Note: Must run the step 1 script first
to smooth the daily input data.

Inputs: 
    concThresh = concentration threshold - value between 0 and 100
    n = the moving average size (e.g., 5)
    maxmo, minmo = months in which to search for maximum and minimum SIC values
    ymax = last year wih valid data - must have data through the maxmo months
        of the following year!
    experiment = experiment for data - e.g., "ssp585", "historical"
    
Outputs: A netcdf file with fields for first/last retreat/advance day, the 
    length of the open-water period, and the date and value of the sea ice 
    minimum. Retreat and advance days are recorded as a "DOY", with 1 being
    the 1st day of January. If there is no retreat or advance, either 0 or 365
    is used as default values for the open-water period, but the retreat and
    advance are both assigned NA.
    
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
concThresh = 15 # A number between 0 and 100 for the concentration threshold
n = 5 # Moving Average Size


### Time Variables ###
ymax = 2013  # Last year of FULL data, meaning you need to have data from 
## the maxmo months of subsequent year, too
maxmo = [1,4] # months in which the sea ice maximum may occur
minmo = [8,10] # months in which the sea ice minimum may occur

### Path Variables ###
experiment = 'historical'
path = "/Volumes/Troilus"
inpath = path+"/CMIP6/SmoothedMA"+str(n)+"/"+experiment
outpath = path+"/CMIP6/AdvanceRetreat/"+experiment+"/C"+str(concThresh)

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")
# Time Set Up
mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
days = ["01","02","03","04","05","06","07","08","09","10","11","12","13",\
    "14","15","16","17","18","19","20","21","22","23","24","25","26","27",\
    "28","29","30","31"]

# Load file names
files = os.listdir(inpath)
files = [f for f in files if ((f.startswith('.') == 0) & (f.startswith("AWI") == 0))]

cmplfiles = os.listdir(outpath)
cmplfiles = [f for f in cmplfiles if f.startswith('.') == 0]

newfiles = [f for f in files if f not in cmplfiles]
newfiles.sort()

# Run main code
for f in newfiles[:]:
    ncfiles = os.listdir(inpath+"/"+f)
    ncfiles = [f for f in ncfiles if f.endswith('.nc')]
    ncfiles.sort()
    
    lrdL, frdL, ladL, fadL, minL, mndL, opcL, opL = [], [], [], [], [], [], [], []
    yearsnc = np.array([])

    for fi in range(len(ncfiles)):
        print("Starting " + ncfiles[fi])
        ncf = nc.Dataset(inpath+"/"+f+"/"+ncfiles[fi])
        
        # Determine Start Year & Leap Years
        times = ncf.variables['time'][:]
        reftime = [int(ncf.variables['time'].units[11:15]),1,1,0,0,0]
        
        if (times.shape[0]%365 == 0) or ("CESM2" in f) or (times.shape[0]%360 == 0) :
            lyb = 0
        else:
            lyb = 1
        
        # Force the first day to be day 1
        day1 = md.timeAdd(reftime,[0,0,times[0],0,0,0],lyb)[2]
        if day1 != 1:
            times = times - (day1-1)
        
        starttime = float(times[0])
        endtime = float(times[-1])
        
        # Establish Years
        years = np.arange( md.timeAdd(reftime, [0,0,starttime,0,0,0] ,lyb)[0] , md.timeAdd(reftime, [0,0,endtime,0,0,0] ,lyb)[0]+1 )
        
        # Limit to allowed years
        years = years[years <= ymax]
        
        # Store years used for export
        yearsnc = np.concatenate((yearsnc,years))
        
        # Load Data
        for y in years:
            # Subset for current year
            t1 = md.daysBetweenDates(reftime,[y,maxmo[0],1],lyb)
            t2 = md.daysBetweenDates(reftime,[y+1,maxmo[1]+1,1],lyb)
            tis = np.where((times >= t1) & (times < t2) )
            
            arr = ncf['siconc'][tis[0],:,:]
            
            # Subset for next year in cases where next year is in a different file
            if y == np.max(years) and fi+1 < len(ncfiles):
                ncf2 = nc.Dataset(inpath+"/"+f+"/"+ncfiles[fi+1])
                reftime2 = [int(ncf2.variables['time'].units[11:15]),1,1,0,0,0]
                times2 = ncf2.variables['time'][:] + md.daysBetweenDates(reftime,reftime2,lyb)
                
                arr2 = ncf2['siconc'][np.where((times2 >= t1) & (times2 < t2) )[0],:,:]
                arr = np.concatenate((arr, arr2), 0)
                
                ncf2.close()
                del times2, reftime2, arr2 # remove unneeded variables
        
            # Calculate minimum & maximum value for year
            maxi0 = md.daysBetweenDates([y,1,1],[y,maxmo[0],1], lyb)
            maxi1 = md.daysBetweenDates([y,1,1],[y,maxmo[1]+1,1], lyb)
            
            mini0 = md.daysBetweenDates([y,1,1],[y,minmo[0],1], lyb)
            mini1 = md.daysBetweenDates([y,1,1],[y,minmo[1]+1,1], lyb)
            
            maxi2 = md.daysBetweenDates([y,1,1],[y+1,maxmo[0],1], lyb)
            maxi3 = md.daysBetweenDates([y,1,1],[y+1,maxmo[1]+1,1], lyb)
            
            Maxes1 = np.amax(arr[maxi0:maxi1],0)
            Mins = np.amin(arr[mini0:mini1],0)
            Maxes2 = np.amax(arr[maxi2:maxi3],0)
            
            ### Prep Outputs ###
            mndArr, mxdArr = np.zeros_like(Mins)*np.nan, np.zeros_like(Mins)*np.nan
            frdArr, lrdArr = np.zeros_like(Mins)*np.nan, np.zeros_like(Mins)*np.nan
            fadArr, ladArr = np.zeros_like(Mins)*np.nan, np.zeros_like(Mins)*np.nan
            opArr, opcArr = np.zeros_like(Mins)*np.nan, np.zeros_like(Mins)*np.nan
            
            ### Calculate Retreat and Advance Events ###
            validcells = np.where( (np.isnan(Mins) == 0) & (np.isnan(Maxes1) == 0) & (np.isnan(Maxes2) == 0) )
            for i in range(len(validcells[0])):
                r,c = validcells[0][i], validcells[1][i] # Assign row and column
                
                # Calculate index for minimum & maximum
                MaxesI1 = int(np.median(np.where(arr[maxi0:maxi1,r,c] == Maxes1[r,c]))) # Gives first occurrence of maximum if multiples present
                MinsI = mini0 + int(np.median(np.where(arr[mini0:mini1,r,c] == Mins[r,c]))) # Gives first occurrence of minimum if multiples present
                MaxesI2 = maxi2 + int(np.median(np.where(arr[maxi2:maxi3,r,c] == Maxes2[r,c]))) # Gives first occurrence of maximum if multiples present
                
                # Store Minimum Day
                mxdArr[r,c] = MaxesI1 + 1
                mndArr[r,c] = MinsI + 1
                mxd2rr = MaxesI2 + 1
                
                # If it's always above the concentration threshold...
                if Mins[r,c] >= concThresh: 
                    opArr[r,c], opcArr[r,c] = 0, 0
                    
                # If it's never above the concentration threshold... 
                elif (Maxes1[r,c] < concThresh) & (Maxes2[r,c] < concThresh):
                    opArr[r,c], opcArr[r,c] = 365, 365
                    
                # Otherwise...
                else:
                    above = np.where(arr[:,r,c] >= concThresh)[0] # Indices above concentration
                    below = np.where(arr[:,r,c] < concThresh)[0] # Indices below concentration
                    
                    # First Retreat Day
                    # First index after Maxes1 and before/on Mins for which concentration is below threshold
                    try:
                        frdArr[r,c] = below[np.where((below <= MinsI) & (below > MaxesI1))][0] + 1
                    except:
                        frdArr[r,c] = np.nan
                    
                    # Last Retreat Day
                    # Last index after Maxes1 and before/on Mins for which concentration is below threshold
                    try:
                        lrdArr[r,c] = above[np.where((above < MinsI) & (above >= MaxesI1))][-1] + 1
                    except:
                        lrdArr[r,c] = np.nan
        
                    # First Advance Day
                    # First index after Mins and before/on Maxes2 for which concentration is above threshold
                    try: 
                        fadArr[r,c] = above[np.where((above > MinsI) & (above <= MaxesI2))][0] + 1
                    except:
                        fadArr[r,c] = np.nan
                    
                    # Last Advance Day
                    # Last index after Mins anbd before/on Maxes2 for which concentration is below threshold
                    try:
                        ladArr[r,c] = below[np.where((below >= MinsI) & (below < MaxesI2))][-1] + 1
                    except:
                        ladArr[r,c] = np.nan
                    
                    # Open Water Periods
                    if (Maxes1[r,c] < concThresh)  & (Maxes2[r,c] >= concThresh): # When it starts below threshold but ends above
                        opArr[r,c] =  ladArr[r,c] - mxdArr[r,c] #np.min([365, ladArr[r,c] - MaxesI1])
                        opcArr[r,c] = fadArr[r,c] - mxdArr[r,c] #np.min([365, fadArr[r,c] - MaxesI1])
                        
                        lrdArr[r,c] = np.nan
                        frdArr[r,c] = np.nan
                        
                    elif (Maxes1[r,c] >= concThresh)  & (Maxes2[r,c] < concThresh): # When it starts above threshold but ends below
                        opArr[r,c] =  mxd2rr - frdArr[r,c] #np.min([365, MaxesI2 - frdArr[r,c]])
                        opcArr[r,c] = mxd2rr - lrdArr[r,c] #np.min([365, MaxesI2 - lrdArr[r,c]])
                        
                        fadArr[r,c] = np.nan
                        ladArr[r,c] = np.nan
                        
                    else: # Simple Case
                        opArr[r,c] =  ladArr[r,c] - frdArr[r,c] #np.min([365, ladArr[r,c] - frdArr[r,c]])
                        opcArr[r,c] =  fadArr[r,c] - lrdArr[r,c] #np.min([365, fadArr[r,c] - lrdArr[r,c]])
                    
                    del above, below
            
            # Append to lists
            minL.append(Mins)
            mndL.append(mndArr)
            lrdL.append(lrdArr)
            ladL.append(ladArr)
            frdL.append(frdArr)
            fadL.append(fadArr)
            opcL.append(opcArr)
            opL.append(opArr)
            
            # Remove objects from memory
            del MinsI, validcells, MaxesI1, MaxesI2, mxdArr, mxd2rr, mndArr
            del frdArr, lrdArr, fadArr, ladArr, opArr, opcArr, Mins, Maxes1, Maxes2
    
    ### Write Outputs ###
    outName = "siphenologyC"+str(concThresh)+"_"+f+"_"+str(int(yearsnc[0]))+"-"+str(int(yearsnc[-1]))+".nc"
    try:
        os.chdir(outpath+"/"+f)
    except:
        os.mkdir(outpath+"/"+f)
        os.chdir(outpath+"/"+f)
    
    ncf1 = nc.Dataset(outpath+"/"+f+"/"+outName, 'w', format='NETCDF4')
    ncf1.createDimension('y', arr.shape[1])
    ncf1.createDimension('x', arr.shape[2])
    ncf1.createDimension('time', len(minL))
    
    yNC = ncf1.createVariable('y', np.float32, ('y',))
    xNC = ncf1.createVariable('x', np.float32, ('x',))
    tNC = ncf1.createVariable('time', np.float32, ('time',))
    
    minNC = ncf1.createVariable('min', np.float32, ('time','y','x',))
    mndNC = ncf1.createVariable('mnd', np.float32, ('time','y','x',))
    lrdNC = ncf1.createVariable('lrd', np.float32, ('time','y','x',))
    ladNC = ncf1.createVariable('lad', np.float32, ('time','y','x',))
    frdNC = ncf1.createVariable('frd', np.float32, ('time','y','x',))
    fadNC = ncf1.createVariable('fad', np.float32, ('time','y','x',))
    opcNC = ncf1.createVariable('opc', np.float32, ('time','y','x',))
    opNC = ncf1.createVariable('op', np.float32, ('time','y','x',))
    
    latNC = ncf1.createVariable('lat', np.float32, ('y','x',))
    lonNC = ncf1.createVariable('lon', np.float32, ('y','x',))
    
    ncf1.description = '''Phenology of Sea Ice Concentration with 5-Day Moving
    Average based on Concentration Threshold of ''' + str(concThresh) + '%'
    ncf1.source = 'netCDF4 python module'
    tNC.units = 'years'
    latNC.units = 'degrees north'
    lonNC.units = 'degrees east'      
    minNC.units = 'percentage'
    mndNC.units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
    lrdNC.units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
    frdNC.units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
    ladNC.units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
    fadNC.units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
    opcNC.units = 'days'
    opNC.units = 'days'
   
    tNC[:] = yearsnc
    latNC[:] = ncf.variables['lat'][:]
    lonNC[:] = ncf.variables['lon'][:]
    minNC[:] = np.array(minL)
    mndNC[:] = np.array(mndL)
    lrdNC[:] = np.array(lrdL)
    ladNC[:] = np.array(ladL)
    frdNC[:] = np.array(frdL)
    fadNC[:] = np.array(fadL)
    opcNC[:] = np.array(opcL)
    opNC[:] = np.array(opL)
    
    ncf1.close()
    ncf.close()
    
    now = clock()
    print(" -- Completed " + f + ", Elapsed Time: " + str(now-start))
    start = clock()
    
print("Complete.")