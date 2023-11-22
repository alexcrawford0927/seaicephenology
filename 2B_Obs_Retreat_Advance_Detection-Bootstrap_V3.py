'''*********************************************
Authors: Alex Crawford
Date Created: 6/9/15
Date Modified: 1/28/19; 8/20/19 edited for Python 3;
6/4/20 switched to movingAverage2
6/11/20 simplified the maximum day storage to be an integer, not an array,
since the array was never saved
7/7/20 adapted for Bootstrap method and converted output to netcdf instead of geotiff
1/21/22 added x and y coordinates using nsidc.org/data/polar-stereo/ps_grids.html

Purpose: To calculate the day of the year on which sea ice retreats and
advances in each grid cell of each year.

Inputs:
    ct = concentration threshold - value between 0 and 1
    n = the moving average size (e.g., 5)
    maxmo, minmo = months in which to search for maximum and minimum SIC values

Outputs: One netcdf file for all years
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
import pandas as pd
import netCDF4 as nc
import numpy as np
import CycloneModule_13_2 as md

'''*******************************************
Declare Variables
*******************************************'''
print("Declaring Variables")

### Input Variables ###
ct = 0.80 # A number between 0 and 1 for the concentration threshold
n = 5 # Moving Average Size (# observations total, including the central
# value -- should always be an odd number)
ver = "Bootstrap"

### Time Variables ###
ymin, ymax = 1979, 2021 # Years of interest
maxmo = [1,4] # months in which the sea ice maximum may occur
minmo = [8,10] # months in which the sea ice minimum may occur

### Path Variables ###
path = "/Volumes/Miranda/SeaIce/"+ver
inpath = path+"/Daily"
inpath2 = path+"/DailyFiller"
mapath = path+"/SmoothedMA"+str(n)
outpath = path+"/AdvanceRetreat2"
suppath = "/Volumes/Miranda/Projections"

#### File Variables ###
maskMin = 2
maskMax = 19

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")
# Load Lat and Long Rasters
projnc = nc.Dataset(suppath+"/psn_projection.nc")
lons = projnc['lon'][:].data
lats = projnc['lat'][:].data
areas = projnc['area'][:].data
mask = projnc['reg'][:].data

# Time Set Up
years = range(ymin,ymax+1)

noretreat = np.zeros(len(years))
files = os.listdir(inpath) + os.listdir(inpath2)
files = [f for f in files if f.startswith('bt')]
files.sort()

# Prep outputs
fadlist, frdlist, ladlist, lrdlist, opclist, oplist, mndlist, minlist = [], [], [], [], [], [], [], []

for y in years:
    print(" Year: "+str(y))
    dateref = [y,1,0]

    ### Prep Inputs ###
    # Load Current Year
    try: # If a smoothed SIC array already exists, open it
        arr1 = pd.read_pickle(mapath+"/SIC_SmoothedMA"+str(n)+"_"+str(y)+".pkl")

    except: # If not, you must create it anew first
        # Load Annual Files
        fileyList = [f for f in files if ( (int(f[3:7]) == y) )]
        fileyList_2 = [f for f in files if ( (int(f[3:7]) == y+1 and int(f[7:9]) == 1 and int(f[9:11]) <= n) )]
        fileyList_0 = [f for f in files if ( (int(f[3:7]) == y-1 and int(f[7:9]) == 12 and int(f[9:11]) > 31-n) )]
        # Load Arrays into a single 3-D array
        arrList = []
        for f in fileyList_0 + fileyList + fileyList_2: # Skip the last file because it's only referenced if no advance occurs
            # Reclassify array and mask out anything that's not in the regions of the mask
            try:
                arr0 = np.fromfile(inpath+'/'+f,dtype='uint16')[:].reshape((448,304))/1000
            except:
                arr0 = np.fromfile(inpath2+'/'+f,dtype='uint16')[:].reshape((448,304))/1000

            arrList.append(  np.where((mask <= maskMax) & (mask >= maskMin) & (arr0 <= 1), arr0, np.nan ) )

        # Smooth array
        if len(fileyList_2) == 0:
            ei = len(fileyList_0)+len(fileyList)
        else:
            ei = -len(fileyList_2)

        arr1 = np.apply_along_axis(md.movingAverage2,0,np.array(arrList),n)[len(fileyList_0):ei,:,:]
        pd.to_pickle(arr1,mapath+"/SIC_SmoothedMA"+str(n)+"_"+str(y)+".pkl")

    # Load Following Year
    try:
        arr2 = pd.read_pickle(mapath+"/SIC_SmoothedMA"+str(n)+"_"+str(y+1)+".pkl")

    except: # If not, you must create it anew first
        print("--- Applying moving average... (this may take a minute or two)")
        # Load Annual Files
        fileyList = [f for f in files if ( (int(f[3:7]) == y+1) )]
        fileyList_2 = [f for f in files if ( (int(f[3:7]) == y+2 and int(f[7:9]) == 1 and int(f[9:11]) <= n) )]
        fileyList_0 = [f for f in files if ( (int(f[3:7]) == y and int(f[7:9]) == 12 and int(f[9:11]) > 31-n) )]
        # Load Arrays into a single 3-D array
        arrList = []
        for f in fileyList_0 + fileyList + fileyList_2: # Skip the last file because it's only referenced if no advance occurs
            # Reclassify array and mask out anything that's not in the regions of the mask
            try:
                arr0 = np.fromfile(inpath+'/'+f,dtype='uint16')[:].reshape((448,304))/1000
            except:
                arr0 = np.fromfile(inpath2+'/'+f,dtype='uint16')[:].reshape((448,304))/1000

            arrList.append(  np.where((mask <= maskMax) & (mask >= maskMin) & (arr0 <= 1), arr0, np.nan ) )

        # Smooth array
        if len(fileyList_2) == 0:
            ei = len(fileyList_0)+len(fileyList)
        else:
            ei = -len(fileyList_2)

        arr2 = np.apply_along_axis(md.movingAverage2,0,np.array(arrList),n)[len(fileyList_0):ei,:,:]
        pd.to_pickle(arr2,mapath+"/SIC_SmoothedMA"+str(n)+"_"+str(y+1)+".pkl")

    # Load file list
    fileList = [f for f in files if ( ( (int(f[3:7]) == y) or (int(f[3:7]) == y+1 ) ))]
    doys = np.array( [md.daysBetweenDates([y,1,1],[int(f[3:7]),int(f[7:9]),int(f[9:11])]) for f in fileList] )

    arr = np.concatenate( (arr1, arr2), axis=0 )

    # Fill in Pole Hole
    for i in range(arr.shape[0]):
        holeR, holeC = np.where((lats >= 83) & (mask == 15) & (np.isnan(arr[i,:,:])))
        fillR, fillC = np.where((lats >= np.min(lats[holeR,holeC])-1) & (np.isfinite(arr[i,:,:])))
        arr[i,holeR,holeC] = np.nanmean(arr[i,fillR,fillC])

    # Calculate minimum & maximum value for year
    maxi0 = np.where( doys >= md.daysBetweenDates([y,1,1],[y,maxmo[0],1]) )[0][0]
    maxi1 = np.where( doys <= md.daysBetweenDates([y,1,1],[y,maxmo[1]+1,1]) )[0][-1]

    mini0 = np.where( doys >= md.daysBetweenDates([y,1,1],[y,minmo[0],1]) )[0][0]
    mini1 = np.where( doys <= md.daysBetweenDates([y,1,1],[y,minmo[1]+1,1]) )[0][-1]

    maxi2 = np.where( doys >= md.daysBetweenDates([y,1,1],[y+1,maxmo[0],1]) )[0][0]
    maxi3 = np.where( doys <= md.daysBetweenDates([y,1,1],[y+1,maxmo[1]+1,1]) )[0][-1]

    Maxes1 = np.amax(arr[maxi0:maxi1],0)
    Mins = np.amin(arr[mini0:mini1],0)
    Maxes2 = np.amax(arr[maxi2:maxi3],0)

    print("--- Calculating Retreat/Advance...")
    CT = str(int(ct*100))
    ### Prep Outputs ###
    mndArr = np.zeros_like(lats)-99
    frdArr, lrdArr = np.zeros_like(lats)-99, np.zeros_like(lats)-99
    fadArr, ladArr = np.zeros_like(lats)-99, np.zeros_like(lats)-99
    opArr, opcArr = np.zeros_like(lats)-99, np.zeros_like(lats)-99

    ### Calculate Retreat and Advance Events ###
    validcells = np.where( (np.isnan(Mins) == 0) & (np.isnan(Maxes1) == 0) & (np.isnan(Maxes2) == 0) )
    for i in range(len(validcells[0])):
        r,c = validcells[0][i], validcells[1][i] # Assign row and column

        # Calculate index for minimum & maximum
        MaxesI1 = int(np.median(np.where(arr[maxi0:maxi1,r,c] == Maxes1[r,c]))) # Gives first occurrence of minimum if multiples present
        MinsI = mini0 + int(np.median(np.where(arr[mini0:mini1,r,c] == Mins[r,c]))) # Gives first occurrence of minimum if multiples present
        MaxesI2 = maxi2 + int(np.median(np.where(arr[maxi2:maxi3,r,c] == Maxes2[r,c]))) # Gives first occurrence of minimum if multiples present

        # Store Minimum Day
        mndArr[r,c] = md.daysBetweenDates(dateref,[int(fileList[MinsI][3:7]),int(fileList[MinsI][7:9]),int(fileList[MinsI][9:11]),0,0,0])
        mxd1rr = md.daysBetweenDates(dateref,[int(fileList[MaxesI1][3:7]),int(fileList[MaxesI1][7:9]),int(fileList[MaxesI1][9:11]),0,0,0])
        mxd2rr = md.daysBetweenDates(dateref,[int(fileList[MaxesI2][3:7]),int(fileList[MaxesI2][7:9]),int(fileList[MaxesI2][9:11]),0,0,0])

        # If it's always above the concentration threshold...
        if Mins[r,c] >= ct:
            opArr[r,c], opcArr[r,c] = 0, 0

            noretreat[y-min(years)] = noretreat[y-min(years)]+1

        # If it's never above the concentration threshold...
        elif (Maxes1[r,c] < ct) & (Maxes2[r,c] < ct):
            opArr[r,c], opcArr[r,c] = 365, 365

        # Otherwise...
        else:
            above = np.where(arr[:,r,c] >= ct)[0] # Indices above concentration
            below = np.where(arr[:,r,c] < ct)[0] # Indices below concentration

            # First Retreat Day
            # First index after Maxes1 and before/on Mins for which concentration is below threshold
            try:
                fri = below[np.where((below <= MinsI) & (below > MaxesI1))][0]
                frdArr[r,c] = md.daysBetweenDates(dateref,[int(fileList[fri][3:7]),int(fileList[fri][7:9]),int(fileList[fri][9:11]),0,0,0])
            except:
                frdArr[r,c] = -99

            # Last Retreat Day
            # Last index after Maxes1 and before/on Mins for which concentration is below threshold
            try:
                lri = above[np.where((above < MinsI) & (above >= MaxesI1))][-1]
                lrdArr[r,c] = md.daysBetweenDates(dateref,[int(fileList[lri][3:7]),int(fileList[lri][7:9]),int(fileList[lri][9:11])+1,0,0,0])
            except:
                lrdArr[r,c] = -99

            # First Advance Day
            # First index after Mins and before/on Maxes2 for which concentration is above threshold
            try:
                fai = above[np.where((above > MinsI) & (above <= MaxesI2))][0]
                fadArr[r,c] = md.daysBetweenDates(dateref,[int(fileList[fai][3:7]),int(fileList[fai][7:9]),int(fileList[fai][9:11]),0,0,0])
            except:
                fadArr[r,c] = -99

            # Last Advance Day
            # Last index after Mins anbd before/on Maxes2 for which concentration is below threshold
            try:
                lai = below[np.where((below >= MinsI) & (below < MaxesI2))][-1]
                ladArr[r,c] = md.daysBetweenDates(dateref,[int(fileList[lai][3:7]),int(fileList[lai][7:9]),int(fileList[lai][9:11])+1,0,0,0])
            except:
                ladArr[r,c] = -99

            # Open Water Periods
            if (Maxes1[r,c] < ct)  & (Maxes2[r,c] >= ct): # When it starts below threshold but ends above
                opArr[r,c] =  ladArr[r,c] - mxd1rr #np.min([365, ladArr[r,c] - MaxesI1])
                opcArr[r,c] = fadArr[r,c] - mxd1rr #np.min([365, fadArr[r,c] - MaxesI1])

                lrdArr[r,c] = -99
                frdArr[r,c] = -99

            elif (Maxes1[r,c] >= ct)  & (Maxes2[r,c] < ct): # When it starts above threshold but ends below
                opArr[r,c] =  mxd2rr - frdArr[r,c] #np.min([365, MaxesI2 - frdArr[r,c]])
                opcArr[r,c] = mxd2rr - lrdArr[r,c] #np.min([365, MaxesI2 - lrdArr[r,c]])

                fadArr[r,c] = -99
                ladArr[r,c] = -99

            else: # Simple Case
                opArr[r,c] =  ladArr[r,c] - frdArr[r,c] #np.min([365, ladArr[r,c] - frdArr[r,c]])
                opcArr[r,c] =  fadArr[r,c] - lrdArr[r,c] #np.min([365, fadArr[r,c] - lrdArr[r,c]])

    ### Store Outputs ###
    minlist.append(Mins)
    mndlist.append(mndArr)
    frdlist.append(frdArr)
    fadlist.append(fadArr)
    lrdlist.append(lrdArr)
    ladlist.append(ladArr)
    opclist.append(opcArr)
    oplist.append(opArr)

    # Remove objects from memory
    del fai, lai, fri, lri, frdArr, lrdArr, fadArr, ladArr, opArr, opcArr
    del MinsI, validcells, above, below, MaxesI1, MaxesI2, mxd1rr, mxd2rr, mndArr

    del Maxes1, Maxes2, Mins, doys
    del maxi0, maxi1, maxi2, maxi3, mini0, mini1, arr, arr1, arr2

    # Print elapsed time
    print('No Retreat Cells:'+ str(noretreat[y-min(years)]) +'; Elapsed time:',round(clock()-start,2),'seconds')

# Write to File
ncf = nc.Dataset(outpath+"/C"+CT+"/siphenologyC"+CT+"_"+ver+"_"+str(ymin)+"-"+str(ymax)+".nc",'w', format='NETCDF4')
ncf.createDimension('y', lats.shape[0])
ncf.createDimension('x', lats.shape[1])
ncf.createDimension('time', ymax-ymin+1)

yNC = ncf.createVariable('y', np.float32, ('y',))
xNC = ncf.createVariable('x', np.float32, ('x',))
tNC = ncf.createVariable('time', np.int32, ('time',))

latNC = ncf.createVariable('lat', np.float32, ('y','x',))
lonNC = ncf.createVariable('lon', np.float32, ('y','x',))
areaNC = ncf.createVariable('area', np.float32, ('y','x',))
regNC = ncf.createVariable('reg',np.int16,('y','x',))

minNC = ncf.createVariable('min', np.float32, ('time','y','x',))
mndNC = ncf.createVariable('mnd', np.int16, ('time','y','x',))
lrdNC = ncf.createVariable('lrd', np.int16, ('time','y','x',))
ladNC = ncf.createVariable('lad', np.int16, ('time','y','x',))
frdNC = ncf.createVariable('frd', np.int16, ('time','y','x',))
fadNC = ncf.createVariable('fad', np.int16, ('time','y','x',))
opcNC = ncf.createVariable('opc', np.int16, ('time','y','x',))
opNC = ncf.createVariable('op', np.int16, ('time','y','x',))

ncf.description = '''Phenology of Sea Ice Concentration with 5-Day Moving
Average based on Concentration Threshold of ''' + CT + '%. Note: -99 = N/A.'
ncf.projection = '''Polar stereographic north (EPSG = 3413);
central longitude = 45 deg N, central laitude = 90 deg N,
standard parallel = 70 deg N'''
tNC.units = 'years'
xNC.units = 'm'
yNC.units = 'm'
latNC.units = 'degrees north'
lonNC.units = 'degrees east'
areaNC.units = 'sq. km'
regNC.units = 'none'
minNC.units = 'percentage'
mndNC.units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
lrdNC.units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
frdNC.units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
ladNC.units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
fadNC.units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
opcNC.units = 'days'
opNC.units = 'days'

tNC.long_name = 'year'
latNC.long_name = 'Latitude'
lonNC.long_name = 'Longitude'
areaNC.long_name = 'Grid cell area'
regNC.long_name = '''Region: 0 = Other Ocean, 2 = Sea of Okhotsk, 3 = Bering Sea,
    4 = Hudson Bay, 5 = Gulf of St. Lawrence, 6 = Labrador Sea, 7 = Greenland Sea,
    8 = Barents Sea, 9 = Kara Sea, 10 = Laptev Sea, 11 = East Siberian Sea,
    12 = Chukchi Sea, 13 = Beaufort Sea, 14 = CAA, 15 = CAO, 16 = Baffin Bay, 21 = Land'''
minNC.long_name = 'Minimum daily SIC for calendar year'
mndNC.long_name = 'Day of SIC minimum'
lrdNC.long_name = 'Last Retreat Day'
frdNC.long_name = 'First Retreat Day'
ladNC.long_name = 'Last Advance Day'
fadNC.long_name = 'First Advance Day'
opcNC.long_name = 'Continuous Open Period (FAD - LRD)'
opNC.long_name = 'Open Period (LAD - FRD)'

tNC[:] = np.array(list(range(ymin,ymax+1)))
xNC[:] = (np.arange(-3850,3750,25)+12.5)*1000
yNC[:] = (np.arange(5850,-5350,-25)-12.5)*1000
latNC[:] = lats
lonNC[:] = lons
areaNC[:] = areas
regNC[:] = mask
minNC[:] = np.array(minlist)
mndNC[:] = np.array(mndlist)
lrdNC[:] = np.array(lrdlist)
ladNC[:] = np.array(ladlist)
frdNC[:] = np.array(frdlist)
fadNC[:] = np.array(fadlist)
opcNC[:] = np.array(opclist)
opNC[:] = np.array(oplist)

ncf.close()

print("Complete.")