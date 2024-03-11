'''*********************************************
Authors: Alex Crawford
Date Created: 6/9/15
Date Modified: 16 May 2022 --> Adapted from Satellite Grids


Purpose: To calculate the day of the year on which sea ice retreats and
advances in each grid cell of each year.

Inputs:
    ct = concentration threshold - value between 0 and 1
    n = the moving average size (e.g., 5)
    maxmo, minmo = months in which to search for maximum and minimum SIC values
    ncvar = the variable name for sea ice concentration
    i = the indeices i file names that are necessary for identifying the date

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
ver = "PIOMAS"

### Time Variables ###
ymin, ymax = 2020, 2021 # Years of interest
maxmo = [1,4] # months in which the sea ice maximum may occur
minmo = [8,10] # months in which the sea ice minimum may occur
lys, dpy = 0, 365

### Path Variables ###
path = "/Volumes/Miranda/SeaIce/"+ver
inpath = path+"/SIC_Daily"
mapath = path+"/SIC_Daily_SmoothedMA"+str(n)
outpath = path+"/AdvanceRetreat2"
suppath = "/Volumes/Miranda/Projections/PIOMAS_projection.nc"


#### File Variables ###
valMin = 0.00005
valMax = 1

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")
# Set up projection info
pnc = nc.Dataset(suppath)
lats, lons = pnc['lat'][:].data, pnc['lon'][:].data

# Time Set Up
years = range(ymin,ymax+1)
noretreat = np.zeros(len(years))

# Prep outputs
fadlist, frdlist, ladlist, lrdlist, opclist, oplist, mndlist, minlist = [], [], [], [], [], [], [], []

for y in years:
    print(" Year: "+str(y))

    ### Prep Inputs ###
    # Load Current Year
    try: # If a smoothed SIC array already exists, open it
        arr1 = pd.read_pickle(mapath+"/SIC_SmoothedMA"+str(n)+"_"+str(y)+".pkl")

    except: # If not, you must create it anew first

        # Load current year
        sic = np.fromfile(inpath+"/aiday.H"+str(y),dtype='float32').reshape(365,120,360)
        heff = np.fromfile(path+"/Thickness_Daily/hiday.H"+str(y),dtype='float32').reshape(365,120,360)

        # Load prior and subsequent years, if they exist, and concatenate to current year
        try:
            sic0 = np.fromfile(inpath+"/aiday.H"+str(y-1),dtype='float32').reshape(365,120,360)[-(n//2):,:,:]
            heff0 = np.fromfile(path+"/Thickness_Daily/hiday.H"+str(y-1),dtype='float32').reshape(365,120,360)[-(n//2):,:,:]

            sic = np.concatenate((sic0,sic),axis=0)
            heff = np.concatenate((heff0,heff),axis=0)

            # Note index of Jan 1 of current year
            si = n//2

        except:
            si = 0

            print("--- No prior year")

        try:
            sic2 = np.fromfile(inpath+"/aiday.H"+str(y+1),dtype='float32').reshape(365,120,360)[:n//2,:,:]
            heff2 = np.fromfile(path+"/Thickness_Daily/hiday.H"+str(y+1),dtype='float32').reshape(365,120,360)[:n//2,:,:]

            sic = np.concatenate((sic,sic2),axis=0)
            heff = np.concatenate((heff,heff2),axis=0)

            # Note index of Jan 1 of subsequent year
            ei = -(n//2)

        except:
            ei = sic.shape[0]

            print("--- No subsequent year")

        # Keep all data between 0 and 1
        sic[sic < valMin] = 0
        sic[sic > valMax] = 1

        # Mask out land cells (for some reason, sic = 1 for land, but heff = 0)
        sic[(heff == 0) & (sic == 1)] = np.nan

        # Smooth array
        arr1 = np.apply_along_axis(md.movingAverage,0,sic,n)[si:ei,:,:]
        pd.to_pickle(arr1,mapath+"/SIC_SmoothedMA"+str(n)+"_"+str(y)+".pkl")

    # Load Following Year
    try: # If a smoothed SIC array already exists, open it
       arr2 = pd.read_pickle(mapath+"/SIC_SmoothedMA"+str(n)+"_"+str(y+1)+".pkl")

    except: # If not, you must create it anew first

        # Load following year
        sic = np.fromfile(inpath+"/aiday.H"+str(y+1),dtype='float32').reshape(365,120,360)
        heff = np.fromfile(path+"/Thickness_Daily/hiday.H"+str(y+1),dtype='float32').reshape(365,120,360)

        # Load prior and subsequent years, if they exist, and concatenate to current year
        try:
            sic0 = np.fromfile(inpath+"/aiday.H"+str(y),dtype='float32').reshape(365,120,360)[-(n//2):,:,:]
            heff0 = np.fromfile(path+"/Thickness_Daily/hiday.H"+str(y),dtype='float32').reshape(365,120,360)[-(n//2):,:,:]

            sic = np.concatenate((sic0,sic),axis=0)
            heff = np.concatenate((heff0,heff),axis=0)

            # Note index of Jan 1 of current year
            si = n//2

        except:
            si = 0

            print("--- No prior year")

        try:
            sic2 = np.fromfile(inpath+"/aiday.H"+str(y+2),dtype='float32').reshape(365,120,360)[:n//2,:,:]
            heff2 = np.fromfile(path+"/Thickness_Daily/hiday.H"+str(y+2),dtype='float32').reshape(365,120,360)[:n//2,:,:]

            sic = np.concatenate((sic,sic2),axis=0)
            heff = np.concatenate((heff,heff2),axis=0)

            # Note index of Jan 1 of subsequent year
            ei = -(n//2)

        except:
            ei = sic.shape[0]

            print("--- No subsequent year")

        # Keep all data between 0 and 1
        sic[sic < valMin] = 0
        sic[sic > valMax] = 1

        # Mask out land cells (for some reason, sic = 1 for land, but heff = 0)
        sic[(heff == 0) & (sic == 1)] = np.nan

        # Smooth array
        arr2 = np.apply_along_axis(md.movingAverage,0,sic,n)[si:ei,:,:]
        pd.to_pickle(arr2,mapath+"/SIC_SmoothedMA"+str(n)+"_"+str(y+1)+".pkl")

    # Combine current and following year
    arr = np.concatenate( (arr1, arr2), axis=0 )

    # Record days used by PIOMAS (if the first year is a leap year, skip December 31 -- i.e., doy 366)
    if md.leapyearBoolean([y])[0] == 0:
        doys = np.arange(720)+1
    else:
        doys = np.array(list(range(1,366))+list(range(367,722)))

    # Calculate minimum & maximum value for year
    maxi0 = np.where( doys >= md.daysBetweenDates([y,1,1],[y,maxmo[0],1], lys, dpy) )[0][0]
    maxi1 = np.where( doys <= md.daysBetweenDates([y,1,1],[y,maxmo[1]+1,1], lys, dpy)  )[0][-1]

    mini0 = np.where( doys >= md.daysBetweenDates([y,1,1],[y,minmo[0],1], lys, dpy) )[0][0]
    mini1 = np.where( doys <= md.daysBetweenDates([y,1,1],[y,minmo[1]+1,1], lys, dpy) )[0][-1]

    maxi2 = np.where( doys >= md.daysBetweenDates([y,1,1],[y+1,maxmo[0],1], lys, dpy) )[0][0]
    maxi3 = np.where( doys <= md.daysBetweenDates([y,1,1],[y+1,maxmo[1]+1,1], lys, dpy) )[0][-1]

    Maxes1 = np.amax(arr[maxi0:maxi1],0)
    Mins = np.amin(arr[mini0:mini1],0)
    Maxes2 = np.amax(arr[maxi2:maxi3],0)

    print("--- Calculating Retreat/Advance...")
    CT = str(int(ct*100))
    ### Prep Outputs ###
    mndArr = np.zeros_like(lats)*np.nan
    frdArr, lrdArr = np.zeros_like(lats)*np.nan, np.zeros_like(lats)*np.nan
    fadArr, ladArr = np.zeros_like(lats)*np.nan, np.zeros_like(lats)*np.nan
    opArr, opcArr = np.zeros_like(lats)*np.nan, np.zeros_like(lats)*np.nan

    ### Calculate Retreat and Advance Events ###
    validcells = np.where( (np.isnan(Mins) == 0) & (np.isnan(Maxes1) == 0) & (np.isnan(Maxes2) == 0) )
    for j in range(len(validcells[0])):
        r,c = validcells[0][j], validcells[1][j] # Assign row and column

        # Calculate index for minimum & maximum
        MaxesI1 = int(np.median(np.where(arr[maxi0:maxi1,r,c] == Maxes1[r,c]))) # Gives first occurrence of maximum if multiples present
        MinsI = mini0 + int(np.median(np.where(arr[mini0:mini1,r,c] == Mins[r,c]))) # Gives first occurrence of minimum if multiples present
        MaxesI2 = maxi2 + int(np.median(np.where(arr[maxi2:maxi3,r,c] == Maxes2[r,c]))) # Gives first occurrence of maximum if multiples present

        # Store Minimum Day & and (temporarily) maximum days
        mndArr[r,c] = doys[MinsI]
        mxd1rr = doys[MaxesI1]
        mxd2rr = doys[MaxesI2]

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
                frdArr[r,c] = doys[fri]
            except:
                frdArr[r,c] = np.nan

            # Last Retreat Day
            # Last index after Maxes1 and before/on Mins for which concentration is below threshold
            try:
                lri = above[np.where((above < MinsI) & (above >= MaxesI1))][-1]
                lrdArr[r,c] = doys[lri]
            except:
                lrdArr[r,c] = np.nan

            # First Advance Day
            # First index after Mins and before/on Maxes2 for which concentration is above threshold
            try:
                fai = above[np.where((above > MinsI) & (above <= MaxesI2))][0]
                fadArr[r,c] = doys[fai]
            except:
                fadArr[r,c] = np.nan

            # Last Advance Day
            # Last index after Mins and before/on Maxes2 for which concentration is below threshold
            try:
                lai = below[np.where((below >= MinsI) & (below < MaxesI2))][-1]
                ladArr[r,c] = doys[lai]
            except:
                ladArr[r,c] = np.nan

            # Open Water Periods
            if (Maxes1[r,c] < ct)  & (Maxes2[r,c] >= ct): # When it starts below threshold but ends above
                opArr[r,c] =  ladArr[r,c] - mxd1rr #np.min([365, ladArr[r,c] - MaxesI1])
                opcArr[r,c] = fadArr[r,c] - mxd1rr #np.min([365, fadArr[r,c] - MaxesI1])

                lrdArr[r,c] = np.nan
                frdArr[r,c] = np.nan

            elif (Maxes1[r,c] >= ct)  & (Maxes2[r,c] < ct): # When it starts above threshold but ends below
                opArr[r,c] =  mxd2rr - frdArr[r,c] #np.min([365, MaxesI2 - frdArr[r,c]])
                opcArr[r,c] = mxd2rr - lrdArr[r,c] #np.min([365, MaxesI2 - lrdArr[r,c]])

                fadArr[r,c] = np.nan
                ladArr[r,c] = np.nan

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
try:
    os.chdir(outpath+"/C"+CT)
except:
    os.mkdir(outpath+"/C"+CT)
    os.chdir(outpath+"/C"+CT)

ncf = nc.Dataset(outpath+"/C"+CT+"/siphenologyC"+CT+"_"+ver+"_"+str(ymin)+"-"+str(ymax)+".nc",'w', format='NETCDF4')
ncf.createDimension('y', lats.shape[0])
ncf.createDimension('x', lats.shape[1])
ncf.createDimension('time', ymax-ymin+1)

yNC = ncf.createVariable('y', np.float32, ('y',))
xNC = ncf.createVariable('x', np.float32, ('x',))
tNC = ncf.createVariable('time', np.int32, ('time',))

latNC = ncf.createVariable('lat', np.float32, ('y','x',))
lonNC = ncf.createVariable('lon', np.float32, ('y','x',))

minNC = ncf.createVariable('min', np.float32, ('time','y','x',))
mndNC = ncf.createVariable('mnd', np.int16, ('time','y','x',))
lrdNC = ncf.createVariable('lrd', np.int16, ('time','y','x',))
ladNC = ncf.createVariable('lad', np.int16, ('time','y','x',))
frdNC = ncf.createVariable('frd', np.int16, ('time','y','x',))
fadNC = ncf.createVariable('fad', np.int16, ('time','y','x',))
opcNC = ncf.createVariable('opc', np.int16, ('time','y','x',))
opNC = ncf.createVariable('op', np.int16, ('time','y','x',))

ncf.description = '''Phenology of Sea Ice Concentration with 5-Day Moving
Average based on Concentration Threshold of ''' + CT + '%. Note: Negative Values = N/A'
ncf.source = 'netCDF4 python module'
tNC.units = 'years'
xNC.units = 'm'
yNC.units = 'm'
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

tNC[:] = np.array(list(range(ymin,ymax+1)))
xNC[:] = np.arange(lats.shape[1])
yNC[:] = np.arange(lats.shape[0])
latNC[:] = lats
lonNC[:] = lons
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