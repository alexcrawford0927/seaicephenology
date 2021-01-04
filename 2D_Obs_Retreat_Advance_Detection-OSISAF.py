'''*********************************************
Authors: Alex Crawford
Date Created: 6/9/15
Date Modified: 1/28/19; 8/20/19 edited for Python 3; 
6/4/20 switched to movingAverage2
6/11/20 simplified the maximum day storage to be an integer, not an array, 
since the array was never saved
7/7/20 adapted for OSISAF method and converted output to netcdf instead of geotiff

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
import MERRA_Module as md

'''*******************************************
Declare Variables
*******************************************'''
print("Declaring Variables")

### Input Variables ###
ct = 0.80 # A number between 0 and 1 for the concentration threshold
n = 5 # Moving Average Size (# observations total, including the central
# value -- should always be an odd number)
ver = "OSISAF"
ncvar = 'ice_conc'
i = [-15, -11, -11, -9, -9, -7] # Indexes for subsetting year, month, day

### Time Variables ###
ymin, ymax = 1979, 2013 # Years of interest
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
maskMin = 0
maskMax = 100

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")
# Set up projection info
ncprj = nc.Dataset(inpath+"/1979/ice_conc_nh_ease2-250_cdr-v2p0_197901021200.nc")
lats = ncprj['lat'][:]
lons = ncprj['lon'][:]

# Time Set Up
years = range(ymin,ymax+1)
mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
days = ["01","02","03","04","05","06","07","08","09","10","11","12","13",\
    "14","15","16","17","18","19","20","21","22","23","24","25","26","27",\
    "28","29","30","31"]
    
noretreat = np.zeros(len(years))

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
        fileyList = [f for f in os.listdir(inpath+"/"+str(y)) if (f.startswith('.') == 0)]
        fileyList_2 = [f for f in os.listdir(inpath+"/"+str(y+1)) if (f.startswith('.') == 0 and int(f[i[2]:i[3]]) == 1 and int(f[i[4]:i[5]]) <= n)]
        fileyList_0 = [f for f in os.listdir(inpath+"/"+str(y-1))  if ( f.startswith('.') == 0 and int(f[i[2]:i[3]]) == 12 and int(f[i[4]:i[5]]) > 31-n)]
        
        # Load Filled Files
        fileyListA = [f for f in os.listdir(inpath2+"/"+str(y)) if (f.startswith('.') == 0)]
        fileyList_2A = [f for f in os.listdir(inpath2+"/"+str(y+1)) if (f.startswith('.') == 0 and int(f[i[2]:i[3]]) == 1 and int(f[i[4]:i[5]]) <= n)]
        fileyList_0A = [f for f in os.listdir(inpath2+"/"+str(y-1))  if ( f.startswith('.') == 0 and int(f[i[2]:i[3]]) == 12 and int(f[i[4]:i[5]]) > 31-n)]
        
        # Combine
        filelist = fileyList_0 + fileyList_0A + fileyList + fileyListA + fileyList_2 + fileyList_2A
        filelist.sort()
        
        # Load Arrays into a single 3-D array
        arrList = []
        for f in filelist:
            # Reclassify array and mask out anything that's not in the regions of the mask
            try:
                nc0 = nc.Dataset(inpath+'/'+f[i[0]:i[1]]+'/'+f,'r')
            except:
                nc0 = nc.Dataset(inpath2+'/'+f[i[0]:i[1]]+'/'+f,'r')
            
            arr0 = nc0[ncvar][0,:,:].data*0.01
            arrList.append(  np.where((arr0 <= maskMax) & (arr0 >= maskMin), arr0, np.nan ) )

        # Smooth array
        if len(fileyList_2+fileyList_2A) == 0:
            ei = len(fileyList_0)+len(fileyList_0A)+len(fileyList)+len(fileyListA)
        else:
            ei = -len(fileyList_2+fileyList_2A)
            
        arr1 = np.apply_along_axis(md.movingAverage2,0,np.array(arrList),n)[len(fileyList_0+fileyList_0A):ei,:,:]
        pd.to_pickle(arr1,mapath+"/SIC_SmoothedMA"+str(n)+"_"+str(y)+".pkl")
    
    # Load Following Year
    try:
        arr2 = pd.read_pickle(mapath+"/SIC_SmoothedMA"+str(n)+"_"+str(y+1)+".pkl")
 
    except: # If not, you must create it anew first
        print("--- Applying moving average... (this may take a minute or two)")
        # Load Annual Files
        fileyList = [f for f in os.listdir(inpath+"/"+str(y+1)) if (f.startswith('.') == 0)]
        fileyList_2 = [f for f in os.listdir(inpath+"/"+str(y+2)) if (f.startswith('.') == 0 and int(f[i[2]:i[3]]) == 1 and int(f[i[4]:i[5]]) <= n)]
        fileyList_0 = [f for f in os.listdir(inpath+"/"+str(y))  if ( f.startswith('.') == 0 and int(f[i[2]:i[3]]) == 12 and int(f[i[4]:i[5]]) > 31-n)]
        
        # Load Filled Files
        fileyListA = [f for f in os.listdir(inpath2+"/"+str(y+1)) if (f.startswith('.') == 0)]
        fileyList_2A = [f for f in os.listdir(inpath2+"/"+str(y+2)) if (f.startswith('.') == 0 and int(f[i[2]:i[3]]) == 1 and int(f[i[4]:i[5]]) <= n)]
        fileyList_0A = [f for f in os.listdir(inpath2+"/"+str(y))  if ( f.startswith('.') == 0 and int(f[i[2]:i[3]]) == 12 and int(f[i[4]:i[5]]) > 31-n)]
        
        # Combine
        filelist = fileyList_0 + fileyList_0A + fileyList + fileyListA + fileyList_2 + fileyList_2A
        filelist.sort()
        
        # Load Arrays into a single 3-D array
        arrList = []
        for f in filelist: # Skip the last file because it's only referenced if no advance occurs
            # Reclassify array and mask out anything that's not in the regions of the mask
            try:
                nc0 = nc.Dataset(inpath+'/'+f[i[0]:i[1]]+'/'+f,'r')
            except:
                nc0 = nc.Dataset(inpath2+'/'+f[i[0]:i[1]]+'/'+f,'r')
            
            arr0 = nc0[ncvar][0,:,:]*0.01
            arrList.append(  np.where((arr0 <= maskMax) & (arr0 >= maskMin) & (arr0 <= 1), arr0, np.nan ) )

        # Smooth array
        if len(fileyList_2+fileyList_2A) == 0:
            ei = len(fileyList_0)+len(fileyList_0A)+len(fileyList)+len(fileyListA)
        else:
            ei = -len(fileyList_2+fileyList_2A)
            
        arr2 = np.apply_along_axis(md.movingAverage2,0,np.array(arrList),n)[len(fileyList_0+fileyList_0A):ei,:,:]
        pd.to_pickle(arr2,mapath+"/SIC_SmoothedMA"+str(n)+"_"+str(y+1)+".pkl")
    
    # Load file list 
    fileList = os.listdir(inpath+"/"+str(y)) + os.listdir(inpath+"/"+str(y+1)) + os.listdir(inpath2+"/"+str(y)) + os.listdir(inpath2+"/"+str(y+1))
    fileList = [f for f in fileList if ( f.startswith('.') == 0  and ( (int(f[i[0]:i[1]]) == y) or (int(f[i[0]:i[1]]) == y+1 ) ) )]
    fileList.sort()
    doys = np.array( [md.daysBetweenDates([y,1,1],[int(f[i[0]:i[1]]),int(f[i[2]:i[3]]),int(f[i[4]:i[5]])]) for f in fileList] )
    
    arr = np.concatenate( (arr1, arr2), axis=0 )
    
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
    mndArr = np.zeros_like(lats)*np.nan
    frdArr, lrdArr = np.zeros_like(lats)*np.nan, np.zeros_like(lats)*np.nan
    fadArr, ladArr = np.zeros_like(lats)*np.nan, np.zeros_like(lats)*np.nan
    opArr, opcArr = np.zeros_like(lats)*np.nan, np.zeros_like(lats)*np.nan
    
    ### Calculate Retreat and Advance Events ###
    validcells = np.where( (np.isnan(Mins) == 0) & (np.isnan(Maxes1) == 0) & (np.isnan(Maxes2) == 0) )
    for j in range(len(validcells[0])):
        r,c = validcells[0][j], validcells[1][j] # Assign row and column
        
        # Calculate index for minimum & maximum
        MaxesI1 = int(np.median(np.where(arr[maxi0:maxi1,r,c] == Maxes1[r,c]))) # Gives first occurrence of minimum if multiples present
        MinsI = mini0 + int(np.median(np.where(arr[mini0:mini1,r,c] == Mins[r,c]))) # Gives first occurrence of minimum if multiples present
        MaxesI2 = maxi2 + int(np.median(np.where(arr[maxi2:maxi3,r,c] == Maxes2[r,c]))) # Gives first occurrence of minimum if multiples present
        
        # Store Minimum Day & and (temporarily) maximum days
        mndArr[r,c] = md.daysBetweenDates(dateref,[int(fileList[MinsI][i[0]:i[1]]),int(fileList[MinsI][i[2]:i[3]]),int(fileList[MinsI][i[4]:i[5]]),0,0,0])
        mxd1rr = md.daysBetweenDates(dateref,[int(fileList[MaxesI1][i[0]:i[1]]),int(fileList[MaxesI1][i[2]:i[3]]),int(fileList[MaxesI1][i[4]:i[5]]),0,0,0])
        mxd2rr = md.daysBetweenDates(dateref,[int(fileList[MaxesI2][i[0]:i[1]]),int(fileList[MaxesI2][i[2]:i[3]]),int(fileList[MaxesI2][i[4]:i[5]]),0,0,0])
        
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
                frdArr[r,c] = md.daysBetweenDates(dateref,[int(fileList[fri][i[0]:i[1]]),int(fileList[fri][i[2]:i[3]]),int(fileList[fri][i[4]:i[5]]),0,0,0])
            except:
                frdArr[r,c] = np.nan
            
            # Last Retreat Day
            # Last index after Maxes1 and before/on Mins for which concentration is below threshold
            try:
                lri = above[np.where((above < MinsI) & (above >= MaxesI1))][-1]
                lrdArr[r,c] = md.daysBetweenDates(dateref,[int(fileList[lri][i[0]:i[1]]),int(fileList[lri][i[2]:i[3]]),int(fileList[lri][i[4]:i[5]])+1,0,0,0])
            except:
                lrdArr[r,c] = np.nan

            # First Advance Day
            # First index after Mins and before/on Maxes2 for which concentration is above threshold
            try:
                fai = above[np.where((above > MinsI) & (above <= MaxesI2))][0]
                fadArr[r,c] = md.daysBetweenDates(dateref,[int(fileList[fai][i[0]:i[1]]),int(fileList[fai][i[2]:i[3]]),int(fileList[fai][i[4]:i[5]]),0,0,0])
            except:
                fadArr[r,c] = np.nan
            
            # Last Advance Day
            # Last index after Mins anbd before/on Maxes2 for which concentration is below threshold
            try:
                lai = below[np.where((below >= MinsI) & (below < MaxesI2))][-1]
                ladArr[r,c] = md.daysBetweenDates(dateref,[int(fileList[lai][i[0]:i[1]]),int(fileList[lai][i[2]:i[3]]),int(fileList[lai][i[4]:i[5]])+1,0,0,0])
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
Average based on Concentration Threshold of ''' + CT + '%. Note: -99 = N/A'
ncf.source = 'netCDF4 python module'
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

tNC[:] = np.array(list(range(ymin,ymax+1)))
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