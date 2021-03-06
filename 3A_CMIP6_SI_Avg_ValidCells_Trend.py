'''*********************************************
Authors: Alex Crawford
Date Created: 18 May 2020
Date Modified: 21 May 2020
Purpose: To calculate the sea ice phenology averages, trends, and number 
of valid years for each grid cell of the input data.

Inputs: 
    concThresh = concentration threshold - value between 0 and 100
    variables = variables of interest - a list of strings
    minvy = minimum number of valid years for trend analysis
    ymin, ymax = years of interest (inclusive)

Outputs: A netcdf file of the average and trends in various sea ice
    phenology variables. Also included are the number of valid years for the
    averaging and trend time periods. Finally, a CSV file is exported with the 
    annual total number of valid cells (experiencing both retreat and advance).
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
from scipy import stats
import pandas as pd

'''*******************************************
Declare Variables
*******************************************'''
print("Declaring Variables")

### Input Variables ###
concThresh = 15 # A number between 0 and 100 for the concentration threshold
experiment = 'historical'

variables = ['lrd','fad','frd','lad','op','opc']

### Time Variables ###
ymin, ymax = 1981, 2010

minvy = 10 # minimum number of valid years

### Path Variables ###
path = "/Volumes/Troilus"
inpath = path+"/CMIP6/AdvanceRetreat/"+experiment+"/C"+str(concThresh)

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

# Identify models that have not been processed yet
models = os.listdir(inpath)
models = [f for f in models if f.startswith('.') == 0]
models.sort()
newmodels = [f for f in models if sum([ff[14:-12]+str(ymin)+'-'+str(ymax)+'.nc' in ff for ff in os.listdir(inpath+"/"+f)]) == 0]

for mod in newmodels[:]:
    
    ###### PREP STEPS ######
    # Load files
    os.chdir(inpath+"/"+mod)
    files = os.listdir()
    files = [f for f in files if f.startswith('siphenologyC'+str(concThresh)) == 1]
    
    ncf = nc.Dataset(files[0])
    
    # Load Times
    nctime = ncf['time'][:].data
    
    # Initiate New NetCDF Output File
    ncanom = nc.Dataset(files[0][:14]+'_avg_and_trend_N'+str(minvy)+files[0][14:-12]+str(ymin)+'-'+str(ymax)+'.nc', 'w', format='NETCDF4')
    ncanom.description = '''Anomalies & Trends of Phenology of Sea Ice Concentration 
    with 5-Day Moving Average based on Concentration Threshold of ''' + str(concThresh) + '''%
    with respect to ''' + str(ymin) + ' to ' + str(ymax)
    ncanom.source = 'netCDF4 python module'
    
    # Create Dimensions
    ncanom.createDimension('y', ncf['y'].shape[0])
    ncanom.createDimension('x', ncf['x'].shape[0])
    
    yNC = ncanom.createVariable('y', np.float32, ('y',))
    xNC = ncanom.createVariable('x', np.float32, ('x',))
    
    latNC = ncanom.createVariable('lat', np.float32, ('y','x',))
    latNC.units = 'degrees N'
    latNC[:] = ncf.variables['lat'][:]
    
    lonNC = ncanom.createVariable('lon', np.float32, ('y','x',))
    lonNC.units = 'degrees E'
    lonNC[:] = ncf.variables['lon'][:]

    ###### NUMBER OF VALID YEARS #####
    nvy = np.isfinite(ncf['fad'][:].data-ncf['lrd'][:].data) # Valid Gridcells
    nvyts = [np.sum(nvy[i,:,:]) for i in range(nvy.shape[0])]
    pdf = pd.DataFrame(columns=("Year","ValidCells"))
    pdf['Year'] = pd.Series(ncf['time'][:])
    pdf['ValidCells'] = pd.Series(nvyts[:])
    pdf.to_csv('ValidCells_C'+str(concThresh)+files[0][14:]+'.csv',index=False)
    
    # Number of Valid Years for Averages
    itime = np.in1d(nctime,np.arange(ymin,ymax+1))
    nvyNC = ncanom.createVariable('validyearcount_avg',np.int16,('y','x',))
    nvyNC.units = 'number of years from ' + str(ymin) + ' to ' + str(ymax) +'''
        with both sea ice retreat and advance recorded'''
    nvyNC[:] = np.apply_along_axis(np.sum,0,nvy[itime,:,:])
    
    # Number of Valid Years for Trends
    ittime = np.in1d(nctime,np.arange(ymin,ymax+1))
    nty = np.apply_along_axis(np.sum,0,nvy[ittime,:,:])
    irows, icols = np.where(nty >= minvy)
    
    # For Each Variable
    for v in variables:
        
        ###### AVG #####
        val = ncf[v][:].data # Extract Value
        avg = np.apply_along_axis(np.nanmean,0,val[itime,:,:]) # Calculate average
        
        # Assign to New NetCDF File
        vavgNC = ncanom.createVariable(v+'_avg',np.float,('y','x',))
        vavgNC.units = 'average ' + v + ' (in days) for period ' + str(ymin) + '-' + str(ymax)
        vavgNC[:] = avg
        
        ###### TRENDS #####
        b, a, r, p, e = np.zeros_like(nty)*np.nan, np.zeros_like(nty)*np.nan, np.zeros_like(nty)*np.nan, np.ones_like(nty)*1.0, np.zeros_like(nty)*np.nan
        
        for i in range(irows.shape[0]):
            ri, ci = irows[i], icols[i]
            
            tmask = np.isfinite(val[ittime,ri,ci])
            yvals = val[ittime,ri,ci][tmask]
            xvals = np.arange(ymax-ymin+1)[tmask]
            b[ri,ci], a[ri,ci], r[ri,ci], p[ri,ci], e[ri,ci] = stats.linregress(xvals,yvals)

        # Assign to New NetCDF File
        vbNC = ncanom.createVariable(v+'_trend',np.float32,('y','x',))
        vbNC.units = 'trend in ' + v + ' (in days/yr) with respect to period ' + str(ymin) + '-' + str(ymax)
        vbNC[:] = b    
        
        # Assign to New NetCDF File
        vaNC = ncanom.createVariable(v+'_yint',np.float32,('y','x',))
        vaNC.units = 'y-intercept in ' + v + ' (in days) with respect to period ' + str(ymin) + '-' + str(ymax)
        vaNC[:] = a
        
        # Assign to New NetCDF File
        vr2NC = ncanom.createVariable(v+'_rsq',np.float32,('y','x',))
        vr2NC.units = 'r-squared for trend in ' + v + ' with respect to period ' + str(ymin) + '-' + str(ymax)
        vr2NC[:] = np.square(r)
        
        # Assign to New NetCDF File
        vpNC = ncanom.createVariable(v+'_pvalue',np.float32,('y','x',))
        vpNC.units = 'p-value for trend in ' + v + ' with respect to period ' + str(ymin) + '-' + str(ymax)
        vpNC[:] = p
        
        # Assign to New NetCDF File
        veNC = ncanom.createVariable(v+'_stderr',np.float32,('y','x',))
        veNC.units = 'standard error for trend in ' + v + ' with respect to period ' + str(ymin) + '-' + str(ymax)
        veNC[:] = np.square(r)
        
    # Write New NetCDF File
    ncanom.close()
    
    now = clock()
    print(" -- Completed " + mod + ", Elapsed Time: " + str(now-start))
    start = clock()