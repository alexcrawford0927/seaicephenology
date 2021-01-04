'''*********************************************
Authors: Alex Crawford
Date Created: 8 Jul 2020
Date Modified: 12 Jul 2020

Purpose: To create files that fill in the gaps between sea ice data.

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
v = 0

vers = ["OSISAF"] # version of sea ice data
ncvars = ['ice_conc']
ii = [[-15,-11,-9,-7]] # Indices for subsetting the file names

path = "/Volumes/Miranda/SeaIce/"+vers[v]
inpath = path+"/Daily"
outpath = path+"/DailyFiller"

lyb = 1
reftime = [1978,1,0,12,0,0]
years = list(range(1979,2018+1))
mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
days = ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15",\
        "16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"]

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")
# Store correct version-specific variables
i = ii[v]
ncvar = ncvars[v]

for y in years:
    Y = str(y)
    
    # Create output folders
    try:
        os.chdir(outpath+"/"+Y)
    except:
        os.mkdir(outpath+"/"+Y)
        os.chdir(outpath+"/"+Y)
    
    # Load valid files
    files = os.listdir(inpath+"/"+Y)
    files = [f for f in files if (f.startswith('.') == 0)]
    files.sort()
    
    if y > min(years): # Add last known day from prior year
        files0 = os.listdir(inpath+"/"+str(y-1))
        files0 = [f for f in files0 if (f.startswith('.') == 0)]
        files0.sort()
        files = [files0[-1]] + files
    if y < max(years): # Add first known day from subsequent year
        files2 = os.listdir(inpath+"/"+str(y+1))
        files2 = [f for f in files2 if (f.startswith('.') == 0)]
        files2.sort()
        files =  files + [files2[0]]
    
    start = [int(files[0][i[0]:i[1]]),int(files[0][i[1]:i[2]]),int(files[0][i[2]:i[3]]),12,0,0]
    end = [int(files[-1][i[0]:i[1]]),int(files[-1][i[1]:i[2]]),int(files[-1][i[2]:i[3]]),12,0,0]
    
    filePr = files[0]
    
    files2 = files[1:]
    t = md.timeAdd(start,[0,0,1],lyb)
    
    while t != end:
        
        # If this file already exists...
        if str(t[0])+mons[t[1]-1]+days[t[2]-1] in files2[0]:
            filePr = files2[0] # it becomes the new "previous" file
            files2 = files2[1:] # remove it from consideration
            t = md.timeAdd(t,[0,0,1]) # advance by 1 day
        
        # If this file doesn't exist but is in the urrent year
        elif t[0] == y: # Take a weighted average of the "previous" file and the "next" file (a linear interpolation)
            # Identify times
            t0 = [int(filePr[i[0]:i[1]]),int(filePr[i[1]:i[2]]),int(filePr[i[2]:i[3]]),0,0,0]
            t2 = [int(files2[0][i[0]:i[1]]),int(files2[0][i[1]:i[2]]),int(files2[0][i[2]:i[3]]),0,0,0]
            
            # Generate inverse-distance weights (so weight for t0 is based on distance of t2-t) -- must add to 1
            w2 = md.daysBetweenDates(t0,t)/md.daysBetweenDates(t0,t2)
            w0 = md.daysBetweenDates(t,t2)/md.daysBetweenDates(t0,t2)

            try:
                nc0 = nc.Dataset(inpath+'/'+Y+'/'+filePr,'r')
            except:
                nc0 = nc.Dataset(inpath+'/'+str(y-1)+'/'+filePr,'r')
            try:
                nc2 = nc.Dataset(inpath+'/'+Y+'/'+files2[0],'r')
            except:
                nc2 = nc.Dataset(inpath+'/'+str(y+1)+'/'+files2[0],'r')
            
            arr0 = nc0[ncvar][:]
            arr2 = nc2[ncvar][:]
            
            arr = arr0*w0 + arr2*w2
            
            # Save as a netcdf file with the same formatting as inputs
            nc1 = nc.Dataset(files2[0][:i[0]]+str(t[0])+mons[t[1]-1]+days[t[2]-1]+files2[0][i[3]:],'w',format='NETCDF4')
            
            nc1.createDimension('time',1)
            nc1.createDimension('yc',432)
            nc1.createDimension('xc',432)
            
            nctime = nc1.createVariable('time',np.float64,('time'))
            nctime.units = nc0['time'].units
            nctime[:] = md.daysBetweenDates(reftime,t)*86400 # convert time to seconds
            
            ncyc = nc1.createVariable('xc',np.float64,('xc'))
            ncyc.units = nc0['xc'].units
            ncyc[:] = nc0['xc'][:]
            
            ncxc = nc1.createVariable('yc',np.float64,('yc'))
            ncxc.units = nc0['yc'].units
            ncxc[:] = nc0['yc'][:]
            
            ncsic = nc1.createVariable(ncvar,np.int32,('time','yc','xc'))
            ncsic.units = nc0[ncvar].units
            ncsic.scale_factor = 0.01
            ncsic[:] = arr
            
            nc1.close()
            
            print("Saved "+str(t[0])+mons[t[1]-1]+days[t[2]-1])
            # Advance by 1 day
            t = md.timeAdd(t,[0,0,1])
        
        else:
            # Advance by 1 day
            t = md.timeAdd(t,[0,0,1])
            