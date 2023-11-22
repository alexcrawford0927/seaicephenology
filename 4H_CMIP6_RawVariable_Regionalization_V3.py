"""
Author: Alex Crawford
Date Created: 3 Aug 2020
Date Modified: 29 Aug 2022
Purpose: Creates a pandas data frame of various regionalized CMIP6 variables for
each year/month.
"""
'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import numpy as np
import pandas as pd
import os
import xarray as xr
import CycloneModule_12_4 as md

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

'''#####################
 Declare Variables
#####################'''

### Input Variables ###
j = 3
regs = [ [40,41,42] , [43,44,45,46,47,48,49] , [40,41,42] , [0]+list(range(2,16+1)) , [-1] ]
nanmax = 9999
regvar = ['regHB3', 'regHB2' ,'regHB' , 'reg' , 'reg']
ncvar = 'sisnthick' #'tas' # 'sisnthick' #'friver' # 'tos' #     #'zos' #  'uas' # 'sos' # 'vas' # 'uas' #
tres = 'f_mon' #'t_Amon' #
domain = 'SeaIce' #'Atmosphere' # 'SeaIce' # 'Ocean' #

experiment = 'ssp585' # 'historical' # 

### Path Variables ###
ver = "CMIP6"

path = "/Volumes/Cassandra/"+ver+"/" # '/project/6061839/crawfora/'+ver+'/'  # '/project/6061839/crawfora/'  #  '/Users/acrawfora/Downloads/' #
inpath = path+"/data/e_"+experiment+"/v_"+ncvar
mskpath =  path+"/RegionMasks/"+domain # "/Volumes/Cassandra/"+ver+"/RegionMasks/"+domain #
outpath = path+"/RegionalStats/"+ncvar # "/Volumes/Cassandra/"+ver+"/RegionalStats/"+ncvar

modstoskip = ['CAS-ESM2-0'] #
# modstokeep = ['BCC-CSM2-MR','AWI-CM-1-1-MR'] # ['EC-Earth3'] # ['AWI-ESM-1-1-LR'] # ['IPSL-CM6A-LR'] # ['MPI-ESM1-2-LR'] # ['EC-Earth3','IPSL-CM6A-LR','MPI-ESM1-2-LR']  # ['EC-Earth3-Veg','GFDL-CM4'] #
mods1d = ['ICON-ESM-LR','AWI-CM-1-1-MR','AWI-ESM-1-1-LR']

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

if np.max(regs[j]) < 0:
    regvar2 = 'global'
else:
    regvar2 = regvar[j]

models = md.listdir(inpath)
models = [model for model in models if model.split('_')[1] not in modstoskip]
# models = [model for model in models if model.split('_')[1] in modstokeep]
# models = [model for model in models if model.startswith(modstokeep[0])]
models.sort()

for model in models[53:]:
    mod = model.split("s_")[1].split("_e_")[0]

    completedfiles = md.listdir(outpath)
    outname = ver+"_"+ncvar+"_Regional-"+regvar2+"_"+experiment+"_"+mod+"_"+model.split("vl_")[-1].split('_')[0]+".csv"

    if outname not in completedfiles:
        print("Starting "+model)
        # Establish file names
        files = os.listdir(inpath+"/"+model)
        files = [f for f in files if f.startswith('.') == 0 and f.endswith('.nc')]
        files.sort()

        # Load phenology & Regions & Areas
        ncf = xr.open_dataset(inpath+"/"+model+"/"+files[0],decode_times=False)

        # Load lats, areas, & regions
        regnc = xr.open_dataset(mskpath+"/Regions_"+mod+".nc")
        regarr = xr.open_dataset(mskpath+"/Regions_"+mod+".nc")[regvar[j]][:].data
        areaarr = regnc['area'][:].data
        latname = [key for key in ncf.variables.keys() if key.startswith('lat') or key.startswith('nav_lat')][0]
        lats = ncf[latname][:].data
        lats = np.where((lats > 90) | (lats < -90), np.nan, lats)

        if (len(lats.shape) == 1) and (model.split("_")[1] not in mods1d):
            lonname = [key for key in ncf.variables.keys() if key.startswith('lon') or key.startswith('nav_lon')][0]
            lons = ncf[lonname][:].data
            lons = np.where((lons < -180) | (lons > 360), np.nan, lons)
            lons, lats = np.meshgrid( lons,lats )

        # Load times
        times = ncf['time'][:].data
        if len(files) > 1:
            for fi in range(1,len(files)):
                times = np.concatenate( (times, xr.open_dataset(inpath+"/"+model+"/"+files[fi], decode_times=False)['time'][:].data),axis=0)

        # Identify calendar #
        if 'mon' in tres:
            if (times[-1] - times[11])%360 == 0:
                dpy, lys = 360, 0
            elif (times[-1] - times[11])%365 == 0:
                dpy, lys = 365, 0
            else:
                dpy, lys = 365, 1

        else: # for daily or hourly...
            if times.shape[0]%360 == 0:
                dpy, lys = 360, 0
            elif times.shape[0]%365 == 0:
                dpy, lys = 365, 0
            else:
                dpy, lys = 365, 1

        del times


        # Create empty lists
        nlist, modellist, memberlist, familylist, ylist, mlist = [[] for i in range(6)]
        reglist, area, pern, arravg = [[] for i in range(4)]

        for f in files:
            # print(f)
            ncf = xr.open_dataset(inpath+"/"+model+"/"+f, decode_times=False)
            times = ncf['time'][:].data.astype(int)
            refyear = int(ncf['time'].units.split('days since ')[-1][0:4])

            if 'mon' in tres:
                if int(f.split('-')[-2][-2:]) > 1:
                    ymin = int(f.split('-')[-2][-6:-2])+1
                else:
                    ymin = int(f.split('-')[-2][-6:-2])

                ymax = int(f.split('-')[-1][0:4])

            else:
                ymin = int(f.split('-')[-2][-8:-4])
                if int(f.split('-')[-1][4:8]) >= 1200:
                    ymax = int(f.split('-')[-1][0:4])
                else:
                    ymax = int(f.split('-')[-1][0:4])-1

            for y in range(ymin,ymax+1):
                for m in range(12):
                    sday, eday = md.daysBetweenDates([refyear,1,1,0,0,0],[y,m+1,1,0,0,0],lys,dpy), md.daysBetweenDates([refyear,1,1,0,0,0],md.timeAdd([y,m+1,1,0,0,0],[0,1,0,0,0,0],lys,dpy),lys,dpy)
                    ti = np.where((times >= sday) & (times < eday))[0]

                    for reg in regs[j]: # For each region
                        if len(areaarr.shape) == 2:
                            if reg > 0:
                                rows, cols = np.where(regarr == reg)
                            elif reg == 0:
                                rows, cols = np.where( (regarr >= regs[j][1]) & (regarr <= regs[j][-1]) )
                            else:
                                rows, cols = np.where( (regarr >= -1*np.inf) & (regarr <= np.inf) )

                            # Extract data for each year/month
                            arr = ncf[ncvar][ti,:,:].data[:,rows,cols]
                            arr[np.abs(arr) >= nanmax] = np.nan

                            # Calculate number of cells and area of cells
                            n = np.sum(np.isfinite(arr)) # number of valid observations
                            pern.append( n / (rows.shape[0]*ti.shape[0]) ) # percentage of possible observations that are valid
                            area.append( np.nansum(areaarr[rows,cols]) ) # area of valid observations

                            # Save parameters -- aggregating with area-weighted grid cells ###
                            arravg.append( np.nansum( arr*areaarr[rows,cols] ) / np.nansum( np.isfinite(arr)*areaarr[rows,cols] ) )

                        else:
                            if reg > 0:
                                rows = np.where(regarr == reg)[0]
                            elif reg == 0:
                                rows = np.where( (regarr >= regs[j][1]) & (regarr <= regs[j][-1]) )[0]
                            else:
                                rows = np.where( (regarr >= -1*np.inf) & (regarr <= np.inf) )[0]

                            # Extract data for each year/month
                            arr = ncf[ncvar][ti,rows].data
                            arr[np.abs(arr) >= nanmax] = np.nan

                            # Calculate number of cells and area of cells
                            n = np.sum(np.isfinite(arr)) # number of valid observations
                            pern.append( n / (rows.shape[0]*ti.shape[0]) ) # percentage of possible observations that are valid
                            area.append( np.nansum(areaarr[rows]) ) # area of valid observations

                            # Save parameters -- aggregating with area-weighted grid cells ###
                            arravg.append( np.nansum( arr*areaarr[rows] ) / np.nansum( np.isfinite(arr)*areaarr[rows] ) )

                        # Add accessories
                        # modellist.append( model )
                        familylist.append( mod )
                        memberlist.append( int(model.split("vl_r")[-1].split("i")[0]) )
                        reglist.append( reg )
                        ylist.append( y )
                        mlist.append( m+1 )
                        nlist.append( n )

        ### Store parameters ###
        pdf = pd.DataFrame({"Family":familylist, "Member":memberlist,
                              "Region":reglist,"Year":ylist,"Month":mlist,"Area":area,
                              "NCells":nlist,"PerCells":pern,ncvar:arravg})


        pdf.to_csv(outpath+"/"+outname,index=False)
