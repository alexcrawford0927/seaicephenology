"""
Author: Alex Crawford
Date Created: 3 Aug 2020
Date Modified: 14 Mar 2023

Purpose: Creates a pandas data frame of various regionalized sea ice phenology
parameters for each year and model.
"""
'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import numpy as np
import pandas as pd
import os
import netCDF4 as nc

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


'''#####################
 Declare Variables
#####################'''

### Input Variables ###
thresh = '10cm'
regs =    [40,41,42] # [43,44,45,46,47,48,49] #[0]+list(range(2,16+1)) # [40,41,42] #
nanmin = 0
regvar =  'regHB3' # 'regHB2' # 'reg' #  'regHB' # 'regHB3' #
minlat = 45
experiment = 'historical' #'ssp585' #
sspymin = 2015 # irrevelvant for non-ssp experiments

### Path Variables ###
ver = "CMIP6"

path = "/Volumes/Cassandra/CMIP6"
inpath = path+"/SeaIce/ThicknessPhenology/"+experiment
mskpath = path+"/RegionMasks/SeaIce"
areapath = path+"/GridAreas/SeaIce"
outpath = path+"/RegionalStats/ThicknessPhenology/"+regvar

modstoskip = ['AWI-ESM-1-1-LR-PSN_historical_r1i1p1f1_gn', 'AWI-CM-1-1-MR-PSN_historical_r1i1p1f1_gn', 'ICON-ESM-LR-PSN_historical_r1i1p1f1_gn']
# modstokeep = ['MPI-ESM1-2-LR_historical_r'+str(i)+'i1p1f1_gn' for i in range(11,12)] #  ['IPSL-CM6A-LR_historical_r33i1p1f1_gn'] #

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

completed = os.listdir(outpath)

models = os.listdir(inpath+"/"+thresh)
models = [model for model in models if (model.startswith('.') == 0) and (model not in modstoskip)]
# models = [model for model in models if (model in modstokeep)]
models.sort()

for model in models:
    # Establish file names
    outname = ver+"_ThicknessPhenology_"+thresh+"_Regional-"+regvar+"_"+experiment+"_"+model.split("_")[1]+"_"+model.split("_")[5]+".csv"

    if outname not in completed:
        print("Starting "+model)

        # Create empty lists
        nlist, modellist, memberlist, reglist, familylist, ylist = [[] for i in range(6)]
        area, pern, lrd_avg, cip_avg, fad_avg, opc_avg = [[] for i in range(6)]

        # Load phenology & Regions & Areas
        files = os.listdir(inpath+"/"+thresh+"/"+model)
        infile = [file for file in files if file.startswith(model) and 'ValidCells' not in file and '_avg_' not in file][0]

        ncf = nc.Dataset(inpath+"/"+thresh+"/"+model+"/"+infile)
        times = ncf['time'][:].data

        latarr = ncf['lat'][:].data

        regarr = nc.Dataset(mskpath+"/Regions_"+model.split("_")[1]+".nc")[regvar][:].data
        areaarr = nc.Dataset(mskpath+"/Regions_"+model.split("_")[1]+".nc")['area'][:].data

        if regarr.shape != latarr.shape:
            reglat = nc.Dataset(mskpath+"/Regions_"+model.split("_")[1]+".nc")['lat'][:].data
            irows = np.unique(np.where((np.isfinite(reglat) == 1) & (reglat >= minlat))[0])

            if len(regarr.shape) == 2:
                regarr = regarr[irows,:]
                areaarr = areaarr[irows,:]
            else:
                regarr = regarr[irows]
                areaarr = areaarr[irows]

        for y in times: # For each year
            Y = str(y)

            # Read in Sea Ice Files
            lrd = ncf['lrd'][np.where(times == y)[0][0],:].data.astype(float)
            fad = ncf['fad'][np.where(times == y)[0][0],:].data.astype(float)
            cip = ncf['cip'][np.where(times == y)[0][0],:].data.astype(float)

            # Set NaNs
            lrd[lrd < nanmin] = np.nan
            fad[fad < nanmin] = np.nan
            cip[cip < nanmin] = np.nan

            # Identify Finite Values of LRD and FAD and cip
            lrdbool = np.isfinite(lrd)
            fadbool = np.isfinite(fad)
            cipbool = np.isfinite(cip)

            for reg in regs: # For each region
                if len(areaarr.shape) == 2:
                    if reg != 0:
                        rows, cols = np.where(regarr == reg)
                    else:
                        rows, cols = np.where( (regarr >= regs[1]) & (regarr <= regs[-1]) )

                    # Calculate number of cells and area of cells
                    n = np.sum(np.isfinite(cip[rows,cols])) # number of valid cip observations
                    pern.append( np.sum(np.isfinite(lrd[rows,cols]))/n ) # percentage of valid cip obs that are also valid lrd obs
                    area.append( np.nansum(areaarr[rows,cols]) ) # area of valid observations

                    ### Save parameters -- aggregating with area-weighted grid cells ###
                    # The Average LRD
                    lrd_avg.append( np.nansum(lrd[rows,cols]*areaarr[rows,cols]) / np.nansum( lrdbool[rows,cols]*areaarr[rows,cols] ) )

                    # The Average FAD
                    fad_avg.append( np.nansum(fad[rows,cols]*areaarr[rows,cols]) / np.nansum( fadbool[rows,cols]*areaarr[rows,cols] ) )

                    # The average continuous open water period
                    cip_avg.append( np.nansum( cip[rows,cols]*areaarr[rows,cols] ) / np.nansum( cipbool[rows,cols]*areaarr[rows,cols] ) )
                    opc_avg.append( fad_avg[-1] - (lrd_avg[-1]-365) )

                else:
                    if reg != 0:
                        rows = np.where(regarr == reg)[0]
                    else:
                        rows = np.where( (regarr >= regs[1]) & (regarr <= regs[-1]) )[0]

                    # Calculate number of cells and area of cells
                    n = np.sum(np.isfinite(cip[rows])) # number of valid cip observations
                    pern.append( np.sum(np.isfinite(lrd[rows]))/n ) # percentage of valid cip obs that are also valid lrd obs
                    area.append( np.nansum(areaarr[rows]) ) # area of valid observations

                    ### Save parameters -- aggregating with area-weighted grid cells ###
                    # The Average LRD
                    lrd_avg.append( np.nansum(lrd[rows]*areaarr[rows]) / np.nansum( lrdbool[rows]*areaarr[rows] ) )

                    # The Average FAD
                    fad_avg.append( np.nansum(fad[rows]*areaarr[rows]) / np.nansum( fadbool[rows]*areaarr[rows] ) )

                    # The average continuous open water period
                    cip_avg.append( np.nansum( cip[rows]*areaarr[rows] ) / np.nansum( cipbool[rows]*areaarr[rows] ) )
                    opc_avg.append( fad_avg[-1] - (lrd_avg[-1]-365) )

                # Add accessories
                modellist.append( model )
                familylist.append( model.split("_")[1] )
                memberlist.append( int(model.split("_vl_r")[-1].split("i")[0]) )
                reglist.append( reg )
                ylist.append( y )
                nlist.append( n )

        ### Store parameters ###
        pdf = pd.DataFrame({"Model":modellist,"Family":familylist, "Member":memberlist,
                              "Region":reglist,"Year":ylist,"Area":area,
                              "NCells":nlist,"PerCells":pern,"OPCavg":opc_avg,
                              "LRDavg":lrd_avg, "FADavg":fad_avg,"CIPavg":cip_avg})

        pdf['OPCavg'] = np.where( (np.isfinite(np.array(pdf['OPCavg'])) == 0) & (np.array(pdf['CIPavg']) >= 0), 365-np.array(pdf['CIPavg']), np.array(pdf['OPCavg']) )

        if experiment.startswith('ssp'):
            pdf = pdf.loc[pdf['Year'] >= sspymin]

        pdf.to_csv(outpath+"/"+outname,index=False)
