"""
Author: Alex Crawford
Date Created: 3 Aug 2020
Date Modified: 3 Aug 2020
Purpose: Creates a pandas data frame of various regionalized sea ice phenology 
parameters for each year and model. 
"""
'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import numpy as np
import pandas as pd
# from osgeo import gdal, gdalnumeric
import os
import netCDF4 as nc

np.seterr(all='ignore')

'''#####################
 Declare Variables
#####################'''

### Input Variables ###
cts = [15,30,50,80] # A number between 0 and 100 for the concentration threshold
reglist = [0]+list(range(2,16+1))
ymin, ymax = 2015, 2099
nanmin = 0
minN = 10

experiment='ssp126'

# Plotting Variables
bb = [35,-45,35,135] # [ll_lat, ll_lon, ur_lat, ur_lon]
lon_0 = 0 # Central Meridian
lat_0 = 90 # Latitude of Origin
levs = list(range(8,21,1)) #[0,100,200,400,800,1600,3200]

### Path Variables ###
ver = "CMIP6"

path = "/Volumes/Troilus/CMIP6"
inpath = "/Volumes/Prospero/CMIP6/AdvanceRetreat/"+experiment
mskpath = path+"/RegionMasks"
areapath = path+"/GridAreas/SeaIce"
outpath = path+"/RegionalStats" 

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

# Load Regional File & Thresholds
threshs = pd.read_csv(path+"/Regional_DOY_Thresholds25_C50_V2.csv")

for ct in cts:
    models = os.listdir(inpath+"/C"+str(ct))
    models = [model for model in models if model.startswith('.') == 0]
    models.sort()
    
    # Initiate Data Frame
    pdf = pd.DataFrame()
    
    for model in models[:]:
        print("Starting "+model)
        # Establish file names
        files = os.listdir(inpath+"/C"+str(ct)+"/"+model)
        infile = [file for file in files if file.startswith("siphenologyC"+str(ct)+"_"+model)][0]
        #amlyfile = [file for file in files if file.startswith("siphenologyC"+str(ct)+"_anomalies_and_trendsN"+str(minN)+"_"+model)][0]
        
        # Load phenology & Regions & Areas
        regarr = pd.read_pickle(mskpath+"/Regions_"+model.split("_"+experiment)[0]+".pkl")
        areaarr = pd.read_pickle(areapath+"/Areas_"+model.split("_"+experiment)[0]+".pkl")
        ncf = nc.Dataset(inpath+"/C"+str(ct)+"/"+model+"/"+infile)
        #ncfa = nc.Dataset(inpath+"/C"+str(ct)+"/"+model+"/"+amlyfile)
        times = ncf['time'][:].data
        
        for y in range(ymin, ymax+1): # For each year
            Y = str(y)
            
            # Read in Sea Ice Files
            lrd = ncf['lrd'][np.where(times == y)[0][0],:,:].data.astype(float)
            fad = ncf['fad'][np.where(times == y)[0][0],:,:].data.astype(float)
            opc = ncf['opc'][np.where(times == y)[0][0],:,:].data.astype(float)
            #opca = ncfa['opc_amly'][np.where(times == y)[0][0],:,:].data.astype(float)
            
            # Set NaNs
            lrd[lrd < nanmin] = np.nan
            fad[fad < nanmin] = np.nan
            opc[opc < nanmin] = np.nan
            #opca[np.isnan(opc)] = np.nan
            
            for reg in reglist: # For each region
                if reg != 0:
                    rows, cols = np.where(regarr == reg)
                else:
                    rows, cols = np.where( (regarr >= reglist[1]) & (regarr <= reglist[-1]) )
                
                lrdth = int(threshs[threshs['Region'] == reg]['LRDThresh'])
                fadth = int(threshs[threshs['Region'] == reg]['FADThresh'])
                opcth = int(threshs[threshs['Region'] == reg]['OPCThresh'])
                
                # Use Continuous Open Water Period as benchmark
                n = np.sum(np.isfinite(opc[rows,cols])) # number of valid observations
                area = np.nansum(areaarr[rows,cols]) # area of valid observations
                
                ### Save parameters -- aggregating with area-weighted grid cells ###
                # The % of grid cells experiencing retreat before X OR open conditions year-round
                lrd_bool = (lrd[rows,cols] < lrdth) | ( (opc[rows,cols] >= 365) & (np.isnan(lrd[rows,cols])))
                lrd_lt = np.nansum( lrd_bool*areaarr[rows,cols] ) / area * 100
                        
                # The % of grid cells experiencing advance after X OR open conditions year-round
                fad_bool = (fad[rows,cols] >  fadth) | ( (opc[rows,cols] >= 365) & (np.isnan(fad[rows,cols])))
                fad_gt = np.nansum( fad_bool*areaarr[rows,cols] ) / area * 100
        
                # The % of grid cells experiencing at least X weeks of continuous open water
                opc_gt = np.nansum( (opc[rows,cols] > opcth)*areaarr[rows,cols] ) / area * 100
                
                # The average continuous open water period
                opc_avg = np.nansum( opc[rows,cols]*areaarr[rows,cols] ) / area
                
                # The average anomaly of the continuous open water period
                #opca_avg = np.nansum( opca[rows,cols]*areaarr[rows,cols] ) / area
                
                ### Store parameters ###
                prow = pd.DataFrame([{"Model":model,"Family":model.split("_"+experiment)[0],
                                      "Member":model.split("_"+experiment+"_r")[-1].split("i")[0],
                                      "Region":reg,"Year":y,"Area":area,
                                      "NCells":n,"OPCavg":opc_avg,"OPCgt":opc_gt,
                                      "LRDlt":lrd_lt, "FADgt":fad_gt},])
                
                pdf = pdf.append(prow,ignore_index=True,sort=False)

    pdf.to_csv(outpath+"/"+ver+"_Regionalized_"+experiment+"_C"+str(ct)+"_V6.csv",index=False)
