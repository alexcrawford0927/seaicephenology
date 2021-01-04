"""
Author: Alex Crawford
Date Created: 3 Jun 2020
Date Modified: 11 Jun 2020
Purpose: Creates a pandas data frame of various regionalized sea ice phenology 
parameters for each year of an observational dataset. 
"""
'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import numpy as np
import pandas as pd
import netCDF4 as nc

'''#####################
 Declare Variables
#####################'''

### Input Variables ###
cts = [15,30,50,80] # A number between 0 and 100 for the concentration threshold
reglist = [0]+list(range(2,16+1))
ymin, ymax = 1979, 2013
nanmin = 0 # If inputs are lower than this value, convert to NA 
spres = 25 # in km

### Mask Variables ###
lonmin, lonmax = -180, 360
latmin, latmax = -90,90

# Plotting Variables
bb = [35,-45,35,135] # [ll_lat, ll_lon, ur_lat, ur_lon]
lon_0 = 0 # Central Meridian
lat_0 = 90 # Latitude of Origin
levs = list(range(8,21,1)) #[0,100,200,400,800,1600,3200]

### Path Variables ###
ver = "OSISAF" # Observational Dataset

path = "/Volumes/Miranda"
inpath = path+"/SeaIce/"+ver+"/AdvanceRetreat2/"
prjpath = path+"/Projections"
mskpath = path+"/SeaIce/"+ver
outpath = inpath #"/Volumes/Troilus/CMIP6/RegionalStats" 

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

# Load Regional File & Thresholds
regarr = pd.read_pickle(mskpath+"/Regions_"+ver+".pkl")
threshs = pd.read_csv(path+"/SeaIce/Regional_DOY_Thresholds25_C50_V2.csv")

if ver in ["NASATeam","Bootstrap"]:
    areaarr = nc.Dataset(prjpath+'/psn_projection.nc')['area'][:]
    
    for ct in cts:
        # Load phenology
        ncf = nc.Dataset(inpath+"/C"+str(ct)+"/siphenologyC"+str(ct)+"_"+ver+"_"+str(ymin)+"-"+str(ymax)+".nc")
        times = ncf['time'][:].data
        
        # Initiate Data Frame
        pdf = pd.DataFrame()
        
        for y in range(ymin, ymax+1): # For each year
            Y = str(y)
            
            # Read in Sea Ice Files
            lrd = ncf['lrd'][np.where(times == y)[0][0],:,:].data.astype(float)
            fad = ncf['fad'][np.where(times == y)[0][0],:,:].data.astype(float)
            opc = ncf['opc'][np.where(times == y)[0][0],:,:].data.astype(float)
            
            # Set NaNs
            lrd[lrd < nanmin] = np.nan
            fad[fad < nanmin] = np.nan
            opc[opc < nanmin] = np.nan
            
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
                
                ### Store parameters ###
                prow = pd.DataFrame([{"Model":ver,"Region":reg,"Year":y,"Area":area,
                                      "NCells":n,"OPCavg":opc_avg,"OPCgt":opc_gt,
                                      "LRDlt":lrd_lt, "FADgt":fad_gt},])
                
                pdf = pdf.append(prow,ignore_index=True,sort=False)
        
        pdf.to_csv(outpath+"/C"+str(ct)+"/"+ver+"_Regionalized_Historical_C"+str(ct)+"_V6.csv",index=False)

else:
    for ct in cts:
        # Load phenology
        ncf = nc.Dataset(inpath+"/C"+str(ct)+"/siphenologyC"+str(ct)+"_"+ver+"_"+str(ymin)+"-"+str(ymax)+".nc")
        times = ncf['time'][:].data
        
        # Initiate Data Frame
        pdf = pd.DataFrame()
        
        for y in range(ymin, ymax+1): # For each year
            Y = str(y)
            
            # Read in Sea Ice Files
            lrd = ncf['lrd'][np.where(times == y)[0][0],:,:].data.astype(float)
            fad = ncf['fad'][np.where(times == y)[0][0],:,:].data.astype(float)
            opc = ncf['opc'][np.where(times == y)[0][0],:,:].data.astype(float)
            
            # Set NaNs
            lrd[lrd < nanmin] = np.nan
            fad[fad < nanmin] = np.nan
            opc[opc < nanmin] = np.nan
            
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
                
                ### Save parameters -- aggregating with area-weighted grid cells ###
                # The % of grid cells experiencing retreat before X OR open conditions year-round
                lrd_bool = (lrd[rows,cols] < lrdth) | ( (opc[rows,cols] >= 365) & (np.isnan(lrd[rows,cols])))
                lrd_lt = np.nansum( lrd_bool ) / n * 100
                        
                # The % of grid cells experiencing advance after X OR open conditions year-round
                fad_bool = (fad[rows,cols] >  fadth) | ( (opc[rows,cols] >= 365) & (np.isnan(fad[rows,cols])))
                fad_gt = np.nansum( fad_bool ) / n * 100
        
                # The % of grid cells experiencing at least X weeks of continuous open water
                opc_gt = np.nansum( (opc[rows,cols] > opcth)) / n * 100
                
                # The average continuous open water period
                opc_avg = np.nansum( opc[rows,cols] ) / n

                ### Store parameters ###
                prow = pd.DataFrame([{"Model":ver,"Region":reg,"Year":y,"Area":n*spres*spres,
                                      "NCells":n,"OPCavg":opc_avg,"OPCgt":opc_gt,
                                      "LRDlt":lrd_lt, "FADgt":fad_gt},])
                
                pdf = pdf.append(prow,ignore_index=True,sort=False)
        
        pdf.to_csv(outpath+"/C"+str(ct)+"/"+ver+"_Regionalized_Historical_C"+str(ct)+"_V6.csv",index=False)
