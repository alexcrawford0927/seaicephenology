'''*********************************************
Authors: Alex Crawford
Date Created: 20 Oct 2020
Date Modified: 20 Oct 2020
Purpose: Combines all observational temperature records into a single 
observational average. Since not all observatuonal records extend back to
1850, the linear regression is run between the full records and partial records
for the overlap periods (1880-2019) to estimate the values for the partial 
records during the missing period (1850-1879).
*********************************************'''

'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")
import numpy as np
import pandas as pd
from scipy import stats

'''*******************************************
Declare Variables
*******************************************'''
print("Declaring Variables")

path = "/Volumes/Troilus/SurfaceTemperature/"

ymin, ymax = 1850, 2019 # Target Range of Years
bmin, bmax = 1880, 1900 # Common Baseline Period
imin, imax = 1850, 1900 # Target Baseline Period

'''*******************************************
Main Analysis
*******************************************'''

# Read in files
noaa = pd.read_csv(path+'/NOAAGlobalTemp/NOAAGlobalTemp_Annual_ts_1880-2019.csv')
best = pd.read_csv(path+'/BEST/BEST_Annual_ts_1850-2019.csv')
hadcrut = pd.read_csv(path+'/HadCRUT/HadCRUT_Annual_ts_1850-2019.csv')
gistemp = pd.read_csv(path+'/GISTEMP/GISTEMP_Annual_ts_1880-2019.csv')

# Initialize data frame
avg = pd.DataFrame()
avg['Year'] = range(ymin,ymax+1)

#### Standardize the Baseline Period for All Datasets ####
noaa['tglobal'] = noaa['tglobal'] - noaa.loc[ (noaa['Year'] >= bmin) & (noaa['Year'] <= bmax) , 'tglobal' ].mean()
noaa['tregion'] = noaa['tregion'] - noaa.loc[ (noaa['Year'] >= bmin) & (noaa['Year'] <= bmax) , 'tregion' ].mean()

best['tglobal'] = best['tglobal'] - best.loc[ (best['Year'] >= bmin) & (best['Year'] <= bmax) , 'tglobal' ].mean()
best['tregion'] = best['tregion'] - best.loc[ (best['Year'] >= bmin) & (best['Year'] <= bmax) , 'tregion' ].mean()

hadcrut['tglobal'] = hadcrut['tglobal'] - hadcrut.loc[ (hadcrut['Year'] >= bmin) & (hadcrut['Year'] <= bmax) , 'tglobal' ].mean()
hadcrut['tregion'] = hadcrut['tregion'] - hadcrut.loc[ (hadcrut['Year'] >= bmin) & (hadcrut['Year'] <= bmax) , 'tregion' ].mean()

gistemp['tglobal'] = gistemp['tglobal'] - gistemp.loc[ (gistemp['Year'] >= bmin) & (gistemp['Year'] <= bmax) , 'tglobal' ].mean()
gistemp['tregion'] = gistemp['tregion'] - gistemp.loc[ (gistemp['Year'] >= bmin) & (gistemp['Year'] <= bmax) , 'tregion' ].mean()

#### Find relationship between Avg(HadCRUT,BEST) & Avg(GISTEMP, NOAA) ####

# Create data frame for NOAA/GISTEMP
noaagis = pd.DataFrame()
noaagis['Year'] = range(ymin,ymax+1)
noaagis['tglobal'] = np.nan
noaagis['tregion'] = np.nan
noaagis.loc[noaagis['Year'] >= bmin,'tglobal'] = np.array((noaa['tglobal'] + gistemp['tglobal'])/2)
noaagis.loc[noaagis['Year'] >= bmin,'tregion'] = np.array((noaa['tregion'] + gistemp['tregion'])/2)

# Create data frame for HadCRUT/BEST
besthad = pd.DataFrame()
besthad['Year'] = range(ymin,ymax+1)
besthad['tglobal'] = ( np.array(best.loc[best['Year'] >= ymin,'tglobal']) + np.array(hadcrut['tglobal'])) / 2
besthad['tregion'] = ( np.array(best.loc[best['Year'] >= ymin,'tregion']) + np.array(hadcrut['tregion'])) / 2

# Calculate statistics
statsglobal = stats.linregress( np.array(besthad.loc[besthad['Year'] >= bmin,'tglobal']) , np.array(noaagis.loc[besthad['Year'] >= bmin,'tglobal']) )
statsregion = stats.linregress( np.array(besthad.loc[besthad['Year'] >= bmin,'tregion']) , np.array(noaagis.loc[besthad['Year'] >= bmin,'tregion']) )

# Fill in missing data for noaagis
noaagis.loc[noaagis['Year'] < bmin,'tglobal'] = np.array(besthad.loc[besthad['Year'] < bmin,'tglobal']) * statsglobal[0] + statsglobal[1]
noaagis.loc[noaagis['Year'] < bmin,'tregion'] = np.array(besthad.loc[besthad['Year'] < bmin,'tregion']) * statsregion[0] + statsregion[1]

#### Create a Final Temperature Record with the same baseline as IPCC ####

# Average Everything!
tavg = pd.DataFrame()
tavg['Year'] = range(ymin,ymax+1)
tavg['tglobal'] = np.array((noaagis['tglobal'] + besthad['tglobal'])/2)
tavg['tregion'] = np.array((noaagis['tregion'] + besthad['tregion'])/2)

# Adjust Baseline Period
tavg['tglobal'] = tavg['tglobal'] - tavg.loc[ (tavg['Year'] >= imin) & (tavg['Year'] <= imax), 'tglobal'].mean()
tavg['tregion'] = tavg['tregion'] - tavg.loc[ (tavg['Year'] >= imin) & (tavg['Year'] <= imax), 'tregion'].mean()

# Save to File
tavg.to_csv(path+"/ObsAvg_Annual_ts_1850-2019.csv",index=False)
