"""
Date Created: 4 Apr 2022
Date Modified: 17 Mar 2023
Author: Alex Crawford

Purpose: Calculates the regional average, trend, temperature sensitivity, and error for
sea ice thickness phenology variables. Error is taken as the maximum internal variability 
from a set of single-model ensembles.
 
--> Modified in Mar 2023 to work with historical + ssp instead of one experiment at a time
"""

'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import numpy as np
import scipy.stats as stats
import pandas as pd
import CycloneModule_12_4 as md

'''#####################
 Declare Variables
#####################'''

### Input Variables ###
thres = '10cm' # threshold
minN = 15
regvar = 'reg'
V = 'V9'

aymin, aymax = 1979, 2021

experiment1='historical'
experiment2='ssp585'

varlist = ['OPCavg','CIPavg','LRDavg','FADavg']

### Path Variables ###
ver = "CMIP6"

path = "/Volumes/Cassandra/CMIP6"
inpath = path+"/RegionalStats/ThicknessPhenology/"+regvar
obpath = '/Volumes/Miranda/SeaIce'
obpath2 = '/Volumes/Prospero/'
tobpath = '/Volumes/Theseus/SurfaceTemperature/BEST/BEST_Annual_ts_1850-2022.csv'
# tobpath2 = '/Volumes/Theseus/NCEP-NCAR/RegionalStats/tas/NCEP-NCAR_tasRegionalized_global_1979-2014.csv'
outpath = path+"/RegionalStats/"+V

# Models To Keep
# modstokeep = ["ACCESS-CM2","AWI-CM-1-1-MR","BCC-CSM2-MR",
#              "CESM2","CESM2-FV2","CESM2-WACCM",
#               'CMCC-CM2-SR5','CMCC-ESM2',"CNRM-CM6-1","CNRM-ESM2-1",
#               "EC-Earth3-CC","IPSL-CM6A-LR","KIOST-ESM","MIROC-ES2L","MIROC6",
#               "MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM"]

modsforerror = ["MPI-ESM1-2-LR"]

# modstokeep = ["ACCESS-CM2","AWI-ESM-1-1-LR","AWI-CM-1-1-MR","BCC-CSM2-MR","BCC-ESM1",
#               "CanESM5","CESM2","CESM2-FV2","CESM2-WACCM","CESM2-WACCM-FV2",
#               'CMCC-CM2-HR4','CMCC-CM2-SR5','CMCC-ESM2',"CNRM-CM6-1","CNRM-CM6-1-HR","CNRM-ESM2-1","EC-Earth3",
#               "EC-Earth3-AerChem","EC-Earth3-CC","EC-Earth3-Veg","EC-Earth3-Veg-LR","GFDL-CM4","ICON-ESM-LR",
#               "IPSL-CM5A2-INCA","IPSL-CM6A-LR-INCA","IPSL-CM6A-LR","KIOST-ESM","MIROC-ES2L","MIROC6","MPI-ESM1-2-HAM",
#               "MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM","SAM0-UNICON","UKESM1-0-LL"]

# modsforerror = ["EC-Earth3","IPSL-CM6A-LR","MPI-ESM1-2-LR"]

obsmods = ['PIOMAS'] # 'ERA5'
obspaths = [obpath] # obpath2

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

######################
# STEP 0 - LOAD OBS DATA #
######################

# Load Main Data
odf = pd.DataFrame()
for o in range(len(obsmods)):
    odf1 = pd.read_csv(obspaths[o]+'/'+obsmods[o]+'/ThicknessPhenology/'+thres+'/'+obsmods[o]+'_Regionalized_'+regvar+'_Historical_'+thres+'.csv')
    odf = odf.append(odf1,ignore_index=True)

odf = odf.loc[ (odf.Year >= aymin) & (odf.Year <= aymax) ]
regs = np.unique(odf['Region'])

# Load Surface Temperature
tof = pd.read_csv(tobpath)
tof = tof.loc[ (tof.Year >= aymin) & (tof.Year <= aymax) ]

# tof2 = pd.read_csv(tobpath2)
# tof2 = tof2.loc[ (tof2.Year >= aymin) & (tof2.Year <= aymax) ]
# tof2 = pd.DataFrame({'Year':np.arange(aymin,aymax+1),'tglobal':np.array(tof2.groupby('Year').mean()['tas'])})

##################################
# STEP 1 - CMIP6 AVERAGE & TREND #
##################################
print('Step 1')

# Load CMIP6 Data And Subset
filelist = md.listdir(inpath)
filelist = [f for f in filelist if f.split('_')[3].split('Regional-')[-1] == regvar and f.split('_')[2] == thres]
# filelist = [f for f in filelist if f.split('_')[3].split('Regional-')[-1] == regvar and f.split('_')[-2] in modstokeep and f.split('_')[2] == thres]

# Identify models that exist in most experiments
runs1 = [f for f in filelist if f.split('_')[-3] == experiment1]
runs2 = [f for f in filelist if f.split('_')[-3] == experiment2]
models1 = np.unique([f.split('_')[-2]+'_'+f.split('_')[-1].split('.csv')[0] for f in runs1])
models2 = np.unique([f.split('_')[-2]+'_'+f.split('_')[-1].split('.csv')[0] for f in runs2])
models = np.intersect1d(models1,models2)

# Load Temperature Data
tdf1 = pd.read_csv(path+"/"+ver+"_Annual_tas_"+experiment1+".csv")
tdf2 = pd.read_csv(path+"/"+ver+"_Annual_tas_"+experiment2+".csv")
tdf = tdf1.append(tdf2, ignore_index=True)
tdf = tdf.loc[ (tdf.Year >= aymin) & (tdf.Year <= aymax) & (tdf.Year != 2014) ]
tdf = tdf.sort_values(by=['Family','Member','Year'])

pdfavg, pdfsd = pd.DataFrame(), pd.DataFrame()
pdftrendlists = [[[] for i in range(5)] for v in range(len(varlist))  ]
pdftsenslists = [[[] for i in range(5)] for v in range(len(varlist))  ]
families = []
for m in models:

    pdf = pd.read_csv(inpath+"/"+[f for f in filelist if f.split('_')[-2]+'_'+f.split('_')[-1].split('.csv')[0] == m and f.split('_')[-3] == experiment1][0])
    pdf = pdf.append( pd.read_csv(inpath+"/"+[f for f in filelist if f.split('_')[-2]+'_'+f.split('_')[-1].split('.csv')[0] == m and f.split('_')[-3] == experiment2][0]), ignore_index=True)
    pdf = pdf[(pdf.Year >= aymin) & (pdf.Year <= aymax) ]
    pdf = pdf.sort_values(by=['Family','Member','Year'])

    pdfavg = pdfavg.append( pdf[['Model','Region','Family','Member']+varlist].groupby(by=['Family','Member','Region']).mean() )
    pdfsd = pdfsd.append( pdf[['Model','Region','Family','Member']+varlist].groupby(by=['Family','Member','Region']).std() )
    families.append( m.split('_')[-2] )

    for reg in regs:
        mdf = pdf[(pdf['Region'] == reg)]

        for v, var in enumerate(varlist):
            # Trend wrt Year
            lm = stats.linregress( mdf['Year'] , mdf[var] )

            for i in range(5):
                pdftrendlists[v][i].append( lm[i] )

            # Trend wrt Temperature
            lm = stats.linregress( np.array(tdf[(tdf['Family'] == m.split('_')[-2]) & (tdf['Member'] == mdf.iloc[0]['Member'])]['tglobal']) , np.array(mdf[var]) )

            for i in range(5):
                pdftsenslists[v][i].append( lm[i] )

# Reset indices to be columns
pdfavg = pdfavg.reset_index()
pdfsd = pdfsd.reset_index()

# Build overall pdf
pdftrend = pdfavg.iloc[:,0:3]
pdftrend['Family'] = np.repeat(families,len(regs))

for v, var in enumerate(varlist):
    pdftrend[var+'_Avg'] = np.array( pdfavg[var] )
    pdftrend[var+'_SD'] = np.array( pdfsd[var] )

    pdftrend[var+'_Trend'] = np.array(pdftrendlists[v][0])
    pdftrend[var+'_Trend_yint'] = np.array(pdftrendlists[v][1])
    pdftrend[var+'_Trend_r2'] = np.array(pdftrendlists[v][2])**2
    pdftrend[var+'_Trend_p'] = np.array(pdftrendlists[v][3])
    pdftrend[var+'_Trend_se'] = np.array(pdftrendlists[v][4])

    pdftrend[var+'_TSens'] = np.array(pdftsenslists[v][0])
    pdftrend[var+'_TSens_yint'] = np.array(pdftsenslists[v][1])
    pdftrend[var+'_TSens_r2'] = np.array(pdftsenslists[v][2])**2
    pdftrend[var+'_TSens_p'] = np.array(pdftsenslists[v][3])
    pdftrend[var+'_TSens_se'] = np.array(pdftsenslists[v][4])

del pdfavg, pdfsd, pdftrendlists, pdftsenslists

#####################################
# STEP 2 - OBSERVED AVERAGE & TREND #
#####################################
print('Step 2')

odfavg = odf.groupby(by=['Model','Region']).mean()
odfsd = odf.groupby(by=['Model','Region']).std()
odfavg = odfavg.reset_index()
odfsd = odfsd.reset_index()

odftrend = pd.DataFrame()
odftrend['Model'] = np.repeat(np.unique(odfavg['Model']),len(regs))
odftrend['Member'] = 1
odftrend['Region'] = np.repeat( [np.unique(odfavg['Region'])], len(obsmods), axis=0 ).flatten()

for var in varlist:
    odftrend[var+'_Avg'] = np.array( odfavg[var] )
    odftrend[var+'_SD'] = np.array( odfsd[var] )

for run in obsmods:
    for reg in regs:
        mdf = odf[(odf['Region'] == reg) & (odf['Model'] == run)]

        for var in varlist:
            # Trend wrt Year
            lm = stats.linregress( mdf['Year'] , mdf[var] )
            odftrend.loc[(odftrend['Region'] == reg) & (odftrend['Model'] == run),var+'_Trend'] = lm[0]
            odftrend.loc[(odftrend['Region'] == reg) & (odftrend['Model'] == run),var+'_Trend_yint'] = lm[1]
            odftrend.loc[(odftrend['Region'] == reg) & (odftrend['Model'] == run),var+'_Trend_r2'] = lm[2]
            odftrend.loc[(odftrend['Region'] == reg) & (odftrend['Model'] == run),var+'_Trend_p'] = lm[3]
            odftrend.loc[(odftrend['Region'] == reg) & (odftrend['Model'] == run),var+'_Trend_se'] = lm[4]

            # Trend wrt Temperature
            lm = stats.linregress( np.array(tof['tglobal']) , np.array(mdf[var]) )
            odftrend.loc[(odftrend['Region'] == reg) & (odftrend['Model'] == run),var+'_TSens'] = lm[0]
            odftrend.loc[(odftrend['Region'] == reg) & (odftrend['Model'] == run),var+'_TSens_yint'] = lm[1]
            odftrend.loc[(odftrend['Region'] == reg) & (odftrend['Model'] == run),var+'_TSens_r2'] = lm[2]
            odftrend.loc[(odftrend['Region'] == reg) & (odftrend['Model'] == run),var+'_TSens_p'] = lm[3]
            odftrend.loc[(odftrend['Region'] == reg) & (odftrend['Model'] == run),var+'_TSens_se'] = lm[4]

del odfavg, odfsd

##############################
# STEP 3 - MULTI-MODEL MEANS #
##############################
print('Step 3')

pdfmmm = pdftrend.loc[((pdftrend['Member'] == 1) & (pdftrend['Family'] != 'CESM2')) | ((pdftrend['Member'] == 4) & (pdftrend['Family'] == 'CESM2')),['Region']+[var+'_Avg' for var in varlist]+[var+"_Trend" for var in varlist]+[var+"_TSens" for var in varlist]].groupby('Region').mean()
pdfmmm = pdfmmm.reset_index()

odfmmm = odftrend.loc[np.in1d(odftrend['Model'],['Bootstrap','NASATeam','OSISAF'])][['Region']+[var+'_Avg' for var in varlist]+[var+"_Trend" for var in varlist]+[var+"_TSens" for var in varlist]].groupby('Region').mean()
odfmmm = odfmmm.reset_index()

#################################
# STEP 4 - INTERNAL VARIABILITY #
#################################
print('Step 4')

pdfse, pdfsemx = pdfmmm+0, pdfmmm+0
pdfsetemp = pdftrend.loc[np.in1d(pdftrend['Family'],modsforerror)]

# Identify the minimum numbers of members in the models chosen for internal variability
nmods = []
for model in modsforerror:
    nmods.append( np.unique(pdftrend.loc[pdftrend['Family'] == model,'Member']).shape[0] )
minn = min(nmods) - 1

# Calculate internal variability within each model
for model in modsforerror:
    mdf = pdftrend.loc[(pdftrend['Family'] == model)]

    # Identify maximum member for this model
    mems = np.unique(mdf['Member'])
    mems.sort()

    # If there are at least three runs for this model, calculate a standard deviation
    for reg in regs:
        for var in varlist:
            pdfsetemp.loc[(pdfsetemp['Family'] == model) & (pdfsetemp['Region'] == reg),var+'_Avg'] = md.sd_unbiased( mdf.loc[(mdf['Region'] == reg) & (mdf['Member'] <= mems[minn]),var+"_Avg"] )
            pdfsetemp.loc[(pdfsetemp['Family'] == model) & (pdfsetemp['Region'] == reg),var+'_Trend'] = md.sd_unbiased( mdf.loc[(mdf['Region'] == reg) & (mdf['Member'] <= mems[minn]),var+"_Trend"] )
            pdfsetemp.loc[(pdfsetemp['Family'] == model) & (pdfsetemp['Region'] == reg),var+'_TSens'] = md.sd_unbiased( mdf.loc[(mdf['Region'] == reg) & (mdf['Member'] <= mems[minn]),var+"_TSens"] )

# Calculate overall error in multi-model means
for reg in regs:
    for var in varlist:

        # CMIP6 Models - Mean
        pdfse.loc[pdfse['Region'] == reg,var+"_Avg"] = pdfsetemp.loc[pdfsetemp['Region'] == reg,var+"_Avg"].mean()
        pdfse.loc[pdfse['Region'] == reg,var+"_Trend"] = pdfsetemp.loc[pdfsetemp['Region'] == reg,var+"_Trend"].mean()
        pdfse.loc[pdfse['Region'] == reg,var+"_TSens"] = pdfsetemp.loc[pdfsetemp['Region'] == reg,var+"_TSens"].mean()

        # CMIP6 Models - Max
        pdfsemx.loc[pdfse['Region'] == reg,var+"_Avg"] = pdfsetemp.loc[pdfsetemp['Region'] == reg,var+"_Avg"].max()
        pdfsemx.loc[pdfse['Region'] == reg,var+"_Trend"] = pdfsetemp.loc[pdfsetemp['Region'] == reg,var+"_Trend"].max()
        pdfsemx.loc[pdfse['Region'] == reg,var+"_TSens"] = pdfsetemp.loc[pdfsetemp['Region'] == reg,var+"_TSens"].max()


##########################
# STEP 5 - WRITE TO FILE #
##########################
print('Step 5')

odftrend.to_csv(outpath+"/ThicknessAvgTrend"+str(aymin)+"-"+str(aymax)+"/ObsAvgTrend_"+regvar+"_"+experiment1+"-"+experiment2+"_"+thres+"_"+V+".csv",index=False)
odfmmm.to_csv(outpath+"/ThicknessAvgTrend"+str(aymin)+"-"+str(aymax)+"/ObsAvgTrend_MultiModelMean_"+regvar+"_"+experiment1+"-"+experiment2+"_"+thres+"_"+V+".csv",index=False)

pdftrend.to_csv(outpath+"/ThicknessAvgTrend"+str(aymin)+"-"+str(aymax)+"/CMIP6AvgTrend_"+regvar+"_"+experiment1+"-"+experiment2+"_"+thres+"_"+V+".csv",index=False)
# pdfmmm.to_csv(outpath+"/ThicknessAvgTrend"+str(aymin)+"-"+str(aymax)+"/CMIP6AvgTrend_MultiModelMedian_"+regvar+"_"+experiment1+"-"+experiment2+"_"+thres+"_"+V+".csv",index=False)
pdfse.to_csv(outpath+"/ThicknessAvgTrend"+str(aymin)+"-"+str(aymax)+"/CMIP6AvgTrend_MultiModelStdDev_"+regvar+"_"+experiment1+"-"+experiment2+"_"+thres+"_"+V+".csv",index=False)
pdfsemx.to_csv(outpath+"/ThicknessAvgTrend"+str(aymin)+"-"+str(aymax)+"/CMIP6AvgTrend_MultiModelMaxStdDev_"+regvar+"_"+experiment1+"-"+experiment2+"_"+thres+"_"+V+".csv",index=False)

pdfmmm.to_csv(outpath+"/ThicknessAvgTrend"+str(aymin)+"-"+str(aymax)+"/CMIP6AvgTrend_MultiModelMean_"+regvar+"_"+experiment1+"-"+experiment2+"_"+thres+"_"+V+".csv",index=False)


print('Complete')
