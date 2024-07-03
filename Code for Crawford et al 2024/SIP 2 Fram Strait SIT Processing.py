#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 03 Apr 2024
Modified: 03 Jul 2024

Author: Alex Crawford

Purpose: Process Fram Strait Mooring Data in the style of Sumata et al. (2023)

Sumata, H., L. de Steur, D. V. Divine, M. A. Granskog, and S. Gerland, 2023: 
Regime shift in Arctic Ocean sea ice thickness. Nature, 615, 443–449, 
https://doi.org/10.1038/s41586-022-05686-x.
"""

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import lognorm
import matplotlib.pyplot as plt
import datetime, sys, netCDF4, os, csv, pickle
import calendar, subprocess
from dateutil.relativedelta import relativedelta
import pandas as pd
import uls_util

frampath = '/Volumes/Miranda/SeaIce/Ice_thickness_distribution_Fram_Strait/data_draft_PDF_range_0.0-8.0_bin0.1m'
outpath = '/Volumes/Miranda/SeaIce/Ice_thickness_distribution_Fram_Strait/FramStrait_SIT_Summary.csv'
# outpath2 = '/Volumes/Miranda/SeaIce/Ice_thickness_distribution_Fram_Strait/FramStrait_SIT_Summary-PercentilesOnly.csv'

files = os.listdir(frampath)

latlist, lonlist, smpintlist, timeperclist, timestlist, timeedlist, nsamplist = [[] for i in range(7)]
meanlist, mode2list, modefittedlist, medianlist, p75list, p90list = [[] for i in range(6)]

for file in files:
    
    with open(frampath+"/"+file, 'rb') as f:
        pdf = pickle.load(f) # pdf is uls_util.IceDraftPdf object
    
    # Store Attributes for Entire File
    attributes = [att for att in dir(pdf) if att.startswith('__') == 0]
    latlist.append( pdf.lsat ), lonlist.append( pdf.lon ), nsamplist.append( pdf.num_samples )
    smpintlist.append( pdf.smpl_intval ), timeperclist.append( pdf.temporal_data_coverage )
    timestlist.append( pdf.sampling_period[0] ), timeedlist.append( pdf.sampling_period[1] )
    
    # Load the probability density function for the file as a pandas dataframe 
    df = pd.DataFrame({'xbin':pdf.x_bin, 'ycount':pdf.y_count})

    # Identify the mean, median
    meanlist.append( np.sum(df['xbin'] * df['ycount']) / np.sum(df['ycount']) )
    
    allvals = np.concatenate([np.repeat(df['xbin'][i],df['ycount'][i]) for i in range(len(df))])
    medianlist.append( np.median(allvals) )
    
    p75list.append(np.percentile(allvals,75))
    p90list.append(np.percentile(allvals,90))
    
    # Identify the second mode
    df['ycount-min5'] = df['ycount'].rolling(5,center=True).min() # Finding the minimum value for every set of five consecutive bins
    
    # Start the search for a maximum after the first local minimum
    # If there is no local minimum, start at beginning
    try:
        firstmin = np.where(df['ycount'] == df['ycount-min5'])[0][0]
    except:
        firstmin = 0
    
    # Find the maximum after this minimum threshold
    mode2 = np.where(df['ycount'].values == np.nanmax(df['ycount'].values[firstmin:]))[0][0] # max freq bin following the first minimum
    mode2list.append(df['xbin'].values[mode2])
    
    # Find the maximum of a log-normal distrubtion fit to the data
    valstofit = np.concatenate([np.repeat(df['xbin'][i],df['ycount'][i]) for i in range(firstmin,len(df))])
    xvals = df['xbin'].values[firstmin:(len(df))]

    param = lognorm.fit(valstofit,loc=df['xbin'].values[firstmin])
    valsfitted = lognorm.pdf(xvals, param[0], loc=param[1], scale=param[2])
    modefitted = np.where(valsfitted == np.max(valsfitted))[0][0] + firstmin
    modefittedlist.append(df['xbin'].values[modefitted])
    
    
# Save File
output = pd.DataFrame({'lat':latlist, 'lon':lonlist, 'time_start':timestlist,
                        'time_end':timeedlist, 'sample_count':nsamplist,
                        'valid_percent':timeperclist, 'sample_interal':smpintlist,
                        'avg_sit':meanlist,'median_sit':medianlist, 'p75_sit':p75list,
                        'mode2_sit':mode2list,'modefitted_sit':modefittedlist})

output.to_csv(outpath,index=False)

# output2 = pd.DataFrame({'lat':latlist, 'lon':lonlist, 'time_start':timestlist,
#                        'time_end':timeedlist, 'sample_count':nsamplist,
#                        'valid_percent':timeperclist, 'sample_interal':smpintlist,
#                        'median_sit':medianlist, 'p75_sit':p75list, 'p90_sit':p90list})

# output2.to_csv(outpath2,index=False)

### Figure for testing ###
# fig = plt.figure(figsize=(6,6))
# # plt.hist(allvals,bins=len(df))
# plt.hist(valstofit, bins=len(xvals), density=True)
# plt.plot( xvals, valsfitted )