#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 01 May 2024
Modified: 03 Jul 2024

Author: Alex Crawford

Purpose: Identify pauses in autumn SIE expansion based on a combination of
CDR, Bootstrap, and NASA Team products

Note: Using the Goddard versions of Bootstrap and NASA Team
Note: Earlier versions also used Sea Ice Index or a percentile definition for
"pause". Some steps have been left but commented out in case of future reversion
to that method.
"""

import pandas as pd
import numpy as np
import CycloneModule_13_2 as md

ymin, ymax = 1979, 2022
dmin, dmax = 276, 428 # 276, 369
inpath = '/home/alex/OneDrive/CMIP6 OWP/Sea Ice Extent Expansion Pauses'

# inpath = '/Users/acrawfora/Library/CloudStorage/OneDrive-UniversityofManitoba/CMIP6 OWP/Sea Ice Extent Expansion Pauses/'

siedfall = pd.read_csv(inpath+"/SIE_Daily_ArcticWide.csv")
siedfall = siedfall.loc[:,[col for col in siedfall.columns if col not in ['SeaIceIndex']]]
# siidf = pd.read_csv(inpath+"/SIE_Daily_ArcticWide_sii.csv")
# siedfall = pd.merge(siedfall,siidf.loc[:,['JulianDay','SeaIceIndex','sie_sii','sii_delsie_5day']],how='inner',on='JulianDay')
siedfall['SYear'] = np.where(siedfall['month'] > 8, siedfall['year'], siedfall['year']-1)
siedfall = siedfall.loc[(siedfall['SYear'] >= ymin) & (siedfall['SYear'] <= ymax)]
siedfall['DOY'] = [md.daysBetweenDates([siedfall.loc[i,'SYear'],1,0],[siedfall.loc[i,'year'],siedfall.loc[i,'month'],siedfall.loc[i,'day']]) for i in siedfall.index]

cdrall = siedfall.loc[(siedfall['cdr_delsie_5day'] != 0) & (siedfall['DOY'] > dmin) & ( (siedfall['DOY'] < dmax) )] 
ntall = siedfall.loc[(siedfall['nt_delsie_5day'] != 0) & (siedfall['DOY'] > dmin) & ( (siedfall['DOY'] < dmax) )] 
btall = siedfall.loc[(siedfall['bt_delsie_5day'] != 0) & (siedfall['DOY'] > dmin) & ( (siedfall['DOY'] < dmax) )] 

# Find Low percentile growth
# lp = 5
# btlow = btall.loc[btall['bt_delsie_5day'] < np.percentile(btall['bt_delsie_5day'], lp)]
# ntlow = ntall.loc[ntall['nt_delsie_5day'] < np.percentile(ntall['nt_delsie_5day'], lp)]
# siilow = siiall.loc[siiall['sii_delsie_5day'] < np.percentile(siiall['sii_delsie_5day'], lp)]
# cdrlow = cdrall.loc[cdrall['cdr_delsie_5day'] < np.percentile(cdrall['cdr_delsie_5day'], lp)]

btlow = btall.loc[btall['bt_delsie_5day'] < 0.1]
ntlow = ntall.loc[ntall['nt_delsie_5day'] < 0.1]
# siilow = siiall.loc[siiall['sii_delsie_5day'] < 0.1]
cdrlow = cdrall.loc[cdrall['cdr_delsie_5day'] < 0.1]

# Find Negative Growth
btneg = btall.loc[btall['bt_delsie_5day'] < 0]
ntneg = ntall.loc[ntall['nt_delsie_5day'] < 0]
# siineg = siiall.loc[siiall['sii_delsie_5day'] < 0]
cdrneg = cdrall.loc[cdrall['cdr_delsie_5day'] < 0]

# Find intersection of low and/or negative growth amongst multiple methods
siedfall['LowVal'] = 0
siedfall['NegVal'] = 0

for jd in siedfall['JulianDay']:  
    siedfall.loc[siedfall['JulianDay'] == jd,'LowVal'] = int(jd in btlow['JulianDay'].values) + int(jd in ntlow['JulianDay'].values) + int(jd in cdrlow['JulianDay'].values)
    siedfall.loc[siedfall['JulianDay'] == jd,'NegVal'] = int(jd in btneg['JulianDay'].values) + int(jd in ntneg['JulianDay'].values) + int(jd in cdrneg['JulianDay'].values)

pauseprelim = siedfall.loc[(siedfall['NegVal'] >= 1)]
pausefinal = siedfall.loc[(siedfall['LowVal'] >= 3) & (siedfall['NegVal'] >= 1)]

# Write to file
pauseprelim.to_csv(inpath+"/Pauses_combined_anyneg_bytime_"+str(dmin)+"-"+str(dmax)+".csv",index=False)
pausefinal.to_csv(inpath+"/Pauses_combined_final_bytime_"+str(dmin)+"-"+str(dmax)+".csv",index=False)
