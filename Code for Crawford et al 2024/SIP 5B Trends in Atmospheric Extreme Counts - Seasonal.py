#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 24 May 2024
Modified: 03 Jul 2024

Author: Alex Crawford

Purpose: Calculate trends in atmospheric variables around the time of SIE expansion pauses
"""

import pandas as pd
import CycloneModule_13_2 as md
import numpy as np
from scipy.stats import mannwhitneyu

inpath = '/home/alex/OneDrive/CMIP6 OWP/Sea Ice Extent Expansion Pauses/ERA5_and_SST_Avgs_forbboxes_forevents_V2.csv'
outpath = '/home/alex/OneDrive/CMIP6 OWP/Sea Ice Extent Expansion Pauses/ERA5_and_SST_Avgs_forbboxes_seasonal_extremecounts.csv'

y2 = 2005 # first year of second period

vlist = ['sst','t','msl','seb','u10','v10']

pdf = pd.read_csv(inpath)
# If the value is exactly 0, don't trust it -- set to NaN
for var in vlist:
    pdf[var] = np.where(pdf[var] == 0,np.nan, pdf[var])
pdf['syear'] = np.where(pdf['month'] < 3,pdf['year']-1,pdf['year'])
pdf = pdf.loc[np.in1d(pdf['enddate'],['2016-11-19','2008-12-17'])]
pdf90 = pdf.loc[:,['bbox']+vlist].groupby(by=['bbox']).quantile(0.9).reset_index()
pdf10 = pdf.loc[:,['bbox']+vlist].groupby(by=['bbox']).quantile(0.1).reset_index()

bboxlist, varlist, p90list, p10list, countlist, c90aftlist, c90beflist, c90plist, c10aftlist, c10beflist, c10plist = [[] for i in range(11)]

for reg in np.unique(pdf['bbox']):
    df = pdf.loc[(pdf['bbox'] == reg)]
    
    for var in vlist:
        
        # determine if 6-day periods are above p90 or below p10
        p90 = pdf90.loc[(pdf90['bbox'] == reg), var].values[0]
        p10 = pdf10.loc[(pdf10['bbox'] == reg), var].values[0]
       
        df[var+'_gtp90'] = df[var] > p90
        df[var+'_ltp10'] = df[var] < p10
        
        # count number of extremes before and after user-defined break point
        count = np.sum(np.isfinite(df[var]))
        c90bef = np.sum(df.loc[df['syear'] < y2,var+'_gtp90'])
        c90aft = np.sum(df.loc[df['syear'] >= y2,var+'_gtp90'])
        c10bef = np.sum(df.loc[df['syear'] < y2,var+"_ltp10"])
        c10aft = np.sum(df.loc[df['syear'] >= y2,var+"_ltp10"])
        
        # Bootstrap a p-value
        naft = len(df.loc[df['syear'] >= y2,])
        
        c90afttest, c10afttest = [], []
        for i in range(1000):
            c90afttest.append( np.random.choice(df[var+'_gtp90'],size=naft, replace=False).sum() )  
            c10afttest.append( np.random.choice(df[var+'_ltp10'],size=naft, replace=False).sum() )
            
        c90p = np.mean( c90aft > np.array(c90afttest) ) * 100
        c10p = np.mean( c10aft > np.array(c10afttest) ) * 100
        
        
        bboxlist.append(reg), varlist.append(var) 
        countlist.append(count), p90list.append(p90), p10list.append(p10)
        c90aftlist.append(c90aft), c90beflist.append(c90bef), c90plist.append(c90p)
        c10aftlist.append(c10aft), c10beflist.append(c10bef), c10plist.append(c10p)
        
tdf = pd.DataFrame({'bbox':bboxlist,'var':varlist,
                    'ntotal':countlist,'p90':p90list, 'p10':p10list,
                    'n90aft':c90aftlist,'n90bef':c90beflist,'n90p':c90plist,
                    'n10aft':c10aftlist,'n10bef':c10beflist,'n10p':c10plist})

tdf.to_csv(outpath,index=False)
