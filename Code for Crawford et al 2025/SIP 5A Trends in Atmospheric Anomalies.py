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

inpath = '/home/alex/OneDrive/CMIP6 OWP/Sea Ice Extent Expansion Pauses/ERA5_and_SST_Avgs_forbboxes_forevents_V2.csv'
outpath = '/home/alex/OneDrive/CMIP6 OWP/Sea Ice Extent Expansion Pauses/ERA5_and_SST_Avgs_forbboxes_forevents_trends_V2.csv'

pdf = pd.read_csv(inpath)

bboxlist, varlist, datelist, trendlist, intlist, r2list, plist = [[] for i in range(7)]

for reg in np.unique(pdf['bbox']):
    for date in np.unique(pdf['enddate']):
        df = pdf.loc[(pdf['enddate'] == date) & (pdf['bbox'] == reg) & (pdf['shift'] == 0)]
        
        for var in ['sst','t','msl','seb','u10','v10']:
            x = df['year'].values
            y = df[var].values
            mask = np.isfinite(x+y)
            x, y = x[mask], y[mask]
    
            ts = md.theilsenslope(x, y)
            
            bboxlist.append(reg), datelist.append(date)
            varlist.append(var), trendlist.append(ts[0])
            intlist.append(ts[1]), r2list.append(ts[2]), plist.append(ts[3])
            
tdf = pd.DataFrame({'enddate':datelist,'bbox':bboxlist,'var':varlist,'trend':trendlist,
                    'yint':intlist, 'r2':r2list, 'pvalue':plist})
tdf.to_csv(outpath,index=False)
