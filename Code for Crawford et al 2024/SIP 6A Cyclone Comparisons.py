#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 21 May 2024
Modified: 03 Jul 2024

Author: Alex Crawford

Purpose: Manual identification of extratropical cyclones from the CNECT database
that impacted sea ice extent pause events.

Note: Only useful in combination with iterative creation of maps of cyclone tracks
(Script SIP 6B)
"""
import numpy as np
import pandas as pd

bbox = '45_-90_55_-70' #'65_-30_80_30' # '50_130_65_-160' # '65_10_85_70' #  
pdf = pd.read_csv('/Volumes/Cressida/CycloneTracking/tracking13_2E5R/BBox10/AggregationCyclone/BBox10_AvgInBBox_CycloneTracksEvents_13_2E5R_1940_2023_bbox'+bbox+'.csv')

pdf2 = pdf.loc[(pdf['year'] == 2010) & (pdf['month'] == 12)]
tids = pdf2.tid.values
print(tids)

tdf = pd.DataFrame({'tid':[409,484]})

minp, maxpgrad = [], []
for tid in tdf['tid']:
    minp.append( np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 1)),'minp'].values < pdf2.loc[pdf2['tid'] == tid,'minp'].values[0]) )
    maxpgrad.append( np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 1)),'maxpgrad'].values < pdf2.loc[pdf2['tid'] == tid,'maxpgrad'].values[0]) )

tdf['minp'] = np.array(minp)*100
tdf['maxpgrad'] = np.array(maxpgrad)*100


###############
'''
16 Dec 1980, Nordic: 395,429*,436*,501
26 Dec 1982, Nordic: 613, 702
15 Dec 1983, Pacific: 255, 301, 349, 360, 381
25 Nov 1984, Nordic: 612, 665, 706
24 Dec 1990, Nordic: 624, 661, 673*, 696*, 740, 1*
26 Dec 1999, Nordic: 510, 551... Pacific: 569
23 Dec 2004, Nordic: 530, 671, 715... Pacific: 608, 659, 660, 672, 722
11 Dec 2008, Nordic: 257, 291, 359, 419, 443... Pacific: 289, 295, 367*, 403, 412
14 Dec 2010*, Eastern Canada: 409, 484
02 Nov 2013, Nordic: 29, 31, 45, 112 ... Pacific: 60
14 Dec 2013, Nordic: 352, 359, 360, 362
26 Dec 2015, Nordic: 589, 621, 673, 8
13 Nov 2016, Nordic: 220, 281, 302
'''

trs = pd.read_pickle('/Volumes/Cressida/CycloneTracking/tracking13_2E5R/BBox10/CycloneTracks/1990/BBox10cyclonetracks199012.pkl')
tr = [tr for tr in trs if tr.tid == 673][0]
np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 2)),'minp'].values < tr.data.p_cent.min()/100)*100
np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 2)),'maxpgrad'].values < tr.data.p_grad.max()/100*1000*1000)*100

trs = pd.read_pickle('/Volumes/Cressida/CycloneTracking/tracking13_2E5R/BBox10/CycloneTracks/1984/BBox10cyclonetracks198411.pkl')
tr = [tr for tr in trs if tr.tid == 675][0]
np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 1)),'minp'].values < tr.data.p_cent.min()/100)*100
np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 1)),'maxpgrad'].values < tr.data.p_grad.max()/100*1000*1000)*100


trs = pd.read_pickle('/Volumes/Cressida/CycloneTracking/tracking13_2E5R/BBox10/CycloneTracks/2015/BBox10cyclonetracks201512.pkl')
tr = [tr for tr in trs if tr.tid == 621][0]
np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 1)),'minp'].values < tr.data.p_cent.min()/100)*100
np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 1)),'maxpgrad'].values < tr.data.p_grad.max()/100*1000*1000)*100

trs = pd.read_pickle('/Volumes/Cressida/CycloneTracking/tracking13_2E5R/BBox10/CycloneTracks/2015/BBox10cyclonetracks201512.pkl')
tr = [tr for tr in trs if tr.tid == 673][0]
np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 1)),'minp'].values < tr.data.p_cent.min()/100)*100
np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 1)),'maxpgrad'].values < tr.data.p_grad.max()/100*1000*1000)*100

trs = pd.read_pickle('/Volumes/Cressida/CycloneTracking/tracking13_2E5R/BBox10/CycloneTracks/2016/BBox10cyclonetracks201601.pkl')
tr = [tr for tr in trs if tr.tid == 8][0]
np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 1)),'minp'].values < tr.data.p_cent.min()/100)*100
np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 1)),'maxpgrad'].values < tr.data.p_grad.max()/100*1000*1000)*100

trs = pd.read_pickle('/Volumes/Cressida/CycloneTracking/tracking13_2E5R/BBox10/CycloneTracks/1999/BBox10cyclonetracks199912.pkl')
tr = [tr for tr in trs if tr.tid == 569][0]
np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 1)),'minp'].values < tr.data.p_cent.min()/100)*100
np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 1)),'maxpgrad'].values < tr.data.p_grad.max()/100*1000*1000)*100

trs = pd.read_pickle('/Volumes/Cressida/CycloneTracking/tracking13_2E5R/BBox10/CycloneTracks/1999/BBox10cyclonetracks199912.pkl')
tr = [tr for tr in trs if tr.tid == 510][0]
np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 1)),'minp'].values < tr.data.p_cent.min()/100)*100
np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 1)),'maxpgrad'].values < tr.data.p_grad.max()/100*1000*1000)*100


trs = pd.read_pickle('/Volumes/Cressida/CycloneTracking/tracking13_2E5R/BBox10/CycloneTracks/1991/BBox10cyclonetracks199101.pkl')
tr = [tr for tr in trs if tr.tid == 1][0]
np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 1)),'minp'].values < tr.data.p_cent.min()/100)*100
np.mean(pdf.loc[(pdf['year'] > 1978) & ((pdf['month'] > 9) | (pdf['month'] < 1)),'maxpgrad'].values < tr.data.p_grad.max()/100*1000*1000)*100
