'''
Author: Alex Crawford
Date Created: 9 Dec 2020
Date modified: 9Dec 2020
Purpose: Identify the percentage of the pan-Arctic sea ice domain 
that is covered by each region
'''

import netCDF4 as nc
import numpy as np
import pandas as pd

proj = nc.Dataset('/Volumes/Miranda/Projections/psn_projection.nc')

reg = proj['reg'][:].data
area = proj['area'][:].data

totalarea = np.sum(area[ (reg >= 2) & (reg <= 16)])


pdf = pd.DataFrame()
pdf['Region'] = [0] + list(range(2,17))
pdf['Area'] = 0

for rg in pdf['Region']:
    if rg == 0:
        pdf.loc[pdf['Region'] == rg,'Area'] = np.sum(area[ (reg >= 2) & (reg <= 16)])
    else:
        pdf.loc[pdf['Region'] == rg,'Area'] = np.sum(area[ reg == rg])

pdf['RelArea'] = pdf['Area'] / totalarea * 100
