#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 15 May 2024
Modified: 03 Jul 2024

Purpose: Calculate the SEB from net LW (downward), net SW (downward), 
surface LH flux (downward), and surface SH flux (downward) from ERA5

Note: The description of the net LW flux describes it as upward, but the convention
are clearly positive downward, for the net LW is nearly always a negative value.

"""
import xarray as xr
import CycloneModule_13_2 as md
import numpy as np

path1 = '/media/alex/Datapool/ERA5/TurbulentEnergyFluxes'
path3 = '/media/alex/Datapool/ERA5/Radiation_SW_Surface'
path4 = '/media/alex/Datapool/ERA5/Radiation_LW_Surface'
outpath = '/media/alex/Datapool/ERA5/SEB'

var1 ='slhf'
var2 = 'sshf'
var3 = 'ssr'
var4 = 'str'

# Identify Files shared by All Variables
files1 = md.listdir(path1)
files3 = md.listdir(path3)
files4 = md.listdir(path4)

times1 = [f.split('_')[-1][:-3] for f in files1]
times3 = [f.split('_')[-1][:-3] for f in files3]
times4 = [f.split('_')[-1][:-3] for f in files4]

times = np.intersect1d(times1, np.intersect1d(times3,times4))

files1 = [f for f in files1 if f.split('_')[-1][:-3] in times]
files3 = [f for f in files3 if f.split('_')[-1][:-3] in times]
files4 = [f for f in files4 if f.split('_')[-1][:-3] in times]

for t in list(range(len(times)))[-2:-1]:
    xr1 = xr.open_dataset(path1+"/"+files1[t])
    xr3 = xr.open_dataset(path3+"/"+files3[t])
    xr4 = xr.open_dataset(path4+"/"+files4[t])

    if t == 0:
        seb = xr3[var3].data[::3,:,:][3:,:,:] + xr1[var1].data + xr1[var2].data + xr4[var4].data[::3,:,:][3:,:,:]
    else:
        seb = xr3[var3].data[::3,:,:] + xr1[var1].data + xr1[var2].data + xr4[var4].data[::3,:,:]
    
    xr1['seb'] = (('time','latitude','longitude'), seb)
    
    xro = xr1.drop_vars((var1,var2))
    xro.to_netcdf(outpath+"/ERA5_SEB_3h_"+times[t]+".nc")
    
        





