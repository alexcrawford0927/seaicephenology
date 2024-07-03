#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 5 Apr 2024
Modified: 3 Jul 2024
Author: Alex Crawford

Purpose:
Calculate ERA5 and SST anomalies for multiple periods spanning about 1 month,
centered on each SIE pause event for bounding boxes of interest

"""

import cftime
import pandas as pd
import xarray as xr
# import xesmf as xe
import numpy as np
import scipy.stats as stats
# from metpy.calc import divergence
# from metpy.units import units
import CycloneModule_13_2 as md

'''*******************************************
Declare Variables
*******************************************'''
path = '/media/alex/Datapool'
sstpath = path+'/NOAA OISST/'
eventpath = '/home/alex/OneDrive/CMIP6 OWP/Sea Ice Extent Expansion Pauses/Pauses_combined_final_bytime_276-400.csv'

epath1, epath2 = path+'/ERA5', path+'/ERA5'
e5filename = ['SLP','HorizontalWind_10m','Temperature_925','Radiation_LW_Surface','TurbulentEnergyFluxes','SEB']
era5path = [epath1,epath1,epath1,epath2,epath1,epath2]
vvalue = [['msl'],['u10','v10'],['t'],['str'],['sshf','slhf'],['seb']]
tstep = ['Hourly','3h','3h','Hourly','3h','3h']

ymin, ymax = 1979, 2023
shifts = np.array([-12, -6, 0, 6, 12])

# bboxname = ['Hudson Bay','Foxe Basin','East Greenland Sea','Barents & Kara Seas','Chukchi Sea','Bering Sea','Sea of Okhotsk']
# bbox = [[55,65,-85,-75],[64,70,-85,-70],[70,80,-15,10],[72,83,20,80],[67,75,-180,-159],[56,64,-180,-159],[50,65,130,165]]

bboxname = ['Baffin Bay']
bbox = [[65,80,-90,-50]] # [latmin, latmax, lonmin, lonmax]

eraprjpath = path+'/ERA5/DailyClimatology_1981-2010/ERA5_SLP_DailyClimatology_1981-2010.nc'
psnprjpath = path+'/Projections/psn_projection.nc'
sstprjpath = sstpath+"/OISST_DailyClimatology_1981-2010.nc"

outpath = '/home/alex/OneDrive/CMIP6 OWP/Sea Ice Extent Expansion Pauses/'

'''*******************************************
Main Analysis
*******************************************'''

# Load Events
pdf = pd.read_csv(eventpath)
pdf = pdf.loc[(pdf['Use'] == 1) & (pdf['year'] > 1978)]
pdf = pdf.reset_index()

# Load projections
era5prj = xr.open_dataset(eraprjpath)
sstprj = xr.open_dataset(sstprjpath)
sstprj['lon2'] = (('lon'),np.where(sstprj['lon'] > 180, sstprj['lon'] - 360, sstprj['lon']))

era5coslats = np.cos(era5prj['latitude']*np.pi/180)
sstcoslats = np.cos(sstprj['lat']*np.pi/180)


# Create masks
era5mask, sstmask = [], []

for b in bbox:
    
   era5mask.append( (era5prj['latitude'] >= b[0]) & (era5prj['latitude'] <= b[1]) & (era5prj['longitude'] >= b[2]) & (era5prj['longitude'] <= b[3]) )
   sstmask.append( (sstprj['lat'] >= b[0]) & (sstprj['lat'] <= b[1]) & (sstprj['lon2'] >= b[2]) & (sstprj['lon2'] <= b[3]) )
    
enddates = [list(pdf.loc[e,['year','month','day']].astype(int)) for e in range(len(pdf))]
enddatesstr = np.array([str(int(ed[0]))+"-"+md.dd[int(ed[1])-1]+"-"+md.dd[int(ed[2])-1] for ed in enddates])


# ### Initiate Output DataFrame ###
alldates = np.array([[md.timeAdd(enddate,[0,0,shift]) for shift in shifts] for enddate in enddates])
years = np.arange(ymin,ymax+1)
ny = ymax-ymin+1
nb = len(bbox)
nd = len(pdf)
ns = len(shifts)

odf = pd.DataFrame({'enddate': np.repeat(enddatesstr,ny*nb*ns),
                    'shift': np.tile(np.repeat(shifts,ny*nb),nd),
                    'year': np.tile(np.repeat(years,nb),nd*ns),
                    'month': np.tile(np.repeat(alldates[:,:,1].flatten(),ny*nb),1),
                    'day': np.tile(np.repeat(alldates[:,:,2].flatten(),ny*nb),1),
                    'bbox': np.tile(bboxname,nd*ns*ny)})
odf['sst'] = np.nan
for val in np.concatenate(vvalue):
    odf[val] = np.nan

print("Load SST")

### SST ###
sstarrs = [ xr.open_dataset(sstpath+"/Daily/sst.day.mean.1981.nc") ]
for y in years[years > 1981]:
    print('--- '+str(y))

    # Load the netcdf files for the current year and prior year
    sstarrs.append( xr.open_dataset(sstpath+"/Daily/sst.day.mean."+str(y)+".nc") )
    sstarr = xr.concat(sstarrs,dim='time')
    sstarrs = sstarrs[1:]
    
    # For each event
    for e in range(len(enddates)):
        ydate = [y,enddates[e][1], enddates[e][2]]
        
        for shift in shifts:
            sd = md.timeAdd(ydate,[0,0,shift-5])
            ed = md.timeAdd(ydate,[0,0,shift])
            
            sdstr = str(sd[0])+"-"+md.dd[sd[1]-1]+"-"+md.dd[sd[2]-1]
            edstr = str(ed[0])+"-"+md.dd[ed[1]-1]+"-"+md.dd[ed[2]-1]
            
            # Temporal average
            sst = sstarr.sel(time=slice(sdstr,edstr)).mean(dim='time')
    
            for b in range(len(bbox)):
                # Spatial average
                sstval = np.nansum(sst['sst'] * sstmask[b] * sstcoslats) / np.sum(np.isfinite(sst['sst']) * sstmask[b] * sstcoslats).values
                
                # Append to Data Frame
                odf.loc[(odf['year'] == y) & (odf['enddate'] == enddatesstr[e]) & (odf['shift'] == shift) & (odf['bbox'] == bboxname[b]),'sst'] = sstval

### Write to File ###
odf.to_csv(outpath+"/ERA5_and_SST_Avgs_forbboxes_forevents_V3-Oct-Jan.csv", index=False)
del sstarrs, sstarr, sstval

odf = pd.read_csv(outpath+"/ERA5_and_SST_Avgs_forbboxes_forevents_V3-Oct-Jan.csv")

### ERA5 ###
for vi in range(len(vvalue)):
    print(vvalue[vi])
    
    for y in years:
        print('--- '+str(y))
        
        # Load only Oct this year through Jan next year
        if y == ymax:
            e5arr = xr.concat([ xr.open_dataset(era5path[vi]+"/"+e5filename[vi]+"/ERA5_"+e5filename[vi]+"_"+tstep[vi]+"_"+str(y)+"10.nc"),
                       xr.open_dataset(era5path[vi]+"/"+e5filename[vi]+"/ERA5_"+e5filename[vi]+"_"+tstep[vi]+"_"+str(y)+"11.nc"),
                       xr.open_dataset(era5path[vi]+"/"+e5filename[vi]+"/ERA5_"+e5filename[vi]+"_"+tstep[vi]+"_"+str(y)+"12.nc") ], dim='time')
        else:
            e5arr = xr.concat([ xr.open_dataset(era5path[vi]+"/"+e5filename[vi]+"/ERA5_"+e5filename[vi]+"_"+tstep[vi]+"_"+str(y)+"10.nc"),
                       xr.open_dataset(era5path[vi]+"/"+e5filename[vi]+"/ERA5_"+e5filename[vi]+"_"+tstep[vi]+"_"+str(y)+"11.nc"),
                       xr.open_dataset(era5path[vi]+"/"+e5filename[vi]+"/ERA5_"+e5filename[vi]+"_"+tstep[vi]+"_"+str(y)+"12.nc"),
                       xr.open_dataset(era5path[vi]+"/"+e5filename[vi]+"/ERA5_"+e5filename[vi]+"_"+tstep[vi]+"_"+str(y+1)+"01.nc") ], dim='time')
            
        # For each event
        for e in range(len(enddates)):
            if enddates[e][1] == 1:
                ydate = [y+1,enddates[e][1], enddates[e][2]]
            else:
                ydate = [y,enddates[e][1], enddates[e][2]]
            
            for shift in shifts: 
                sd = md.timeAdd(ydate,[0,0,shift-5])
                ed = md.timeAdd(ydate,[0,0,shift])
                
                sdstr = str(sd[0])+"-"+md.dd[sd[1]-1]+"-"+md.dd[sd[2]-1]
                edstr = str(ed[0])+"-"+md.dd[ed[1]-1]+"-"+md.dd[ed[2]-1]
                
                # Temporal average
                e5vals = e5arr.sel(time=slice(sdstr,edstr)).mean(dim='time')
        
                for b in range(len(bbox)):
                    
                    for vj in range(len(vvalue[vi])):
                        # Spatial average
                        e5val = np.sum(e5vals[vvalue[vi][vj]] * era5mask[b] * era5coslats).values / np.sum( era5mask[b] * era5coslats ).values
                    
                        # Append to Data Frame
                        odf.loc[(odf['year'] == ydate[0]) & (odf['enddate'] == enddatesstr[e]) & (odf['shift'] == shift) & (odf['bbox'] == bboxname[b]),vvalue[vi][vj]] = e5val

    ### Write to File ###
    odf.to_csv(outpath+"/ERA5_and_SST_Avgs_forbboxes_forevents_V3-Oct-Jan.csv", index=False)
    