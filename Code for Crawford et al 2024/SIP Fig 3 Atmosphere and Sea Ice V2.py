#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 05 Apr 2024
Modified: 03 Jul 2024

Author: Alex Crawford

Purpose: Make a Plot of the SLP and wind, SLP and wind amlies, 2-m temperature & amly
SEB amly, SST and the SIC growth/loss, and SST & SIC amly period of each event
"""

import cftime
import pandas as pd
import xarray as xr
import xesmf as xe
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from metpy.calc import divergence
from metpy.units import units
import CycloneModule_13_2 as md
from matplotlib.colors import ListedColormap

'''*******************************************
Declare Variables
*******************************************'''

path = '/media/alex/Datapool'
sicpath = path+'/SeaIce/G02202_V4/Daily'
siapath = path+'/SeaIce/Age'
simpath = path+'/SeaIce/Motion'
sstpath = path+'/NOAA OISST/'
eventpath = '/home/alex/OneDrive/CMIP6 OWP/Sea Ice Extent Expansion Pauses/Pauses_combined_final_bytime_276-400.csv'

epath1, epath2 = path+'/ERA5', path+'/ERA5'
e5filename = ['SLP','HorizontalWind_10m','Temperature_925','SEB']
era5path = [epath1,epath1,epath1,epath1]
vvalue = ['msl',['u10','v10'],'t','seb']
tstep = ['Hourly','3h','3h','3h']

breaks = [ np.arange(96800,105000,400), '', np.arange(-30,16,5),  np.arange(-600,510, 60)*3600]
breakstick = [ np.arange(96800,105000,800), '', np.arange(-30,16,5),  np.arange(-600,510, 60)*3600]
breakslab = [ np.arange(968,1050,8), '', np.arange(-30,16,5),  np.arange(-600,510, 60)]
colorbarlabel = ['Sea-Level Pressure (hPa)', '', '925 hPa Temperature (°C)', 'Net Surface LW Radiation (W m$^{-2}$)','Surface Turbulent Heat Flux (W m$^{-2}$']
cmap = [ plt.cm.PuOr, '', plt.cm.Reds, plt.cm.Oranges, plt.cm.Blues]

breaksamly = [ md.divbreaks(1200,100), '', md.divbreaks(12,2), md.divbreaks(240,40)*3600]
breakstickamly = [ md.divbreaks(1200,200), '',  md.divbreaks(12,2), md.divbreaks(240,40)*3600]
breakslabamly = [ md.divbreaks(12,2), '', md.divbreaks(12,2), md.divbreaks(240,40)]
colorbarlabelamly = ['SLP Anomaly (hPa)', '10 m s$^{-1}$', '925 hPa Temperature Amly (°C)', 'Net Surface Energy Balance Amly (W m$^{-2}$)']
cmapamly = [ plt.cm.RdBu_r, plt.cm.RdBu_r, plt.cm.RdBu_r, plt.cm.RdBu_r]

cmapSIE = ListedColormap(["#C4E9FF", "magenta", "darkgreen", "white", "0.8"])
cmapSIA = ListedColormap([(46/255,70/255,226/255),(124/255,211/255,224/255),(68/255,148/255,45/255),(255/255,241/255,81/255),(225/255,58/255,46/255),'0.8','0.5'])


extent = [-3800000,3700000,-5300000,5800000]

prj = ccrs.NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)
eraprjpath = path+'/ERA5/DailyClimatology_1981-2010/ERA5_SLP_DailyClimatology_1981-2010.nc'
psnprjpath = path+'/Projections/psn_projection.nc'
siaprjpath = siapath+"/iceage_nh_12.5km_20110101_20111231_v4.1.nc"
simprjpath = simpath+"/icemotion_daily_nh_25km_20110101_20111231_v4.1.nc"
sstprjpath = sstpath+"/OISST_DailyClimatology_1981-2010.nc"

outpath = '/home/alex/OneDrive/CMIP6 OWP/Sea Ice Extent Expansion Pauses/Figures/AtmosphereAnomalies_V2'

'''*******************************************
Main Analysis
*******************************************'''

# Load Events
pdf = pd.read_csv(eventpath)
pdf = pdf.loc[(pdf['Use'] == 1) & (pdf['SYear'] > 1980)]
pdf = pdf.reset_index()
# Load projections
era5prj = xr.open_dataset(eraprjpath)
psnprj = xr.open_dataset(psnprjpath)
siaprj = xr.open_dataset(siaprjpath)
simprj = xr.open_dataset(simprjpath)
sstprj = xr.open_dataset(sstprjpath)

regridder = xe.Regridder(era5prj, psnprj, 'bilinear')
siatopsn = xe.Regridder(siaprj,psnprj, method='nearest_s2d')
simtopsn = xe.Regridder(simprj,psnprj, method='bilinear', ignore_degenerate=True)
ssttopsn = xe.Regridder(sstprj,psnprj, method='bilinear', periodic=True)

outpngs = md.listdir(outpath)

for e in range(len(pdf)):
    
    ### Dates ###
    enddate = list(pdf.loc[e,['year','month','day']].astype(int)) + [0,0,0]
    startdate = md.timeAdd(enddate,[0,0,-6,0,0,0])
    

    cenddate = md.timeAdd(enddate,[0,0,15,0,0,0])
    cstartdate = md.timeAdd(startdate,[0,0,-15,0,0,0])
    
    doyend = md.daysBetweenDates([1,1,0],[1,cenddate[1],cenddate[2]])
    doystart = md.daysBetweenDates([1,1,0],[1,cstartdate[1],cstartdate[2]])
    if doyend < doystart:
        doyend2 = 365
        doystart2 = 1
        
    # Dates as Strings
    enddatestr = str(enddate[0])+"-"+md.dd[enddate[1]-1]+"-"+md.dd[enddate[2]-1]
    startdatestr = str(startdate[0])+"-"+md.dd[startdate[1]-1]+"-"+md.dd[startdate[2]-1]
    
    outpng = "Fig6_AtmosphericSetup_"+enddatestr+"_V2.png"
    
    if outpng in outpngs:
        print(enddatestr + " already complete -- skipping")

    else:
        print(enddatestr)
        print("Load Sea Ice")
        
        ### SEA ICE VARIABLES ###
        # Concentration
        sicnc1 = xr.open_dataset(sicpath+"/"+str(startdate[0])+"/"+[f for f in md.listdir(sicpath+"/"+str(startdate[0])) if ''.join(startdatestr.split('-')) in f][0])
        sicarr1 = sicnc1['cdr_seaice_conc'][0,:,:].data
        sicarr1[ (sicarr1 > 1) | (psnprj['reg0780'] > 20) | (psnprj['reg0780'] < 0)] = np.nan
        
        sicnc2 = xr.open_dataset(sicpath+"/"+str(enddate[0])+"/"+[f for f in md.listdir(sicpath+"/"+str(enddate[0])) if ''.join(enddatestr.split('-')) in f][0])
        sicarr2 = sicnc2['cdr_seaice_conc'][0,:,:].data
        sicarr2[ (sicarr2 > 1) | (psnprj['reg0780'] > 20) | (psnprj['reg0780'] < 0)] = np.nan
        
        siclist = []
        date = startdate + []
        while date != md.timeAdd(enddate,[0,0,1]):
            sicnc = xr.open_dataset(sicpath+"/"+str(date[0])+"/"+[f for f in md.listdir(sicpath+"/"+str(date[0])) if str(date[0])+md.dd[date[1]-1]+md.dd[date[2]-1] in f][0])
            arr = sicnc['cdr_seaice_conc'][0,:,:].data
            arr[(arr > 1) | (psnprj['reg0780'] > 20) | (psnprj['reg0780'] < 0)] = np.nan
            siclist.append(arr)
            
            date = md.timeAdd(date,[0,0,1])
        sicarr = np.array(siclist)
        
        sie = np.where( (sicarr1 >= 0.15) & (sicarr2 >= 0.15), 3, np.where(sicarr2 >= 0.15, 2, np.where(sicarr1 >= 0.15, 1, np.where(np.isfinite(sicarr1*sicarr2), 0, 5))))
    
        # Age
        try:
            sianc = xr.open_dataset(siapath+"/iceage_nh_12.5km_"+str(startdate[0])+"0101_"+str(startdate[0])+"1231_v4.1.nc")
            siaarr = sianc['age_of_sea_ice'].data[int(np.round(md.daysBetweenDates([startdate[0],1,1],startdate)/7,0)),:,:].astype(float)
            siaarr[(siaarr > 4) & (siaarr < 20)] = 5
            siaarr[siaarr == 20] = 6
            siaarr[siaarr == 21] = 7
            siaarr[(siaarr < 1)] = np.nan
            sia = siatopsn(siaarr)
        except:
            sia = np.where(sie == 5, 6, np.where(sie == 0, np.nan, 7) )
            
        # Convergence
        simnc = xr.open_dataset(simpath+"/icemotion_daily_nh_25km_"+str(startdate[0])+"0101_"+str(startdate[0])+"1231_v4.1.nc")  
        simnc = simnc.sel(time=slice(startdatestr,enddatestr))

        if enddate[0] != startdate[0]:
            simnc2 = xr.open_dataset(simpath+"/icemotion_daily_nh_25km_"+str(enddate[0])+"0101_"+str(enddate[0])+"1231_v4.1.nc")  
            simnc2 = simnc2.sel(time=slice(startdatestr,enddatestr))
            simnc = xr.concat((simnc,simnc2),dim='time')
        
        # Use metpy to calculate divergence
        div = divergence(simnc['u']*units('cm/s'), simnc['v']*units('cm/s'), dx=2500000*units('cm'), dy=-2500000*units('cm')) * 86400*units('s/day')
        
        # Remove units
        div = xr.DataArray(div.values, coords={'time':div.time.values,
                           'y':div.y.values, 'x':div.x.values},
                           dims=['time','y','x'])
        
        # Reproject and caluclate area convergence
        vconv = -1 * np.mean(simtopsn(div)*units('1/day') * sicarr * psnprj['area']*units('km')*units('km'), axis=0) # -1 converts divergence to convergence
    
        # Rotate from x y to u v (for plotting)
        inu = np.mean(simnc['u'],axis=0)*np.cos(-1*np.pi*simnc['longitude']/180) + np.mean(simnc['v'],axis=0)*np.sin(-1*np.pi*simnc['longitude']/180)
        inv = np.mean(simnc['v'],axis=0)*np.cos(-1*np.pi*simnc['longitude']/180) - np.mean(simnc['u'],axis=0)*np.sin(-1*np.pi*simnc['longitude']/180)
       
        # Regrid to psn grid (for plotting - step 1)
        ## Also masks landmasses to be NaN
        uarr, varr = simtopsn(inu), simtopsn(inv)
       
        # Rotate from u v to psn x y (for plotting - step 2)
        # uirot, virot = md.rotateCoordsAroundOrigin(uarr, varr, (psnprj['lon']-(-45))*np.pi/-180)
        uirot, virot = md.rotateCoordsAroundOrigin(simtopsn(np.mean(simnc['u'],axis=0)), simtopsn(np.mean(simnc['v'],axis=0)), (45)*np.pi/-180)
    
        print("Load ERA5")
        
        ### ERA5 Variables ###
        values, clims = [], []
        # Load  Values
        for i in range(len(era5path)):
            ncf1 = xr.open_dataset(era5path[i]+"/"+e5filename[i]+"/ERA5_"+e5filename[i]+"_"+tstep[i]+"_"+str(enddate[0])+md.dd[enddate[1]-1]+".nc")
            
            if startdate[1] != enddate[1]:
                ncf1a = xr.open_dataset(era5path[i]+"/"+e5filename[i]+"/ERA5_"+e5filename[i]+"_"+tstep[i]+"_"+str(startdate[0])+md.dd[startdate[1]-1]+".nc")
                ncf1 = xr.concat((ncf1a,ncf1),dim='time')
                del ncf1a
                
            values.append( ncf1.sel(time=slice(startdatestr,enddatestr)).mean(dim='time') )
            
            #Load SLP Climatology
            ccf = xr.open_dataset(era5path[i]+"/DailyClimatology_1981-2010/ERA5_"+e5filename[i]+"_DailyClimatology_1981-2010.nc")
            
            if doyend > doystart:
                ccf1 = ccf.sel(dayofyear=slice(doystart,doyend))
            else:
                ccf1 = ccf.sel(dayofyear=slice(doystart,doyend2))
                ccf2 = ccf.sel(dayofyear=slice(doystart2,doyend))
                ccf1 = xr.concat((ccf1,ccf2),dim='dayofyear')
                del ccf2
            
            clims.append( ccf1.mean(dim='dayofyear') )
            del ccf, ccf1
        
        print("Load SST")
            
        ### SST ###
        sst = xr.open_dataset(sstpath+"/Daily/sst.day.mean."+str(startdate[0])+".nc")
        if startdate[0] != enddate[0]: 
            sst1 = xr.open_dataset(sstpath+"/Daily/sst.day.mean."+str(enddate[0])+".nc")
            sst = xr.concat((sst,sst1),dim='time')
            del sst1
        sst = sst.sel(time=slice(startdatestr,enddatestr)).mean(dim='time')
        
        csst = xr.open_dataset(sstpath+"/OISST_DailyClimatology_1981-2010.nc")
        
        if doyend > doystart:
            csst1 = csst.sel(dayofyear=slice(doystart,doyend))
        else:
            csst1 = csst.sel(dayofyear=slice(doystart,doyend2))
            csst2 = csst.sel(dayofyear=slice(doystart2,doyend))
            csst1 = xr.concat((csst1,csst2),dim='dayofyear')
            del csst2
        
        csst = csst1.mean(dim='dayofyear')
        del csst1    
        
        print('Start Making Figure')
        
        ### Start Figure ###
        fig = plt.figure(figsize=(7.5,5.25))
        axs = [ fig.add_subplot(2,3,a+1, projection=prj) for a in range(6)]
            
        # Plot One: Sea Ice Age + SST
        a = 1
        if len(np.unique(sia[np.isfinite(sia)])) == 2:
            pc0 = axs[a].pcolormesh(psnprj['x'], psnprj['y'], sia, cmap=ListedColormap(['0.8','0.5']), transform=prj)
            cb = fig.colorbar(pc0,ax=axs[a],orientation='vertical', ticks=[6.25,6.75])
            cb.ax.set_yticklabels(['Land', 'Age\nUnknown'], fontsize=7)
        else:
            pc0 = axs[a].pcolormesh(psnprj['x'], psnprj['y'], sia, cmap=cmapSIA, transform=prj)
            cb = fig.colorbar(pc0,ax=axs[a],orientation='vertical', ticks=np.arange(1.42,6.9,0.86))
            cb.ax.set_yticklabels(['≤ 1', '1 - 2', '2 - 3', '3 - 4', '> 4', 'Land', 'Mask'], fontsize=7)
           
        c0 = axs[a].contour(psnprj['x'][::4], psnprj['y'][::4], ssttopsn(sst)['sst'].data[::4,::4], levels=np.arange(0,20,4), transform=prj, colors='k', linewidths=0.5, )
        axs[a].clabel(c0, inline=True, fontsize=5, fmt='%2d')
        axs[a].set_title('b. Sea Ice Age (years)\n& Sea-Surface Temperature', fontsize=7)
    
        # Plot Two: Sea Ice Growth/Loss + SST Amly
        a = 0
        pc1 = axs[a].pcolormesh(psnprj['x'], psnprj['y'], sie, cmap=cmapSIE, transform=prj)
        cb = fig.colorbar(pc1,ax=axs[a],orientation='vertical', ticks=np.arange(0.5,5.5,1))
        cb.ax.set_yticklabels(['Open\nWater', 'Ice\nLoss', 'Ice\nGain', 'Ice\nCover', 'Land'], fontsize=7)
        axs[a].set_title('a. Sea Ice Tendency\n' + startdatestr + " to " + enddatestr, fontsize=7)
    
        # Plot Six: Sea Ice Convergence and Divergence & Motion
        a = 2
        
        plotteramly = ssttopsn( sst - csst )   
        cf0 = axs[a].contourf(psnprj['x'], psnprj['y'], plotteramly['sst'].data, levels=md.divbreaks(5,1), cmap=plt.cm.RdBu_r, transform=prj, extend='both')
        cf1 = axs[a].contourf(psnprj['x'], psnprj['y'], vconv, levels=md.divbreaks(30,5), cmap=plt.cm.PRGn_r, transform=prj, extend='both')
    
        axs[a].add_feature(cfeature.LAND, zorder=6, facecolor='0.8')
        axs[a].set_extent(extent, prj)
        axs[a].contour(psnprj['x'], psnprj['y'], sicarr1, levels=[0.15], transform=prj, colors='magenta', linewidths=1.2)
    
        vec1 = axs[a].quiver(psnprj['x'], psnprj['y'], uirot.data, virot.data, transform=prj, regrid_shape=26, width=0.005, scale=150, color='k')
    
        cb = fig.colorbar(cf1,ax=axs[a],orientation='vertical', ticks=md.divbreaks(30,5))
        cb.ax.set_yticklabels(md.divbreaks(30,5), fontsize=7)
        cb.set_label('Sea Ice Area Convergence (km$^{2}$ per day)', fontsize=7)
        axs[a].set_title('c. Sea Ice Motion & Convergence\nand Sea Surface Temperature Anomaly',fontsize=7)
    
        # Plot Three: SLP & 10-m Wind
        i, a = 0, 5
        plotter = regridder(values[i])['msl']
        cf2 = axs[a].contourf(psnprj['x'],psnprj['y'], plotter, levels=breaks[i], cmap=cmap[i], extend='both', transform=prj)
        axs[a].contour(psnprj['x'], psnprj['y'], sicarr1, levels=[0.15], transform=prj, colors='yellow', linewidths=1.2)
        axs[a].add_feature(cfeature.COASTLINE, zorder=6)
        axs[a].set_extent(extent, prj)
        cb = fig.colorbar(cf2,ax=axs[a],orientation='vertical', ticks=breakstick[i])
        cb.ax.set_yticklabels(breakslab[i], fontsize=7)
        cb.set_label(colorbarlabel[i], fontsize=7)
        axs[a].set_title('f. Sea-Level Pressure & 10-m Wind', fontsize=7)
        
        i = 1
        uvpsn = regridder(values[i])
        urot, vrot = md.rotateCoordsAroundOrigin(uvpsn['u10'], uvpsn['v10'], (psnprj['lon']-(-45))*np.pi/-180)
        vec2 = axs[a].quiver(psnprj['x'], psnprj['y'], urot.data, vrot.data, transform=prj, regrid_shape=21, width=0.005, scale=150)
    
        # Plot Four: 925-hPa Temperature & Amly
        i, a = 2, 3 # index, axis
        plotter = regridder(values[i])[vvalue[i]] - 273.15
        plotteramly = plotter - (regridder(clims[i])[vvalue[i]]-273.15)
        # cf2 = axs[a].contourf(psnprj['x'],psnprj['y'], plotter, levels=breaks[i], cmap=cmap[i], extend='both', transform=prj)
        # axs[a].contour(psnprj['x'][::2],psnprj['y'][::2], plotteramly[::2,::2], levels=breaksamly[i], extend='both', transform=prj, colors='0.8', linewidths=0.5)
        cf3 = axs[a].contourf(psnprj['x'],psnprj['y'], plotteramly, levels=breaksamly[i], cmap=cmapamly[i], extend='both', transform=prj)
        c3 = axs[a].contour(psnprj['x'][::2],psnprj['y'][::2], plotter[::2,::2], levels=breaks[i], extend='both', transform=prj, colors='0.3', linewidths=0.6)
        axs[a].contour(psnprj['x'], psnprj['y'], sicarr1, levels=[0.15], transform=prj, colors='yellow', linewidths=1.2)
        axs[a].clabel(c3, inline=True, fontsize=5, fmt='%2d')
        axs[a].set_title('d. 925-hPa Temperature & Anomaly', fontsize=7)
    
        axs[a].add_feature(cfeature.COASTLINE, zorder=6)
        axs[a].set_extent(extent, prj)
        # cb = fig.colorbar(cf2,ax=axs[a],orientation='vertical', ticks=breakstick[i])
        # cb.ax.set_yticklabels(breakslab[i], fontsize=7)
        # cb.set_label(colorbarlabel[i], fontsize=8)   
        cb = fig.colorbar(cf3,ax=axs[a],orientation='vertical', ticks=breakstickamly[i])
        cb.ax.set_yticklabels(breakslabamly[i], fontsize=7)
        cb.set_label(colorbarlabelamly[i], fontsize=7)   
        
        # Plot Five: Net SEB + Amly
        i, a = 3, 4
        plotter = regridder(values[i])[vvalue[i]]
        plotteramly = plotter - regridder(clims[i])[vvalue[i]]
        cf4 = axs[a].contourf(psnprj['x'],psnprj['y'], plotteramly, levels=breaksamly[i], cmap=cmapamly[i], extend='both', transform=prj)
        axs[a].add_feature(cfeature.COASTLINE, zorder=6)
        axs[a].set_extent(extent, prj)
        cb = fig.colorbar(cf4,ax=axs[a],orientation='vertical', ticks=breakstickamly[i])
        cb.ax.set_yticklabels(breakslabamly[i], fontsize=7)
        cb.set_label(colorbarlabelamly[i], fontsize=7)   
        axs[a].set_title('e. Net Surface Energy Balance & Anomaly', fontsize=7)
            
        axs[a].contour(psnprj['x'][::4],psnprj['y'][::4], plotter[::4,::4], levels=breaks[i], extend='both', transform=prj, colors='0.3', linewidths=0.6)  
        axs[a].contour(psnprj['x'], psnprj['y'], sicarr1, levels=[0.15], transform=prj, colors='yellow', linewidths=1)
    
        plt.tight_layout(rect=[0,0,0.95,1])
        
        # Quivers
        axs[a].quiverkey(vec1, 0.65, 0.95, 10, r'10 $\frac{cm}{s}$', coordinates='axes', labelpos='E', fontproperties={'size':7})
        axs[a].quiverkey(vec2, 1.5, 0, 10, r'10 $\frac{m}{s}$', coordinates='axes', labelpos='E', fontproperties={'size':7})
    
        # Extra Color Bar
        cbax = fig.add_axes([0.93, 0.51, 0.012, 0.42])
        cb0 = fig.colorbar(cf0,cbax,orientation='vertical', ticks=md.divbreaks(5,1))
        cb0.ax.set_yticklabels(md.divbreaks(5,1), fontsize=7)
        cb0.set_label('SST Anomaly (°C)', fontsize=7)
        
        plt.savefig(outpath+"/"+outpng, dpi=300)
    
        