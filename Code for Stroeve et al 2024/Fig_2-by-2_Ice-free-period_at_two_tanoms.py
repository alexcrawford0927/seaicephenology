#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 09:10:25 2023

@author: acrawfora
"""

import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import CycloneModule_13_2 as md
import matplotlib.patheffects as pe
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

exp = 'historical-ssp585'
YY = '1979-2021'
tcoords = [2,4]
xlabel = ['Bias-Adjusted Average','Weighted Average']
v = 0
varlist = ['opc','lrd','fad','cip']
varname = ['Continuous Ice-Free Period (Days)', 'Retreat Day', 'Advance Day', 'Continuous Ice-Covered Period']
varname2 = ['IceFreePeriod','RetreatDay','AdvanceDay','IceCoveredPeriod']

inpath = '/Volumes/Cassandra/CMIP6/SeaIce/ThicknessPhenology'
path1 = inpath+'/CMIP6_bias-corrected-by-t2mamly_'+exp+'_baselineyears_'+YY+'_multimodelmean.nc'
path2 = inpath+'/CMIP6_weightedaverage_tanombins_'+exp+'_baselineyears_'+YY+'.nc'
# inpath+'/CMIP6_multimodelmean_historical_climatology_'+YY+'.nc'
ease2path = "/Volumes/Miranda/Projections/EASE2_N0_25km_Projection.nc"
outpath = "/Users/acrawfora/OneDrive - University of Manitoba/CMIP6 OWP/Hudson Bay Figures/ThicknessPhenologyFigures/"+varname2[v]+"Comparison_T"+'-'.join([str(t) for t in tcoords])+"_"+"-".join(xlabel)+".png"

# Plotting Variables
prjdata = ccrs.LambertAzimuthalEqualArea(central_longitude=0, central_latitude=90)
prjplot = ccrs.LambertAzimuthalEqualArea(central_longitude=-80, central_latitude=90)
bbox = [-1000000,500000,-2700000,-4400000] # left, right, top, bottom

levs = [ np.arange(165,315+15,15), [135,142,152,159,166,173,182,189,196], [335,342,349,356,366,373,380,387,397]]
# levs = [ np.arange(90,240+15,15) , [152,159,166,173,182,189,196,203,213],] # 
levsdiff = np.arange(14,49+1,7) # md.divbreaks(49,7)
ticks = [levs[0], levs[1][::2], levs[2][::2]] #np.arange(90,180+15,15)
ticksdiff = np.arange(14,49+1,7) # md.divbreaks(42,14)
ticklabels = [ levs[0], ['May 15','Jun 1','Jun 15','Jul 1','Jul 15'], ['Dec 1','Dec 15','Jan 1','Jan 15','Feb 1'] ]
# ticklabels = [ levs[0], ['Jun 1','Jun 15','Jul 1','Jul 15','Aug 1'], ['Dec 1','Dec 15','Jan 1','Jan 15','Feb 1'] ]

lines = [(63,-91,63,-88),(63,-88,56,-88),(60,-88,60,-77)]

##### MAIN ANALYSIS #####
ease2 = xr.open_dataset(ease2path)
ds1 = xr.open_dataset(path1)
ds2 = xr.open_dataset(path2)

fig = plt.figure(figsize=(7.5,8))

a = 1
for t in range(len(tcoords)):
    
    if varlist[v] == 'opc':
        arrs = [ 365 - ds1.sel(tanom=tcoords[t])['cip'].data,
                365 - ds2.sel(tanom=tcoords[t])['cip'].data]
    elif varlist[v] == 'lrd':
        arrs = [ ds1.sel(tanom=tcoords[t])[varlist[v]].data,
                ds2.sel(tanom=tcoords[t])[varlist[v]].data - 365]
    else:
        arrs = [ds1.sel(tanom=tcoords[t])[varlist[v]].data,
        ds2.sel(tanom=tcoords[t])[varlist[v]].data]
        
    # Band air to adjsut opcarr2 to have the same extent as opcarr1
    arrs[1][np.isnan(arrs[0])] = np.nan
    
    ax = [ fig.add_subplot( 2, 2, i+a, projection=prjplot) for i in range(2) ]
    
    for i in range(len(ax)):
        ax[i].add_feature(cfeature.COASTLINE)
        ax[i].set_extent(bbox, crs=prjplot)
            
        ax[i].contourf(ds1['x'], ds1['y'], np.where(np.isfinite(arrs[i]), np.nan, 1), [0,1,2,3,4], cmap=plt.cm.gray_r,transform=prjdata)
        cf1 = ax[i].contourf(ds1['x'], ds1['y'], arrs[i], levs[v], cmap=plt.cm.viridis, transform=prjdata, extend='both')
    
        gdl = ax[i].gridlines(xlocs=range(-180,180,15),ylocs=range(0,90,5),linestyle='dashed',linewidth=0.5,draw_labels=True)
        gdl.top_labels = False
        gdl.right_labels = False
        gdl.xformatter = LongitudeFormatter()
        gdl.yformatter = LatitudeFormatter()
        gdl.xlabel_style = {'size':8}
        gdl.ylabel_style = {'size':8}
    
        if i > 2:
            gdl.bottom_labels = False  
            
        txt = ax[i].annotate(md.abc[a+i-1],xy=(0.02,0.98),xycoords='axes fraction',ha='left',va='top',weight='bold')
        txt.set_path_effects([pe.withStroke(linewidth=3, foreground='w')])
    
        if a > 2:
            gdl.bottom_labels = False 
            
        # Annotate
        if a == 1:
            ax[i].annotate(xlabel[i], xy=(0.50, 1.05), xycoords='axes fraction',va='bottom',ha='center',weight='bold', size=9)
        
        # Add Study Area Lines
        for l in range(len(lines)):
            if lines[l][0] == lines[l][2]:
                ax[i].plot( np.arange(np.min([lines[l][1],lines[l][3]]),np.max([lines[l][1],lines[l][3]])+0.5,0.5), np.repeat(lines[l][0],np.ceil(1+(abs(lines[l][3]-lines[l][1]))/0.5)), color='w', transform=ccrs.Geodetic())
            else:
                ax[i].plot(np.repeat(lines[l][1],np.ceil(1+(abs(lines[l][2]-lines[l][0]))/0.5)), np.arange(np.min([lines[l][2],lines[l][0]]),np.max([lines[l][2],lines[l][0]])+0.5,0.5), color='w', transform=ccrs.Geodetic())
        
        # Labels
        txt = [ ax[i].annotate('SHB', xy=(0.50,0.50), xycoords='axes fraction',va='center',ha='center',size=8),
               ax[i].annotate('WHB', xy=(0.25,0.65), xycoords='axes fraction',va='center',ha='center',size=8)]
        [ t.set_path_effects([pe.withStroke(linewidth=2, foreground='w')]) for t in txt ]
        

    a += 2

    ax[0].annotate("Warming Level: " + str(tcoords[t])+"°C", xy=(-0.2, 0.50), rotation=90, 
                 xycoords='axes fraction',va='center',ha='center',size=9, weight='bold')
    


plt.tight_layout(rect=[0.03,0.1,1,0.99]) # left, bottom, right, top

### Color Bars ###
cbar_ax1 = fig.add_axes([0.25, 0.07, 0.50, 0.015]) # xloc, yloc, xsize, ysize
cbar1 = fig.colorbar(cf1,cbar_ax1,orientation='horizontal',ticks=ticks[v])
cbar1.set_label(varname[v],size=8)
cbar1.set_ticklabels(ticklabels[v])
cbar1.ax.tick_params(labelsize=8)

# plt.savefig(outpath, dpi=300)


