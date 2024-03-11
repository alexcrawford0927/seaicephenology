"""
Author: Alex Crawford
Date Created: 16 Oct 2020
Date Modified: 16 Oct 2020
Purpose: Plot maps of the day on which a variable exceeds some value
variable = OPC, LRD, FAD
value = (90, 180, 270), (May 1, June 1, July 1, Aug 1), (Oct 31, Nov 30, Dec 31, Jan 31, Feb 28)
"""

'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import numpy as np
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


np.seterr(all='ignore')

'''#####################
 Declare Variables
#####################'''

experiment1 = 'historical'
experiments2 = ['ssp126','ssp245','ssp585']
var = 'opc'
varName = "Open Water Period"
thress = [90,180,270]
ct = 80

outpath1 = "/Volumes/Troilus/CMIP6/RegionalStats/ExceedanceDates/"
outpath2 = "/Volumes/Troilus/CMIP6/RegionalStats/ExceedanceTemp/"
figpath = "/Volumes/Troilus/CMIP6/Figures/ExceedanceMaps"

# Inputs for reprojection
breaks1 = np.arange(1980,2100,10)
breaks2 = np.arange(0.5,6,0.5)

abc = ['a','b','c','d','e','f','g','h','i','j','k','l']
bb2 = [30,-45,30,135]
xsize, ysize = 50000, -50000 # in meters
nx, ny = int(180*(100000/xsize)), int(180*(100000/xsize)) # number of grid cells
lon_0 = 0 # Central Meridian (which longitude is at the 6 o'clock position)
lat_0 = 90 # Reference Latitude (center of projection)
box = dict(boxstyle='square,pad=0.3',fc='0.95',ec='0.95',lw=0.5) # white rectangle with white outline
box2 = dict(boxstyle='square,pad=0.3',fc='0.85',ec='0.85',lw=0.5) # white rectangle with white outline

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

# Load output projection data
projnc = nc.Dataset("/Volumes/Troilus/Projections/EASE2_N0_"+str(int(xsize/1000))+"km_Projection.nc")
outlons = projnc['lon'][:]
outlats = projnc['lat'][:]

fig = plt.figure(figsize=(8,9.5))

i=0
# BY DATE
for experiment2 in experiments2:
    for thres in thress:
        ### Find Multi-Model Mean ###
        ncf = nc.Dataset(outpath1+"/C"+str(ct)+"/ExceedanceDates_C"+str(ct)+"_"+var+str(thres)+"_"+experiment2+".nc", 'r')
        mmm = np.apply_along_axis(np.nanmean,0,ncf[var][:].data)
        mmm[(outlats < 50) & (outlons > 25) & (outlons < 60)] = np.nan

        ### Plot Multi-Model Mean ###
        ax = fig.add_subplot(4,3,i+1)
        mp = Basemap(projection='laea',lat_0=lat_0,lon_0=lon_0,\
            llcrnrlat=bb2[0], llcrnrlon=bb2[1],urcrnrlat=bb2[2], urcrnrlon=bb2[3],resolution='c')
        
        cf1 = mp.contourf(outlons,outlats,mmm,breaks1,latlon=1,extend='both')
        mp.drawmeridians(np.arange(0,360,30),color="0.4",linewidth=0.8)
        mp.drawparallels(np.arange(0,90,10),color="0.4",linewidth=0.8)
        mp.fillcontinents(color='0.95',lake_color='0.95')
        mp.drawmapboundary(fill_color='0.95' )
        
        ax.annotate(abc[i],xy=(0.03,0.97), xycoords='axes fraction', va='top', ha= 'left', size=9, weight='bold',bbox=box)
        i += 1

# BY TEMPERATURE
for thres in thress:        
    ### Find Multi-Model Mean ###
    ncf = nc.Dataset(outpath2+"/C"+str(ct)+"/ExceedanceTemp_C"+str(ct)+"_"+var+str(thres)+"_"+experiment2+".nc", 'r')
    mmm = np.apply_along_axis(np.nanmean,0,ncf[var][:].data)
    mmm[(outlats < 50) & (outlons > 25) & (outlons < 60)] = np.nan
    
    ### Plot Multi-Model Mean ###
    ax = fig.add_subplot(4,3,i+1)
    mp = Basemap(projection='laea',lat_0=lat_0,lon_0=lon_0,\
        llcrnrlat=bb2[0], llcrnrlon=bb2[1],urcrnrlat=bb2[2], urcrnrlon=bb2[3],resolution='c')
    
    cf2 = mp.contourf(outlons,outlats,mmm,breaks2,latlon=1,cmap=plt.cm.magma,extend='both')
    mp.drawmeridians(np.arange(0,360,30),color="0.4",linewidth=0.8)
    mp.drawparallels(np.arange(0,90,10),color="0.4",linewidth=0.8)
    mp.fillcontinents(color='0.85',lake_color='0.85')
    mp.drawmapboundary(fill_color='0.85' )
    
    ax.annotate(abc[i],xy=(0.03,0.97), xycoords='axes fraction', va='top', ha= 'left', size=9, weight='bold',bbox=box2)
    # txt = plt.text(8900000, 100000, str(thress[i]) + " days",fontsize=9, ha='right',va='bottom',color='k')
    # txt.set_path_effects([PathEffects.withStroke(linewidth=3, foreground='0.95')])
    i += 1

plt.tight_layout(rect=(0.03,0,0.86,0.98),pad=0.8,w_pad=0.6) # [left, bottom, right, top]
plt.annotate('',xy=(0.02,0.25),xytext=(0.855,0.25),xycoords="figure fraction",arrowprops=dict(arrowstyle='-',ls='dashed'))
plt.annotate('',xy=(.855,0.42),xytext=(0.855,0.248),xycoords="figure fraction",arrowprops=dict(arrowstyle='-',ls='dashed'))
plt.annotate('',xy=(.854,0.42),xytext=(0.985,0.42),xycoords="figure fraction",arrowprops=dict(arrowstyle='-',ls='dashed'))

cbar_ax1 = fig.add_axes([0.87, 0.55, 0.025, 0.34])
cbar_ax2 = fig.add_axes([0.87, 0.05, 0.025, 0.34])
cbar1 = fig.colorbar(cf1,cbar_ax1,orientation='vertical')
cbar2 = fig.colorbar(cf2,cbar_ax2,orientation='vertical')
cbar1.set_label("First Year of Long-Term Exceedance",fontsize=9)
cbar2.set_label("Exceedance Global Temperature Anomaly (°C)",fontsize=9)

plt.annotate("SSP126",(0.02,0.87),rotation='vertical',va='center',fontsize=9,xycoords="figure fraction")
plt.annotate("SSP245",(0.02,0.62),rotation='vertical',va='center',fontsize=9,xycoords="figure fraction")
plt.annotate("SSP585",(0.02,0.37),rotation='vertical',va='center',fontsize=9,xycoords="figure fraction")
plt.annotate("Temperature Sensitivity",(0.03,0.14),rotation='vertical',va='center',fontsize=9,xycoords="figure fraction")

plt.annotate("Days: "+str(thress[0]),(0.17,0.98),rotation='horizontal',ha='center',fontsize=9,xycoords="figure fraction")
plt.annotate("Days: "+str(thress[1]),(0.45,0.98),rotation='horizontal',ha='center',fontsize=9,xycoords="figure fraction")
plt.annotate("Days: "+str(thress[2]),(0.72,0.98),rotation='horizontal',ha='center',fontsize=9,xycoords="figure fraction")

plt.savefig(figpath+"/ExceedanceDateTemp_CT"+str(ct)+"_AllSSPs_"+str(thress[0])+"-"+str(thress[1])+"-"+str(thress[2])+"_pallette3.eps",dpi=300)
