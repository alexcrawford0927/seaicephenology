"""
Author: Alex Crawford
Date Created: 3 Jun 2020
Date Modified: 3 Jun 2020
Purpose: Creates a numpy array with categorical values for each sea ice region.
"""

'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pandas as pd
import netCDF4 as nc
import os

'''#####################
 Declare Variables
#####################'''

lonmin, lonmax = -360, 360
latmin, latmax = -90,90

var = 'siconc'
frequency = 'day'
tableid = 'SIday' # may be identical to frequency, but not always
experiment = 'ssp126'

# Plotting Variables
bb = [35,-45,35,135] # [ll_lat, ll_lon, ur_lat, ur_lon]
lon_0 = 0 # Central Meridian
lat_0 = 90 # Latitude of Origin
levs = list(range(2,10,1)) #[0,100,200,400,800,1600,3200]

### Path Variables ###
inpath = "/Volumes/Prospero/CMIP6/SmoothedMA5/"+experiment
outpath = "/Volumes/Troilus/CMIP6/RegionMasks"

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

##################
# Create regions #
##################
# Load file names
files = os.listdir(inpath)
files = [f for f in files if f.startswith('.') == 0]
files.sort()

models = np.unique([f.split("_"+experiment)[0] for f in files])

cmplmods = os.listdir(outpath)
cmplmods = [f[8:-4] for f in cmplmods if f.startswith('.') == 0]
newmods = [f for f in models if f not in cmplmods]

for mod in newmods:
    # 1) Load lat, lon, and a file to use as the region base (e.g., land/sea mask)
    modfiles = [f for f in files if mod+"_"+experiment in f]
    os.chdir(inpath+"/"+modfiles[0])
    ncf = nc.Dataset([f for f in os.listdir() if f.startswith(var)][0])   
    
    lons = ncf['lon'][:].data
    lats = ncf['lat'][:].data
    
    # Make sure all lats and lons are -90 to 90 and -180 to 180
    lons[(lons < lonmin) | (lons > lonmax)] = np.nan
    lons = np.where(lons < -180, 360+lons, np.where(lons > 180, lons-360, lons))
    lats[(lats < latmin) | (lats > latmax)] = np.nan
    
    try: # Try with three dimensions first (t * y * x)
        regs = np.where(np.isfinite(ncf['siconc'][0,:,:].data), 1, 0)
    except: # If that doesn't work try 2-D (t * x)
        regs = np.where(np.isfinite(ncf['siconc'][0,:].data), 1, 0)
    
    # 2) Use lats and longs to assign region numbers
    
    # Okhotsk - 2
    regs[(lons > 126.75) & (lons <= 162.1) & (lats > 55.7) & (lats <= 64) & (regs == 1)] = 2
    regs[(lons > 161.15) & (lons <= 164.9) & (lats > 60.2) & (lats <= 64) & (regs == 1)] = 2
    regs[(lons > 126.75) & (lons <= 161.8) & (lats > 54.7) & (lats <= 55.7) & (regs == 1)] = 2
    regs[(lons > 126.75) & (lons <= 156.6) & (lats > 50.9) & (lats <= 54.7) & (regs == 1)] = 2
    regs[(lons > 126.75) & (lons <= 155.6) & (lats > 49.9) & (lats <= 50.9) & (regs == 1)] = 2
    regs[(lons > 126.75) & (lons <= 154.1) & (lats > 48.9) & (lats <= 49.9) & (regs == 1)] = 2
    regs[(lons > 126.75) & (lons <= 152.6) & (lats > 47.9) & (lats <= 48.9) & (regs == 1)] = 2
    regs[(lons > 126.75) & (lons <= 151.1) & (lats > 46.9) & (lats <= 47.9) & (regs == 1)] = 2
    regs[(lons > 126.75) & (lons <= 149.6) & (lats > 45.9) & (lats <= 46.9) & (regs == 1)] = 2
    regs[(lons > 126.75) & (lons <= 147.6) & (lats > 44.9) & (lats <= 45.9) & (regs == 1)] = 2
    regs[(lons > 126.75) & (lons <= 145.6) & (lats > 43.5) & (lats <= 44.9) & (regs == 1)] = 2
    regs[(lons > 126.75) & (lons <= 140.8) & (lats > 40) & (lats <= 43.9) & (regs == 1)] = 2
    
    # CAA - 14
    regs[(lons > -123) & (lons <= -60) & (lats > 66.6) & (lats <= 82) & (regs == 1)] = 14
    
    # Labrador - 6
    regs[(lons > -69.7) & (lons <= -45) & (lats > 45.0) & (lats <= 66.5) & (regs > 0)] = 6
    
    # CAO - 15
    regs[(lons > -55) & (lons <= 180) & (lats > 75) & (lats <= 90) & (regs > 0)] = 15
    regs[(lons > -181) & (lons <= -90) & (lats > 80) & (lats <= 90) & (regs > 0)] = 15
    regs[(lons > -90) & (lons <= -55) & (lats > 82) & (lats <= 90) & (regs > 0)] = 15
    regs[(lons > -114) & (lons <= -96.4) & (lats > 79) & (lats <= 90) & (regs > 0)] = 15
    regs[(lons > -123) & (lons <= -114) & (lats > 77.5) & (lats <= 90) & (regs > 0)] = 15
    regs[(lons > -181) & (lons <= -123) & (lats > 75) & (lats <= 90) & (regs > 0)] = 15
    
    # Beaufort - 13
    regs[(lons > -156.6) & (lons <= -123) & (lats > 65.5) & (lats <= 80) & (regs == 1)] = 13
    
    # Barents - 8
    regs[(lons > 16.5) & (lons <= 30) & (lats > 66.7) & (lats <= 80) & (regs > 0)] = 8
    regs[(lons > 30) & (lons <= 103.1) & (lats > 60) & (lats <= 80) & (regs > 0)] = 8
    
    # Kara - 9
    regs[(lons > 67.8) & (lons <= 103.1) & (lats > 60) & (lats <= 80) & (regs > 0)] = 9
    regs[(lons > 61.0) & (lons <= 67.9) & (lats > 60) & (lats <= 75.6) & (regs > 0)] = 9
    regs[(lons > 55.1) & (lons <= 67.9) & (lats > 70.8) & (lats <= 74.1) & (regs > 0)] = 9
    regs[(lons > 57.1) & (lons <= 67.9) & (lats > 74.1) & (lats <= 74.9) & (regs > 0)] = 9
    regs[(lons > 60.0) & (lons <= 67.9) & (lats > 74.9) & (lats <= 75.6) & (regs > 0)] = 9
    regs[(lons > 65.5) & (lons <= 67.9) & (lats > 75.6) & (lats <= 76.1) & (regs > 0)] = 9
    regs[(lons > 57.6) & (lons <= 67.9) & (lats > 70.5) & (lats <= 72.0) & (regs > 0)] = 9
    regs[(lons > 58.8) & (lons <= 67.9) & (lats > 70.1) & (lats <= 72.0) & (regs > 0)] = 9
    regs[(lons > 60.2) & (lons <= 67.9) & (lats > 69.7) & (lats <= 72.0) & (regs > 0)] = 9
    
    # Laptev - 10
    regs[(lons > 96.9) & (lons <= 103.1) & (lats > 79.1) & (lats <= 80) & (regs > 0)] = 10
    regs[(lons > 100.7) & (lons <= 103.1) & (lats > 78.3) & (lats <= 80) & (regs > 0)] = 10
    regs[(lons > 103.1) & (lons <= 142.5) & (lats > 60) & (lats <= 80) & (regs > 0)] = 10
    
    # East Siberian - 11
    regs[(lons > 142.5) & (lons <= 146) & (lats > 66) & (lats <= 79) & (regs > 0)] = 11
    regs[(lons > 146) & (lons <= 150) & (lats > 66) & (lats <= 78.5) & (regs > 0)] = 11
    regs[(lons > 150) & (lons <= 155) & (lats > 66) & (lats <= 78) & (regs > 0)] = 11
    regs[(lons > 155) & (lons <= 160) & (lats > 66) & (lats <= 77.5) & (regs > 0)] = 11
    regs[(lons > 160) & (lons <= 165) & (lats > 66) & (lats <= 77) & (regs > 0)] = 11
    regs[(lons > 165) & (lons <= 170) & (lats > 66) & (lats <= 76.5) & (regs > 0)] = 11
    regs[(lons > 170) & (lons <= 175) & (lats > 66) & (lats <= 76) & (regs > 0)] = 11
    regs[(lons > 175) & (lons <= 180) & (lats > 66) & (lats <= 75.5) & (regs > 0)] = 11
    
    # Chukchi - 12
    regs[(lons > -181) & (lons <= -156.6) & (lats > 66) & (lats <= 70) & (regs > 0)] = 12
    regs[(lons > -181) & (lons <= -156.6) & (lats > 70) & (lats <= 75) & (regs > 0)] = 12
    
    # Bering - 3
    regs[(lons > -181) & (lons <= -177) & (lats > 66) & (lats <= 67) & (regs > 0)] = 3
    regs[((lons > 161.7) | (lons <= -157.1)) & (lats > 57) & (lats <= 66) & (regs == 1)] = 3
    regs[((lons > 161.7) | (lons <= -159.1)) & (lats > 56.3) & (lats <= 57) & (regs == 1)] = 3
    regs[((lons > 163.1) | (lons <= -160.5)) & (lats > 55.7) & (lats <= 56.3) & (regs == 1)] = 3
    regs[((lons > 165.8) | (lons <= -161.9)) & (lats > 55.3) & (lats <= 55.7) & (regs == 1)] = 3
    regs[((lons > 167.7) | (lons <= -163.6)) & (lats > 54.7) & (lats <= 55.3) & (regs == 1)] = 3
    regs[((lons > 169.0) | (lons <= -166.3)) & (lats > 54.0) & (lats <= 54.7) & (regs == 1)] = 3
    regs[((lons > 171.9) | (lons <= -168.0)) & (lats > 53.4) & (lats <= 54.0) & (regs == 1)] = 3
    regs[((lons > 172.9) | (lons <= -169.0)) & (lats > 52.9) & (lats <= 53.4) & (regs == 1)] = 3
    regs[((lons > 175.4) | (lons <= -170.7)) & (lats > 52.6) & (lats <= 52.9) & (regs == 1)] = 3
    regs[((lons > 176.5) | (lons <= -172.2)) & (lats > 52.2) & (lats <= 52.6) & (regs == 1)] = 3
    regs[((lons > 177.5) | (lons <= -176.0)) & (lats > 52.0) & (lats <= 52.6) & (regs == 1)] = 3
    
    # Hudson - 4
    regs[(lons > -96.3) & (lons <= -64.8) & (lats > 50.6) & (lats <= 61.5) & (regs > 0)] = 4
    regs[(lons > -96.3) & (lons <= -66.5) & (lats > 61.5) & (lats <= 62.3) & (regs > 0)] = 4
    regs[(lons > -94.0) & (lons <= -69.4) & (lats > 62.3) & (lats <= 66.7) & (regs > 0)] = 4
    regs[(lons > -84.4) & (lons <= -69.4) & (lats > 66.7) & (lats <= 69.6) & (regs > 0)] = 4
    regs[(lons > -83.4) & (lons <= -76.0) & (lats > 69.6) & (lats <= 71.1) & (regs > 0)] = 4
    
    # Gulf of St Law - 5
    regs[(lons > -72.2) & (lons <= -56.2) & (lats > 50.7) & (lats <= 52.1) & (regs > 0)] = 5
    regs[(lons > -72.2) & (lons <= -57.1) & (lats > 48.0) & (lats <= 50.7) & (regs > 0)] = 5
    regs[(lons > -65.0) & (lons <= -53.0) & (lats > 45.0) & (lats <= 48.0) & (regs > 0)] = 5
    
    # Greenland Sea - 7
    regs[(lons > -45) & (lons <= 16.5) & (lats > 71) & (lats <= 80) & (regs > 0)] = 7
    regs[(lons > -45) & (lons <= 10) & (lats > 70.5) & (lats <= 71) & (regs > 0)] = 7
    regs[(lons > -45) & (lons <= 6) & (lats > 70) & (lats <= 70.5) & (regs > 0)] = 7
    regs[(lons > -45) & (lons <= 4) & (lats > 69) & (lats <= 70) & (regs > 0)] = 7
    regs[(lons > -45) & (lons <= 2) & (lats > 68.5) & (lats <= 69) & (regs > 0)] = 7
    regs[(lons > -45) & (lons <= 0) & (lats > 67.5) & (lats <= 68.5) & (regs > 0)] = 7
    regs[(lons > -45) & (lons <= -3) & (lats > 66.5) & (lats <= 67.5) & (regs > 0)] = 7
    regs[(lons > -45) & (lons <= -4) & (lats > 65) & (lats <= 66.5) & (regs > 0)] = 7
    regs[(lons > -45) & (lons <= -7) & (lats > 64) & (lats <= 65) & (regs > 0)] = 7
    regs[(lons > -45) & (lons <= -10) & (lats > 63) & (lats <= 64) & (regs > 0)] = 7
    regs[(lons > -45) & (lons <= -13) & (lats > 62) & (lats <= 63) & (regs > 0)] = 7
    regs[(lons > -45) & (lons <= -16) & (lats > 61) & (lats <= 62) & (regs > 0)] = 7
    regs[(lons > -45) & (lons <= -20) & (lats > 60.5) & (lats <= 61) & (regs > 0)] = 7
    regs[(lons > -45) & (lons <= -25) & (lats > 60) & (lats <= 60.5) & (regs > 0)] = 7
    regs[(lons > -45) & (lons <= -29) & (lats > 59.5) & (lats <= 60) & (regs > 0)] = 7
    regs[(lons > -45) & (lons <= -32) & (lats > 58.5) & (lats <= 59.5) & (regs > 0)] = 7
    regs[(lons > -45) & (lons <= -35) & (lats > 58) & (lats <= 58.5) & (regs > 0)] = 7
    
    # Baffin - 16
    regs[(lons > -71) & (lons <= -50) & (lats > 66.5) & (lats <= 77) & (regs > 0)] = 16
    regs[(lons > -79.3) & (lons <= -71) & (lats > 70.1) & (lats <= 77) & (regs > 0)] = 16
    
    pd.to_pickle(regs,outpath+"/Regions_"+mod+".pkl")
    
    ##########
    # Plotting
    
    # Set up figure
    # fig = plt.figure(figsize=(8.5,8.5))
    # ax1 = fig.add_subplot(1,1,1)
    # plt.title("Regions",fontsize=11)
    
    # # Set up map
    # mp = Basemap(projection='laea',lat_0=lat_0,lon_0=lon_0,\
    # llcrnrlat=bb[0], llcrnrlon=bb[1],urcrnrlat=bb[2], urcrnrlon=bb[3],resolution='l')
    
    # # Plot contour
    # cf2 = mp.contourf(lons,lats,regs,levs,latlon=True,cmap=plt.cm.Pastel2,extend="max")
    
    # # Add map accessories
    # mp.fillcontinents(color='0.8',lake_color='0.8')
    # mp.drawmapboundary(fill_color='1')
    # mp.drawmeridians(np.arange(0,360,30),color="0.4",linewidth=0.8)
    # mp.drawparallels(np.arange(0,90,15),color="0.4",linewidth=0.8)
    
    # # Add color bar
    # cbar_ax = fig.add_axes([0.87, 0.25, 0.02, 0.5])
    # cbar1 = fig.colorbar(cf2,cbar_ax,orientation='vertical',ticks=levs)
    # cbar1.ax.set_yticklabels(levs) 
    # cbar1.set_label('m/s')