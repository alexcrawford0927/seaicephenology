'''*********************************************
Author: Alex Crawford
Date Created: 14 May 2019
Date Modified: 5 Jan 2022
Purpose: Plots all tracks in a given timespan as points & lines
*********************************************'''

'''**********
Import necessary modules for the analysis
**********'''
print("Loading modules and licenses.")
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection
from matplotlib.path import Path
from copy import deepcopy
import numpy as np
import pandas as pd
import netCDF4 as nc
import CycloneModule_13_2 as md

'''**********
Set Environement & Define Variables
**********'''
# Plotting Parameters
bbox = [-45,135,10,40] # ll lon, ur lon, ll lat, lr lat
lon_0 = 0 # 0 means prime meridian is at 6 o'clock in image
lat_0 = 90
plotprj = ccrs.LambertAzimuthalEqualArea(lon_0, lat_0)  

minls = 0
mintl = 0
minlat = 0

# Path Variables
vers = "13_2E5R" # Version of the Tracking Algorithm
typ = "System"
subset = "BBox10"
path = "/Volumes/Cressida" # "/media/alex/Datapool" #
inpath = path+"/CycloneTracking/tracking"+vers+"/"+subset+"/"+typ+"Tracks"
figpath = inpath+"/Figures/"

# Time Parameters
starttime = [2022,4,14,0,0,0] # [2010,12,17,0,0,0]
endtime = [2022,4,20,0,0,0] # [2010,12,21,0,0,0]
redtime = [2022,4,17,12,0,0] #[2010,12,19,0,0,0] # Time to plot with a special color
monthstep = [0,1,0,0,0,0]
reftime = [1900,1,1,0,0,0]

months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]

colors = ['purple','orange','pink','0.5','maroon','teal','tan','brown','k','r','b','g']

'''**********
Main Analysis
**********'''
endmonth = md.timeAdd(endtime,[0,1,0,0,0,0])
endmonth[2:] = [1,0,0,0]
etm = md.daysBetweenDates(reftime, endmonth)

et = md.daysBetweenDates(reftime, endtime)
st = md.daysBetweenDates(reftime, starttime)


try:
    spres = pd.read_pickle(path+"/CycloneTracking/tracking"+vers+"/cycloneparams.pkl")['spres']
except:
    spres = 100

# Load Projection
dataprj = ccrs.Geodetic()

# Set up figure
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)

# set up map projection
ax = fig.add_subplot(1,1,1, projection = plotprj ) # (nrows, ncols)
ax.set_extent(bbox, ccrs.Geodetic())
ax.gridlines(xlocs=range(-180,180,45),ylocs=range(0,90,15),linestyle='dashed',linewidth=0.5)
ax.coastlines('110m')  

time = deepcopy(starttime)
time[2:] = [1,0,0,0]
t0 = md.daysBetweenDates(reftime,redtime)
t = md.daysBetweenDates(reftime,time)
while t <= etm:
    Y = str(time[0])
    M = mons[time[1]-1]
    print(Y+M)
    
    col = colors[time[1]-1]
    
    # Load Data
    ct = pd.read_pickle(inpath+"/"+Y+"/"+subset+typ.lower()+"tracks"+Y+M+".pkl")
    
    # Subset Data
    ct = [c for c in ct if ((c.lifespan() >= minls) & (c.trackLength() >= mintl)) & (list(c.data.time)[0] <= et) & (list(c.data.time)[-1] >= st)]
    
    # Plot Data
    for c in ct:
        times = np.array(c.data["time"])
        
        # Plot each point as a blue dot and genesis point as a black dot
        lons = np.array(c.data["lon"])
        lats = np.array(c.data["lat"])
        ax.plot(lons,lats,marker='o',markersize=1,linewidth=1, transform = dataprj)
        ax.plot(lons[0],lats[0],color='k',marker='o',markersize=3,linewidth=1, transform = dataprj)
        
        # Plot location of each storm at time of interest
        x, y = lons[np.where( (lats >= minlat) &  (times == t0))], lats[np.where( (lats >= minlat) &  (times == t0))]
        ax.plot(x,y,color='r',marker='o',markersize=3, transform = dataprj)
        
        if len(lats) > 0:

            if typ == 'System':
                ax.text(lons[0],lats[0]-1,c.sid, transform = dataprj)
            else:
                ax.text(lons[0],lats[0]-1,c.tid, transform = dataprj)

    time = md.timeAdd(time,monthstep)
    t = md.daysBetweenDates(reftime,time)

plt.title(typ+" Tracks " + str(starttime[3]) + " UTC "+str(starttime[2])+" "+months[starttime[1]-1]+" "+str(starttime[0])+" - "+str(endtime[3])+" UTC "+str(endtime[2])+" "+months[endtime[1]-1]+" "+str(endtime[0])+", "+vers+"\n\
          \n mintl = "+str(mintl)+" gridcells, minls = "+str(minls)+" day")


# import matplotlib.pyplot as plt

# lons = range(50)
# lats = range(25, 75)
# vals = range(50,100)

# plt.scatter(lons, lats, c=vals, cmap='afmhot')
# plt.colorbar()
# plt.show()

# c182 = [c for c in ct if c.tid == 182][0]
# c201 = [c for c in ct if c.tid == 201][0]
# c244 = [c for c in ct if c.tid == 244][0]

# c182.data.time
# c201.data.lat
# c201.data.long
# c244.data.lat
# c244.data.long



# md.timeAdd([1900,1,1,0,0,0],[0,0,np.array(c182.data.time)[-1],0,0,0])
