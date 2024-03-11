"""
Date Created: 9 Mar 2023
Date Modified: 9 Mar 2023
Author: Alex Crawford
Purpose: Create a set of maps to compare the difference in ice-free period
between two periods using three different measurements for "ice-free period"

Also includes annotation for two regional averages

"""
import numpy as np
import netCDF4 as nc
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import xesmf as xe
import CycloneModule_13_2 as md

'''*******************************************
Declare Variables
*******************************************'''
print("Declaring Variables")

ymin, ymax = 1979, 2021 # for input file
ymin1, ymax1 = 1980, 1989 # earlier period
ymin2, ymax2 = 2012, 2021 # later period
bymin, bymax = 1850, 1900 # baseline for temperature

regvar = 'regHB3'
regs = [41,42]
nanmin = 0
nmin = ymax1-ymin1+1

thres = [['C10','C10'],['C10','C10']] # [['C50','C10'],['C50','C10'],['10cm','10cm']]
mod = ['Bootstrap','NASATeam']
xlabel = [str(ymin1)+"-"+str(ymax1),str(ymin2)+"-"+str(ymax2),'Difference']
ylabel = ['Ice-Free Period: Bootstrap','Ice-Free Period: NASA Team', 'Surface Air Temperature: BEST']

inpaths1 = ['/Volumes/Miranda/SeaIce/', '/Volumes/Miranda/SeaIce/']
inpaths2 = ['AdvanceRetreat2', 'AdvanceRetreat2']
prjpath = ['/Volumes/Miranda/Projections/psn_projection.nc',
           '/Volumes/Miranda/Projections/psn_projection.nc']
tpath = "/Volumes/Theseus/SurfaceTemperature/BEST/BEST_Land_and_Ocean_LatLong1.nc"
tprjpath = "/Volumes/Theseus/Projections/BEST_Projection.nc"
plotprjpath = '/Volumes/Miranda/Projections/EASE2_N0_25km_Projection.nc'
outpath = "/Users/acrawfora/OneDrive - University of Manitoba/CMIP6 OWP/Hudson Bay Figures/ThicknessPhenologyFigures/Temperature_and_Ice-FreePeriodComparison_"+str(ymin1)+"-"+str(ymax1)+"_vs_"+str(ymin2)+"-"+str(ymax2)+"_"+"-".join(ylabel)+"_V3.png"

# Plotting Variables
prjdata = ccrs.LambertAzimuthalEqualArea(central_longitude=0, central_latitude=90)
prjplot = ccrs.LambertAzimuthalEqualArea(central_longitude=-80, central_latitude=90)
bbox = [-1000000,500000,-2700000,-4400000] # left, right, top, bottom
lines = [(63,-91,63,-88),(63,-88,56,-88),(60,-88,60,-77)]

levs = np.arange(90,180+15,15)
levsdiff = np.arange(14,49+1,7) # md.divbreaks(49,7)
ticks = np.arange(90,180+15,15)
ticksdiff = np.arange(14,49+1,7) # md.divbreaks(42,14)

tlevs = np.arange(-0.4,2.0,0.2)
tlevsdiff = np.arange(0.6,2.0,0.2)

'''*******************************************
Pre-Processing
*******************************************'''
print('Pre-Processing')
YY = str(ymin)+"-"+str(ymax)
years = np.arange(ymin,ymax+1)

# Load plot projection
plotprjnc = nc.Dataset(plotprjpath)
xarr, yarr = plotprjnc['x'][:].data, plotprjnc['y'][:].data

tprjnc = nc.Dataset(tprjpath)

# Prep lists for outputs
avg1, avg2, pval  = [[[],[]] for i in range(3)]
opc1, opc2, pvals = [[] for i in range(3)]
for i in range(len(thres)):
    ### Projection Info
    prjnc = nc.Dataset(prjpath[i])
    regarr = prjnc[regvar][:].data
    areaarr = prjnc['area'][:].data

    lonname = [key for key in prjnc.variables.keys() if key.startswith('lon') or key.startswith('nav_lon')][0]
    latname = [key for key in prjnc.variables.keys() if key.startswith('lat') or key.startswith('nav_lat')][0]
    inlons, inlats = prjnc[lonname][:].data, prjnc[latname][:].data

    regridder = xe.Regridder({'lon':inlons,'lat':inlats},
                              {'lon':plotprjnc['lon'][:].data,'lat':plotprjnc['lat'][:].data},'bilinear')

    #### Load Data
    retnc = nc.Dataset(inpaths1[i]+"/"+mod[i]+"/"+inpaths2[i]+"/"+thres[i][0]+"/siphenology"+thres[i][0]+"_"+mod[i]+"_"+YY+".nc")
    advnc = nc.Dataset(inpaths1[i]+"/"+mod[i]+"/"+inpaths2[i]+"/"+thres[i][1]+"/siphenology"+thres[i][1]+"_"+mod[i]+"_"+YY+".nc")

    retarr = retnc['lrd'][:].data.astype(float)
    advarr = advnc['fad'][:].data.astype(float)

    retarr[retarr < nanmin] = np.nan
    advarr[advarr < nanmin] = np.nan

    # Adjust year of retreat if using thickness phenology
    if 'cm' in thres[i][0]:
        retarr = np.concatenate( (np.zeros_like(retarr[:1,:])*np.nan,retarr[:-1,:]), axis=0 ) - 365

    # Calculate ice-free period
    opcarr = advarr - retarr


    ### Take spatial averages
    opc_reg1 = np.nansum(areaarr[regarr == regs[0]] * opcarr[:,regarr == regs[0]],axis=1) / np.sum(np.isfinite(opcarr[:,regarr == regs[0]])*areaarr[regarr == regs[0]],axis=1)
    opc_reg2 = np.nansum(areaarr[regarr == regs[1]] * opcarr[:,regarr == regs[1]],axis=1) / np.sum(np.isfinite(opcarr[:,regarr == regs[1]])*areaarr[regarr == regs[1]],axis=1)

    # Take temporal averages of the spatial averages (and compare)
    avg1[0].append( opc_reg1[np.where((years >= ymin1) & (years <= ymax1))].mean() )
    avg2[0].append( opc_reg1[np.where((years >= ymin2) & (years <= ymax2))].mean() )
    pval[0].append( stats.mannwhitneyu(opc_reg1[np.where((years >= ymin1) & (years <= ymax1))], opc_reg1[np.where((years >= ymin2) & (years <= ymax2))])[1] )

    avg1[1].append( opc_reg2[np.where((years >= ymin1) & (years <= ymax1))].mean() )
    avg2[1].append( opc_reg2[np.where((years >= ymin2) & (years <= ymax2))].mean() )
    pval[1].append( stats.mannwhitneyu(opc_reg2[np.where((years >= ymin1) & (years <= ymax1))], opc_reg2[np.where((years >= ymin2) & (years <= ymax2))])[1] )
    r2 = np.nanmean( opcarr[np.where((years >= ymin2) & (years <= ymax2))[0],:,:], axis=0)

    # Identify p-values from a test of difference
    pvalarr = np.ones_like(regarr).astype(np.float64)
    
    opc1arr = opcarr[np.where((years >= ymin1) & (years <= ymax1))]
    opc2arr = opcarr[np.where((years >= ymin2) & (years <= ymax2))]
    nopc1 = np.isfinite(opc1arr).sum(axis=0)
    nopc2 = np.isfinite(opc2arr).sum(axis=0)
    
    rows, cols = np.where((nopc1 >= nmin) & (nopc2 >= nmin))
    for j in range(len(rows)):
        r, c = rows[j], cols[j]

        pvalarr[r,c] = stats.mannwhitneyu( opc1arr[:,r,c], opc2arr[:,r,c] )[1]

    # Reproject temporal averages (by gridcell) for plotting
    pvals.append(regridder(pvalarr))
    opc1.append(regridder(np.nanmean(opc1arr,axis=0)))
    opc2.append(regridder(np.nanmean(opc2arr,axis=0)))

### LOAD TEMPERATURE DATA ###
tnc = nc.Dataset(tpath)
tlats, tlons = tnc['latitude'][:].data, tnc['longitude'][:].data
tyears = tnc['time'][::12].data.astype(int)

# Convert monthly to annual
tann = tnc['temperature'][:,np.where(tlats > 40)[0],np.where((tlons < -15) & (tlons > -150))[0]].data

tann = tann.reshape((int(tann.shape[0]/12),12,tann.shape[-2],tann.shape[-1])).mean(axis=1)
# # baseline average
# tbase = np.nanmean(tann[(tyears >= bymin) & (tyears <= bymax)],axis=0) 
# # Temperature anomaly
# tann_amly = tann - tbase

tlons2, tlats2 = np.meshgrid(tlons[(tlons < -15) & (tlons > -150)],tlats[(tlats > 40)])

# Calculate temporal averages
tavg1 = np.nanmean(tann[(tyears >= ymin1) & (tyears <= ymax1)],axis=0)
tavg2 = np.nanmean(tann[(tyears >= ymin2) & (tyears <= ymax2)],axis=0)
tpval = stats.mannwhitneyu(tann[(tyears >= ymin1) & (tyears <= ymax1)], tann[(tyears >= ymin2) & (tyears <= ymax2)])[1]

# Calculate spatial averages
treg = tprjnc[regvar][np.where(tlats > 40)[0],np.where((tlons < -15) & (tlons > -150))[0]].data
coslats = np.cos(tlats2*np.pi/180)

tann_reg1 = np.nansum( (coslats*tann)[:,treg == regs[0]], axis=1 ) / np.sum( (coslats*np.isfinite(tann))[:,treg == regs[0]], axis=1 )
tann_reg2 = np.nansum( (coslats*tann)[:,treg == regs[1]], axis=1 ) / np.sum( (coslats*np.isfinite(tann))[:,treg == regs[1]], axis=1 )

avg1[0].append( tann_reg1[np.where((tyears >= ymin1) & (tyears <= ymax1))].mean() )
avg2[0].append( tann_reg1[np.where((tyears >= ymin2) & (tyears <= ymax2))].mean() )
pval[0].append( stats.mannwhitneyu(tann_reg1[np.where((tyears >= ymin1) & (tyears <= ymax1))], tann_reg1[np.where((tyears >= ymin2) & (tyears <= ymax2))])[1])

avg1[1].append( tann_reg2[np.where((tyears >= ymin1) & (tyears <= ymax1))].mean() )
avg2[1].append( tann_reg2[np.where((tyears >= ymin2) & (tyears <= ymax2))].mean() )
pval[1].append( stats.mannwhitneyu(tann_reg2[np.where((tyears >= ymin1) & (tyears <= ymax1))], tann_reg2[np.where((tyears >= ymin2) & (tyears <= ymax2))])[1])

pval = [['*' if p < 0.05 else '' for p in pv] for pv in pval]

'''*******************************************
Plotting
*******************************************'''
print('Plotting')

fig = plt.figure(figsize=(7.2,8.5))

### Temperature Maps ###

# PERIOD 1 #
ax = fig.add_subplot(3,3,1, projection=prjplot)
ax.add_feature(cfeature.COASTLINE)
ax.set_extent(bbox, crs=prjplot)
ct = ax.contourf(tlons2, tlats2, tavg1, tlevs, cmap=plt.cm.inferno, transform=ccrs.PlateCarree(), extend='both')

gdl = ax.gridlines(xlocs=range(-180,180,15),ylocs=range(0,90,5),linestyle='dashed',linewidth=0.5,draw_labels=True)
gdl.top_labels = False
gdl.right_labels = False
gdl.bottom_labels = False
gdl.xformatter = LongitudeFormatter()
gdl.yformatter = LatitudeFormatter()
gdl.xlabel_style = {'size':8}
gdl.ylabel_style = {'size':8}

# Add Study Area Lines
for l in range(len(lines)):
    if lines[l][0] == lines[l][2]:
        ax.plot( np.arange(np.min([lines[l][1],lines[l][3]]),np.max([lines[l][1],lines[l][3]])+0.5,0.5), np.repeat(lines[l][0],np.ceil(1+(abs(lines[l][3]-lines[l][1]))/0.5)), color='k', transform=ccrs.Geodetic())
    else:
        ax.plot(np.repeat(lines[l][1],np.ceil(1+(abs(lines[l][2]-lines[l][0]))/0.5)), np.arange(np.min([lines[l][2],lines[l][0]]),np.max([lines[l][2],lines[l][0]])+0.5,0.5), color='k', transform=ccrs.Geodetic())

# Annotate
ax.annotate(xlabel[0], xy=(0.50, 1.02), xycoords='axes fraction',va='bottom',ha='center', size=9, weight='bold')
txt = [ ax.annotate('SHB', xy=(0.50,0.50), xycoords='axes fraction',va='center',color='white', ha='center'),
       ax.annotate('WHB', xy=(0.23,0.63), xycoords='axes fraction',va='center',color='white', ha='center')]
# [ t.set_path_effects([pe.withStroke(linewidth=3, foreground='w')]) for t in txt ]


ax.annotate(ylabel[-1],xy=(-0.30,0.50),xycoords='axes fraction',rotation=90, va='center',ha='right', size=9, weight='bold')

txt = ax.annotate(md.abc[0],xy=(0.02,0.98),xycoords='axes fraction',ha='left',va='top',weight='bold')
txt.set_path_effects([pe.withStroke(linewidth=3, foreground='w')])

txt = ax.annotate('SHB = ' + str((np.round(avg1[0][-1],1))) +" °C\nWHB = " + str((np.round(avg1[1][-1],1))) + " °C",
            xy=(0.02,0.02),xycoords='axes fraction', va='bottom', ha='left', size=8, color='white')
# txt.set_path_effects([pe.withStroke(linewidth=2, foreground='w')])

# PERIOD 2 #
ax = fig.add_subplot(3,3,2, projection=prjplot)
ax.add_feature(cfeature.COASTLINE)
ax.set_extent(bbox, crs=prjplot)
ax.contourf(tlons2, tlats2, tavg2, tlevs, cmap=plt.cm.inferno, transform=ccrs.PlateCarree(), extend='both')

ax.gridlines(xlocs=range(-180,180,15),ylocs=range(0,90,5),linestyle='dashed',linewidth=0.5,draw_labels=False)

# Add Study Area Lines
for l in range(len(lines)):
    if lines[l][0] == lines[l][2]:
        ax.plot( np.arange(np.min([lines[l][1],lines[l][3]]),np.max([lines[l][1],lines[l][3]])+0.5,0.5), np.repeat(lines[l][0],np.ceil(1+(abs(lines[l][3]-lines[l][1]))/0.5)), color='k', transform=ccrs.Geodetic())
    else:
        ax.plot(np.repeat(lines[l][1],np.ceil(1+(abs(lines[l][2]-lines[l][0]))/0.5)), np.arange(np.min([lines[l][2],lines[l][0]]),np.max([lines[l][2],lines[l][0]])+0.5,0.5), color='k', transform=ccrs.Geodetic())

ax.annotate(xlabel[1], xy=(0.50, 1.02), xycoords='axes fraction',va='bottom',ha='center', size=9, weight='bold')
txt = ax.annotate(md.abc[1],xy=(0.02,0.98),xycoords='axes fraction',ha='left',va='top',weight='bold')
txt.set_path_effects([pe.withStroke(linewidth=3, foreground='w')])

txt = ax.annotate('SHB = ' + str((np.round(avg2[0][-1],1))) +" °C\nWHB = " + str((np.round(avg2[1][-1],1))) + " °C",
            xy=(0.02,0.02),xycoords='axes fraction', va='bottom', ha='left', size=8, color='white')
# txt.set_path_effects([pe.withStroke(linewidth=2, foreground='w')])

# DIFFERENCE #
ax = fig.add_subplot(3,3,3, projection=prjplot)
ax.add_feature(cfeature.COASTLINE)
ax.set_extent(bbox, crs=prjplot)
ctd = ax.contourf(tlons2, tlats2, tavg2-tavg1, tlevsdiff, cmap=plt.cm.Reds, transform=ccrs.PlateCarree(), extend='both')
ax.contourf(tlons2, tlats2, np.where( tpval < 0.05, 0, 1), [0,0.05,1], colors="none", hatches=['...',''], transform=ccrs.PlateCarree())

ax.gridlines(xlocs=range(-180,180,15),ylocs=range(0,90,5),linestyle='dashed',linewidth=0.5,draw_labels=False)

# Add Study Area Lines
for l in range(len(lines)):
    if lines[l][0] == lines[l][2]:
        ax.plot( np.arange(np.min([lines[l][1],lines[l][3]]),np.max([lines[l][1],lines[l][3]])+0.5,0.5), np.repeat(lines[l][0],np.ceil(1+(abs(lines[l][3]-lines[l][1]))/0.5)), color='k', transform=ccrs.Geodetic())
    else:
        ax.plot(np.repeat(lines[l][1],np.ceil(1+(abs(lines[l][2]-lines[l][0]))/0.5)), np.arange(np.min([lines[l][2],lines[l][0]]),np.max([lines[l][2],lines[l][0]])+0.5,0.5), color='k', transform=ccrs.Geodetic())

ax.annotate(xlabel[2], xy=(0.50, 1.02), xycoords='axes fraction',va='bottom',ha='center', size=9, weight='bold')
txt = ax.annotate(md.abc[2],xy=(0.02,0.98),xycoords='axes fraction',ha='left',va='top',weight='bold')
txt.set_path_effects([pe.withStroke(linewidth=3, foreground='w')])

txt=ax.annotate('SHB = ' + str((np.round(avg2[0][-1]-avg1[0][-1],1))) +" °C" + pval[0][-1] + "\nWHB = " + str((np.round(avg2[1][-1]-avg1[1][-1],1))) + " °C" + pval[1][-1],
            xy=(0.02,0.02),xycoords='axes fraction', va='bottom', ha='left', size=8)
# txt.set_path_effects([pe.withStroke(linewidth=2, foreground='w')])

### Sea Ice Maps ###
for i in range(len(mod)):

    ### PERIOD 1 ###
    ax = fig.add_subplot(3,3,3*(i+1)+1, projection=prjplot)
    ax.add_feature(cfeature.COASTLINE)
    ax.set_extent(bbox, crs=prjplot)
    ax.contourf(xarr,yarr, np.where(np.isfinite(opc1[i]), np.nan, 1), [0,1,2,3,4], cmap=plt.cm.gray_r,transform=prjdata)
    cf = ax.contourf(xarr, yarr, opc1[i], levs, cmap=plt.cm.viridis, transform=prjdata, extend='both')

    gdl = ax.gridlines(xlocs=range(-180,180,15),ylocs=range(0,90,5),linestyle='dashed',linewidth=0.5,draw_labels=True)
    gdl.top_labels = False
    gdl.right_labels = False
    gdl.xformatter = LongitudeFormatter()
    gdl.yformatter = LatitudeFormatter()
    gdl.xlabel_style = {'size':8}
    gdl.ylabel_style = {'size':8}

    if i == 0:
        gdl.bottom_labels = False

    # Add Study Area Lines
    for l in range(len(lines)):
        if lines[l][0] == lines[l][2]:
            ax.plot( np.arange(np.min([lines[l][1],lines[l][3]]),np.max([lines[l][1],lines[l][3]])+0.5,0.5), np.repeat(lines[l][0],np.ceil(1+(abs(lines[l][3]-lines[l][1]))/0.5)), color='k', transform=ccrs.Geodetic())
        else:
            ax.plot(np.repeat(lines[l][1],np.ceil(1+(abs(lines[l][2]-lines[l][0]))/0.5)), np.arange(np.min([lines[l][2],lines[l][0]]),np.max([lines[l][2],lines[l][0]])+0.5,0.5), color='k', transform=ccrs.Geodetic())

    ax.annotate(ylabel[i],xy=(-0.30,0.50),xycoords='axes fraction',rotation=90, va='center',ha='right', size=9, weight='bold')

    txt = ax.annotate(md.abc[3*(i+1)],xy=(0.02,0.98),xycoords='axes fraction',ha='left',va='top',weight='bold')
    txt.set_path_effects([pe.withStroke(linewidth=3, foreground='w')])

    ax.annotate('SHB = ' + str(int(np.round(avg1[0][i],0))) +" days\nWHB = " + str(int(np.round(avg1[1][i],0))) + " days",xy=(0.02,0.02),xycoords='axes fraction', va='bottom', ha='left', size=8)

    ### PERIOD 2 ###
    ax = fig.add_subplot(3,3,3*(i+1)+2, projection=prjplot)
    gdl = ax.gridlines(xlocs=range(-180,180,15),ylocs=range(0,90,5),linestyle='dashed',linewidth=0.5,draw_labels=True)
    gdl.top_labels, gdl.right_labels, gdl.left_labels = False, False, False
    gdl.xformatter = LongitudeFormatter()
    gdl.yformatter = LatitudeFormatter()
    gdl.xlabel_style = {'size':8}
    gdl.ylabel_style = {'size':8}

    if i == 0:
        gdl.bottom_labels = False

    ax.add_feature(cfeature.COASTLINE)
    ax.set_extent(bbox, crs=prjplot)
    ax.contourf(xarr,yarr, np.where(np.isfinite(opc2[i]), np.nan, 1), [0,1,2,3,4], cmap=plt.cm.gray_r,transform=prjdata)
    ax.contourf(xarr, yarr, opc2[i], levs, cmap=plt.cm.viridis, transform=prjdata, extend='both')

    # Add Study Area Lines
    for l in range(len(lines)):
        if lines[l][0] == lines[l][2]:
            ax.plot( np.arange(np.min([lines[l][1],lines[l][3]]),np.max([lines[l][1],lines[l][3]])+0.5,0.5), np.repeat(lines[l][0],np.ceil(1+(abs(lines[l][3]-lines[l][1]))/0.5)), color='k', transform=ccrs.Geodetic())
        else:
            ax.plot(np.repeat(lines[l][1],np.ceil(1+(abs(lines[l][2]-lines[l][0]))/0.5)), np.arange(np.min([lines[l][2],lines[l][0]]),np.max([lines[l][2],lines[l][0]])+0.5,0.5), color='k', transform=ccrs.Geodetic())

    # Annotate
    txt = ax.annotate(md.abc[3*(i+1)+1],xy=(0.02,0.98),xycoords='axes fraction',ha='left',va='top',weight='bold')
    txt.set_path_effects([pe.withStroke(linewidth=3, foreground='w')])

    ax.annotate('SHB = ' + str(int(np.round(avg2[0][i],0))) +" days\nWHB = " + str(int(np.round(avg2[1][i],0))) + " days", xy=(0.02,0.02),xycoords='axes fraction', va='bottom', ha='left', size=8)

    ### DIFFERENCE ###
    ax = fig.add_subplot(3,3,3*(i+1)+3, projection=prjplot)
    gdl = ax.gridlines(xlocs=range(-180,180,15),ylocs=range(0,90,5),linestyle='dashed',linewidth=0.5,draw_labels=True)
    gdl.top_labels, gdl.right_labels, gdl.left_labels = False, False, False
    gdl.xformatter = LongitudeFormatter()
    gdl.yformatter = LatitudeFormatter()
    gdl.xlabel_style = {'size':8}
    gdl.ylabel_style = {'size':8}

    if i == 0:
        gdl.bottom_labels = False

    ax.add_feature(cfeature.COASTLINE)
    ax.set_extent(bbox, crs=prjplot)
    ax.contourf(xarr,yarr, np.where(np.isfinite(opc2[i]-opc1[i]), np.nan, 1), [0,1,2,3,4], cmap=plt.cm.gray_r,transform=prjdata)
    
    # plotarr = np.where( pvals[i] < 0.05, opc2[i]-opc1[i], np.nan)
    plotarr = opc2[i]-opc1[i]
    cd = ax.contourf(xarr, yarr, plotarr, levsdiff, cmap=plt.cm.Reds, transform=prjdata, extend='both')
    ax.contourf(xarr, yarr, np.where( pvals[i] < 0.05, 0, 1), [0,0.05,1], colors="none", hatches=['...',''], transform=prjdata)

    # Add Study Area Lines
    for l in range(len(lines)):
        if lines[l][0] == lines[l][2]:
            ax.plot( np.arange(np.min([lines[l][1],lines[l][3]]),np.max([lines[l][1],lines[l][3]])+0.5,0.5), np.repeat(lines[l][0],np.ceil(1+(abs(lines[l][3]-lines[l][1]))/0.5)), color='k', transform=ccrs.Geodetic())
        else:
            ax.plot(np.repeat(lines[l][1],np.ceil(1+(abs(lines[l][2]-lines[l][0]))/0.5)), np.arange(np.min([lines[l][2],lines[l][0]]),np.max([lines[l][2],lines[l][0]])+0.5,0.5), color='k', transform=ccrs.Geodetic())

    # Annotate
    txt = ax.annotate(md.abc[3*(i+1)+2],xy=(0.02,0.98),xycoords='axes fraction',ha='left',va='top',weight='bold')
    txt.set_path_effects([pe.withStroke(linewidth=3, foreground='w')])

    ax.annotate('SHB ='+str(int(np.round(avg2[0][i]-avg1[0][i],0))) + " days" + pval[0][i]+"\nWHB = "+str(int(np.round(avg2[1][i]-avg1[1][i],0))) + " days" + pval[1][i],xy=(0.02,0.02),xycoords='axes fraction', va='bottom', ha='left', size=8)


### Accessories ###
plt.tight_layout(rect=[0.03,0,1.0,1.0],h_pad=1.1)

### Color Bars ###
cbar_ax1 = fig.add_axes([0.21, 0.355, 0.4, 0.01]) # xloc, yloc, xsize, ysize
cbar1 = fig.colorbar(cf,cbar_ax1,orientation='horizontal',ticks=ticks)
cbar1.set_label("Continuous Ice-Free Period (Days)",size=7)
cbar1.ax.tick_params(labelsize=7)

cbar_ax2 = fig.add_axes([0.69, 0.355, 0.28, 0.01]) # xloc, yloc, xsize, ysize
cbar2 = fig.colorbar(cd,cbar_ax2,orientation='horizontal',ticks=ticksdiff)
cbar2.set_label("Difference (Days)",size=7)
cbar2.ax.tick_params(labelsize=7)

cbar_ax3 = fig.add_axes([0.21, 0.675, 0.4, 0.01])
cbar3 = fig.colorbar(ct, cbar_ax3, orientation='horizontal')
cbar3.set_label("Avg. Surface Air Temperature (°C)",size=7)
cbar3.ax.tick_params(labelsize=7)

cbar_ax4 = fig.add_axes([0.69, 0.675, 0.28, 0.01])
cbar4 = fig.colorbar(ctd, cbar_ax4, orientation='horizontal')
cbar4.set_label("Temperature Difference (°C)",size=7)
cbar4.ax.tick_params(labelsize=7)

plt.annotate('', xy=(0.19,0.31),xytext=(0.19, 0.375),xycoords='figure fraction',
             textcoords='figure fraction',arrowprops=dict(arrowstyle='<->',facecolor='k'))
plt.annotate('', xy=(0.625,0.31),xytext=(0.625, 0.375),xycoords='figure fraction',
             textcoords='figure fraction',arrowprops=dict(arrowstyle='<->',facecolor='k'))
plt.annotate('', xy=(0.98,0.31),xytext=(0.98, 0.375),xycoords='figure fraction',
             textcoords='figure fraction',arrowprops=dict(arrowstyle='<->',facecolor='k'))

plt.annotate('', xy=(0.19,0.69),xytext=(0.19, 0.65),xycoords='figure fraction',
             textcoords='figure fraction',arrowprops=dict(arrowstyle='->',facecolor='k'))
plt.annotate('', xy=(0.625,0.69),xytext=(0.625, 0.65),xycoords='figure fraction',
             textcoords='figure fraction',arrowprops=dict(arrowstyle='->',facecolor='k'))
plt.annotate('', xy=(0.98,0.69),xytext=(0.98, 0.65),xycoords='figure fraction',
             textcoords='figure fraction',arrowprops=dict(arrowstyle='->',facecolor='k'))

# plt.savefig(outpath,dpi=300)
