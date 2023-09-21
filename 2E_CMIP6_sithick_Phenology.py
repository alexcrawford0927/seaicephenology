"""
Created: 31 Jan 2023
Modified: 2 Feb 2023
Author: Alex Crawford

Purpose: Calculates a "ice-thickness season" using sea ice thickness (note, not
effective thickness, just raw thickness); sea ice concentraton is needed
to distinguish whether a NAN in thickness means ocean with no  sea ice or land
"""

'''*******************************************
Load Modules
*******************************************'''
import os
import pandas as pd
import netCDF4 as nc
import numpy as np
import CycloneModule_12_4 as md

'''*******************************************
Declare Variables
*******************************************'''

### File Variables ###
experiment = 'historical' # 'ssp585' #
nch = 'sithick'
nci = 'siconc'

### Process Variables ###
minh = 0.10 # minimum sea ice thickness (in meters)
n = 5 # Moving Average Size (# observations total, including the central
# value -- should always be an odd number)
nanmax = 9999

dateref = [1850,1,1,0,0,0]

years = np.arange(1850,2013+1) # np.arange(2015,2099+1)
minmo = [8,10] # Range of months to search for minimum
maxmo = [2,6] # Range of months to search for maximum

### Path Variables ###
mod = "CMIP6"

path = '/project/6061839/crawfora' # '/Users/acrawfora/Documents/' #   "/Volumes/Troilus"

try:
    os.chdir(path)
except:
    path = "/Volumes/Cassandra"
    os.chdir(path)

inpathh = path+"/"+mod+"/data/e_"+experiment+"/v_"+nch
inpathi = path+"/"+mod+"/data/e_"+experiment+"/v_"+nci
outpath = path+"/"+mod+"/SeaIce/ThicknessPhenology/"+experiment+"/"+str(int(minh*100))+"cm"
suppath = path+"/"+mod+"/RegionMasks/SeaIce"
calendarpath = path+"/CMIP6/calendars.csv"

modstoskip = ['AWI-CM-1-1-MR','AWI-ESM-1-1-LR','ICON-ESM-LR',"KIOST-ESM"] # ['EC-Earth3','MPI-ESM1-2-LR','IPSL-CM6A-LR']
modstokeep = ["ACCESS-CM2","AWI-CM-1-1-MR","BCC-CSM2-MR",
              "CESM2","CESM2-FV2","CESM2-WACCM",
              'CMCC-CM2-SR5','CMCC-ESM2',"CNRM-CM6-1","CNRM-ESM2-1",
              "EC-Earth3-CC","IPSL-CM6A-LR","KIOST-ESM","MIROC-ES2L","MIROC6",
              "MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM"]

# modstokeep = ["NESM3"]


'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

# Load calendar
calendar = pd.read_csv(calendarpath)

# Load file names
hfiles = md.listdir(inpathh)
ifiles = md.listdir(inpathi)
hfiles = [f for f in hfiles if f.split('_')[1] not in modstoskip and f.split('_')[1] in modstokeep]
imods = [f.split("_v_"+nci)[0] for f in ifiles]
hmods = [f.split("_v_"+nch)[0] for f in hfiles]
matchedfiles = np.intersect1d(imods, hmods)

cmplfiles = md.listdir(outpath)

newfiles = [f for f in matchedfiles if f not in cmplfiles]
newfiles.sort()

# Run main code
for f in newfiles[:]:
    print(f)
    # Calendar
    dpy, lys = int(calendar[calendar['Model']==f.split('_')[1]]['dpy']), int(calendar[calendar['Model']==f.split('_')[1]]['lyb'])

    # Load files
    nchfiles = md.listdir(inpathh+"/"+f+"_v_"+nch)
    nchfiles = [f for f in nchfiles if f.endswith('.nc')]
    nchfiles.sort()

    ncifiles = md.listdir(inpathi+"/"+f+"_v_"+nci)
    ncifiles = [f for f in ncifiles if f.endswith('.nc')]
    ncifiles.sort()

    # Identify valid times for files
    c6hstarts = [ md.daysBetweenDates(dateref,[int(nf.split('_')[-1].split('-')[0][0:4]),int(nf.split('_')[-1].split('-')[0][4:6]),int(nf.split('_')[-1].split('-')[0][6:8]),0,0,0],lys,dpy) for nf in nchfiles]
    c6hends = [ md.daysBetweenDates(dateref,[int(nf.split('_')[-1].split('-')[1][0:4]),int(nf.split('_')[-1].split('-')[1][4:6]),int(nf.split('_')[-1].split('-')[1][6:8]),0,0,0],lys,dpy) for nf in nchfiles]

    c6istarts = [ md.daysBetweenDates(dateref,[int(nf.split('_')[-1].split('-')[0][0:4]),int(nf.split('_')[-1].split('-')[0][4:6]),int(nf.split('_')[-1].split('-')[0][6:8]),0,0,0],lys,dpy) for nf in ncifiles]
    c6iends = [ md.daysBetweenDates(dateref,[int(nf.split('_')[-1].split('-')[1][0:4]),int(nf.split('_')[-1].split('-')[1][4:6]),int(nf.split('_')[-1].split('-')[1][6:8]),0,0,0],lys,dpy) for nf in ncifiles]

    # Set up projection info
    pnc = nc.Dataset(suppath+"/Regions_"+f.split('_')[1]+".nc")
    lats, lons = pnc['lat'][:].data, pnc['lon'][:].data

    # Prep outputs
    fadlist, frdlist, ladlist, lrdlist, chdayslist, hdayslist, mxhlist, mxdlist = [], [], [], [], [], [], [], []

    for y in years:
        print("   Year: "+str(y))
        doys = np.arange(md.daysBetweenDates([y,1,1,0,0,0],[y+2,1,1,0,0,0],lys,dpy)).astype(int)+1
        dpy_y = dpy + md.leapyearBoolean([y+1])[0]*lys

        # Identify files with current & subseqent year
        timestart, timeend = md.daysBetweenDates(dateref,[y,1,1,0,0,0],lys,dpy), md.daysBetweenDates(dateref,[y+1,12,31,0,0,0],lys,dpy)

        fileshi = [i for i in range(len(c6hends)) if (c6hends[i] >= timestart) and (c6hstarts[i] <= timeend)]
        filesii = [i for i in range(len(c6iends)) if (c6iends[i] >= timestart) and (c6istarts[i] <= timeend)]

        # Prep outputs for this year
        mxdArr, mxh = np.zeros_like(lats)*np.nan, np.zeros_like(lats)*np.nan
        frdArr, lrdArr = np.zeros_like(lats)*np.nan, np.zeros_like(lats)*np.nan
        fadArr, ladArr = np.zeros_like(lats)*np.nan, np.zeros_like(lats)*np.nan
        chdaysArr, hdaysArr = np.zeros_like(lats)*np.nan, np.zeros_like(lats)*np.nan

        # Load the file(s) and store sithick in memory
        sithick = np.empty((0,lats.shape[-2],lats.shape[-1]))
        for fi in fileshi:
            ncf = nc.Dataset(inpathh+"/"+f+"_v_"+nch+"/"+nchfiles[fi])
            nctimes = np.arange(c6hstarts[fi],c6hends[fi]+1)
            ti = np.where( (nctimes >= timestart) & (nctimes <= timeend))[0]
            sithick = np.concatenate( (sithick, ncf[nch][ti,:,:].data), axis=0)

        # Load the file(s) and store siconc in memory
        siconc = np.empty((0,lats.shape[-2],lats.shape[-1]))
        for fi in filesii:
            ncf = nc.Dataset(inpathi+"/"+f+"_v_"+nci+"/"+ncifiles[fi])
            nctimes = np.arange(c6istarts[fi],c6iends[fi]+1)
            ti = np.where( (nctimes >= timestart) & (nctimes <= timeend))[0]
            siconc = np.concatenate( (siconc, ncf[nci][ti,:,:].data), axis=0)

        if sithick.shape[0] > 550: # Only perform processing if there's more than 1.5 years (should be 720-731 days)
            # If there's no sea ice, thickness should be 0....
            # but no need to do all the calculations if there's never sea ice - so make that NaN
            sithick[np.abs(siconc) >= nanmax] = np.nan
            sithick[np.abs(sithick) >= nanmax] = 0

            # Idetify maximum sea ice thickness for the given year
            # Subset data to be Sep 15 to Sep 15 (since PIOMAS has no leap years, this is kind of Sep 16 during leap years)
            mini0 = np.where( doys >= md.daysBetweenDates([y,1,1],[y,minmo[0],1]) )[0][0]
            mini1 = np.where( doys <= md.daysBetweenDates([y,1,1],[y,minmo[1]+1,1]) )[0][-1]

            maxi0 = np.where( doys >= md.daysBetweenDates([y,1,1],[y+1,maxmo[0],1]) )[0][0]
            maxi1 = np.where( doys <= md.daysBetweenDates([y,1,1],[y+1,maxmo[1]+1,1]) )[0][-1]

            mini2 = np.where( doys >= md.daysBetweenDates([y,1,1],[y+1,minmo[0],1]) )[0][0]
            mini3 = np.where( doys <= md.daysBetweenDates([y,1,1],[y+1,minmo[1]+1,1]) )[0][-1]

            mxh = np.amax(sithick[maxi0:maxi1,:,:],0) # Limit to the 180 days in the middle of the dataset
            mnh1 = np.amin(sithick[mini0:mini1,:,:],0)
            mnh2 = np.amin(sithick[mini2:mini3,:,:],0)

            ### Calculate Retreat and Advance Events ###
            validcells = np.where( (mxh > 0) )
            for j in range(len(validcells[0])):
                r,c = validcells[0][j], validcells[1][j] # Assign row and column

                # Identify when the maximum thickness occurs (median day if it happens more than once)
                mnd1i = mini0 + int(np.median(np.where(sithick[mini0:mini1,r,c] == mnh1[r,c])))
                mxdi = maxi0 + int(np.median(np.where(sithick[maxi0:maxi1,r,c] == mxh[r,c])))
                mnd2i = mini2 + int(np.median(np.where(sithick[mini2:mini3,r,c] == mnh2[r,c])))

                # Store Minimum Day & and (temporarily) maximum days
                mnd1rr = doys[mnd1i]
                mxdArr[r,c] = doys[mxdi]
                mnd2rr = doys[mnd2i]

                # If it's always below the thickness threshold...
                if mxh[r,c] < minh:
                    chdaysArr[r,c], hdaysArr[r,c] = 0, 0

                # If it's always above the thickness threshold...
                elif (mnh1[r,c] >= minh) & (mnh2[r,c] >= minh):
                    chdaysArr[r,c], hdaysArr[r,c] = dpy_y, dpy_y

                else:
                    above = np.where(sithick[:,r,c] >= minh)[0] # Indices above thickness
                    below = np.where(sithick[:,r,c] < minh)[0] # Indices below thickness

                    # First Retreat Day
                    # First index after max day that is below threshold
                    try:
                        fri = below[(below > mxdi) & (below <= mnd2i)][0]
                        frdArr[r,c] = doys[fri]
                    except:
                        frdArr[r,c] = np.nan

                    # First Advance Day
                    # First index above threshold
                    try:
                        fai = above[(above > mnd1i) & (above <= mxdi)][0]
                        fadArr[r,c] = doys[fai]
                    except:
                        fadArr[r,c] = np.nan

                    # Last Retreat Day
                    # Day after last day above threshold
                    try:
                        lri = above[(above < mnd2i) & (above >= mxdi)][-1]+1
                        lrdArr[r,c] = doys[lri]
                    except:
                        lrdArr[r,c] = np.nan

                    # Last Advance Day
                    # Day after last day below threshold before the max day
                    try:
                        lai = below[(below < mxdi) & (below >= mnd1i)][-1]+1
                        ladArr[r,c] = doys[lai]
                    except:
                        ladArr[r,c] = np.nan

                    # Ice-Cover Period
                    if (mnh1[r,c] < minh) and (mnh2[r,c] >= minh): # When it starts below threshold but ends above
                        chdaysArr[r,c] = mnd2rr - ladArr[r,c] + 1
                        hdaysArr[r,c] = mnd2rr - fadArr[r,c] + 1

                    elif (mnh1[r,c] >= minh) and (mnh2[r,c] < minh): # When it starts above threshold but ends below
                        chdaysArr[r,c] = frdArr[r,c] - mnd1rr
                        hdaysArr[r,c] = lrdArr[r,c] - mnd1rr
                    else:
                        chdaysArr[r,c] = frdArr[r,c] - ladArr[r,c]
                        hdaysArr[r,c] = lrdArr[r,c] - fadArr[r,c]

            ### Set Ocean Cells to 0 days if they never have sea ice ###
            chdaysArr[np.isnan(chdaysArr) & np.isfinite(mxh)] = 0
            hdaysArr[np.isnan(hdaysArr) & np.isfinite(mxh)] = 0

            # Remove objects from memory
            del validcells, above, below, mnh1, mnh2, doys, mnd1rr, mnd2rr, fai, lai, fri, lri

        ### Store Outputs ###
        mxdlist.append(mxdArr)
        mxhlist.append(mxh)
        frdlist.append(frdArr)
        fadlist.append(fadArr)
        lrdlist.append(lrdArr)
        ladlist.append(ladArr)
        chdayslist.append(chdaysArr)
        hdayslist.append(hdaysArr)

        # Remove objects from memory
        del frdArr, lrdArr, fadArr, ladArr, chdaysArr, hdaysArr, mxdArr, mxh


    print("Write to file")
    outName = f+"_v_siphenology"+str(int(minh*100))+"cm_"+str(min(years))+"-"+str(max(years))+".nc"
    try:
        os.chdir(outpath+"/"+f)
    except:
        os.mkdir(outpath+"/"+f)
        os.chdir(outpath+"/"+f)

    ncf = nc.Dataset(outpath+"/"+f+"/"+outName,'w', format='NETCDF4')
    ncf.createDimension('y', lats.shape[0])
    ncf.createDimension('x', lats.shape[1])
    ncf.createDimension('time', len(years))

    ncf.createVariable('y', np.float32, ('y',))
    ncf.createVariable('x', np.float32, ('x',))
    ncf.createVariable('time', np.int32, ('time',))

    ncf.createVariable('lat', np.float32, ('y','x',))
    ncf.createVariable('lon', np.float32, ('y','x',))

    ncf.createVariable('mxh', np.float32, ('time','y','x',))
    ncf.createVariable('mxd', np.int16, ('time','y','x',))
    ncf.createVariable('lrd', np.float32, ('time','y','x',))
    ncf.createVariable('lad', np.float32, ('time','y','x',))
    ncf.createVariable('frd', np.float32, ('time','y','x',))
    ncf.createVariable('fad', np.float32, ('time','y','x',))
    ncf.createVariable('cip', np.float32, ('time','y','x',))
    ncf.createVariable('ip', np.float32, ('time','y','x',))

    ncf.description = '''Phenology of Sea Ice Thickness based on Threshold
    of ''' + str(int(minh*100)) + '''cm. Note: Negative Values = N/A'''
    ncf.source = 'netCDF4 python module'
    ncf['time'].units = 'years'
    ncf['x'].units = 'm'
    ncf['y'].units = 'm'
    ncf['lat'].units = 'degrees north'
    ncf['lon'].units = 'degrees east'
    ncf['mxh'].units = 'm'
    ncf['mxd'].units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
    ncf['lrd'].units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
    ncf['frd'].units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
    ncf['lad'].units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
    ncf['fad'].units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
    ncf['ip'].units = 'days'
    ncf['cip'].units = 'days'

    ncf['time'][:] = years
    ncf['x'][:] = np.arange(lats.shape[1])
    ncf['y'][:] = np.arange(lats.shape[0])
    ncf['lat'][:] = lats
    ncf['lon'][:] = lons
    ncf['mxh'][:] = np.array(mxhlist)
    ncf['mxd'][:] = np.array(mxdlist)
    ncf['lrd'][:] = np.array(lrdlist)
    ncf['lad'][:] = np.array(ladlist)
    ncf['frd'][:] = np.array(frdlist)
    ncf['fad'][:] = np.array(fadlist)
    ncf['ip'][:] = np.array(hdayslist)
    ncf['cip'][:] = np.array(chdayslist)

    ncf.close()

print("Complete.")
