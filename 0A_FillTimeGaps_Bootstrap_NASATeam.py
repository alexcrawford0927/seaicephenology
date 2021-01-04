'''*********************************************
Authors: Alex Crawford
Date Created: 8 Jul 2020
Date Modified: 8 Jul 2020

Purpose: To create files that fill in the gaps between sea ice data.

*********************************************'''
# Import clock:
from time import perf_counter as clock
# Start script stopwatch. The clock starts running when time is imported
start = clock()

'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import os
import numpy as np
import MERRA_Module as md

'''*******************************************
Declare Variables
*******************************************'''
print("Declaring Variables")
v = 1

vers = ["NASATeam","Bootstrap"] # version of sea ice data
dtypes = ['uint8','uint16'] # binary data types
nptypes = [np.uint8,np.uint16] # numpy data types
heads = [300,0] # Header length of bin files
ii = [[3,7,9,11],[3,7,9,11]] # Indices for subsetting the file names


path = "/Volumes/Miranda/SeaIce/"+vers[v]
inpath = path+"/Daily"
outpath = path+"/DailyFiller"


lyb = 1

mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
days = ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15",\
        "16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"]

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")
# Load valid files
files = os.listdir(inpath)
files = [f for f in files if (f.startswith('.') == 0)]
files.sort()

# Store correct version-specific variables
dtype = dtypes[v]
head = heads[v]
nptype= nptypes[v]
i = ii[v]

start = [int(files[0][i[0]:i[1]]),int(files[0][i[1]:i[2]]),int(files[0][i[2]:i[3]]),0,0,0]
end = [int(files[-1][i[0]:i[1]]),int(files[-1][i[1]:i[2]]),int(files[-1][i[2]:i[3]]),0,0,0]

filePr = files[0]

files2 = files[1:]
t = md.timeAdd(start,[0,0,1],lyb)

while t != end:
    
    # If this file already exists...
    if str(t[0])+mons[t[1]-1]+days[t[2]-1] in files2[0]:
        filePr = files2[0] # it becomes the new "previous" file
        files2 = files2[1:] # remove it from consideration
        t = md.timeAdd(t,[0,0,1]) # advance by 1 day
    
    # If this file doesn't exist...
    else:
        # Take the average of the "previous" file and the "next" file
        arrPr = np.fromfile(inpath+'/'+filePr,dtype=dtype)[head:].reshape((448,304)).astype(int)
        arrNx = np.fromfile(inpath+'/'+files2[0],dtype=dtype)[head:].reshape((448,304)).astype(int)
        
        arr = (arrPr + arrNx)/2
        
        # Save as a binary file with the same format
        outfile = files2[0][:i[0]]+str(t[0])+mons[t[1]-1]+days[t[2]-1]+files2[0][i[3]:]
        arr.astype(nptype).tofile(outpath+"/"+outfile)
        
        print("Saved "+str(t[0])+mons[t[1]-1]+days[t[2]-1])
        # Advance by 1 day
        t = md.timeAdd(t,[0,0,1])