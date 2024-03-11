# Purpose: Calculates the multi-model mean and internal variability of the sensitivity 
# of sea ice variables to a unit change (1°C) in regional or global temperature.

# Load Modules
library(stringr)

sd_unbiased <- function(data){
  # Calculates an unbiased estimate of standard deviation using the method
  # described by Ben W. Bolch in "More on unbiased estimation of the standard deviation", 
  # The American Statistician, 22(3), p. 27 (1968). Degrees of freedom are set to n-1 by default.
  # data = a list or array of values in the sample
  # For n > 100,  don't use this -- just use the normal standard deviation with Bessel correction (dof = n-1)
  
  # Initial standard deviation calculation for sample (dof = n-1)
  s = sd(data)
  
  # Identify sample size
  n = length(data)
  
  # Calculate the correction using dof = n-1
  c4 = sqrt(2/(n-1)) * gamma(n/2) / gamma((n-1)/2)
  
  # Apply correction
  return(s/c4)
}

#################################
########## PREP WORK ############
#################################
# Define Variables
ct = 15
ver = 6
vars = c("OPCavg","LRDlt","FADgt")
tvar = 'tregion' # 'tregion' or 'tglobal'
path = '/Volumes/Troilus/CMIP6'
cmip6path = paste0(path,"/RegionalStats")
ssmipath = paste0(path,"/RegionalStats/SSMI")

modstokeep = c("ACCESS-CM2","BCC-CSM2-MR","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-CM6-1-HR",
               "CNRM-ESM2-1","EC-Earth3","EC-Earth3-Veg","GFDL-CM4","IPSL-CM6A-LR","MIROC-ES2L","MIROC6",
               "MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM","UKESM1-0-LL")

tymin = 1979
tymax = 2013

bmin = 1880
bmax = 1900

# Read in Sea Ice Data
thresh = read.csv(paste0(path,"/Regional_DOY_Thresholds25_C50_V2.csv"))
boot <- read.csv(paste0(ssmipath,"/Bootstrap_Regionalized_Historical_C",ct,"_V",ver,".csv"))
osi <- read.csv(paste0(ssmipath,"/OSISAF_Regionalized_Historical_C",ct,"_V",ver,".csv"))
nasa <- read.csv(paste0(ssmipath,"/NASATeam_Regionalized_Historical_C",ct,"_V",ver,".csv"))
cmip6 <- read.csv(paste0(cmip6path,"/CMIP6/CMIP6_Regionalized_Historical_C",ct,"_V",ver,".csv"))

# Read in Temperature Data
tbest = read.csv(paste0('/Volumes/Troilus/SurfaceTemperature/BEST/BEST_Annual_ts_1750-2019.csv'))
tgis = read.csv(paste0('/Volumes/Troilus/SurfaceTemperature/GISTEMP/GISTEMP_Annual_ts_1880-2019.csv'))
tnoaa = read.csv(paste0('/Volumes/Troilus/SurfaceTemperature/NOAAGlobalTemp/NOAAGlobalTemp_Annual_ts_1880-2019.csv'))
thad = read.csv(paste0('/Volumes/Troilus/SurfaceTemperature/HadCRUT/HadCRUT_Annual_ts_1850-2019.csv'))
tcmip6 = read.csv(paste0(path,'/CMIP6_Annual_tas_historical.csv'))

###############
# PROCESSING DATA
################
# Convert model names to characters
cmip6$Family <- as.character(cmip6$Family)
tcmip6$Family <- as.character(tcmip6$Family)

# Re-name "Models" to account for differences
cmip6$Model <- with(cmip6, paste0(Family,'--',Member))
tcmip6$Model <- with(tcmip6, paste0(Family,'--',Member))

# Subset to chosen models #
cmip6 = cmip6[cmip6$Family %in% modstokeep,]
tcmip6 = tcmip6[tcmip6$Family %in% modstokeep,]

# Store Unique Values
allmods = unique(cmip6$Model)
tmods = unique(tcmip6$Model)
mods = c()
for(mod in allmods){if(mod %in% tmods){mods <- c(mods,mod)}}

regs = unique(cmip6$Region)
years = unique(cmip6$Year)
obsyears = unique(osi$Year)
modfams <- unique(cmip6$Family) 

# Align Temperature data so that each record has the same baseline
tbest[,tvar] <- tbest[,tvar] - mean(tbest[tbest$Year >= bmin & tbest$Year <= bmax,tvar],na.rm=T)
tnoaa[,tvar] <- tnoaa[,tvar] - mean(tnoaa[tnoaa$Year >= bmin & tnoaa$Year <= bmax,tvar],na.rm=T)
tgis[,tvar]<- tgis[,tvar] - mean(tgis[tgis$Year >= bmin & tgis$Year <= bmax,tvar],na.rm=T)
thad[,tvar] <- thad[,tvar] - mean(thad[thad$Year >= bmin & thad$Year <= bmax,tvar],na.rm=T)

# Combine Temperature and Sea Ice Observations into Single Dataframes
tbest$Model <- "BEST"
tnoaa$Model <- "NOAAGlobalTemp"
tgis$Model <- "GISTEMP"
thad$Model <- "HADCRUT"

temperature <- rbind(tbest,tnoaa,tgis,thad) 
seaice <- rbind(osi,nasa,boot)

# Identify first member of each model
cmip6$First <- 0
tcmip6$First <- 0
for (mod in modfams){
  cmip6[cmip6$Family == mod & cmip6$Member == min(cmip6[cmip6$Family == mod,]$Member),]$First <- 1
  tcmip6[tcmip6$Family == mod & tcmip6$Member == min(tcmip6[tcmip6$Family == mod,]$Member),]$First <- 1
}

###############
# ANALYZE DATA 
################

#### Observational Temperature Sensitivity (1979-2013) ####
### Sensitivity for each ice/temperature model combination ###
obssen <- data.frame(TempModel=rep(c("BEST","NOAAGlobalTemp","GISTEMP","HADCRUT"),rep(3*length(regs),4)),
                     IceModel=rep(rep(c("Bootstrap","NASATeam","OSISAF"),rep(length(regs),3)),4),
                     Region= rep(regs,12),OPCavg=NA,LRDlt=NA,FADgt=NA)
for (tempmod in unique(temperature$Model)){
  x <- temperature[temperature$Model == tempmod & temperature$Year >= tymin & temperature$Year <= tymax,tvar]
  for (icemod in unique(seaice$Model)){
    y <- seaice[seaice$Model == icemod & seaice$Year >= tymin & seaice$Year <= tymax,]
    for (reg in regs){
      for (var in vars){
        obssen[obssen$TempModel == tempmod & obssen$IceModel == icemod & obssen$Region == reg,var] <- lm(y[y$Region == reg,var]~x)$coefficients[2]
}}}}

# boxplot(OPCavg~Region,data=obssen)
# boxplot(LRDlt~Region,data=obssen)
# boxplot(FADgt~Region,data=obssen)

#### Average Sensitivity & Spread of Sensitivity ####
obssenAvg <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)
obssenSD <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)

for (reg in regs){
  for (var in vars){
    obssenAvg[obssenAvg$Region == reg,var] <- mean(obssen[obssen$Region == reg,var])
    obssenSD[obssenSD$Region == reg,var] <- (max(obssen[obssen$Region == reg,var]) - min(obssen[obssen$Region == reg,var]))/2
}}

#### CMIP6 Temperature Sensitivity (1979-2013) ####
### Sensitivity for each model ###
c6senMem <- data.frame(Model=rep(mods,rep(length(regs),length(mods))),Family=NA,Member=NA,Region=rep(regs,length(mods)),OPCavg=NA,LRDlt=NA,FADgt=NA)
for (mod in mods){
  c6senMem[c6senMem$Model == mod,"Member"] <- cmip6[cmip6$Model == mod,"Member"][1]
  c6senMem[c6senMem$Model == mod,"Family"] <- cmip6[cmip6$Model == mod,"Family"][1]
  x <- tcmip6[tcmip6$Model == mod & tcmip6$Year >= tymin & tcmip6$Year <= tymax,tvar]
  y <- cmip6[cmip6$Model == mod & cmip6$Year >= tymin & cmip6$Year <= tymax,]
  for (reg in regs){
    for (var in vars){
      c6senMem[c6senMem$Region == reg & c6senMem$Model == mod,var] <- lm(y[y$Region == reg,var]~x)$coefficients[2]
}}}

c6senMem$First <- 0
for (mod in modfams){
  c6senMem[c6senMem$Family == mod & c6senMem$Member == min(c6senMem[c6senMem$Family == mod,]$Member),]$First <- 1
}

#### Average Sensitivity (1st Member Only) & Std Dev of Sensitivity ####
c6senAvg <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)
c6senSD <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)
c6senSDS <- data.frame(Family=rep(modfams,rep(length(regs),length(modfams))),Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)

for(reg in regs){
  for (var in vars){
    for (fam in modfams){
      if (length(unique(c6senMem[c6senMem$Family == fam,]$Member)) > 2){
        c6senSDS[c6senSDS$Family == fam & c6senSDS$Region == reg,var] <- sd_unbiased(c6senMem[c6senMem$Family == fam & c6senMem$Region == reg,var])
      }}
    c6senSD[c6senSD$Region == reg,var] <- mean(c6senSDS[c6senSDS$Region == reg & is.finite(c6senSDS[,var]),var])
    c6senAvg[c6senAvg$Region == reg,var] <- mean(c6senMem[c6senMem$Region == reg & c6senMem$First == 1,var])
  }}

#### Overall Uncertainty and Maxes & Mins for Plotting ####
### Overall Uncertainty
obssenErr <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)
for (reg in regs){
  for (var in c("OPCavg","LRDlt","FADgt")){
    obssenErr[obssenErr$Region == reg,var] <- 2*sqrt((obssenSD[obssenSD$Region == reg,var]**2) + (c6senSD[c6senSD$Region == reg,var]**2) )
  }}

### Min & Max Plotting Points for gray shading
c6senMinMax <- data.frame(Region=regs,Family=rep(modfams,rep(length(regs),length(modfams))),
                          OPCavgmax=NA,LRDltmax=NA,FADgtmax=NA,OPCavgmin=NA,LRDltmin=NA,FADgtmin=NA)
for (fam in modfams){
  if (length(unique(cmip6[cmip6$Family == fam,]$Member)) > 2){
    for (var in c("OPCavg","LRDlt","FADgt")){
      c6senMinMax[c6senMinMax$Family == fam,paste0(var,"min")] <- obssenAvg[,var] - c6senSDS[c6senSDS$Family == fam,var]
      c6senMinMax[c6senMinMax$Family == fam,paste0(var,"max")] <- obssenAvg[,var] + c6senSDS[c6senSDS$Family == fam,var]
    }}}
c6senMinMax <- c6senMinMax[is.na(c6senMinMax$OPCavgmax) == 0,]

### Min & Max Plotting Points for ± 2 sigma markers
obssenAvgMin <- obssenAvg - obssenErr
obssenAvgMax <- obssenAvg + obssenErr
obssenAvgMin$Region <- obssenAvg$Region
obssenAvgMax$Region <- obssenAvg$Region

# Write Files to Disk
write.csv(c6senMem,paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_",tvar,"_historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(c6senAvg,paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_MultiModelMean_",tvar,"_historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(c6senMinMax,paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_MultiModelMinMax_",tvar,"_historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(c6senSDS,paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_ModelStdDev_",tvar,"_historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(c6senSD,paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_MultiModelStdDev_",tvar,"_historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(obssenAvg,paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/ObsTempSens_MultiModelMean_",tvar,"_historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(obssenSD,paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/ObsTempSens_MultiModelStdDev_",tvar,"_historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(obssenErr,paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/ObsTempSens_MultiModelError_",tvar,"_historical_C",ct,"_V",ver,".csv"),row.names=F)
