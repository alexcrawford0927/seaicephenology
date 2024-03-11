# Purpose: Zip Together Temperature and Sea Ice Data for the Observational Record

#################################
########## PREP WORK ############
#################################
# Define Variables
ct = 15
ver = 6
vars = c("OPCavg","LRDlt","FADgt")
path = '/Volumes/Troilus/CMIP6'
cmip6path = paste0(path,"/RegionalStats")

ymin = 1950
ymax = 2099

# Load Files
siobs <- read.csv(paste0(cmip6path,"/V",ver,"/ByYear/Obs_MultiModelMean_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
tanomobs <- read.csv('/Volumes/Troilus/SurfaceTemperature/ObsAvg_Annual_tas_1850-2019.csv',as.is=T)

###############
# PROCESSING DATA
################
# Subset to chosen models #
regs <- unique(siobs$Region)

### Zip Together Temperature and Sea Ice Data ###
# Observations
siobs$tregionA <- NA
siobs$tglobalA <- NA
for (y in unique(siobs$Year)){
  siobs[siobs$Year == y,"tregionA"] <- tanomobs[tanomobs$Year == y,"tregion"]
  siobs[siobs$Year == y,"tglobalA"] <- tanomobs[tanomobs$Year == y,"tglobal"]
}

write.csv(siobs,paste0(cmip6path,"/V",ver,"/ByTemp/Obs_MultiModelMean_Regionalized_C",ct,"_V",ver,".csv"),row.names=F)
