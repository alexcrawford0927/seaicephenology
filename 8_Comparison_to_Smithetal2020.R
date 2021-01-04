# Purpose: Test the sensitivity of results to the subset of models used and the 
# spatial domain used -- specifically by comparing results to the same domain 
# and/or model subset as Smith et al. (2020) 
# ("Seasonal transition dates can reveal biases in Arctic sea ice simulations"
# The Cryosphere, 2977-2997)

#################################
########## N66 Only #############
#################################
# Define Variables
ct = 15
ver = 6
cmip6path = "/Volumes/Troilus/CMIP6/RegionalStats/byLat"

aymin = 1979
aymax = 2013

##### My Subset of Models #####

obsavgAvg <- read.csv(paste0(cmip6path,"/Avg",aymin,"-",aymax,"/ObsAvg_MultiModelMean_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
opc <- read.csv(paste0(cmip6path,"/Avg",aymin,"-",aymax,"/CMIP6_OPCavg_AvgTrend_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
lrd <- read.csv(paste0(cmip6path,"/Avg",aymin,"-",aymax,"/CMIP6_LRDlt_AvgTrend_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
fad <- read.csv(paste0(cmip6path,"/Avg",aymin,"-",aymax,"/CMIP6_FADgt_AvgTrend_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)

opc <- opc[opc$First == 1,]
fad <- fad[fad$First == 1,]
lrd <- lrd[lrd$First == 1,]

hist(opc$Avg - obsavgAvg$OPCavg)
hist(lrd$Avg - obsavgAvg$LRDlt)
hist(fad$Avg - obsavgAvg$FADgt)

mean(opc$Avg - obsavgAvg$OPCavg)
median(lrd$Avg - obsavgAvg$LRDlt)
median(fad$Avg - obsavgAvg$FADgt)

##### Smith et al. (2020) Subset of Models #####

obsavgAvg <- read.csv(paste0(cmip6path,"/Avg",aymin,"-",aymax," - Smith et al. 2020 Subset/ObsAvg_MultiModelMean_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
opc <- read.csv(paste0(cmip6path,"/Avg",aymin,"-",aymax," - Smith et al. 2020 Subset/CMIP6_OPCavg_AvgTrend_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
lrd <- read.csv(paste0(cmip6path,"/Avg",aymin,"-",aymax," - Smith et al. 2020 Subset/CMIP6_LRDlt_AvgTrend_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
fad <- read.csv(paste0(cmip6path,"/Avg",aymin,"-",aymax," - Smith et al. 2020 Subset/CMIP6_FADgt_AvgTrend_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)

opc <- opc[opc$First == 1,]
fad <- fad[fad$First == 1,]
lrd <- lrd[lrd$First == 1,]

hist(opc$Avg - obsavgAvg$OPCavg)
hist(lrd$Avg - obsavgAvg$LRDlt)
hist(fad$Avg - obsavgAvg$FADgt)

mean(opc$Avg - obsavgAvg$OPCavg)
median(lrd$Avg - obsavgAvg$LRDlt)
median(fad$Avg - obsavgAvg$FADgt)

#################################
########## Pan-Arctic ###########
#################################
cmip6path = "/Volumes/Troilus/CMIP6/RegionalStats/V6"

obsavgAvg <- read.csv(paste0(cmip6path,"/Avg",aymin,"-",aymax,"/ObsAvg_MultiModelMean_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
opc <- read.csv(paste0(cmip6path,"/Avg",aymin,"-",aymax,"/CMIP6_OPCavg_AvgTrend_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
lrd <- read.csv(paste0(cmip6path,"/Avg",aymin,"-",aymax,"/CMIP6_LRDlt_AvgTrend_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
fad <- read.csv(paste0(cmip6path,"/Avg",aymin,"-",aymax,"/CMIP6_FADgt_AvgTrend_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)

opc <- opc[opc$First == 1 & opc$Region == 0,]
fad <- fad[fad$First == 1 & fad$Region == 0,]
lrd <- lrd[lrd$First == 1 & lrd$Region == 0,]

mean(opc$Avg - obsavgAvg[obsavgAvg$Region == 0,]$OPCavg)
median(lrd$Avg - obsavgAvg[obsavgAvg$Region == 0,]$LRDlt)
median(fad$Avg - obsavgAvg[obsavgAvg$Region == 0,]$FADgt)
