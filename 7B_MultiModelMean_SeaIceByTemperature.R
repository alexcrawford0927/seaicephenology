# Purpose: Zipping together temperature and sea ice data for CMIP6 models by a pre-defined temperature interval.

#################################
########## PREP WORK ############
#################################
# Define Variables
ct = 15
ver = 6
e = 1 # Which experiment to use

experiment = c('historical','ssp126','ssp245','ssp585')
tgmin = c(-0.4,0.8,0.8,0.8) # Minimum global temperature
tgmax = c(1.0,6.2,6.2,6.2) # Maximum global temperature
tgint = c(0.2,0.2,0.2,0.2) # Global temperature interval
trmin = c(-2,2,2,2) # Minimum regional temperature
trmax = c(5,16,16,16) # Maximum regional temperature
trint = c(0.25,0.25,0.25,0.25) # Regional temperature interval

vars = c("OPCavg","LRDlt","FADgt")
path = '/Volumes/Troilus/CMIP6'
cmip6path = paste0(path,"/RegionalStats")
figpath = paste0("/Volumes/Troilus/CMIP6/Figures/RegionalStats/C",ct,"_V",ver)

ymin = 1950
ymax = 2099

bmin = 1850
bmax = 1900

modstokeep = c("ACCESS-CM2","BCC-CSM2-MR","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-CM6-1-HR",
               "CNRM-ESM2-1","EC-Earth3","EC-Earth3-Veg","GFDL-CM4","IPSL-CM6A-LR","MIROC-ES2L","MIROC6",
               "MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM","UKESM1-0-LL")

# Load Files
sihist <- read.csv(paste0(cmip6path,"/CMIP6/CMIP6_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
sissp <- read.csv(paste0(cmip6path,"/CMIP6/CMIP6_Regionalized_",experiment[e],"_C",ct,"_V",ver,".csv"),as.is=T)

thist <- read.csv(paste0(path,'/CMIP6_Annual_tas_historical.csv'),as.is=T)
tssp <- read.csv(paste0(path,'/CMIP6_Annual_tas_',experiment[e],'.csv'),as.is=T)

###############
# PROCESSING DATA
################
# Subset to chosen models #
sihist = sihist[sihist$Family %in% modstokeep,]
thist = thist[thist$Family %in% modstokeep,]
sissp = sissp[sissp$Family %in% modstokeep,]
tssp = tssp[tssp$Family %in% modstokeep,]

# Re-name "Models" to account for differences
thist$Model <- with(thist, paste0(Family,'--',Member))
sihist$Model <- with(sihist, paste0(Family,'--',Member))
tssp$Model <- with(tssp, paste0(Family,'--',Member))
sissp$Model <- with(sissp, paste0(Family,'--',Member))

sspmods <- intersect( unique(tssp$Model), unique(sissp$Model) )
histmods <- intersect( unique(thist$Model), unique(sihist$Model) )
regs <- unique(sissp$Region)

###############
# ANALYZING DATA
################
### Calculate the temperature baseline for each model ###
## Only for models that have data in both ssp and hist
tclim <- data.frame(Model=histmods,Family=NA,Member=NA,tregion=NA,tglobal=NA)
for (mod in histmods){
  tclim[tclim$Model == mod,"Family"] <- thist[thist$Model == mod,"Family"][1]
  tclim[tclim$Model == mod,"Member"] <- thist[thist$Model == mod,"Member"][1]
  tclim[tclim$Model == mod,"tregion"] <- mean(thist[thist$Model == mod & thist$Year >= bmin & thist$Year <= bmax,"tregion"])
  tclim[tclim$Model == mod,"tglobal"] <- mean(thist[thist$Model == mod & thist$Year >= bmin & thist$Year <= bmax,"tglobal"])
}
tclim$Model <- as.character(tclim$Model)
tclimmods <- unique(tclim$Model)

# Subset datasets to align
tssp <- tssp[tssp$Model %in% sspmods & tssp$Year >= ymin & tssp$Year <= ymax & tssp$Year != 2014,]
sissp <- sissp[sissp$Model %in% sspmods & sissp$Year >= ymin & sissp$Year <= ymax & sissp$Year != 2014,]

### Calculate temperature anomalies for each model year ###
tssp$tregionA <- NA
tssp$tglobalA <- NA
for (mod in intersect(sspmods,tclimmods)){
    tssp[tssp$Model == mod,"tregionA"] <- tssp[tssp$Model == mod,"tregion"] - tclim[tclim$Model == mod,"tregion"]
    tssp[tssp$Model == mod,"tglobalA"] <- tssp[tssp$Model == mod,"tglobal"] - tclim[tclim$Model == mod,"tglobal"]
}

### Zip Together Temperature and Sea Ice Data ###
sissp$tregionA <- NA
sissp$tglobalA <- NA
for (mod in intersect(sspmods,tclimmods)){
  sissp[sissp$Model == mod,]$tregionA <- rep( tssp[tssp$Model == mod,]$tregionA, rep(length(regs),nrow(tssp[tssp$Model == mod,])) )
  sissp[sissp$Model == mod,]$tglobalA <- rep( tssp[tssp$Model == mod,]$tglobalA, rep(length(regs),nrow(tssp[tssp$Model == mod,])) )
}

# Identify first member of each model
sissp$First <- 0
for (mod in unique(sissp$Family)){
  sissp[sissp$Family == mod & sissp$Member == min(sissp[sissp$Family == mod,]$Member),]$First <- 1
}

write.csv(sissp,paste0(cmip6path,"/V",ver,"/ByTemp/CMIP6_SeaIceByTemp_Regionalized_",experiment[e],"_C",ct,"_V",ver,".csv"),row.names=F)

#######################
#### REGIONAL-ONLY ####
#######################
### Calculate Avg & SD of all sea ice data relative to each $tint$ degrees of warming by region ###
trange = seq(trmin[e],trmax[e],trint[e])
siavgssp <- data.frame(Model="CMIP6 Multi-Model Mean",Region=regs,tanom=rep( trange , rep(length(regs),length( trange )) ),
                      OPCavg=NA,LRDlt=NA,FADgt=NA)
sisdssp <- data.frame(Model="CMIP6 Multi-Model SD",Region=regs,tanom=rep( trange , rep(length(regs),length( trange )) ),
                       OPCavg=NA,LRDlt=NA,FADgt=NA)
for (reg in regs){
  for (t in trange){
    for (var in vars){
      siavgssp[siavgssp$tanom == t & siavgssp$Region == reg,var] <- mean(sissp[sissp$First == 1 & sissp$tregionA >= t & sissp$tregionA < t+trint[e] & sissp$Region == reg,var])
      sisdssp[sisdssp$tanom == t & sisdssp$Region == reg,var] <- sd(sissp[sissp$First == 1 & sissp$tregionA >= t & sissp$tregionA < t+trint[e] & sissp$Region == reg,var])
    }}}

# Create a ribbon for plotting
siminmaxssp <- data.frame(Region=regs,tanom=rep(trange,rep(length(regs),length(trange))),
                          OPCavgmax=NA,LRDltmax=NA,FADgtmax=NA,OPCavgmin=NA,LRDltmin=NA,FADgtmin=NA)
for (var in vars){
  siminmaxssp[,paste0(var,"min")] <- siavgssp[,var] - sisdssp[,var]
  siminmaxssp[,paste0(var,"max")] <- siavgssp[,var] + sisdssp[,var]
  
  if (var == "OPCavg"){siminmaxssp[,paste0(var,"max")] <- ifelse(siminmaxssp[,paste0(var,"max")] > 365, 365, siminmaxssp[,paste0(var,"max")])}
  if (var != "OPCavg"){siminmaxssp[,paste0(var,"max")] <- ifelse(siminmaxssp[,paste0(var,"max")] > 100, 100, siminmaxssp[,paste0(var,"max")])}
  siminmaxssp[,paste0(var,"min")] <- ifelse(siminmaxssp[,paste0(var,"min")] < 0, 0, siminmaxssp[,paste0(var,"min")])
}

write.csv(siavgssp,paste0(cmip6path,"/V",ver,"/ByTemp/CMIP6_MultiModelMean_Regionalized_",experiment[e],"_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(sisdssp,paste0(cmip6path,"/V",ver,"/ByTemp/CMIP6_MultiModelSD_Regionalized_",experiment[e],"_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(siminmaxssp,paste0(cmip6path,"/V",ver,"/ByTemp/CMIP6_MultiModelMinMax_Regionalized_",experiment[e],"_C",ct,"_V",ver,".csv"),row.names=F)

##########################
#### GLOBAL TEMP ONLY ####
##########################
sissp <- read.csv(paste0(cmip6path,"/V",ver,"/ByTemp/CMIP6_SeaIceByTemp_Regionalized_",experiment[e],"_C",ct,"_V",ver,".csv"),as.is = T)

### Calculate Avg & SD of all sea ice data relative to each $tgint$ degrees of warming by region ###
trange = seq(tgmin[e],tgmax[e],tgint[e])

siavgssp <- data.frame(Model="CMIP6 Multi-Model Mean",Region=regs,tanom=rep( trange , rep(length(regs),length( trange )) ),
                       OPCavg=NA,LRDlt=NA,FADgt=NA)
sisdssp <- data.frame(Model="CMIP6 Multi-Model SD",Region=regs,tanom=rep( trange , rep(length(regs),length( trange )) ),
                      OPCavg=NA,LRDlt=NA,FADgt=NA)
for (reg in regs){
  for (t in trange){
    for (var in vars){
      siavgssp[siavgssp$tanom == t & siavgssp$Region == reg,var] <- mean(sissp[sissp$First == 1 & sissp$tglobalA >= t & sissp$tglobalA < t+tgint[e] & sissp$Region == reg,var])
      sisdssp[sisdssp$tanom == t & sisdssp$Region == reg,var] <- sd(sissp[sissp$First == 1 & sissp$tglobalA >= t & sissp$tglobalA < t+tgint[e] & sissp$Region == reg,var])
    }}}

# Create a ribbon for plotting
siminmaxssp <- data.frame(Region=regs,tanom=rep(trange,rep(length(regs),length(trange))),
                          OPCavgmax=NA,LRDltmax=NA,FADgtmax=NA,OPCavgmin=NA,LRDltmin=NA,FADgtmin=NA)
for (var in vars){
  siminmaxssp[,paste0(var,"min")] <- siavgssp[,var] - sisdssp[,var]
  siminmaxssp[,paste0(var,"max")] <- siavgssp[,var] + sisdssp[,var]
  
  if (var == "OPCavg"){siminmaxssp[,paste0(var,"max")] <- ifelse(siminmaxssp[,paste0(var,"max")] > 365, 365, siminmaxssp[,paste0(var,"max")])}
  if (var != "OPCavg"){siminmaxssp[,paste0(var,"max")] <- ifelse(siminmaxssp[,paste0(var,"max")] > 100, 100, siminmaxssp[,paste0(var,"max")])}
  siminmaxssp[,paste0(var,"min")] <- ifelse(siminmaxssp[,paste0(var,"min")] < 0, 0, siminmaxssp[,paste0(var,"min")])
}

write.csv(siavgssp,paste0(cmip6path,"/V",ver,"/ByTemp/CMIP6_MultiModelMean_Globalized_",experiment[e],"_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(sisdssp,paste0(cmip6path,"/V",ver,"/ByTemp/CMIP6_MultiModelSD_Globalized_",experiment[e],"_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(siminmaxssp,paste0(cmip6path,"/V",ver,"/ByTemp/CMIP6_MultiModelMinMax_Globalized_",experiment[e],"_C",ct,"_V",ver,".csv"),row.names=F)
