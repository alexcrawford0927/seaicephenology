# Purpose: Determines the multi-model mean and standard deviation of first ensemble 
# members for each year of a CMIP6 experiment.

#################################
########## PREP WORK ############
#################################
# Define Variables
ct = 15
ver = 6
experiment = 'ssp245'
vars = c("OPCavg","LRDlt","FADgt")
path = '/Volumes/Troilus/CMIP6'
cmip6path = paste0(path,"/RegionalStats")
figpath = paste0("/Volumes/Troilus/CMIP6/Figures/RegionalStats/C",ct,"_V",ver)

tymin = 2015 # 1950 or 2015
tymax = 2099 # 2013 or 2099

modstokeep = c("ACCESS-CM2","BCC-CSM2-MR","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-CM6-1-HR",
               "CNRM-ESM2-1","EC-Earth3","EC-Earth3-Veg","GFDL-CM4","IPSL-CM6A-LR","MIROC-ES2L","MIROC6",
               "MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM","UKESM1-0-LL")

# Load Files
si <- read.csv(paste0(cmip6path,"/CMIP6/CMIP6_Regionalized_",experiment,"_C",ct,"_V",ver,".csv"),as.is=T)

###############
# PROCESSING DATA
################
# Convert model names to characters
si$Model <- as.character(si$Model)
si$Family <- as.character(si$Family)

# Re-name "Models" to account for differences
si$Model <- with(si, paste0(Family,'--',Member))

# Subset to chosen models #
si = si[si$Family %in% modstokeep,]

# Store Unique Values
allmods = unique(si$Model)
regs = unique(si$Region)
years = unique(si$Year)
modfams <- unique(si$Family)

# Identify first member of each model
si$First <- 0
for (mod in modfams){
  si[si$Family == mod & si$Member == min(si[si$Family == mod,]$Member),]$First <- 1
}

# Subset to first ensemble members
si <- si[si$First == 1,]

###############
# CALCULATING DATA
################
# Calculate a Multi-Model Mean from 1st Ensemble Members for the experiment
siavg <- data.frame(Model="CMIP6 Multi-Model Mean",Region=regs,Year=rep(years,rep(length(regs),length(years))),
                    OPCavg=NA,LRDlt=NA,FADgt=NA)
for (reg in regs){
  for (y in years){
    for (var in vars){
      siavg[siavg$Year == y & siavg$Region == reg,var] <- mean(si[si$Year == y & si$Region == reg,var])
    }}}

# Calculate the SD of the multi-model mean of each year
sisd <- data.frame(Model="CMIP6 Multi-Model SD",Region=regs,Year=rep(years,rep(length(regs),length(years))),
                   OPCavg=NA,LRDlt=NA,FADgt=NA)
for (reg in regs){
  for (y in years){
    for (var in vars){
      sisd[sisd$Year == y & sisd$Region == reg,var] <- sd(si[si$Year == y & si$Region == reg,var])
    }}}

# Create Ribbon for Plotting as the CMIP6 Average ± 2*CMIP6 Std. Dev. (1st Ensemble Members Only)
siminmax <- data.frame(Region=regs,Year=rep(years,rep(length(regs),length(years))),
                        OPCavgmax=NA,LRDltmax=NA,FADgtmax=NA,OPCavgmin=NA,LRDltmin=NA,FADgtmin=NA)
for (var in vars){
  siminmax[,paste0(var,"min")] <- siavg[,var] - sisd[,var]
  siminmax[,paste0(var,"max")] <- siavg[,var] + sisd[,var]
  
  if (var == "OPCavg"){siminmax[,paste0(var,"max")] <- ifelse(siminmax[,paste0(var,"max")] > 365, 365, siminmax[,paste0(var,"max")])}
  if (var != "OPCavg"){siminmax[,paste0(var,"max")] <- ifelse(siminmax[,paste0(var,"max")] > 100, 100, siminmax[,paste0(var,"max")])}
  siminmax[,paste0(var,"min")] <- ifelse(siminmax[,paste0(var,"min")] < 0, 0, siminmax[,paste0(var,"min")])
}

# Write to File
write.csv(siavg,paste0(cmip6path,"/V",ver,"/ByYear/CMIP6_MultiModelMean_Regionalized_",experiment,"_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(sisd,paste0(cmip6path,"/V",ver,"/ByYear/CMIP6_MultiModelSD_Regionalized_",experiment,"_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(siminmax,paste0(cmip6path,"/V",ver,"/ByYear/CMIP6_MultiModelMinMax_Regionalized_",experiment,"_C",ct,"_V",ver,".csv"),row.names=F)
