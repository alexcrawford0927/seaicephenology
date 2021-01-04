# Purpose: Calculates the sensitivity of sea ice parameters to a unit change in
# regional or global temperature in CMIP6 emissions scenarios. Includes both
# a slope and y-intercept... calculates multi-model mean and the standard 
# deviation of first ensemble members.

# Define function
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
experiment = 'ssp585'
vars = c("OPCavg","LRDlt","FADgt")
tvar= 'tglobal'
path = '/Volumes/Troilus/CMIP6'
cmip6path = paste0(path,"/RegionalStats")
figpath = paste0("/Volumes/Troilus/CMIP6/Figures/RegionalStats/C",ct,"_V",ver)

tymin = 1850
tymax = 2099

modstokeep = c("ACCESS-CM2","BCC-CSM2-MR","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-CM6-1-HR",
               "CNRM-ESM2-1","EC-Earth3","EC-Earth3-Veg","GFDL-CM4","IPSL-CM6A-LR","MIROC-ES2L","MIROC6",
               "MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM","UKESM1-0-LL")

# Load Files
si <- read.csv(paste0(cmip6path,"/V",ver,"/ByTemp/CMIP6_SeaIceByTemp_Regionalized_",experiment,"_C",ct,"_V",ver,".csv"),as.is=T)

###############
# PROCESSING DATA
################
# Convert model names to characters
si$Model <- as.character(si$Model)
si$Family <- as.character(si$Family)

# Only want finite temperatures
si <- si[is.finite(si[,paste0(tvar,"A")])==1,]

# Re-name "Models" to account for differences
si$Model <- with(si, paste0(Family,'--',Member))

# Subset to chosen models #
si = si[si$Family %in% modstokeep,]

# Store Unique Values
allmods = unique(si$Model)
mods = unique(si$Model)
regs = unique(si$Region)
years = unique(si$Year)
modfams <- unique(si$Family)

# Identify first member of each model
si$First <- 0
for (mod in modfams){
  si[si$Family == mod & si$Member == min(si[si$Family == mod,]$Member),]$First <- 1
}

####################################
########## SENSITIVITY #############
####################################

#### CMIP6 Temperature Sensitivity (Whole Record) ####
### Sensitivity for each model ###
senMem <- data.frame(Model=rep(mods,rep(length(regs),length(mods))),Family=NA,Member=NA,Region=rep(regs,length(mods)),OPCavg=NA,LRDlt=NA,FADgt=NA)
intMem <- data.frame(Model=rep(mods,rep(length(regs),length(mods))),Family=NA,Member=NA,Region=rep(regs,length(mods)),OPCavg=NA,LRDlt=NA,FADgt=NA)
for (mod in mods){
  senMem[senMem$Model == mod,"Member"] <- si[si$Model == mod,"Member"][1]
  senMem[senMem$Model == mod,"Family"] <- si[si$Model == mod,"Family"][1]
  
  intMem[intMem$Model == mod,"Member"] <- si[si$Model == mod,"Member"][1]
  intMem[intMem$Model == mod,"Family"] <- si[si$Model == mod,"Family"][1]
  
  x <- si[si$Model == mod & si$Year >= tymin & si$Year <= tymax & si$Region == 0,paste0(tvar,"A")]
  for (reg in regs){
    y <- si[si$Model == mod & si$Year >= tymin & si$Year <= tymax,] 
    for (var in vars){
      senMem[senMem$Region == reg & senMem$Model == mod,var] <- lm(y[y$Region == reg,var]~x)$coefficients[2]
      intMem[intMem$Region == reg & intMem$Model == mod,var] <- lm(y[y$Region == reg,var]~x)$coefficients[1]
    }}}

# Identify first member of each model
senMem$First <- 0
intMem$First <- 0
for (mod in modfams){
  senMem[senMem$Family == mod & senMem$Member == min(senMem[senMem$Family == mod,]$Member),]$First <- 1
  intMem[intMem$Family == mod & intMem$Member == min(intMem[intMem$Family == mod,]$Member),]$First <- 1
}

#### Average Sensitivity (1st Member Only) & Std Dev of Sensitivity ####
senAvg <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)
senSD <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)
senSDS <- data.frame(Family=rep(modfams,rep(length(regs),length(modfams))),Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)

intAvg <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)
intSD <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)
intSDS <- data.frame(Family=rep(modfams,rep(length(regs),length(modfams))),Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)

for(reg in regs){
  for (var in vars){
    for (fam in modfams){
      if (length(unique(senMem[senMem$Family == fam,]$Member)) > 2){
        senSDS[senSDS$Family == fam & senSDS$Region == reg,var] <- sd_unbiased(senMem[senMem$Family == fam & senMem$Region == reg,var])
        intSDS[intSDS$Family == fam & intSDS$Region == reg,var] <- sd_unbiased(intMem[intMem$Family == fam & intMem$Region == reg,var])
      }}
    senSD[senSD$Region == reg,var] <- mean(senSDS[senSDS$Region == reg & is.finite(senSDS[,var]),var])
    senAvg[senAvg$Region == reg,var] <- mean(senMem[senMem$Region == reg & senMem$First == 1,var])
    
    intSD[intSD$Region == reg,var] <- mean(intSDS[intSDS$Region == reg & is.finite(intSDS[,var]),var])
    intAvg[intAvg$Region == reg,var] <- mean(intMem[intMem$Region == reg & intMem$First == 1,var])
  }}

# Write Files to Disk
write.csv(senMem,paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_",tvar,"_",experiment,"_C",ct,"_V",ver,".csv"),row.names = F)
write.csv(senAvg,paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_MultiModelMean_",tvar,"_",experiment,"_C",ct,"_V",ver,".csv"),row.names = F)
write.csv(senSDS,paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_ModelStdDev_",tvar,"_",experiment,"_C",ct,"_V",ver,".csv"),row.names = F)
write.csv(senSD,paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_MultiModelStdDev_",tvar,"_",experiment,"_C",ct,"_V",ver,".csv"),row.names = F)

write.csv(intMem,paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSensIntercept_",tvar,"_",experiment,"_C",ct,"_V",ver,".csv"),row.names = F)
write.csv(intAvg,paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSensIntercept_MultiModelMean_",tvar,"_",experiment,"_C",ct,"_V",ver,".csv"),row.names = F)
write.csv(intSDS,paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSensIntercept_ModelStdDev_",tvar,"_",experiment,"_C",ct,"_V",ver,".csv"),row.names = F)
write.csv(intSD,paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSensIntercept_MultiModelStdDev_",tvar,"_",experiment,"_C",ct,"_V",ver,".csv"),row.names = F)
