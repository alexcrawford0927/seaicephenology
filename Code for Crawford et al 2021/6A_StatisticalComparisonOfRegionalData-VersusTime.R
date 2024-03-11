# Purpose: Calculates multi-model mean and internal variability of the temporal 
# average and trend in historical experiments and observations.

# Load Modules
library(ggplot2)
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
vars = c("OPCavg","OPCgt","LRDlt","FADgt")
cmip6path = "/Volumes/Troilus/CMIP6/RegionalStats"
ssmipath = "/Volumes/Troilus/CMIP6/RegionalStats/SSMI"

modstokeep = c("ACCESS-CM2","BCC-CSM2-MR","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-CM6-1-HR",
               "CNRM-ESM2-1","EC-Earth3","EC-Earth3-Veg","GFDL-CM4","IPSL-CM6A-LR","MIROC-ES2L","MIROC6",
               "MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM","UKESM1-0-LL")

aymin = 1979
aymax = 2013

tymin = 1979
tymax = 2013

# Read in Data
thresh = read.csv(paste0("/Volumes/Troilus/CMIP6/Regional_DOY_Thresholds25_C50_V2.csv"))
boot <- read.csv(paste0(ssmipath,"/Bootstrap_Regionalized_Historical_C",ct,"_V",ver,".csv"))
osi <- read.csv(paste0(ssmipath,"/OSISAF_Regionalized_Historical_C",ct,"_V",ver,".csv"))
nasa <- read.csv(paste0(ssmipath,"/NASATeam_Regionalized_Historical_C",ct,"_V",ver,".csv"))
cmip6 <- read.csv(paste0(cmip6path,"/CMIP6/CMIP6_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is = T)

###############
# PROCESSING DATA 
################
# Subset to chosen models #
cmip6 = cmip6[cmip6$Family %in% modstokeep,]

# Store Unique Values
cmip6$Model <- as.character(cmip6$Model)
allmods = unique(cmip6$Model)
regs = unique(cmip6$Region)
modfams <- unique(cmip6$Family)

# Identify first member of each model
cmip6$First <- 0
for (mod in modfams){
  cmip6[cmip6$Family == mod & cmip6$Member == min(cmip6[cmip6$Family == mod,]$Member),]$First <- 1
}

####################################
# ANALYZING DATA - AVERAGE & TREND #
####################################

# Step 1: Calculate the Average and Trend for each Ensemble Member
opc <- data.frame(Family=NA,Member=NA,Model=rep(allmods,length(regs)),Region=rep(regs,rep(length(allmods),length(regs))),Avg=0,SD=0,Trend=0,Conf95=0)
lrd <- data.frame(Family=NA,Member=NA,Model=rep(allmods,length(regs)),Region=rep(regs,rep(length(allmods),length(regs))),Avg=0,SD=0,Trend=0,Conf95=0)
fad <- data.frame(Family=NA,Member=NA,Model=rep(allmods,length(regs)),Region=rep(regs,rep(length(allmods),length(regs))),Avg=0,SD=0,Trend=0,Conf95=0)

obsopc <- data.frame(Family="SSMI",Member=rep(c("Bootstrap","NASATeam","OSISAF"),length(regs)),Model=rep(c("Bootstrap","NASATeam","OSISAF"),length(regs)),Region=rep(regs,rep(3,length(regs))),Avg=0,SD=0,Trend=0,Conf95=0)
obslrd <- data.frame(Family="SSMI",Member=rep(c("Bootstrap","NASATeam","OSISAF"),length(regs)),Model=rep(c("Bootstrap","NASATeam","OSISAF"),length(regs)),Region=rep(regs,rep(3,length(regs))),Avg=0,SD=0,Trend=0,Conf95=0)
obsfad <- data.frame(Family="SSMI",Member=rep(c("Bootstrap","NASATeam","OSISAF"),length(regs)),Model=rep(c("Bootstrap","NASATeam","OSISAF"),length(regs)),Region=rep(regs,rep(3,length(regs))),Avg=0,SD=0,Trend=0,Conf95=0)

for (rg in regs){
  # CMIP6 Models
  for (mod in allmods){
  # Calculate Average
  opcavg <- mean(cmip6[cmip6$Model == mod & cmip6$Region == rg & cmip6$Year >= aymin & cmip6$Year <= aymax,"OPCavg"])
  lrdavg <- mean(cmip6[cmip6$Model == mod & cmip6$Region == rg & cmip6$Year >= aymin & cmip6$Year <= aymax,"LRDlt"])
  fadavg <- mean(cmip6[cmip6$Model == mod & cmip6$Region == rg & cmip6$Year >= aymin & cmip6$Year <= aymax,"FADgt"])
  
  # Calculate S.D.
  opcsd <- sd(cmip6[cmip6$Model == mod & cmip6$Region == rg & cmip6$Year >= aymin & cmip6$Year <= aymax,"OPCavg"])
  lrdsd <- sd(cmip6[cmip6$Model == mod & cmip6$Region == rg & cmip6$Year >= aymin & cmip6$Year <= aymax,"LRDlt"])
  fadsd <- sd(cmip6[cmip6$Model == mod & cmip6$Region == rg & cmip6$Year >= aymin & cmip6$Year <= aymax,"FADgt"])
  
  # Calculate Trend
  input <- cmip6[cmip6$Model == mod & cmip6$Region == rg & cmip6$Year >= tymin & cmip6$Year <= tymax,]
  opclmd <- lm(OPCavg~Year,input) 
  lrdlmd <- lm(LRDlt~Year,input) 
  fadlmd <- lm(FADgt~Year,input) 
  
  # Store Results
  opc[opc$Region == rg & opc$Model == mod,]$Member <-  cmip6[cmip6$Model == mod & cmip6$Region == rg,]$Member[1]
  opc[opc$Region == rg & opc$Model == mod,]$Family <-  cmip6[cmip6$Model == mod & cmip6$Region == rg,]$Family[1]
  opc[opc$Region == rg & opc$Model == mod,]$Avg <-  opcavg
  opc[opc$Region == rg & opc$Model == mod,]$SD <-  opcsd
  opc[opc$Region == rg & opc$Model == mod,]$Trend <- opclmd$coefficients[2]
  opc[opc$Region == rg & opc$Model == mod,]$Conf95 <- confint(opclmd)[2,2]-opclmd$coefficients[2]
  
  lrd[lrd$Region == rg & lrd$Model == mod,]$Member <-  cmip6[cmip6$Model == mod & cmip6$Region == rg,]$Member[1]
  lrd[lrd$Region == rg & lrd$Model == mod,]$Family <-  cmip6[cmip6$Model == mod & cmip6$Region == rg,]$Family[1]
  lrd[lrd$Region == rg & lrd$Model == mod,]$Avg <-  lrdavg
  lrd[lrd$Region == rg & lrd$Model == mod,]$SD <-  lrdsd
  lrd[lrd$Region == rg & lrd$Model == mod,]$Trend <- lrdlmd$coefficients[2]
  lrd[lrd$Region == rg & lrd$Model == mod,]$Conf95 <- confint(lrdlmd)[2,2]-lrdlmd$coefficients[2]
  
  fad[fad$Region == rg & fad$Model == mod,]$Member <-  cmip6[cmip6$Model == mod & cmip6$Region == rg,]$Member[1]
  fad[fad$Region == rg & fad$Model == mod,]$Family <-  cmip6[cmip6$Model == mod & cmip6$Region == rg,]$Family[1]
  fad[fad$Region == rg & fad$Model == mod,]$Avg <-  fadavg
  fad[fad$Region == rg & fad$Model == mod,]$SD <-  fadsd
  fad[fad$Region == rg & fad$Model == mod,]$Trend <- fadlmd$coefficients[2]
  fad[fad$Region == rg & fad$Model == mod,]$Conf95 <- confint(fadlmd)[2,2]-fadlmd$coefficients[2]
  }
  
  # Observations
  for (mod in c("Bootstrap","NASATeam","OSISAF")){
    if(mod == "Bootstrap"){df <- boot}
    if(mod == "NASATeam"){df <- nasa}
    if(mod == "OSISAF"){df <- osi}
  # Calculate Average
  opcavg <- mean(df[df$Region == rg & df$Year >= aymin & df$Year <= aymax,"OPCavg"])
  lrdavg <- mean(df[df$Region == rg & df$Year >= aymin & df$Year <= aymax,"LRDlt"])
  fadavg <- mean(df[df$Region == rg & df$Year >= aymin & df$Year <= aymax,"FADgt"])
 
  # Calculate S.D.
  opcsd <- sd(df[df$Region == rg & df$Year >= aymin & df$Year <= aymax,"OPCavg"])
  lrdsd <- sd(df[df$Region == rg & df$Year >= aymin & df$Year <= aymax,"LRDlt"])
  fadsd <- sd(df[df$Region == rg & df$Year >= aymin & df$Year <= aymax,"FADgt"])  
  
  # Calculate Trend
  input <- df[df$Region == rg & df$Year >= tymin & df$Year <= tymax,]
  opclmd <- lm(OPCavg~Year,input) 
  lrdlmd <- lm(LRDlt~Year,input) 
  fadlmd <- lm(FADgt~Year,input)  
  
  # Store Results
  obsopc[obsopc$Region == rg & obsopc$Model == mod,]$Avg <-  opcavg
  obsopc[obsopc$Region == rg & obsopc$Model == mod,]$SD <-  opcsd
  obsopc[obsopc$Region == rg & obsopc$Model == mod,]$Trend <- opclmd$coefficients[2]
  obsopc[obsopc$Region == rg & obsopc$Model == mod,]$Conf95 <- confint(opclmd)[2,2]-opclmd$coefficients[2]
  
  obslrd[obslrd$Region == rg & obslrd$Model == mod,]$Avg <-  lrdavg
  obslrd[obslrd$Region == rg & obslrd$Model == mod,]$SD <-  lrdsd
  obslrd[obslrd$Region == rg & obslrd$Model == mod,]$Trend <- lrdlmd$coefficients[2]
  obslrd[obslrd$Region == rg & obslrd$Model == mod,]$Conf95 <- confint(lrdlmd)[2,2]-lrdlmd$coefficients[2]
  
  obsfad[obsfad$Region == rg & obsfad$Model == mod,]$Avg <-  fadavg
  obsfad[obsfad$Region == rg & obsfad$Model == mod,]$SD <-  fadsd
  obsfad[obsfad$Region == rg & obsfad$Model == mod,]$Trend <- fadlmd$coefficients[2]
  obsfad[obsfad$Region == rg & obsfad$Model == mod,]$Conf95 <- confint(fadlmd)[2,2]-fadlmd$coefficients[2]
  }
}

opc$First <- 0
lrd$First <- 0
fad$First <- 0
for (mod in modfams){
  opc[opc$Family == mod & opc$Member == min(opc[opc$Family == mod,]$Member),]$First <- 1
  lrd[lrd$Family == mod & lrd$Member == min(lrd[lrd$Family == mod,]$Member),]$First <- 1
  fad[fad$Family == mod & fad$Member == min(fad[fad$Family == mod,]$Member),]$First <- 1
  
}

# Step 2: Calculate the multi-model mean for CMIP6 and for Observations
c6avgAvg <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)
c6avgTrd <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)
obsavgAvg <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)
obsavgTrd <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)

for (reg in regs){
    for (var in c("OPCavg","LRDlt","FADgt")){
      if (var == "OPCavg"){
        dfc <- opc[opc$First == 1 & opc$Region == reg,]
        dfo <- obsopc[obsopc$Region == reg,]}
      if (var == "LRDlt"){
        dfc <- lrd[lrd$First == 1 & lrd$Region == reg,]
        dfo <- obslrd[obslrd$Region == reg,]}
      if (var == "FADgt"){
        dfc <- fad[fad$First == 1 & fad$Region == reg,]
        dfo <- obsfad[obsfad$Region == reg,]}
      
      c6avgAvg[c6avgAvg$Region == reg,var] <- mean(dfc[dfc$First == 1,"Avg"])
      c6avgTrd[c6avgTrd$Region == reg,var] <- mean(dfc[dfc$First == 1,"Trend"])
      
      obsavgAvg[obsavgAvg$Region == reg,var] <- mean(dfo$Avg)
      obsavgTrd[obsavgTrd$Region == reg,var] <- mean(dfo$Trend)
}}

# Step 3: Calculate a sigma for each CMIP6 model ensemble and for observations
c6sdsAvg <- data.frame(Region=regs,Family=rep(modfams,rep(length(regs),length(modfams))),OPCavg=NA,LRDlt=NA,FADgt=NA)
c6sdsTrd <- data.frame(Region=regs,Family=rep(modfams,rep(length(regs),length(modfams))),OPCavg=NA,LRDlt=NA,FADgt=NA)
c6sdAvg <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)
c6sdTrd <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)
obssdAvg <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)
obssdTrd <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)

for (reg in regs){
  for (var in c("OPCavg","LRDlt","FADgt")){
    if (var == "OPCavg"){
      dfc <- opc[opc$Region == reg,]
      dfo <- obsopc[obsopc$Region == reg,]}
    if (var == "LRDlt"){
      dfc <- lrd[lrd$Region == reg,]
      dfo <- obslrd[obslrd$Region == reg,]}
    if (var == "FADgt"){
      dfc <- fad[fad$Region == reg,]
      dfo <- obsfad[obsfad$Region == reg,]}
    
    for (fam in modfams){
        if (length(unique(cmip6[cmip6$Family == fam,]$Member)) > 2){
          c6sdsAvg[c6sdsAvg$Family == fam & c6sdsAvg$Region == reg,var] <- sd_unbiased(dfc[dfc$Family == fam,"Avg"])
          c6sdsTrd[c6sdsTrd$Family == fam & c6sdsTrd$Region == reg,var] <- sd_unbiased(dfc[dfc$Family == fam,"Trend"])
        }}
    
      c6sdAvg[c6sdAvg$Region == reg,var] <- mean(c6sdsAvg[c6sdsAvg$Region == reg & is.finite(c6sdsAvg[,var]),var])
      c6sdTrd[c6sdTrd$Region == reg,var] <- mean(c6sdsTrd[c6sdsTrd$Region == reg & is.finite(c6sdsTrd[,var]),var])
      obssdAvg[obssdAvg$Region == reg,var] <- (max(dfo$Avg) - min(dfo$Avg))/2
      obssdTrd[obssdTrd$Region == reg,var] <- (max(dfo$Trend) - min(dfo$Trend))/2
  }}

# Step 4: Calculate an overall observational uncertainty
obserrAvg <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)
obserrTrd <- data.frame(Region=regs,OPCavg=NA,LRDlt=NA,FADgt=NA)

for (reg in regs){
    for (var in c("OPCavg","LRDlt","FADgt")){
      obserrAvg[obserrAvg$Region == reg,var] <- 2*sqrt((obssdAvg[obssdAvg$Region == reg,var]**2) + (c6sdAvg[c6sdAvg$Region == reg,var]**2) )
      obserrTrd[obserrTrd$Region == reg,var] <- 2*sqrt((obssdTrd[obssdTrd$Region == reg,var]**2) + (c6sdTrd[c6sdTrd$Region == reg,var]**2) )
    }}

c6minmaxAvg <- data.frame(Region=regs,Family=rep(modfams,rep(length(regs),length(modfams))),
                        OPCavgmax=NA,LRDltmax=NA,FADgtmax=NA,OPCavgmin=NA,LRDltmin=NA,FADgtmin=NA)
c6minmaxTrd <- data.frame(Region=regs,Family=rep(modfams,rep(length(regs),length(modfams))),
                          OPCavgmax=NA,LRDltmax=NA,FADgtmax=NA,OPCavgmin=NA,LRDltmin=NA,FADgtmin=NA)
for (fam in modfams){
  if (length(unique(cmip6[cmip6$Family == fam,]$Member)) > 2){
    for (var in c("OPCavg","LRDlt","FADgt")){
      c6minmaxAvg[c6minmaxAvg$Family == fam,paste0(var,"min")] <- obsavgAvg[,var] - c6sdsAvg[c6sdsAvg$Family == fam,var]
      c6minmaxAvg[c6minmaxAvg$Family == fam,paste0(var,"max")] <- obsavgAvg[,var] + c6sdsAvg[c6sdsAvg$Family == fam,var]
      c6minmaxTrd[c6minmaxTrd$Family == fam,paste0(var,"min")] <- obsavgTrd[,var] - c6sdsTrd[c6sdsTrd$Family == fam,var]
      c6minmaxTrd[c6minmaxTrd$Family == fam,paste0(var,"max")] <- obsavgTrd[,var] + c6sdsTrd[c6sdsTrd$Family == fam,var]
    }}}
c6minmaxAvg <- c6minmaxAvg[is.na(c6minmaxAvg$OPCavgmax) == 0,]
c6minmaxTrd <- c6minmaxTrd[is.na(c6minmaxTrd$OPCavgmax) == 0,]

# Cap the mins and maxes for the averages
for (r in 1:nrow(c6minmaxAvg)){
  for (c in 3:ncol(c6minmaxAvg)){
    if(c6minmaxAvg[r,c] < 0){c6minmaxAvg[r,c] <- 0}
  }
  for (c in c(4,5,7,8)){
    if(c6minmaxAvg[r,c] > 100){c6minmaxAvg[r,c] <- 100}
  }
  for (c in c(3,6)){
    if(c6minmaxAvg[r,c] > 365){c6minmaxAvg[r,c] <- 365}
  }
} # Cap mins and maxes

# Step 5: Make Error Bars for Plotting (capped at 0 and 365 or 0 and 100)
# Step 5A: It's more complicated for Average...
obsavgAvgmax <- obsavgAvg + obserrAvg
obsavgAvgmin <- obsavgAvg - obserrAvg
obsavgAvgmax$Region <- obsavgAvg$Region
obsavgAvgmin$Region <- obsavgAvg$Region
for (r in 1:nrow(obsavgAvgmin)){
  for (c in 2:ncol(obsavgAvgmin)){
    if(obsavgAvgmin[r,c] < 0){obsavgAvgmin[r,c] <- 0}
    if(obsavgAvgmax[r,c] < 0){obsavgAvgmax[r,c] <- 0}
  }
  for (c in 3:ncol(obsavgAvgmin)){
    if(obsavgAvgmin[r,c] > 100){obsavgAvgmin[r,c] <- 100}
    if(obsavgAvgmax[r,c] > 100){obsavgAvgmax[r,c] <- 100}    
  }
  for (c in 2:2){
    if(obsavgAvgmin[r,c] > 365){obsavgAvgmin[r,c] <- 365}
    if(obsavgAvgmax[r,c] > 365){obsavgAvgmax[r,c] <- 365}
  }
} # Cap mins and maxes

# Step 5B: It's simpler for Trend...
obsavgTrdmax <- obsavgTrd + obserrTrd
obsavgTrdmin <- obsavgTrd - obserrTrd
obsavgTrdmax$Region <- obsavgTrd$Region
obsavgTrdmin$Region <- obsavgTrd$Region

# Write Files to Disk
write.csv(opc,paste0(cmip6path,"/V",ver,"/Avg",aymin,"-",aymax,"/CMIP6_OPCavg_AvgTrend_Regionalized_Historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(lrd,paste0(cmip6path,"/V",ver,"/Avg",aymin,"-",aymax,"/CMIP6_LRDlt_AvgTrend_Regionalized_Historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(fad,paste0(cmip6path,"/V",ver,"/Avg",aymin,"-",aymax,"/CMIP6_FADgt_AvgTrend_Regionalized_Historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(c6avgAvg,paste0(cmip6path,"/V",ver,"/Avg",aymin,"-",aymax,"/CMIP6Avg_MultiModelMean_Regionalized_Historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(c6avgTrd,paste0(cmip6path,"/V",ver,"/Trend",tymin,"-",tymax,"/CMIP6Trend_MultiModelMean_Regionalized_Historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(c6minmaxAvg,paste0(cmip6path,"/V",ver,"/Avg",aymin,"-",aymax,"/CMIP6Avg_MultiModelMinMax_Regionalized_Historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(c6minmaxTrd,paste0(cmip6path,"/V",ver,"/Trend",tymin,"-",tymax,"/CMIP6Trend_MultiModelMinMax_Regionalized_Historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(obsavgTrd,paste0(cmip6path,"/V",ver,"/Trend",tymin,"-",tymax,"/ObsTrend_MultiModelMean_Regionalized_Historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(obserrTrd,paste0(cmip6path,"/V",ver,"/Trend",tymin,"-",tymax,"/ObsTrend_MultiModelError_Regionalized_Historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(obsavgAvg,paste0(cmip6path,"/V",ver,"/Avg",aymin,"-",aymax,"/ObsAvg_MultiModelMean_Regionalized_Historical_C",ct,"_V",ver,".csv"),row.names=F)
write.csv(obserrAvg,paste0(cmip6path,"/V",ver,"/Avg",aymin,"-",aymax,"/ObsAvg_MultiModelError_Regionalized_Historical_C",ct,"_V",ver,".csv"),row.names=F)
