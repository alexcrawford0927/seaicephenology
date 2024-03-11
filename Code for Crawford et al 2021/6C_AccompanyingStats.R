# Purpose: Supporting statistical tests that determine accuracy/bias of particular CMIP6 models used in 6A and 6B scripts.

# Load Modules
library(ggplot2)
library(stringr)
library(gridExtra)

#################################
########## PREP WORK ############
#################################
# Define Variables
ct = 15
ver = 6
aymin = 1979
aymax = 2013

tymin = 1979
tymax = 2013

axistitlesize = 8
axistextsize = 7
xpos = 1

REGS2 = c('All Regions','Okhotsk','Bering','Hudson','Baffin','St. Lawrence','Labrador','Greenland',
          'Barents','Kara','Laptev','E. Siberian','Chukchi','Beaufort',
          'CAA','CAO') # 0 - 7, 8 - 13, 14 - 16

cmip6path = "/Volumes/Troilus/CMIP6/RegionalStats"
figpath = paste0("/Volumes/Troilus/CMIP6/Figures/RegionalStats/C",ct,"_V",ver)

#################################
########## PREP WORK ############
#################################
### Read in Data ###
thresh = read.csv(paste0("/Volumes/Troilus/CMIP6/Regional_DOY_Thresholds25_C50_V2.csv"))
opc <- read.csv(paste0(cmip6path,"/V",ver,"/Avg",aymin,"-",aymax,"/CMIP6_OPCavg_AvgTrend_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
lrd <- read.csv(paste0(cmip6path,"/V",ver,"/Avg",aymin,"-",aymax,"/CMIP6_LRDlt_AvgTrend_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
fad <- read.csv(paste0(cmip6path,"/V",ver,"/Avg",aymin,"-",aymax,"/CMIP6_FADgt_AvgTrend_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)

c6avgAvg <- read.csv(paste0(cmip6path,"/V",ver,"/Avg",aymin,"-",aymax,"/CMIP6Avg_MultiModelMean_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
c6avgTrd <- read.csv(paste0(cmip6path,"/V",ver,"/Trend",tymin,"-",tymax,"/CMIP6Trend_MultiModelMean_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
c6minmaxTrd <- read.csv(paste0(cmip6path,"/V",ver,"/Trend",tymin,"-",tymax,"/CMIP6Trend_MultiModelMinMax_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
c6minmaxAvg <- read.csv(paste0(cmip6path,"/V",ver,"/Avg",aymin,"-",aymax,"/CMIP6Avg_MultiModelMinMax_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
obsavgTrd <- read.csv(paste0(cmip6path,"/V",ver,"/Trend",tymin,"-",tymax,"/ObsTrend_MultiModelMean_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
obserrTrd <- read.csv(paste0(cmip6path,"/V",ver,"/Trend",tymin,"-",tymax,"/ObsTrend_MultiModelError_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
obsavgAvg <- read.csv(paste0(cmip6path,"/V",ver,"/Avg",aymin,"-",aymax,"/ObsAvg_MultiModelMean_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
obserrAvg <- read.csv(paste0(cmip6path,"/V",ver,"/Avg",aymin,"-",aymax,"/ObsAvg_MultiModelError_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)

c6senMem <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_tregion_historical_C",ct,"_V",ver,".csv"),as.is=T)
c6senAvg <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_MultiModelMean_tregion_historical_C",ct,"_V",ver,".csv"),as.is=T)
c6senMinMax <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_MultiModelMinMax_tregion_historical_C",ct,"_V",ver,".csv"),as.is=T)
c6senSDS <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_ModelStdDev_tregion_historical_C",ct,"_V",ver,".csv"),as.is=T)
c6senSD <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_MultiModelStdDev_tregion_historical_C",ct,"_V",ver,".csv"),as.is=T)
obssenAvg <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/ObsTempSens_MultiModelMean_tregion_historical_C",ct,"_V",ver,".csv"),as.is=T)
obssenSD <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/ObsTempSens_MultiModelStdDev_tregion_historical_C",ct,"_V",ver,".csv"),as.is=T)
obssenErr <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/ObsTempSens_MultiModelError_tregion_historical_C",ct,"_V",ver,".csv"),as.is=T)

### Make Error Bars ###
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

### Min & Max Plotting Points for ± 2 sigma markers
obssenAvgMin <- obssenAvg - obssenErr
obssenAvgMax <- obssenAvg + obssenErr
obssenAvgMin$Region <- obssenAvg$Region
obssenAvgMax$Region <- obssenAvg$Region

######### Re-Regionalize ###############
opc$Region2 <- ifelse(opc$Region == 0, 1, ifelse(opc$Region < 5, opc$Region, ifelse(opc$Region == 16, 5, opc$Region+1)))
lrd$Region2 <- ifelse(lrd$Region == 0, 1, ifelse(lrd$Region < 5, lrd$Region, ifelse(lrd$Region == 16, 5, lrd$Region+1)))
fad$Region2 <- ifelse(fad$Region == 0, 1, ifelse(fad$Region < 5, fad$Region, ifelse(fad$Region == 16, 5, fad$Region+1)))

c6minmaxAvg$Region2 <- ifelse(c6minmaxAvg$Region == 0, 1, ifelse(c6minmaxAvg$Region < 5, c6minmaxAvg$Region, ifelse(c6minmaxAvg$Region == 16, 5, c6minmaxAvg$Region+1)))
c6avgAvg$Region2 <- ifelse(c6avgAvg$Region == 0, 1, ifelse(c6avgAvg$Region < 5, c6avgAvg$Region, ifelse(c6avgAvg$Region == 16, 5, c6avgAvg$Region+1)))
obsavgAvgmax$Region2 <- ifelse(obsavgAvgmax$Region == 0, 1, ifelse(obsavgAvgmax$Region < 5, obsavgAvgmax$Region, ifelse(obsavgAvgmax$Region == 16, 5, obsavgAvgmax$Region+1)))
obsavgAvgmin$Region2 <- ifelse(obsavgAvgmin$Region == 0, 1, ifelse(obsavgAvgmin$Region < 5, obsavgAvgmin$Region, ifelse(obsavgAvgmin$Region == 16, 5, obsavgAvgmin$Region+1)))
obsavgAvg$Region2 <- ifelse(obsavgAvg$Region == 0, 1, ifelse(obsavgAvg$Region < 5, obsavgAvg$Region, ifelse(obsavgAvg$Region == 16, 5, obsavgAvg$Region+1)))

c6minmaxTrd$Region2 <- ifelse(c6minmaxTrd$Region == 0, 1, ifelse(c6minmaxTrd$Region < 5, c6minmaxTrd$Region, ifelse(c6minmaxTrd$Region == 16, 5, c6minmaxTrd$Region+1)))
c6avgTrd$Region2 <- ifelse(c6avgTrd$Region == 0, 1, ifelse(c6avgTrd$Region < 5, c6avgTrd$Region, ifelse(c6avgTrd$Region == 16, 5, c6avgTrd$Region+1)))
obsavgTrdmax$Region2 <- ifelse(obsavgTrdmax$Region == 0, 1, ifelse(obsavgTrdmax$Region < 5, obsavgTrdmax$Region, ifelse(obsavgTrdmax$Region == 16, 5, obsavgTrdmax$Region+1)))
obsavgTrdmin$Region2 <- ifelse(obsavgTrdmin$Region == 0, 1, ifelse(obsavgTrdmin$Region < 5, obsavgTrdmin$Region, ifelse(obsavgTrdmin$Region == 16, 5, obsavgTrdmin$Region+1)))
obsavgTrd$Region2 <- ifelse(obsavgTrd$Region == 0, 1, ifelse(obsavgTrd$Region < 5, obsavgTrd$Region, ifelse(obsavgTrd$Region == 16, 5, obsavgTrd$Region+1)))

c6senMinMax$Region2 <- ifelse(c6senMinMax$Region == 0, 1, ifelse(c6senMinMax$Region < 5, c6senMinMax$Region, ifelse(c6senMinMax$Region == 16, 5, c6senMinMax$Region+1)))
c6senMem$Region2 <- ifelse(c6senMem$Region == 0, 1, ifelse(c6senMem$Region < 5, c6senMem$Region, ifelse(c6senMem$Region == 16, 5, c6senMem$Region+1)))
c6senAvg$Region2 <- ifelse(c6senAvg$Region == 0, 1, ifelse(c6senAvg$Region < 5, c6senAvg$Region, ifelse(c6senAvg$Region == 16, 5, c6senAvg$Region+1)))
obssenAvgMax$Region2 <- ifelse(obssenAvgMax$Region == 0, 1, ifelse(obssenAvgMax$Region < 5, obssenAvgMax$Region, ifelse(obssenAvgMax$Region == 16, 5, obssenAvgMax$Region+1)))
obssenAvgMin$Region2 <- ifelse(obssenAvgMin$Region == 0, 1, ifelse(obssenAvgMin$Region < 5, obssenAvgMin$Region, ifelse(obssenAvgMin$Region == 16, 5, obssenAvgMin$Region+1)))
obssenAvg$Region2 <- ifelse(obssenAvg$Region == 0, 1, ifelse(obssenAvg$Region < 5, obssenAvg$Region, ifelse(obssenAvg$Region == 16, 5, obssenAvg$Region+1)))

########################################
####### Which Models are Accurate? #####
########################################
mods <- unique(c6senMem$Model)

# Re-name "Models" to account for differences
opc$Model <- with(opc, paste0(Family,'--',Member))
lrd$Model <- with(lrd, paste0(Family,'--',Member))
fad$Model <- with(fad, paste0(Family,'--',Member))

c6avgMem <- opc[opc$Model %in% mods,c(1:4)]
c6avgMem$OPCavg <- opc[opc$Model %in% mods,"Avg"]
c6avgMem$LRDlt <- lrd[lrd$Model %in% mods,"Avg"]
c6avgMem$FADgt <- fad[fad$Model %in% mods,"Avg"]

accusen <- c6senMem[seq(16,nrow(c6senMem),16),1:3]
accuavg <- c6senMem[seq(16,nrow(c6senMem),16),1:3]
for (var in c("OPCavg","LRDlt","FADgt")){
  for (reg in c(0,2:16)){
    accusen[,paste0(var,reg)] <- c6senMem[c6senMem$Region == reg,var] <= obssenAvgMax[obssenAvgMax$Region == reg,var] & 
                                 c6senMem[c6senMem$Region == reg,var] >= obssenAvgMin[obssenAvgMin$Region == reg,var]
    
    accuavg[,paste0(var,reg)] <- c6avgMem[c6avgMem$Region == reg,var] <= obsavgAvgmax[obsavgAvgmax$Region == reg,var] & 
                                 c6avgMem[c6avgMem$Region == reg,var] >= obsavgAvgmin[obsavgAvgmin$Region == reg,var]
    
  }
}

accusen$OPCTotal <- rowSums(accusen[,4:19])
accusen$LRDTotal <- rowSums(accusen[,20:35])
accusen$FADTotal <- rowSums(accusen[,36:51])
accusen$Total <- rowSums(accusen[,52:54])

accuavg$OPCTotal <- rowSums(accuavg[,4:19])
accuavg$LRDTotal <- rowSums(accuavg[,20:35])
accuavg$FADTotal <- rowSums(accuavg[,36:51])
accuavg$Total <- rowSums(accuavg[,52:54])

hist(accusen$Total)
accusen[accusen$Total >= 32,c(1:3,55)]

hist(accusen$OPCTotal)
hist(accuavg$OPCTotal)
accuavg[accuavg$OPCTotal > 10,c(1:3,52:55)]$Model
accusen[accusen$OPCTotal > 15,c(1:3,52:55)]$Model
accuavg[accuavg$OPCTotal > 10 & accusen$OPCTotal > 10,c(1:3,52:55)]
accusen[accuavg$OPCTotal > 10 & accusen$OPCTotal > 10,c(1:3,52:55)]

####################################
####### T-Tests FAD v. LRD #########
####################################
obssenAvg$FADmLRD <- obssenAvg$FADgt - obssenAvg$LRDlt
c6senMem$FADmLRD <- c6senMem$FADgt - c6senMem$LRDlt
c6senAvg$FADmLRD <- c6senAvg$FADgt - c6senAvg$LRDlt
ggplot() + #geom_linerange(data=c6senMinMax,aes(x=Region2,ymin=FADmLRDmin,ymax=FADmLRDmax),alpha=0.05,size=6) +
  geom_point(data=c6senMem[c6senMem$First == 1,],aes(x=Region2,y=FADmLRD),size=2,shape=95,color='red') +
  geom_point(data=c6senAvg,aes(x=Region2,y=FADmLRD),size=1.5,shape=21,fill='red') +
  #geom_point(data=obssenAvgMax,aes(x=Region2,y=FADmLRD),size=2,shape=4) +
  #geom_point(data=obssenAvgMin,aes(x=Region2,y=FADmLRD),size=2,shape=4) +
  geom_point(data=obssenAvg,aes(x=Region2,y=FADmLRD),size=1.5,shape=21,fill='white') +
  geom_hline(aes(yintercept=0),size=0.25) +
  scale_y_continuous(limits=c(-25,25)) + ylab("Difference in Sensitivity (%/°C)") + 
  scale_x_continuous(breaks=c(1:16),labels=REGS2) + xlab('') + ggtitle(paste0("Difference b/wn Advance & Retreat Sensitivity")) +
  theme_bw() + theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(), plot.title = element_blank(),plot.margin=unit(c(0,5.5,0.5,5.5),'pt'),
                     panel.grid.major.y = element_line(linetype="dashed"), panel.grid.minor.y = element_line(linetype="dashed")) +
  theme(axis.title = element_text(size=axistitlesize), axis.text.y=element_text(size=axistextsize)) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=axistextsize))
