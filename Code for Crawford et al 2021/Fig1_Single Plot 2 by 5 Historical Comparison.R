# Load Modules
library(ggplot2)
library(stringr)
library(gridExtra)

#################################
########## PREP WORK ############
#################################
# Define Variables
ct = 80
ver = 6
aymin = 1979
aymax = 2013

tymin = 1979
tymax = 2013

tvar = 'tregion'

axistitlesize = 8
axistextsize = 7
xpos = 1

REGS2 = c('All Regions','Okhotsk','Bering','Hudson','Baffin','St. Lawrence','Labrador','Greenland',
          'Barents','Kara','Laptev','E. Siberian','Chukchi','Beaufort',
          'CAA','CAO') # 0 - 7, 8 - 13, 14 - 16

modstokeep = c("ACCESS-CM2","BCC-CSM2-MR","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-CM6-1-HR",
               "CNRM-ESM2-1","EC-Earth3","EC-Earth3-Veg","GFDL-CM4","IPSL-CM6A-LR","MIROC-ES2L","MIROC6",
               "MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM","UKESM1-0-LL")

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

c6senMem <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_",tvar,"_historical_C",ct,"_V",ver,".csv"),as.is=T)
c6senAvg <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_MultiModelMean_",tvar,"_historical_C",ct,"_V",ver,".csv"),as.is=T)
c6senMinMax <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_MultiModelMinMax_",tvar,"_historical_C",ct,"_V",ver,".csv"),as.is=T)
c6senSDS <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_ModelStdDev_",tvar,"_historical_C",ct,"_V",ver,".csv"),as.is=T)
c6senSD <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_MultiModelStdDev_",tvar,"_historical_C",ct,"_V",ver,".csv"),as.is=T)
obssenAvg <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/ObsTempSens_MultiModelMean_",tvar,"_historical_C",ct,"_V",ver,".csv"),as.is=T)
obssenSD <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/ObsTempSens_MultiModelStdDev_",tvar,"_historical_C",ct,"_V",ver,".csv"),as.is=T)
obssenErr <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/ObsTempSens_MultiModelError_",tvar,"_historical_C",ct,"_V",ver,".csv"),as.is=T)

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
  
### Labels ###
DOYs = data.frame(
  LabelLong=c("Jan 1","Jan 31","Feb 1","Feb 28","Mar 1","Mar 31","Apr 1","Apr 30","May 1","May 31","Jun 1","Jun 30","Jul 1","Jul 31","Aug 1","Aug 31","Sep 1","Sep 30","Oct 1","Oct 31","Nov 1","Nov 30","Dec 1","Dec 31","Jan 1","Jan 31","Feb 1"),
  LabelShort=c("1/1","1/31","2/1","2/28","3/1","3/31","4/1","4/30","5/1","5/31","6/1","6/30","7/1","7/31","8/1","8/31","9/1","9/30","10/1","10/31","11/1","11/30","12/1","12/31","1/1","1/31","2/1"),
  Label=c("Jan","Jan","Feb","Feb","Mar","Mar","Apr","Apr","May","May","Jun","Jun","Jul","Jul","Aug","Aug","Sep","Sep","Oct","Oct","Nov","Nov","Dec","Dec","Jan","Jan","Feb"),
  DOY=c(1,31,32,59,60,90,91,120,121,151,152,181,182,212,213,243,244,273,274,304,305,334,335,365,366,396,397))
DOYs$Label <- as.character(DOYs$Label)

labels <- data.frame(lrd=NA,fad=NA)
row <- 0
for (thi in 1:nrow(thresh)){
  row <- row + 1
  labels[row,"lrd"] <- DOYs[DOYs$DOY == thresh$LRDThresh[thi],'Label']
  labels[row,"fad"] <- DOYs[DOYs$DOY == thresh$FADThresh[thi],'Label']
}
labels <- rbind(labels[1:4,],labels[16,],labels[5:15,])

###############################
######### PLOTTING ############
###############################

# Average OPC
p1 <- ggplot() + geom_linerange(data=c6minmaxAvg,aes(x=Region2,ymin=OPCavgmin,ymax=OPCavgmax),alpha=0.05,size=6) +
  geom_point(data=opc[opc$First == 1,],aes(x=Region2,y=Avg),size=2,shape=95,color='red') +
  geom_point(data=c6avgAvg,aes(x=Region2,y=OPCavg),size=1.5,shape=21,fill='red') +
  geom_point(data=obsavgAvgmax,aes(x=Region2,y=OPCavg),size=2,shape=4) +
  geom_point(data=obsavgAvgmin,aes(x=Region2,y=OPCavg),size=2,shape=4) +
  geom_point(data=obsavgAvg,aes(x=Region2,y=OPCavg),size=1.5,shape=21,fill='white') +
  scale_y_continuous(limits=c(0,365)) + ylab("Average (days)") + 
  scale_x_continuous(breaks=c(1:16),labels=REGS2) + xlab('') + #ggtitle("Open-Water Period") + 
  theme_bw() + theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(), plot.title = element_blank(), plot.margin=unit(c(5.5,5.5,0.5,5.5),'pt'),
  panel.grid.major.y = element_line(linetype="dashed"), panel.grid.minor.y = element_line(linetype="dashed")) +  
  theme(plot.title = element_text(hjust = 0.5, size=9), axis.title = element_text(size=axistitlesize)) +
  theme(axis.text.x = element_blank(), axis.text.y=element_text(size=axistextsize)) +
  annotate("text",x=xpos,y=360,label="a",size=3,fontface='bold')

# Average LRD %
p2 <- ggplot() + geom_linerange(data=c6minmaxAvg,aes(x=Region2,ymin=LRDltmin,ymax=LRDltmax),alpha=0.05,size=6) +
  geom_point(data=lrd[lrd$First == 1,],aes(x=Region2,y=Avg),size=2,shape=95,color='red') +
  geom_point(data=c6avgAvg,aes(x=Region2,y=LRDlt),size=1.5,shape=21,fill='red') +
  geom_point(data=obsavgAvgmax,aes(x=Region2,y=LRDlt),size=2,shape=4) +
  geom_point(data=obsavgAvgmin,aes(x=Region2,y=LRDlt),size=2,shape=4) +
  geom_point(data=obsavgAvg,aes(x=Region2,y=LRDlt),size=1.5,shape=21,fill='white') +
  scale_y_continuous(limits=c(-5,100),breaks=seq(0,100,20)) + ylab("Average (%)") + 
  scale_x_continuous(breaks=c(1:16),labels=REGS2) + xlab('') + #ggtitle(paste0("% of Area with Retreat Before Given Month")) +
  theme_bw() + theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(), plot.title = element_blank(), plot.margin=unit(c(5.5,5.5,0.5,5.5),'pt'),
  panel.grid.major.y = element_line(linetype="dashed"), panel.grid.minor.y = element_line(linetype="dashed")) +
  theme(plot.title = element_text(hjust = 0.5, size=9),axis.title = element_text(size=axistitlesize)) +
  theme(axis.text.x = element_blank(), axis.text.y=element_text(size=axistextsize)) +
  annotate("text",x=c(1:16),y=rep(-5,16),label=labels$lrd,size=2) +
  annotate("text",x=xpos,y=95,label="d",size=3,fontface='bold')
 
# Average FAD %
p3 <- ggplot() + geom_linerange(data=c6minmaxAvg,aes(x=Region2,ymin=FADgtmin,ymax=FADgtmax),alpha=0.05,size=6) +
  geom_point(data=fad[fad$First == 1,],aes(x=Region2,y=Avg),size=2,shape=95,color='red') +
  geom_point(data=c6avgAvg,aes(x=Region2,y=FADgt),size=1.5,shape=21,fill='red') +
  geom_point(data=obsavgAvgmax,aes(x=Region2,y=FADgt),size=2,shape=4) +
  geom_point(data=obsavgAvgmin,aes(x=Region2,y=FADgt),size=2,shape=4) +
  geom_point(data=obsavgAvg,aes(x=Region2,y=FADgt),size=1.5,shape=21,fill='white') +
  scale_y_continuous(limits=c(-5,100),breaks=seq(0,100,20)) + ylab("Average (%)") +  
  scale_x_continuous(breaks=c(1:16),labels=REGS2) + xlab('') + #ggtitle(paste0("% of Area with Advance After Given Month")) +
  theme_bw() + theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(), plot.title = element_blank(), plot.margin=unit(c(5.5,5.5,0.5,5.5),'pt'),
  panel.grid.major.y = element_line(linetype="dashed"), panel.grid.minor.y = element_line(linetype="dashed")) +
  theme(plot.title = element_text(hjust = 0.5, size=9), axis.title = element_text(size=axistitlesize)) +
  theme(axis.text.x = element_blank(), axis.text.y=element_text(size=axistextsize)) +
  annotate("text",x=c(1:16),y=rep(-5,16),label=labels$fad,size=2) +
  annotate("text",x=xpos,y=95,label="e",size=3,fontface='bold')

# Trend OPC
p4 <- ggplot() + geom_linerange(data=c6minmaxTrd,aes(x=Region2,ymin=OPCavgmin,ymax=OPCavgmax),alpha=0.05,size=6) +
  geom_point(data=opc[opc$First == 1,],aes(x=Region2,y=Trend),size=2,shape=95,color='red') +
  geom_point(data=c6avgTrd,aes(x=Region2,y=OPCavg),size=1.5,shape=21,fill='red') +
  geom_point(data=obsavgTrdmax,aes(x=Region2,y=OPCavg),size=2,shape=4) +
  geom_point(data=obsavgTrdmin,aes(x=Region2,y=OPCavg),size=2,shape=4) +
  geom_point(data=obsavgTrd,aes(x=Region2,y=OPCavg),size=1.5,shape=21,fill='white') +
  geom_hline(aes(yintercept=0),size=0.25) +
  scale_y_continuous(limits=c(-2,4.8)) + ylab(expression("Trend (days  "*yr^-1*")")) + 
  scale_x_continuous(breaks=c(1:16),labels=REGS2) + xlab('') + ggtitle("Open-Water Period") + 
  theme_bw() + theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(), plot.title = element_blank(), plot.margin=unit(c(0,5.5,0.5,5.5),'pt'),
  panel.grid.major.y = element_line(linetype="dashed"), panel.grid.minor.y = element_line(linetype="dashed")) +
  theme(axis.title = element_text(size=axistitlesize), axis.text.y=element_text(size=axistextsize)) +
  theme(axis.text.x = element_blank()) +
  annotate("text",x=xpos,y=4.5,label="b",size=3,fontface='bold')

# Trend LRD %
p5 <- ggplot() + geom_linerange(data=c6minmaxTrd,aes(x=Region2,ymin=LRDltmin,ymax=LRDltmax),alpha=0.05,size=6) +
  geom_point(data=lrd[lrd$First == 1,],aes(x=Region2,y=Trend),size=2,shape=95,color='red') +
  geom_point(data=c6avgTrd,aes(x=Region2,y=LRDlt),size=1.5,shape=21,fill='red') +
  geom_point(data=obsavgTrdmax,aes(x=Region2,y=LRDlt),size=2,shape=4) +
  geom_point(data=obsavgTrdmin,aes(x=Region2,y=LRDlt),size=2,shape=4) +
  geom_point(data=obsavgTrd,aes(x=Region2,y=LRDlt),size=1.5,shape=21,fill='white') +
  geom_hline(aes(yintercept=0),size=0.25) +
  scale_y_continuous(limits=c(-1.6,3.65)) + ylab(expression("Trend (%  "*yr^-1*")")) + 
  scale_x_continuous(breaks=c(1:16),labels=REGS2) + xlab('') + ggtitle(paste0("% of Area with Retreat Before Given Month")) +
  theme_bw() + theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(), plot.title = element_blank(),plot.margin=unit(c(0,5.5,0.5,5.5),'pt'),
  panel.grid.major.y = element_line(linetype="dashed"), panel.grid.minor.y = element_line(linetype="dashed")) +
  theme(axis.title = element_text(size=axistitlesize), axis.text.y=element_text(size=axistextsize)) +
  theme(axis.text.x = element_blank()) +
  annotate("text",x=c(1:16),y=rep(-1.6,16),label=labels$lrd,size=2) +
  annotate("text",x=xpos,y=3.25,label="f",size=3,fontface='bold')

# Trend FAD %
p6 <- ggplot() + geom_linerange(data=c6minmaxTrd,aes(x=Region2,ymin=FADgtmin,ymax=FADgtmax),alpha=0.05,size=6) +
  geom_point(data=fad[fad$First == 1,],aes(x=Region2,y=Trend),size=2,shape=95,color='red') +
  geom_point(data=c6avgTrd,aes(x=Region2,y=FADgt),size=1.5,shape=21,fill='red') +
  geom_point(data=obsavgTrdmax,aes(x=Region2,y=FADgt),size=2,shape=4) +
  geom_point(data=obsavgTrdmin,aes(x=Region2,y=FADgt),size=2,shape=4) +
  geom_point(data=obsavgTrd,aes(x=Region2,y=FADgt),size=1.5,shape=21,fill='white') +
  geom_hline(aes(yintercept=0),size=0.25) +
  scale_y_continuous(limits=c(-1.6,3.65)) + ylab(expression("Trend (%  "*yr^-1*")")) + 
  scale_x_continuous(breaks=c(1:16),labels=REGS2) + xlab('') + ggtitle(paste0("% of Area with Advance After Given Month")) +
  theme_bw() + theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(), plot.title = element_blank(),plot.margin=unit(c(0,5.5,0.5,5.5),'pt'),
  panel.grid.major.y = element_line(linetype="dashed"), panel.grid.minor.y = element_line(linetype="dashed")) +
  theme(axis.title = element_text(size=axistitlesize), axis.text.y=element_text(size=axistextsize)) +
  theme(axis.text.x = element_blank()) +
  annotate("text",x=c(1:16),y=rep(-1.6,16),label=labels$fad,size=2) +
  annotate("text",x=xpos,y=3.25,label="g",size=3,fontface='bold')

# Sensitivity in OPC
p7 <- ggplot() + geom_linerange(data=c6senMinMax,aes(x=Region2,ymin=OPCavgmin,ymax=OPCavgmax),alpha=0.05,size=6) +
  geom_point(data=c6senMem[c6senMem$First == 1,],aes(x=Region2,y=OPCavg),size=2,shape=95,color='red') +
  geom_point(data=c6senAvg,aes(x=Region2,y=OPCavg),size=1.5,shape=21,fill='red') +
  geom_point(data=obssenAvgMax,aes(x=Region2,y=OPCavg),size=2,shape=4) +
  geom_point(data=obssenAvgMin,aes(x=Region2,y=OPCavg),size=2,shape=4) +
  geom_point(data=obssenAvg,aes(x=Region2,y=OPCavg),size=1.5,shape=21,fill='white') +
  geom_hline(aes(yintercept=0),size=0.25) +
  scale_y_continuous(limits=c(-15,55)) + ylab(expression("Sensitivity (days "*degree*C^-1*")")) +
  scale_x_continuous(breaks=c(1:16),labels=REGS2) + xlab('') + ggtitle("Open-Water Period") + 
  theme_bw() + theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(), plot.title = element_blank(),plot.margin=unit(c(0,5.5,0.5,5.5),'pt'),
  panel.grid.major.y = element_line(linetype="dashed"), panel.grid.minor.y = element_line(linetype="dashed")) +
  theme(axis.title = element_text(size=axistitlesize), axis.text.y=element_text(size=axistextsize)) +
  theme(axis.text.x = element_blank()) +
  annotate("text",x=xpos,y=51,label="c",size=3,fontface='bold')

# Sensitivity in LRD
p8 <- ggplot() + geom_linerange(data=c6senMinMax,aes(x=Region2,ymin=LRDltmin,ymax=LRDltmax),alpha=0.05,size=6) +
  geom_point(data=c6senMem[c6senMem$First == 1,],aes(x=Region2,y=LRDlt),size=2,shape=95,color='red') +
  geom_point(data=c6senAvg,aes(x=Region2,y=LRDlt),size=1.5,shape=21,fill='red') +
  geom_point(data=obssenAvgMax,aes(x=Region2,y=LRDlt),size=2,shape=4) +
  geom_point(data=obssenAvgMin,aes(x=Region2,y=LRDlt),size=2,shape=4) +
  geom_point(data=obssenAvg,aes(x=Region2,y=LRDlt),size=1.5,shape=21,fill='white') +
  geom_hline(aes(yintercept=0),size=0.25) +
  scale_y_continuous(limits=c(-9,56)) + ylab(expression("Sensitivity (% "*degree*C^-1*")")) +
  scale_x_continuous(breaks=c(1:16),labels=REGS2) + xlab('') + ggtitle(paste0("% of Area with Retreat Before Given Month")) +
  theme_bw() + theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(), plot.title = element_blank(),plot.margin=unit(c(0,5.5,0.5,5.5),'pt'),
  panel.grid.major.y = element_line(linetype="dashed"), panel.grid.minor.y = element_line(linetype="dashed")) +
  theme(axis.title = element_text(size=axistitlesize), axis.text.y=element_text(size=axistextsize)) +
  theme(axis.text.x = element_blank()) +
  annotate("text",x=c(1:16),y=rep(-9,16),label=labels$lrd,size=2) +
  annotate("text",x=xpos,y=52,label="h",size=3,fontface='bold')

# Sensitivity in FAD
p9 <- ggplot() + geom_linerange(data=c6senMinMax,aes(x=Region2,ymin=FADgtmin,ymax=FADgtmax),alpha=0.05,size=6) +
  geom_point(data=c6senMem[c6senMem$First == 1,],aes(x=Region2,y=FADgt),size=2,shape=95,color='red') +
  geom_point(data=c6senAvg,aes(x=Region2,y=FADgt),size=1.5,shape=21,fill='red') +
  geom_point(data=obssenAvgMax,aes(x=Region2,y=FADgt),size=2,shape=4) +
  geom_point(data=obssenAvgMin,aes(x=Region2,y=FADgt),size=2,shape=4) +
  geom_point(data=obssenAvg,aes(x=Region2,y=FADgt),size=1.5,shape=21,fill='white') +
  geom_hline(aes(yintercept=0),size=0.25) +
  scale_y_continuous(limits=c(-9,56)) + ylab(expression("Sensitivity (% "*degree*C^-1*")")) +
  scale_x_continuous(breaks=c(1:16),labels=REGS2) + xlab('') +
  theme_bw() + theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(), plot.title = element_blank(), plot.margin=unit(c(0,5.5,0.5,5.5),'pt'),
  panel.grid.major.y = element_line(linetype="dashed"), panel.grid.minor.y = element_line(linetype="dashed")) +
  theme(axis.title = element_text(size=axistitlesize), axis.text.y=element_text(size=axistextsize)) +
  theme(axis.text.x = element_blank()) +
  annotate("text",x=c(1:16),y=rep(-9,16),label=labels$fad,size=2) +
  annotate("text",x=xpos,y=52,label="i",size=3,fontface='bold')

ptotal <- grid.arrange(p1,ggplot()+theme_bw(),p4,p7,p2,p3,p5,p6,p8,p9,nrow=5)
ggsave(paste0("SIVars_Historical_Comparison_CT",ct,"_ByRegion_V6.png"),ptotal,'png',figpath,width=7,height=9,dpi=300)
# ggsave(paste0("SIVars_Historical_Comparison_ByRegion_V3.eps"),ptotal,device=cairo_ps,figpath,width=7,height=9)
# 
# ####### T-Tests FAD v. LRD #########
# obssenAvg$FADmLRD <- obssenAvg$FADgt - obssenAvg$LRDlt
# c6senMem$FADmLRD <- c6senMem$FADgt - c6senMem$LRDlt
# c6senAvg$FADmLRD <- c6senAvg$FADgt - c6senAvg$LRDlt
# ggplot() + #geom_linerange(data=c6senMinMax,aes(x=Region2,ymin=FADmLRDmin,ymax=FADmLRDmax),alpha=0.05,size=6) +
#   geom_point(data=c6senMem[c6senMem$First == 1,],aes(x=Region2,y=FADmLRD),size=2,shape=95,color='red') +
#   geom_point(data=c6senAvg,aes(x=Region2,y=FADmLRD),size=1.5,shape=21,fill='red') +
#   #geom_point(data=obssenAvgMax,aes(x=Region2,y=FADmLRD),size=2,shape=4) +
#   #geom_point(data=obssenAvgMin,aes(x=Region2,y=FADmLRD),size=2,shape=4) +
#   geom_point(data=obssenAvg,aes(x=Region2,y=FADmLRD),size=1.5,shape=21,fill='white') +
#   geom_hline(aes(yintercept=0),size=0.25) +
#   scale_y_continuous(limits=c(-25,25)) + ylab("Difference in Sensitivity (%/°C)") + 
#   scale_x_continuous(breaks=c(1:16),labels=REGS2) + xlab('') + ggtitle(paste0("Difference b/wn Advance & Retreat Sensitivity")) +
#   theme_bw() + theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(), plot.title = element_blank(),plot.margin=unit(c(0,5.5,0.5,5.5),'pt'),
#                      panel.grid.major.y = element_line(linetype="dashed"), panel.grid.minor.y = element_line(linetype="dashed")) +
#   theme(axis.title = element_text(size=axistitlesize), axis.text.y=element_text(size=axistextsize)) +
#   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=axistextsize))
