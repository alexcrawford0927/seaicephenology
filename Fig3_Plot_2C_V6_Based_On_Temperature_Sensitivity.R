##### Create a graph for each region that shows the open-water period at 2°C #####

# Load Modules
library(ggplot2)
library(stringr)
library(gridExtra)

#################################
########## PREP WORK ############
#################################
# Define Variables
ct = 15 # Concentration Threshold
ver = 6 # Version
x = 2 # Temperature anomaly
ssp = 585
tymin = 1850
tymax = 2099

vars = c("OPCavg","LRDlt","FADgt")
tvar = 'tglobal'
path = '/Volumes/Troilus/CMIP6'
cmip6path = paste0(path,"/RegionalStats")
figpath = paste0("/Volumes/Troilus/CMIP6/Figures/RegionalStats/C",ct,"_V",ver)

REGS2 = c('All Regions','Okhotsk','Bering','Hudson','Baffin','St. Lawrence','Labrador','Greenland',
          'Barents','Kara','Laptev','E. Siberian','Chukchi','Beaufort',
          'CAA','CAO') # 0 - 7, 8 - 13, 14 - 16

xpos = 0.75
  
DOYs = data.frame(
  LabelLong=c("Jan 1","Jan 31","Feb 1","Feb 28","Mar 1","Mar 31","Apr 1","Apr 30","May 1","May 31","Jun 1","Jun 30","Jul 1","Jul 31","Aug 1","Aug 31","Sep 1","Sep 30","Oct 1","Oct 31","Nov 1","Nov 30","Dec 1","Dec 31","Jan 1","Jan 31"),
  LabelShort=c("1/1","1/31","2/1","2/28","3/1","3/31","4/1","4/30","5/1","5/31","6/1","6/30","7/1","7/31","8/1","8/31","9/1","9/30","10/1","10/31","11/1","11/30","12/1","12/31","1/1","1/31"),
  Label=c("Jan","Jan","Feb","Feb","Mar","Mar","Apr","Apr","May","May","Jun","Jun","Jul","Jul","Aug","Aug","Sep","Sep","Oct","Oct","Nov","Nov","Dec","Dec","Jan","Jan"),
  DOY=c(1,31,32,59,60,90,91,120,121,151,152,181,182,212,213,243,244,273,274,304,305,334,335,365,366,396))
DOYs$Label <- as.character(DOYs$Label)

# Load Files
thresh = read.csv(paste0("/Volumes/Troilus/CMIP6/Regional_DOY_Thresholds25_C50_V2.csv"))
senMemhist <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_",tvar,"_historical_C",ct,"_V",ver,".csv"),as.is=T)
intMemhist <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSensIntercept_",tvar,"_historical_C",ct,"_V",ver,".csv"),as.is=T)

senMemssp <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSens_",tvar,"_ssp",ssp,"_C",ct,"_V",ver,".csv"),as.is=T)
intMemssp <- read.csv(paste0(cmip6path,"/V",ver,"/TempSens",tymin,"-",tymax,"/CMIP6TempSensIntercept_",tvar,"_ssp",ssp,"_C",ct,"_V",ver,".csv"),as.is=T)

regs = unique(senMemssp$Region)

###################################
######### Focus on X°C ############
###################################
# Labels
labels <- data.frame(lrd=NA,fad=NA)
row <- 0
for (thi in 1:nrow(thresh)){
  row <- row + 1
  labels[row,"lrd"] <- DOYs[DOYs$DOY == thresh$LRDThresh[thi],'Label']
  labels[row,"fad"] <- DOYs[DOYs$DOY == thresh$FADThresh[thi],'Label']
}
labels <- rbind(labels[1:4,],labels[16,],labels[5:15,])

# Create versions with only 1st ensemble member
senMemhist2 <- senMemhist[senMemhist$First == 1,]
intMemhist2 <- intMemhist[intMemhist$First == 1,]
senMemssp2 <- senMemssp[senMemssp$First == 1,]
intMemssp2 <- intMemssp[intMemssp$First == 1,]

# Historical: Create version of intercept with 0s instead of negative values
intMemhist0 <- intMemhist2
intMemhist0[intMemhist0$OPCavg < 0,]$OPCavg <- 0
intMemhist0[intMemhist0$LRDlt < 0,]$LRDlt <- 0
intMemhist0[intMemhist0$FADgt < 0,]$FADgt <- 0
intMemhist0[intMemhist0$LRDlt > 100,]$LRDlt <- 100
intMemhist0[intMemhist0$FADgt > 100,]$FADgt <- 100
intMemhist0[intMemhist0$OPCavg > 365,]$OPCavg <- 365
intMemhist0$Anomaly <- "0°C"

# SSP: Create version of intercept with 0s instead of negative values
modfams <- unique(senMemssp2$Family)
tadf <- data.frame(Anomaly=paste0(x,"°C"),Family=rep(modfams,rep(length(regs),length(modfams))), Region=regs, OPCavg=NA, LRDlt=NA, FADgt=NA)
for (mod in modfams){
  for (reg in regs){
    for (var in vars){
      tadf[tadf$Family == mod & tadf$Region == reg,var] <- x*senMemssp2[senMemssp2$Family == mod & senMemssp2$Region == reg,var] + intMemssp2[intMemssp2$Family == mod & intMemssp2$Region == reg,var]
    }}}
tadf[tadf$OPCavg < 0,]$OPCavg <- 0
tadf[tadf$LRDlt < 0,]$LRDlt <- 0
tadf[tadf$FADgt < 0,]$FADgt <- 0
tadf[tadf$LRDlt > 100,]$LRDlt <- 100
tadf[tadf$FADgt > 100,]$FADgt <- 100
tadf[tadf$OPCavg > 365,]$OPCavg <- 365

# PLOTTING
modfams2 <- intersect(unique(intMemhist0$Family),unique(tadf$Family))
plotdf <- rbind(intMemhist0[intMemhist0$Family %in% modfams2,c(-1,-3,-8)],tadf[tadf$Family %in% modfams2,])
plotdf$Region2 <- ifelse(plotdf$Region == 0, 1, ifelse(plotdf$Region < 5, plotdf$Region, ifelse(plotdf$Region == 16, 5, plotdf$Region+1)))

p1 <- ggplot() + geom_boxplot(data=plotdf,aes(y=OPCavg,x=as.factor(Region2),fill=Anomaly)) +
  scale_fill_manual(values=c("darkgray","red")) + scale_x_discrete(labels=REGS2) + 
  xlab('') + ylab('Open-Water Period (Days)') +
  scale_y_continuous(limits=c(0,365)) + theme_bw() +
  theme(panel.grid = element_line(linetype='dashed'),panel.grid.minor.x = element_blank(),
        plot.title = element_blank(), plot.margin=unit(c(5.5,5.5,0,5.5),'pt'),
        axis.text = element_text(size=7), axis.title = element_text(size=7)) +
  theme(legend.position = c(0.92,0.75), legend.text = element_text(size=6), legend.title = element_text(size=7),
        axis.text.x = element_text(angle=45,hjust=1)) +
  annotate('text',x=xpos,y=360,label="a",size=3,fontface='bold')
# ggsave(paste0(figpath,"/OPCavg_",x,"Deg_Comparison_ssp",ssp,".png"),width=width,height=height,dpi=300)

p2 <- ggplot() + geom_boxplot(data=plotdf,aes(y=LRDlt,x=as.factor(Region2),fill=Anomaly)) +
  scale_fill_manual(values=c("darkgray","red")) + scale_x_discrete(labels=REGS2) + 
  xlab('') + ylab("% Area with Retreat Before First Day of Given Month") +
  scale_y_continuous(limits=c(-4,100)) + theme_bw() +
  theme(panel.grid = element_line(linetype='dashed'),panel.grid.minor.x = element_blank(),
        plot.title = element_blank(), plot.margin=unit(c(0,5.5,0.5,5.5),'pt'),
        axis.text = element_text(size=7), axis.title = element_text(size=7)) +
  theme(legend.position = 'none', axis.text.x = element_text(angle=45,hjust=1)) +
  annotate("text",x=c(1:16),y=rep(-4,16),label=labels$lrd,size=2) +
  annotate('text',x=xpos,y=95,label="b",size=3,fontface='bold')
# ggsave(paste0(figpath,"/LRDlt_",x,"Deg_Comparison_ssp",ssp,".png"),width=width,height=height,dpi=300)

p3 <- ggplot() + geom_boxplot(data=plotdf,aes(y=FADgt,x=as.factor(Region2),fill=Anomaly)) +
  scale_fill_manual(values=c("darkgray","red")) + scale_x_discrete(labels=REGS2) + 
  xlab('') + ylab('% Area with Advance After Last Day of Given Month') +
  scale_y_continuous(limits=c(-4,100)) + theme_bw() +
  theme(panel.grid = element_line(linetype='dashed'),panel.grid.minor.x = element_blank(),
        plot.title = element_blank(), plot.margin=unit(c(0,5.5,0.5,5.5),'pt'),
        axis.text = element_text(size=7), axis.title = element_text(size=7)) +
  theme(legend.position = 'none', axis.text.x = element_text(angle=45,hjust=1)) +
  annotate("text",x=c(1:16),y=rep(-4,16),label=labels$fad,size=2) +
  annotate('text',x=xpos,y=95,label="c",size=3,fontface='bold')
# ggsave(paste0(figpath,"/FADgt_",x,"Deg_Comparison_ssp",ssp,".png"),width=width,height=height,dpi=300)

ptotal <- grid.arrange(p1,p2,p3,nrow=3)
ggsave(paste0("SIVars_",x,"Deg_Comparison_CT",ct,"_ssp",ssp,".eps"),ptotal,'eps',figpath,width=5,height=7.5,dpi=300)

###########
## STATS ##
###########
modunique <- intersect(unique(plotdf[plotdf$Anomaly == "0°C","Family"]),unique(plotdf[plotdf$Anomaly == "2°C","Family"]))
pdf <- plotdf[plotdf$Family %in% modunique,]
tadf$LRDltDiff <- plotdf[plotdf$Anomaly == "2°C","LRDlt"] - plotdf[plotdf$Anomaly == "0°C","LRDlt"]
tadf$FADgtDiff <- plotdf[plotdf$Anomaly == "2°C","FADgt"] - plotdf[plotdf$Anomaly == "0°C","FADgt"]
tadf$FADmLRDDiff <- tadf$FADgtDiff - tadf$LRDltDiff
tadf$Region2 <- ifelse(tadf$Region == 0, 1, ifelse(tadf$Region < 5, tadf$Region, ifelse(tadf$Region == 16, 5, tadf$Region+1)))

p4 <- ggplot() + geom_boxplot(data=tadf,aes(y=FADmLRDDiff,x=as.factor(Region2)),fill='gray') +
  theme_bw() + scale_x_discrete(labels=REGS2) + xlab('') + ylab('Advance % Difference - Retreat % Difference') +
  theme(panel.grid = element_line(linetype='dashed'),panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust='0.5',size=10), plot.margin=unit(c(0,5.5,0.5,5.5),'pt'),
        axis.text = element_text(size=7), axis.title = element_text(size=7)) +
  scale_y_continuous(breaks=seq(-60,80,20)) +
  theme(legend.position = 'none', axis.text.x = element_text(angle=45,hjust=1)) +
  geom_hline(aes(yintercept=0)) + ggtitle("Change in Sea Ice Advance Outpaces Change in Sea Ice Retreat\n0°C to 2°C Global Temperature Anomaly") +
  annotate("text",x=1,y=76,label="Trend in Advance is Stronger",size=3,hjust='left') +
  annotate("text",x=1,y=-34,label="Trend in Retreat is Stronger",size=3, hjust='left')

ggsave(paste0("FADmLRDDiff_",x,"Deg_Comparison_ssp",ssp,".eps"),p4,'eps',figpath,width=5,height=3,dpi=300)

tdf <- data.frame(Region=c(0,2:16),OPC0=NA,OPC2=NA,p=NA)
for (reg in tdf$Region){
  x0 <- pdf[pdf$Region == reg & pdf$Anomaly=="0°C","OPCavg"]
  x2 <- pdf[pdf$Region == reg & pdf$Anomaly=="2°C","OPCavg"]
  t <- t.test(x0,x2)
  
  tdf[tdf$Region == reg,"OPC0"] <- mean(x0)
  tdf[tdf$Region == reg,"OPC2"] <- mean(x2)
  tdf[tdf$Region == reg,"p"] <- t$p.value
}
tdf$Diff <- tdf$OPC2 - tdf$OPC0

tdf <- data.frame(Region=c(0,2:16),LRD0=NA,LRD2=NA,p=NA)
for (reg in tdf$Region){
  x0 <- pdf[pdf$Region == reg & pdf$Anomaly=="0°C","LRDlt"]
  x2 <- pdf[pdf$Region == reg & pdf$Anomaly=="2°C","LRDlt"]
  t <- t.test(x0,x2)
  
  tdf[tdf$Region == reg,"LRD0"] <- mean(x0)
  tdf[tdf$Region == reg,"LRD2"] <- mean(x2)
  tdf[tdf$Region == reg,"p"] <- t$p.value
}
tdf$Diff <- tdf$LRD2 - tdf$LRD0

tdf <- data.frame(Region=c(0,2:16),FAD0=NA,FAD2=NA,p=NA)
for (reg in tdf$Region){
  x0 <- pdf[pdf$Region == reg & pdf$Anomaly=="0°C","FADgt"]
  x2 <- pdf[pdf$Region == reg & pdf$Anomaly=="2°C","FADgt"]
  t <- t.test(x0,x2)
  
  tdf[tdf$Region == reg,"FAD0"] <- mean(x0)
  tdf[tdf$Region == reg,"FAD2"] <- mean(x2)
  tdf[tdf$Region == reg,"p"] <- t$p.value
}
tdf$Diff <- tdf$FAD2 - tdf$FAD0
