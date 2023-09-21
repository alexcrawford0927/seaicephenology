# This code is used to plot sea ice parameters versus time for all experiments.

# Load Modules
library(ggplot2)
library(stringr)

#################################
########## PREP WORK ############
#################################
# Define Variables
ct = 80
ver = 6
vars = c("OPCavg","LRDlt","FADgt")
path = '/Volumes/Troilus/CMIP6'
cmip6path = paste0(path,"/RegionalStats")
figpath = paste0("/Volumes/Troilus/CMIP6/Figures/RegionalStats/V",ver,"_C",ct,"_SeaIceByYear_1950-2099")

tymin = 2015
tymax = 2099

regs = c(0,2:16)
REGS = c('All Regions','','Sea of Okhotsk','Bering Sea','Hudson Bay','Gulf of St. Lawrence','Labrador Sea','Greenland Sea',
         'Barents Sea','Kara Sea','Laptev Sea','East Siberian Sea','Chukchi Sea','Beaufort Sea',
         'Canadian Arctic Archipelago','Central Arctic Ocean','Baffin Bay') # 0 - 7, 8 - 13, 14 - 16
REGS2 = c('All Regions','','Okhotsk','Bering','Hudson','St. Lawrence','Labrador','Greenland',
          'Barents','Kara','Laptev','E. Siberian','Chukchi','Beaufort',
          'CAA','CAO','Baffin') # 0 - 7, 8 - 13, 14 - 16

DOYs = data.frame(
  LabelLong=c("Jan 1","Jan 31","Feb 1","Feb 28","Mar 1","Mar 31","Apr 1","Apr 30","May 1","May 31","Jun 1","Jun 30","Jul 1","Jul 31","Aug 1","Aug 31","Sep 1","Sep 30","Oct 1","Oct 31","Nov 1","Nov 30","Dec 1","Dec 31","Jan 1","Jan 31"),
  LabelShort=c("1/1","1/31","2/1","2/28","3/1","3/31","4/1","4/30","5/1","5/31","6/1","6/30","7/1","7/31","8/1","8/31","9/1","9/30","10/1","10/31","11/1","11/30","12/1","12/31","1/1","1/31"),
  Label=c("Jan","Jan","Feb","Feb","Mar","Mar","Apr","Apr","May","May","Jun","Jun","Jul","Jul","Aug","Aug","Sep","Sep","Oct","Oct","Nov","Nov","Dec","Dec","Jan","Jan"),
  DOY=c(1,31,32,59,60,90,91,120,121,151,152,181,182,212,213,243,244,273,274,304,305,334,335,365,366,396))
DOYs$Label <- as.character(DOYs$Label)

# Load Files
thresh = read.csv(paste0(path,"/Regional_DOY_Thresholds25_C50_V2.csv"))
siavgobs <- read.csv(paste0(cmip6path,"/V",ver,"/ByYear/Obs_MultiModelMean_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
siavghist <- read.csv(paste0(cmip6path,"/V",ver,"/ByYear/CMIP6_MultiModelMean_Regionalized_Historical_C",ct,"_V",ver,".csv"),as.is=T)
siavg585 <- read.csv(paste0(cmip6path,"/V",ver,"/ByYear/CMIP6_MultiModelMean_Regionalized_ssp585_C",ct,"_V",ver,".csv"),as.is=T)
siavg245 <- read.csv(paste0(cmip6path,"/V",ver,"/ByYear/CMIP6_MultiModelMean_Regionalized_ssp245_C",ct,"_V",ver,".csv"),as.is=T)
siavg126 <- read.csv(paste0(cmip6path,"/V",ver,"/ByYear/CMIP6_MultiModelMean_Regionalized_ssp126_C",ct,"_V",ver,".csv"),as.is=T)

siminmaxhist <- read.csv(paste0(cmip6path,"/V",ver,"/ByYear/CMIP6_MultiModelMinMax_Regionalized_historical_C",ct,"_V",ver,".csv"),as.is=T)
siminmax585 <- read.csv(paste0(cmip6path,"/V",ver,"/ByYear/CMIP6_MultiModelMinMax_Regionalized_ssp585_C",ct,"_V",ver,".csv"),as.is=T)
siminmax245 <- read.csv(paste0(cmip6path,"/V",ver,"/ByYear/CMIP6_MultiModelMinMax_Regionalized_ssp245_C",ct,"_V",ver,".csv"),as.is=T)
siminmax126 <- read.csv(paste0(cmip6path,"/V",ver,"/ByYear/CMIP6_MultiModelMinMax_Regionalized_ssp126_C",ct,"_V",ver,".csv"),as.is=T)

siavgobs$Reg2 <- ifelse(siavgobs$Region %in% c(2,3,4,16), 1, ifelse(siavgobs$Region %in% c(5,6,7,8), 2, ifelse(siavgobs$Region %in% c(9,10,11,12), 3,4)))
siavgobs$Reg3 <- ifelse(siavgobs$Region %in% c(2,5,9,13), 1, ifelse(siavgobs$Region %in% c(3,6,10,14), 2, ifelse(siavgobs$Region %in% c(4,7,11,15), 3,4)))
siavghist$Reg2 <- ifelse(siavghist$Region %in% c(2,3,4,16), 1, ifelse(siavghist$Region %in% c(5,6,7,8), 2, ifelse(siavghist$Region %in% c(10,11,12,13), 3,4)))
siavghist$Reg3 <- ifelse(siavghist$Region %in% c(2,5,9,13), 1, ifelse(siavghist$Region %in% c(3,6,10,14), 2, ifelse(siavghist$Region %in% c(4,7,11,15), 3,4)))
siavg585$Reg2 <- ifelse(siavg585$Region %in% c(2,3,4,16), 1, ifelse(siavg585$Region %in% c(5,6,7,8), 2, ifelse(siavg585$Region %in% c(10,11,12,13), 3,4)))
siavg585$Reg3 <- ifelse(siavg585$Region %in% c(2,5,9,13), 1, ifelse(siavg585$Region %in% c(3,6,10,14), 2, ifelse(siavg585$Region %in% c(4,7,11,15), 3,4)))
siavg245$Reg2 <- ifelse(siavg245$Region %in% c(2,3,4,16), 1, ifelse(siavg245$Region %in% c(5,6,7,8), 2, ifelse(siavg245$Region %in% c(10,11,12,13), 3,4)))
siavg245$Reg3 <- ifelse(siavg245$Region %in% c(2,5,9,13), 1, ifelse(siavg245$Region %in% c(3,6,10,14), 2, ifelse(siavg245$Region %in% c(4,7,11,15), 3,4)))
siavg126$Reg2 <- ifelse(siavg126$Region %in% c(2,3,4,16), 1, ifelse(siavg126$Region %in% c(5,6,7,8), 2, ifelse(siavg126$Region %in% c(10,11,12,13), 3,4)))
siavg126$Reg3 <- ifelse(siavg126$Region %in% c(2,5,9,13), 1, ifelse(siavg126$Region %in% c(3,6,10,14), 2, ifelse(siavg126$Region %in% c(4,7,11,15), 3,4)))

siminmaxhist$Reg2 <- ifelse(siminmaxhist$Region %in% c(2,3,4,16), 1, ifelse(siminmaxhist$Region %in% c(5,6,7,8), 2, ifelse(siminmaxhist$Region %in% c(10,11,12,13), 3,4)))
siminmaxhist$Reg3 <- ifelse(siminmaxhist$Region %in% c(2,5,9,13), 1, ifelse(siminmaxhist$Region %in% c(3,6,10,14), 2, ifelse(siminmaxhist$Region %in% c(4,7,11,15), 3,4)))
siminmax585$Reg2 <- ifelse(siminmax585$Region %in% c(2,3,4,16), 1, ifelse(siminmax585$Region %in% c(5,6,7,8), 2, ifelse(siminmax585$Region %in% c(10,11,12,13), 3,4)))
siminmax585$Reg3 <- ifelse(siminmax585$Region %in% c(2,5,9,13), 1, ifelse(siminmax585$Region %in% c(3,6,10,14), 2, ifelse(siminmax585$Region %in% c(4,7,11,15), 3,4)))
siminmax245$Reg2 <- ifelse(siminmax245$Region %in% c(2,3,4,16), 1, ifelse(siminmax245$Region %in% c(5,6,7,8), 2, ifelse(siminmax245$Region %in% c(10,11,12,13), 3,4)))
siminmax245$Reg3 <- ifelse(siminmax245$Region %in% c(2,5,9,13), 1, ifelse(siminmax245$Region %in% c(3,6,10,14), 2, ifelse(siminmax245$Region %in% c(4,7,11,15), 3,4)))
siminmax126$Reg2 <- ifelse(siminmax126$Region %in% c(2,3,4,16), 1, ifelse(siminmax126$Region %in% c(5,6,7,8), 2, ifelse(siminmax126$Region %in% c(10,11,12,13), 3,4)))
siminmax126$Reg3 <- ifelse(siminmax126$Region %in% c(2,5,9,13), 1, ifelse(siminmax126$Region %in% c(3,6,10,14), 2, ifelse(siminmax126$Region %in% c(4,7,11,15), 3,4)))

##############################################
########## PLOTTING YEAR-BY-YEAR #############
##############################################
######### OPCavg ###########
# Labels
labels <- data.frame(x=rep(1950,16),y=rep(96,16),lrd=NA,fad=NA)
for (thi in 1:nrow(thresh)){
  labels[thi,"lrd"] <- DOYs[DOYs$DOY == thresh$LRDThresh[thi],'Label']
  labels[thi,"fad"] <- DOYs[DOYs$DOY == thresh$FADThresh[thi],'Label']
}
labels <- rbind(labels[1:4,],labels[16,],labels[5:15,])
labels$Reg2=rep(1:4,rep(4,4))
labels$Reg3=rep(1:4,4)

annot = data.frame(RegLabel = REGS2[c(3:5,17,6:16,1)],Reg2=rep(1:4,rep(4,4)),Reg3=rep(1:4,4),
                   x=rep(1950,16),y=c(rep(15,8),rep(357,8)),hjust=rep(0,16),vjust=rep(1,16))
annot$REGLabelLetter <- c("a)","b)","c)","d)","e)","f)","g)","h)","i)","j)","k)","l)","m)","n)","o)","p)")
annot$REGLabel2 <- annot$REGLabelLetter
for (row in 1:nrow(annot)){
  annot[row,"REGLabel2"] <- paste(annot[row,'REGLabelLetter'],annot[row,'RegLabel'])
}

# legend = data.frame(Reg2=rep(1:4,rep(4,4)),Reg3=rep(1:4,4),
#                     xobs = 2050, yobs = 120, lobs = c(rep('',15),'x Observations'),
#                     xhist = 2050, yhist = 95, lhist = c(rep('',15),'- Historical'),
#                     x126 = 2050, y126 = 70, l126 = c(rep('',15),'- SSP126'),
#                     x245 = 2050, y245 = 45, l245 = c(rep('',15),'- SSP245'),
#                     x585 = 2050, y585 = 20, l585 = c(rep('',15),'- SSP585'))

ggplot() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
panel.grid = element_line(linetype='dashed',size=0.35),strip.text = element_blank(), strip.background = element_blank()) +
geom_ribbon(data=siminmaxhist,aes(ymax=OPCavgmax,ymin=OPCavgmin,x=Year),alpha=0.3,fill='black') +
geom_ribbon(data=siminmax585,aes(ymax=OPCavgmax,ymin=OPCavgmin,x=Year),alpha=0.3,fill='red') +
geom_ribbon(data=siminmax245,aes(ymax=OPCavgmax,ymin=OPCavgmin,x=Year),alpha=0.3,fill='orange') +
geom_ribbon(data=siminmax126,aes(ymax=OPCavgmax,ymin=OPCavgmin,x=Year),alpha=0.3,fill='blue') +

geom_point(data=siavgobs,aes(y=OPCavg,x=Year),size=0.5,shape=4)  + 
geom_line(data=siavghist,aes(y=OPCavg,x=Year),color='#525252',lwd=1) + 
geom_line(data=siavg585,aes(y=OPCavg,x=Year),lwd=1,color='red') + 
geom_line(data=siavg245,aes(y=OPCavg,x=Year),lwd=1,color='orange') + 
geom_line(data=siavg126,aes(y=OPCavg,x=Year),lwd=1,color='blue') + 

facet_grid(rows=vars(Reg2),cols=vars(Reg3)) +
scale_y_continuous(limits=c(0,365),breaks=seq(0,365,90)) + 
scale_x_continuous(breaks=seq(1950,2100,50)) + xlab('') + ylab('Open-Water Period (Days)') +
geom_text(data=annot,aes(x=x,y=y,label=REGLabel2),hjust=0,fontface='bold') +
# geom_text(data=legend,aes(x=xobs,y=yobs,label=lobs),hjust=0,color='black',size=2) +
# geom_text(data=legend,aes(x=xhist,y=yhist,label=lhist),hjust=0,color='#525252',size=2) +
# geom_text(data=legend,aes(x=x126,y=y126,label=l126),hjust=0,color='blue',size=2) +
# geom_text(data=legend,aes(x=x245,y=y245,label=l245),hjust=0,color='orange',size=2) +
# geom_text(data=legend,aes(x=x585,y=y585,label=l585),hjust=0,color='red',size=2)

ggsave(filename = paste0(path,"/Figures/RegionalStats/V",ver,"_SeaIceByYear_C",ct,"_1950-2099/OPCavg_1950-2099_SeaIceByYear_AllRegs_gridlines.png"),width=7.5,height=7.5,dpi=300)

#### NO Gridlines
ggplot() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
                              panel.grid = element_blank(),strip.text = element_blank(), strip.background = element_blank()) +
  geom_ribbon(data=siminmaxhist,aes(ymax=OPCavgmax,ymin=OPCavgmin,x=Year),alpha=0.3,fill='black') +
  geom_ribbon(data=siminmax585,aes(ymax=OPCavgmax,ymin=OPCavgmin,x=Year),alpha=0.3,fill='red') +
  geom_ribbon(data=siminmax245,aes(ymax=OPCavgmax,ymin=OPCavgmin,x=Year),alpha=0.3,fill='orange') +
  geom_ribbon(data=siminmax126,aes(ymax=OPCavgmax,ymin=OPCavgmin,x=Year),alpha=0.3,fill='blue') +
  
  geom_point(data=siavgobs,aes(y=OPCavg,x=Year),size=0.5,shape=4)  + 
  geom_line(data=siavghist,aes(y=OPCavg,x=Year),color='#525252',lwd=1) + 
  geom_line(data=siavg585,aes(y=OPCavg,x=Year),lwd=1,color='red') + 
  geom_line(data=siavg245,aes(y=OPCavg,x=Year),lwd=1,color='orange') + 
  geom_line(data=siavg126,aes(y=OPCavg,x=Year),lwd=1,color='blue') + 
  
  facet_grid(rows=vars(Reg2),cols=vars(Reg3)) +
  scale_y_continuous(limits=c(0,365),breaks=seq(0,365,60)) + 
  scale_x_continuous(breaks=seq(1950,2100,50)) + xlab('') + ylab('Open-Water Period (Days)') +
  geom_text(data=annot,aes(x=x,y=y,label=REGLabel2),hjust=0)

ggsave(filename = paste0(path,"/Figures/RegionalStats/V",ver,"_SeaIceByYear_C",ct,"_1950-2099/OPCavg_1950-2099_SeaIceByYear_AllRegs.png"),width=7.5,height=7.5,dpi=300)

###### LRDlt ########
annot = data.frame(RegLabel = REGS2[c(3:5,17,6:16,1)],Reg2=rep(1:4,rep(4,4)),Reg3=rep(1:4,4),
                   x=rep(2100,16),y=c(rep(5,8),rep(5,8)),hjust=rep(0,16),vjust=rep(1,16))
annot$REGLabelLetter <- c("a)","b)","c)","d)","e)","f)","g)","h)","i)","j)","k)","l)","m)","n)","o)","p)")
annot$REGLabel2 <- annot$REGLabelLetter
for (row in 1:nrow(annot)){
  annot[row,"REGLabel2"] <- paste(annot[row,'REGLabelLetter'],annot[row,'RegLabel'])
}

ggplot() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        panel.grid = element_line(linetype='dashed',size=0.35),strip.text = element_blank(), strip.background = element_blank()) +
  geom_ribbon(data=siminmaxhist,aes(ymax=LRDltmax,ymin=LRDltmin,x=Year),alpha=0.3,fill='black') +
  geom_ribbon(data=siminmax585,aes(ymax=LRDltmax,ymin=LRDltmin,x=Year),alpha=0.3,fill='red') +
  geom_ribbon(data=siminmax245,aes(ymax=LRDltmax,ymin=LRDltmin,x=Year),alpha=0.3,fill='orange') +
  geom_ribbon(data=siminmax126,aes(ymax=LRDltmax,ymin=LRDltmin,x=Year),alpha=0.3,fill='blue') +
  
  geom_point(data=siavgobs,aes(y=LRDlt,x=Year),size=0.5,shape=4)  + 
  geom_line(data=siavghist,aes(y=LRDlt,x=Year),color='#525252',lwd=1) + 
  geom_line(data=siavg585,aes(y=LRDlt,x=Year),lwd=1,color='red') + 
  geom_line(data=siavg245,aes(y=LRDlt,x=Year),lwd=1,color='orange') + 
  geom_line(data=siavg126,aes(y=LRDlt,x=Year),lwd=1,color='blue') + 
  
  facet_grid(rows=vars(Reg2),cols=vars(Reg3)) +
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25)) + 
  scale_x_continuous(breaks=seq(1950,2100,50)) + xlab('') + ylab('% Open Before First of Given Month') +
  geom_text(data=annot,aes(x=x,y=y,label=REGLabel2),hjust=1,vjust=0,fontface='bold',size=3.5) +
  geom_text(data=labels,aes(x=x,y=y,label=lrd),hjust=0,vjust=0,size=3)

ggsave(filename = paste0(path,"/Figures/RegionalStats/V",ver,"_SeaIceByYear_C",ct,"_1950-2099/LRDlt_1950-2099_SeaIceByYear_AllRegs_gridlines.png"),width=7.5,height=7.5,dpi=300)

# No gridlines
ggplot() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
          panel.grid = element_blank(),strip.text = element_blank(), strip.background = element_blank()) +
  geom_ribbon(data=siminmaxhist,aes(ymax=LRDltmax,ymin=LRDltmin,x=Year),alpha=0.3,fill='black') +
  geom_ribbon(data=siminmax585,aes(ymax=LRDltmax,ymin=LRDltmin,x=Year),alpha=0.3,fill='red') +
  geom_ribbon(data=siminmax245,aes(ymax=LRDltmax,ymin=LRDltmin,x=Year),alpha=0.3,fill='orange') +
  geom_ribbon(data=siminmax126,aes(ymax=LRDltmax,ymin=LRDltmin,x=Year),alpha=0.3,fill='blue') +
  
  geom_point(data=siavgobs,aes(y=LRDlt,x=Year),size=0.5,shape=4)  + 
  geom_line(data=siavghist,aes(y=LRDlt,x=Year),color='#525252',lwd=1) + 
  geom_line(data=siavg585,aes(y=LRDlt,x=Year),lwd=1,color='red') + 
  geom_line(data=siavg245,aes(y=LRDlt,x=Year),lwd=1,color='orange') + 
  geom_line(data=siavg126,aes(y=LRDlt,x=Year),lwd=1,color='blue') + 
  
  facet_grid(rows=vars(Reg2),cols=vars(Reg3)) +
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25)) + 
  scale_x_continuous(breaks=seq(1950,2100,50)) + xlab('') + ylab('% Open Before First of Given Month') +
  geom_text(data=annot,aes(x=x,y=y,label=REGLabel2),hjust=1,vjust=0,fontface='bold',size=3.5) +
  geom_text(data=labels,aes(x=x,y=y,label=lrd),hjust=0,vjust=0,size=3)

ggsave(filename = paste0(path,"/Figures/RegionalStats/V",ver,"_SeaIceByYear_C",ct,"_1950-2099/LRDlt_1950-2099_SeaIceByYear_AllRegs.png"),width=7.5,height=7.5,dpi=300)

###### FADgt ########
ggplot() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
                              panel.grid = element_line(linetype='dashed',size=0.35),strip.text = element_blank(), strip.background = element_blank()) +
  geom_ribbon(data=siminmaxhist,aes(ymax=FADgtmax,ymin=FADgtmin,x=Year),alpha=0.3,fill='black') +
  geom_ribbon(data=siminmax585,aes(ymax=FADgtmax,ymin=FADgtmin,x=Year),alpha=0.3,fill='red') +
  geom_ribbon(data=siminmax245,aes(ymax=FADgtmax,ymin=FADgtmin,x=Year),alpha=0.3,fill='orange') +
  geom_ribbon(data=siminmax126,aes(ymax=FADgtmax,ymin=FADgtmin,x=Year),alpha=0.3,fill='blue') +
  
  geom_point(data=siavgobs,aes(y=FADgt,x=Year),size=0.5,shape=4)  + 
  geom_line(data=siavghist,aes(y=FADgt,x=Year),color='#525252',lwd=1) + 
  geom_line(data=siavg585,aes(y=FADgt,x=Year),lwd=1,color='red') + 
  geom_line(data=siavg245,aes(y=FADgt,x=Year),lwd=1,color='orange') + 
  geom_line(data=siavg126,aes(y=FADgt,x=Year),lwd=1,color='blue') + 
  
  facet_grid(rows=vars(Reg2),cols=vars(Reg3)) +
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25)) + 
  scale_x_continuous(breaks=seq(1950,2100,50)) + xlab('') + ylab('% Open After End of Given Month') +
  geom_text(data=annot,aes(x=x,y=y,label=REGLabel2),hjust=1,vjust=0,fontface='bold',size=3.5) +
  geom_text(data=labels,aes(x=x,y=y,label=fad),hjust=0,vjust=0,size=3)

ggsave(filename = paste0(path,"/Figures/RegionalStats/V",ver,"_SeaIceByYear_C",ct,"_1950-2099/FADgt_1950-2099_SeaIceByYear_AllRegs_gridlines.png"),width=7.5,height=7.5,dpi=300)

# No Gridlines
ggplot() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
                              panel.grid = element_blank(),strip.text = element_blank(), strip.background = element_blank()) +
  geom_ribbon(data=siminmaxhist,aes(ymax=FADgtmax,ymin=FADgtmin,x=Year),alpha=0.3,fill='black') +
  geom_ribbon(data=siminmax585,aes(ymax=FADgtmax,ymin=FADgtmin,x=Year),alpha=0.3,fill='red') +
  geom_ribbon(data=siminmax245,aes(ymax=FADgtmax,ymin=FADgtmin,x=Year),alpha=0.3,fill='orange') +
  geom_ribbon(data=siminmax126,aes(ymax=FADgtmax,ymin=FADgtmin,x=Year),alpha=0.3,fill='blue') +
  
  geom_point(data=siavgobs,aes(y=FADgt,x=Year),size=0.5,shape=4)  + 
  geom_line(data=siavghist,aes(y=FADgt,x=Year),color='#525252',lwd=1) + 
  geom_line(data=siavg585,aes(y=FADgt,x=Year),lwd=1,color='red') + 
  geom_line(data=siavg245,aes(y=FADgt,x=Year),lwd=1,color='orange') + 
  geom_line(data=siavg126,aes(y=FADgt,x=Year),lwd=1,color='blue') + 
  
  facet_grid(rows=vars(Reg2),cols=vars(Reg3)) +
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25)) + 
  scale_x_continuous(breaks=seq(1950,2100,50)) + xlab('') + ylab('% Open After End of Given Month') +
  geom_text(data=annot,aes(x=x,y=y,label=REGLabel2),hjust=1,vjust=0,fontface='bold',size=3.5) +
  geom_text(data=labels,aes(x=x,y=y,label=fad),hjust=0,vjust=0,size=3)

ggsave(filename = paste0(path,"/Figures/RegionalStats/V",ver,"_SeaIceByYear_C",ct,"_1950-2099/FADgt_1950-2099_SeaIceByYear_AllRegs.png"),width=7.5,height=7.5,dpi=300)
