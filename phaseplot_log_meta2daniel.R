if (Sys.info()['sysname'] == "Darwin") {
  os.foldername <- "/Volumes"
  source("/Users/daniel/phase/Multiwell.R")
  source("/Users/daniel/phase/general.R")
  source("/Users/daniel/phase/MicroscopeToolBox.R")
  
}else if (Sys.info()['sysname'] == "Linux") {
  os.foldername <- "/media"
  
  source("/data/elevy/70_R_Data/bin/RToolBox/Multiwell.R")
  source("/data/elevy/70_R_Data/bin/RToolBox/general.R")
  source("/data/elevy/70_R_Data/bin/RToolBox/MicroscopeToolBox.R")
}

#### FOLDER NAME + LOCATION
imagefoldername <- "200804_valvar_p1_r5"
pdefname <- "pdef_valvar_phase_only.csv"
folderpath <- paste(os.foldername,"elmicro/danielo",imagefoldername,sep = "/")


design.dan = microscope.get.design(
  F=c(paste(folderpath,"images",sep = "/")),
  D=c("comp"),
  PDEF=c(paste(folderpath,pdefname,sep = "/")),
  FORMAT=384,
  OUT=c("_output"),
  CHANELS=c("GFP","RFP"),  
  MEASURE.COL = "_int_b5",
  DIR.res = paste(folderpath,"results",sep = "/")
)

data.raw = microscope.load.data(design.dan)
data.raw = uscope.process.add.ncells(data = data.raw)
lsos() ## check how much memory the object takes

design.dan = uscope.process.estimate.background(data.raw, design.dan)
data.1    = uscope.process.reorder(data.raw, design=design.dan)
data.1    = uscope.process.remove.first.pic(data.1)
data.1    = uscope.process.remove.background(data.1, design.dan)

uscope.count.cells(data.1)
data.1    = uscope.process.remove.small(data.1, MIN.size=950,MAX.size=2500)
data.1   = uscope.process.BF(data.1)
data.12dan   = uscope.process.remove.BF.outliers(data.1, cutoff=0.9)

data.12dan   = uscope.process.add.ncells(data = data.12dan)


uscope.count.cells(data.12dan)

#save(list=ls(),file=paste(folderpath,"processed_data.RData",sep="/"))
##### create condensed df for log
well = design.dan$PDEF$well

datad = vector("list")
for (w in well) {
  
  datad[[w]] = data.frame (data.12dan$comp[[w]]$cell,
                           data.12dan$comp[[w]]$area,
                              data.12dan$comp[[w]]$pic,
                              data.12dan$comp[[w]]$RFP_int_b5,
                              data.12dan$comp[[w]]$GFP_int_b5,
                              data.12dan$comp[[w]]$f1_inRFP_toRFPmed,
                              data.12dan$comp[[w]]$f1_inGFP_toGFPmed,
                              data.12dan$comp[[w]]$f1_inGFParea,
                           data.12dan$comp[[w]]$inGFPnfoci,
                           data.12dan$comp[[w]]$x,
                           data.12dan$comp[[w]]$y,
                        rep (design.dan$PDEF$Colony.No.[design.dan$PDEF$well==w], times = length(data.12dan$comp[[w]]$cell)),
                        rep (design.dan$PDEF$GFP[design.dan$PDEF$well==w], times = length(data.12dan$comp[[w]]$cell)),
                        rep (design.dan$PDEF$RFP[design.dan$PDEF$well==w], times = length(data.12dan$comp[[w]]$cell)),
                        rep (design.dan$PDEF$well[design.dan$PDEF$well==w], times = length(data.12dan$comp[[w]]$cell))
                        )
  
colnames(datad[[w]]) = c("cell", "cellsize","pic","RFP","GFP","RFP_foci", "GFP_foci", "GFP_foci_size","nfoci","x","y", "rep","xmer","dimer", "well")
  }

tabdan = do.call(rbind, datad)


#save(list=ls(),file=paste(folderpath,"data_w_tabdan.RData",sep="/"))

#### NORMALIZATION
library(scales)
## create liniar model

int.summary = read.csv("/media/elmicro/meta/phase/180905-new-epi-screen/180909-proteinconc/Results-summary.csv")
int.summary = read.csv("/media/elmicro/meta/phase/180905-new-epi-screen/180909-proteinconc/Results-summary.csv")
int.summary = int.summary[2:8,]

int.summary$GFP = int.summary$GFP-100
int.summary$RFP = int.summary$RFP-101






tabdan$GFPraw = tabdan$GFP
tabdan$RFPraw = tabdan$RFP


## The Chromaslide gave a mean intensity of 49930.831 when I did the regression for the protein concentration. 
## Now it gives 20604. 

tabdan$GFP = tabdan$GFP* 2.42


## Chromaslide then: 5900.757 | Chromaslide now: 2875 | 

tabdan$RFP = tabdan$RFP *2.052437

tabdan$GFPconc = 10^(predict(lm_fit.G, tabdan))
tabdan$RFPconc = 10^(predict(lm_fit.R, tabdan))

save(list=ls(),file=paste(folderpath,"data_w_norm.RData",sep="/"))

##### 

library(scales) 
library(ggplot2)

  plotphase.dan= function (g) {
    ggplot ()+ geom_point(data=subset(tabdan, xmer == g  & GFP_foci==0& RFP_foci==0), aes(RFPconc, GFPconc), alpha=0.5)+
      geom_point(data=subset(tabdan, xmer == g  & GFP_foci!=0& RFP_foci!=0), aes(RFPconc, GFPconc), alpha=0.5,color="red")+
       theme(axis.line = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
      ggtitle( paste(as.character(g)))+
      scale_y_continuous(trans=log2_trans(),limits = c(1,50000), breaks= c(0,1,10,100,1000,10000))+
      scale_x_continuous(trans=log2_trans(),limits = c(1,50000), breaks= c(0,1,10,100,1000,10000))+
     geom_abline(slope=1, intercept=(-0):(0), color="red")
    
  }
  
  
  plotphase.dan("v401")

  xmer = unique(design.dan$PDEF$GFP) 
  scatterlistd= vector("list")
  for (x in xmer) {
    scatterlistd[[x]] = plotphase.dan(x) 
  }
  
  pdf(file= paste(folderpath,"Rplots/phase_w_foci_new_and_calibrated.pdf",sep="/"),useDingbats = F,width = 10,height = 10)
  
  
  par(mfrow=c(3,3))
  for (multimer in unique(tabdan$xmer)){
    
    data2plot = subset(tabdan, xmer == multimer)
    plot(log(data2plot$RFP),log(data2plot$GFP), xlab=data2plot$dimer[1],ylab=data2plot$xmer[1], pch=16, xlim=c(0,11),ylim=c(0,11),cex=0.6)
    abline(0,1,col="red")
    
  }
  dev.off()
  
  

  pdf(file= paste(folderpath,"Rplots/phase_w_foci_new_and_calibrated.pdf",sep="/"),useDingbats = F)
  scatterlistd
  dev.off()

  

### get general info and try to plot affinity vs lowest point on diagonal
library(stringi)
library(dplyr)

dan = tabdan %>%
  select(cell, pic, cellsize, RFPconc, GFPconc, nfoci, xmer ) %>%
  group_by(xmer) %>%
  summarise(allfoci=sum(nfoci),
            ncells= length(cell), 
            cellsize = median(cellsize), 
            RFPconc = median(RFPconc), 
            GFPconc = median(GFPconc))


