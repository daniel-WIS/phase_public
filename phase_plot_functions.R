
plotphase = function(...){
  ggplot ()+ geom_point(data=subset(tab, ...), aes(RFP_b5, BFP_b5), size=.7)+
    theme(axis.line = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
    scale_y_continuous(trans=log2_trans(),limits = c(1,50000), breaks= c(0,1,10,100,1000,10000))+
    scale_x_continuous(trans=log2_trans(),limits = c(1,50000), breaks= c(0,1,10,100,1000,10000))+
    geom_abline(slope=1, intercept=(-0):(0), color="grey")+
    #ggtitle(paste("Cells with GFP Intensity >",as.character(GFP_val),"AU"))+
    xlab("2mer (RFP Intensity in AU)")+
    ylab("4mer (BFP Intensity in AU)")
}



plot.overlap.gfp = function(xmer="401", kinase_variant = "K34"){
  set.seed(2)
  #min_points = 1000
  dat_nokinase = subset(tab, tetramer == xmer & kinase == kinase_variant & BFP_foci<=0 & RFP_foci<=0 & GFP_b0 < 100)
  #dat_nokinase = dat_nokinase[sample(x = rownames(dat_nokinase), size = min_points),]
  dat_kinase = subset(tab, tetramer == xmer &  kinase == kinase_variant & BFP_foci<=0 & RFP_foci<=0 & GFP_b0 > 1000)
  #dat_kinase = dat_kinase[sample(x = rownames(dat_kinase), size = min_points),]
  
  ggplot ()+ 
    geom_point(data=dat_kinase, aes(RFP_b5, BFP_b5), alpha=0.7, color="red", size=.2)+
    geom_point(data= dat_nokinase, aes(RFP_b5, BFP_b5), alpha=0.7, size=.2)+
    theme(
      axis.line = element_blank(),
      panel.grid.major = element_line(colour = "grey", size = 0.5),
      panel.grid.minor = element_line(colour = "lightgrey"),
      panel.background = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.margin = margin(r = 50, l = 50)
    )+
    scale_y_continuous(trans=log2_trans(),limits = c(1,10000), breaks= c(0,1,10,100,1000,10000))+
    scale_x_continuous(trans=log2_trans(),limits = c(1,10000), breaks= c(0,1,10,100,1000,10000))+
    geom_abline(slope=1, intercept=(-0):(0), color="grey",alpha=0.7)+
    xlab("2mer (RFP_b5 Intensity in AU)")+
    ylab("4mer (BFP_b5 Intensity in AU)")+
    ggtitle(paste(xmer,kinase_variant,"n_nokinase:",nrow(dat_nokinase),"n_kinase:",nrow(dat_kinase)))
  
}


plotphase.overlap = function(data=tabdan, bottom_layer, top_layer){
  ggplot ()+ geom_point(data=subset(data, xmer == bottom_layer &  BFP_foci!=0 & RFP_foci!=0), aes(RFP_b5, BFP_b5), alpha=0.5, size=0.7)+
    geom_point(data=subset(data, xmer == top_layer &  BFP_foci!=0 & RFP_foci!=0), aes(RFP_b5, BFP_b5),  alpha=0.5, color="red", size=0.7)+
    theme(
      axis.line = element_blank(),
      panel.grid.major = element_line(colour = "grey", size = 0.5),
      panel.grid.minor = element_line(colour = "lightgrey"),
      panel.background = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.margin = margin(r = 50, l = 50)
    )+
    scale_y_continuous(trans=log2_trans(),limits = c(1,50000), breaks= c(0,1,10,100,1000,10000))+
    scale_x_continuous(trans=log2_trans(),limits = c(1,50000), breaks= c(0,1,10,100,1000,10000))+
    geom_abline(slope=1, intercept=(-0):(0), color="grey")+
    xlab("2mer (RFP_b5 Intensity in AU)")+
    ylab("4mer (BFP_b5 Intensity in AU)")+
    ggtitle("Phosphomimetic")
}
