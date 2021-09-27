#library(ggplot2)
#library(scales)
## create liniar model

int.summary = read.csv("/media/elmicro/meta/phase/180905-new-epi-screen/180909-proteinconc/Results-summary.csv")
int.summary = int.summary[2:8,]

int.summary$GFP = int.summary$GFP-100
int.summary$RFP = int.summary$RFP-101


# ggplot()+
# geom_point(data=int.summary, aes(x=GFP, y=GFP.conc), color="green", size=3)+
# geom_point(data=int.summary, aes(x=RFP, y=RFP.conc), color="red", size=3)+
# geom_smooth(method="lm", data=int.summary, aes(x=GFP, y=GFP.conc), color="green", se = F)+
# geom_smooth(method="lm", data=int.summary, aes(x=RFP, y=RFP.conc), color="red", se = F)+
# theme(axis.line = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+ 
# scale_y_continuous(trans=log10_trans(),limits = c(1,20000), breaks = c(1, 10, 100, 1000, 10000))+
# scale_x_continuous(trans=log10_trans(),limits = c(1,20000), breaks = c(1, 10, 100, 1000, 10000))








lm_fit.R = lm( log10(RFP.conc) ~log10(RFP), data= int.summary)
summary(lm_fit.R)


lm_fit.G = lm(log10(GFP.conc)~log10(GFP ), data= int.summary)
summary(lm_fit.G)

#get data of normalizing strain to predict concentrations

datnorm = read.csv("/media/elmicro/meta/phase/180905-new-epi-screen/180909-normstrain/datnorm.csv")

datnorm.l = subset(datnorm, id=="ctrl-log")
datnorm.s = subset(datnorm, id=="ctrl-sat")

datnorm.l $predicted.GFP = 10^(predict(lm_fit.G, datnorm.l))
datnorm.l $predicted.RFP = 10^(predict(lm_fit.R, datnorm.l))

datnorm.s $predicted.GFP = 10^(predict(lm_fit.G, datnorm.s))
datnorm.s $predicted.RFP = 10^(predict(lm_fit.R, datnorm.s))


##################################


  
  
  



#add.concentration = function(data, design){
  
  #for(K in 1:length(data)){
    
    #for(L in 1:length(data[[K]])){
      
      #for( each.ch in design$CHANELS){
        
        #col.of.interest = c("GFP_int_b3","RFP_int_b3")
        
      #  data[[K]][[L]]$GFP.conc = pedict(design.cl$GFP.lm , data[[K]][[L]][,"GFP_int_b3"]) *  design.cl$GFP.ref.ratio
       # below.zero = which(data[[K]][[L]][,each.col] < 0)
      #  if(length(below.zero)>0){
      #    print("Warning, background value may be set too high as values are below it - these are set to 0");
      #    data[[K]][[L]][below.zero,each.col]=0
    #    }
        
      # data[[K]][[L]]$RFP.conc = pedict(design.cl$RFP.lm , data[[K]][[L]][,"RFP_int_b3"]) *  design.cl$RFP.ref.ratio
     #   below.zero = which(data[[K]][[L]][,each.col] < 0)
    #    if(length(below.zero)>0){
     #     print("Warning, background value may be set too high as values are below it - these are set to 0");
    #      data[[K]][[L]][below.zero,each.col]=0
    #    }
        
    #  }
#    }
#  }
#}
#return(data)
#}


