# Source Scripts + Data ---------------------------------------------------

#Check which OS ("Darwin" is MacOS)
if (Sys.info()['sysname'] == "Darwin") {
  os.foldername <- "/Volumes"
  source("/Users/d/phase/Multiwell.R")
  source("/Users/d/phase/general.R")
  source("/Users/d/phase/MicroscopeToolBox.R")
  
}else if (Sys.info()['sysname'] == "Linux") {
  os.foldername <- "/media"
  source("/data/elevy/70_R_Data/bin/RToolBox/Multiwell.R")
  source("/data/elevy/70_R_Data/bin/RToolBox/general.R")
  source("/data/elevy/70_R_Data/bin/RToolBox/MicroscopeToolBox.R")
}

source("/home/d/phase_public_final/postprocess_result_txt.R")
source("/home/d/phase_public_final/valvar_diploid_plot_functions.R")

#!/usr/bin/env Rscript
library("optparse")

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied", call.=FALSE)
}


option_list = list(
  make_option(c("-wd", "--working_dir"), type="character", default=NULL, 
              help="Absolute Path of Working Directory", metavar="character"),
  make_option(c("-nd", "--nd_folder_path"), type="character", default=NULL, 
              help="Absolute Path to Folder that contains .nd file", metavar="character"),
  make_option(c("-p", "--pdef_path"), type="character", default=NULL, 
              help="Absolute Path to Plate Definition File", metavar="character"),
  make_option(c("-r", "--results_txt_path"), type="character", default=NULL, 
              help="Absolute Path to folder that contains result .txt files", metavar="character"),
  make_option(c("-w", "--wells"), type="integer", default=384, 
              help="Absolute Path to folder that contains result .txt files", metavar="number"),
  make_option(c("-c", "--column_to_measure"), type="character", default="_int_b5", 
              help="Absolute Path to folder that contains result .txt files", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default=NULL, 
              help="Absolute Path to Output Folder", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

pdef_head <- colnames(read.csv(opt$pdef_path,sep=",")) #get column names of pdef file
CHANNELS_INPUT <- unlist(lapply(pdef_head,function(x) if(x=="GFP" | x=="RFP" | x=="BFP") x))

design = microscope.get.design(
  F = opt$nd_folder_path,
  D = c("comp"),
  PDEF = opt$pdef_path,
  FORMAT = opt$wells,
  OUT = opt$results_txt_path,
  CHANELS = CHANNELS_INPUT,
  MEASURE.COL = opt$column_to_measure,
  DIR.res = opt$out_dir
)

postprocess_result_txt(design.file = design, min.cell.size = 1000, max.cell.size = 2500, brightfield.cutoff = 0.8)
#do normalization now, but ask for it with if condition


save(list=ls(),file = paste0(data.name,"_ready_to_analyze.RData"))

# Data Manipulation + Cleaning --------------------------------------------
#Raw Image Processing
data.raw = microscope.load.data(design)
data.raw = uscope.process.add.ncells(data = data.raw)
lsos() ## check how much memory the object takes

design = uscope.process.estimate.background(data.raw, design)
data.1    = uscope.process.reorder(data.raw, design=design)
data.1    = uscope.process.remove.first.pic(data.1)
data.1    = uscope.process.remove.background(data.1, design)

uscope.count.cells(data.1)

data.1    = uscope.process.remove.small(data.1, MIN.size=1000,MAX.size=2500)
data.1   = uscope.process.BF(data.1)
data.12dan   = uscope.process.remove.BF.outliers(data.1, cutoff=0.8)

data.12dan   = uscope.process.add.ncells(data = data.12dan)

uscope.count.cells(data.12dan)





#Create condensed df for logs
well = design$PDEF$well

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
                        rep (design$PDEF$tech_replica[design$PDEF$well==w], times = length(data.12dan$comp[[w]]$cell)),
                        rep (design$PDEF$bio_replica[design$PDEF$well==w], times = length(data.12dan$comp[[w]]$cell)),
                        rep (design$PDEF$GFP[design$PDEF$well==w], times = length(data.12dan$comp[[w]]$cell)),
                        rep (design$PDEF$RFP[design$PDEF$well==w], times = length(data.12dan$comp[[w]]$cell)),
                        rep (design$PDEF$well[design$PDEF$well==w], times = length(data.12dan$comp[[w]]$cell))
                        )
  
colnames(datad[[w]]) = c("cell", "cellsize","pic","RFP","GFP","RFP_foci", "GFP_foci", "GFP_foci_size","nfoci","x","y", "tech_replica", "biol_replica", "xmer","dimer", "well")
}

tabdan = do.call(rbind, datad)


# NORMALIZATION -----------------------------------------------------------

tabdan$GFPraw = tabdan$GFP
tabdan$RFPraw = tabdan$RFP

tabdan$GFP = tabdan$GFPraw
tabdan$RFP = tabdan$RFPraw

## The Chromaslide gave a mean intensity of 49930.831 when I did the regression for the protein concentration. 
## Now it gives 38455.519.

tabdan$GFP = tabdan$GFP* 1.298405

## Chromaslide then: 5900.757 | Chromaslide now: 3842.237 | 

tabdan$RFP = tabdan$RFP * 1.535761

tabdan$GFPconc = 10^(predict(lm_fit.G, tabdan))
tabdan$RFPconc = 10^(predict(lm_fit.R, tabdan))
tabdan_no_foci = subset(tabdan, GFP_foci == 0 & RFP_foci == 0)


save(list=ls(),file=paste(folderpath,"data_w_norm.RData",sep="/"))


f.add.binodal <- function(l.in = list(), x_data = "RFPconc", y_data = "GFPconc", diagonal_distance = 1, bandwidth = 0.5, fromto = c(0,10)){
  for (name in names(l.in)) {
    
    l.in[[name]] <- list("df" = l.in[[name]], "x" = log2(l.in[[name]][[x_data]]), "y" = log2(l.in[[name]][[y_data]]))
    
    #TRUE / FALSE vector depending on whether a point is less than "1" from the diagonal
    l.in[[name]][["close_to_diag"]]  <-  abs(l.in[[name]]$x-l.in[[name]]$y) < diagonal_distance
    
    l.in[[name]][["avg_dist"]]  <-  (l.in[[name]]$x+l.in[[name]]$y)/2
    
    #logic: density(avg_dist[close_to_diag], from = , to = , bw = )
    l.in[[name]][["density"]]  <- density(l.in[[name]][["avg_dist"]][l.in[[name]][["close_to_diag"]]],from = .5, to = 8, bw = bandwidth)
    
    l.in[[name]][["peak_and_half"]]  <-  get_peak(l.in[[name]]$density, fromto = c(0,10))
  }
  l.in
}



#https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
#get local maxima, take only the last local maximum (usually the second)

#this function finds all local maxima by checking that the second derivative is negative
#then returns the x and y coordinates of the rightmost? maximum
#input has to be a density
# get_peak ------
get_peak <- function(d, fromto){
  
  maxima_indices <- which(diff(sign(diff(d$y)))==-2)+1
  
  if (missing(fromto)) {
    last_index <- tail(maxima_indices, n=1)
    peak_x <- d$x[last_index]
    peak_y <- d$y[last_index]
    
  } else {
    maxima_x <- d$x[maxima_indices]
    x_in_fromto <- maxima_x[which(maxima_x > fromto[1] & maxima_x < fromto[2])]
    peak_x <- tail(x_in_fromto, n=1)
    peak_y <- d$y[which(d$x == peak_x)]
  }
  
  x_after <- d$x > peak_x
  y2 <- d$y
  y2[!x_after] = -1
  
  #This is how you search which entry in vector x is closest to your.number
  #If an entry is exactly your.number then the result is 0, i.e. the minimum
  #which.min(abs(x-your.number))
  half_index = which.min(abs(y2 - ( peak_y/ 2)))
  
  return_vector <- c("x_peak" = peak_x, "y_peak" = peak_y, 
                     "x_at_half" = d$x[half_index], "y_at_half" = d$y[half_index])
  return(return_vector)
}



# SUBSET no foci ----------------------------------------------------------

dat.no_foci <- subset(tabdan, GFP_foci == 0 & RFP_foci == 0)

#List; Each element is an xmer
l.no.foci <- by(data = dat.no_foci, INDICES = dat.no_foci$xmer, FUN = print)
l.no.foci <- f.add.binodal(l.no.foci)




pdf(file= paste(folderpath,"Rplots/binodal_all_xmers.pdf",sep="/"), useDingbats = F)
par(mfrow=c(2,2))


counter = 1
for(entry in names(l_binodal_multimers)[order(names(l_binodal_multimers))]){
  
  plot(l_binodal_multimers[[entry]]$density,main = "")
  
  x_peak <- l_binodal_multimers[[entry]]$peak_and_half["x_peak"]
  y_peak <- l_binodal_multimers[[entry]]$peak_and_half["y_peak"]
  x_peak_half <- l_binodal_multimers[[entry]]$peak_and_half["x_at_half"]
  y_peak_half <- l_binodal_multimers[[entry]]$peak_and_half["y_at_half"]
  
  points(x = x_peak, y = y_peak, col = "red")
  points(x = x_peak_half, y = y_peak_half, col = "blue")
  
  mtext(paste("x =",round(x_peak,2),"x_half =",round(x_peak_half,2)), side = 3)
  
  plot(l_binodal_multimers[[entry]]$x, 
       l_binodal_multimers[[entry]]$y, 
       main = paste(entry,"/ n =", length(l_binodal_multimers[[entry]]$x)), 
       pch = 20, xlab = "", ylab = "", xlim = c(0,14), ylim = c(0,14)
  )
  
  points(x = l_binodal_multimers[[entry]]$x[ l_binodal_multimers[[entry]]$close_to_diag], y = l_binodal_multimers[[entry]]$y[l_binodal_multimers[[entry]]$close_to_diag], pch = 20, col = 3)
  points(x= x_peak_half, y= x_peak_half, bg="blue", pch=22, cex = 2)
  points(x= x_peak, y= x_peak, bg="red", pch=22, cex = 2)
  counter = counter + 1
}

par(mfrow=c(1,1))
barplot(l_multimer_peaks$x_peak, main = "x_peak", names.arg = rownames(l_multimer_peaks), col = c(rep("orange",times = 8), rep("red", times = 3), rep("purple", times = 2), "green", "blue"), las = 2)
barplot(l_multimer_peaks$x_at_half, main = "x_at_half", names.arg = rownames(l_multimer_peaks), col = c(rep("orange",times = 8), rep("red", times = 3), rep("purple", times = 2), "green", "blue"), las = 2)


dev.off()

#SUBSET 4MERS----------------------------------------------------------

tetramer_names = unique(tabdan$xmer)[grepl("v40",unique(tabdan$xmer))]

#Here we use sum() bc tabdan_no_foci$xmer == "value" outputs a vector of the form c(FALSE,FALSE,FALSE,....,TRUE,TRUE,TRUE,....FALSE,FALSE,FALSE
#create vector with element containing number of rows and then take minimum
min_points_tetramer = min(sapply(tetramer_names,function(x) sum(tabdan_no_foci$xmer == x)))

list_tetramer_no_foci_reduced <- vector(mode="list",length = length(tetramer_names))
set.seed(2)
for (i in 1:length(tetramer_names)){
  dat_temp <- subset(tabdan_no_foci,xmer == tetramer_names[i])
  list_tetramer_no_foci_reduced[[i]] <- dat_temp[sample(x=rownames(dat_temp),size=min_points_tetramer),]
}
names(list_tetramer_no_foci_reduced) <- tetramer_names

dat_tetramers_no_foci_reduced <- do.call(rbind,list_tetramer_no_foci_reduced)


#SUBSET 6MERS and p53-----

# CREATE data frame for p53-E9, v601, v602, v603

dat.hexamers = subset(tabdan, xmer == "v601" | xmer == "v602" | xmer == "v603" | xmer == "p53-E9")

dat.hexamers_no_foci = subset(dat.hexamers, GFP_foci == 0 & RFP_foci == 0)

#Reduce all variants to same amount of cells
#calculate occurrences with table and get the minimum of that
min_points = min(table(dat.hexamers_no_foci$xmer))

# CREATE data frame for p53-E9, v601, v602, v603 - no foci - amount of cells reduced to min_points = 6815
#create empty data frame from an existing one that already contains all the columns you want by removing the rows
dat.hexamers_no_foci_reduced <- dat.hexamers[FALSE,]

for (multimer in unique(dat.hexamers$xmer)) {
  set.seed(2)
  df.temp <- subset(dat.hexamers, xmer == multimer & GFP_foci == 0 & RFP_foci == 0)
  df.temp <- df.temp[sample(x = rownames(df.temp), size = min_points),]
  dat.hexamers_no_foci_reduced <- rbind.data.frame(dat.hexamers_no_foci_reduced, df.temp)
}



# Data to export as csv ---------------------------------------------------
dat.p53_2export = dat.p53_no_foci_reduced[c("RFPconc", "GFPconc")]
dat.601_2export = dat.601_no_foci_reduced[c("RFPconc", "GFPconc")]
dat.602_2export = dat.602_no_foci_reduced[c("RFPconc", "GFPconc")]
dat.603_2export = dat.603_no_foci_reduced[c("RFPconc", "GFPconc")]


# Save Reduced Data -------------------------------------------------------
save(list=ls(),file=paste(folderpath,"data_w_norm_reduced.RData",sep="/"))

#Export data to csv that is: w/o foci, reduced to eq. points, calibrated
write.table(x=dat.p53_2export, file=paste(folderpath,"v400_1olg.csv",sep="/"),sep=",")
write.table(x=dat.601_2export, file=paste(folderpath,"v601_2opa_P1A.csv",sep="/"),sep=",")
write.table(x=dat.602_2export, file=paste(folderpath,"v602_4ntk_C27A.csv",sep="/"),sep=",")
write.table(x=dat.603_2export, file=paste(folderpath,"v603_4ok9_wt.csv",sep="/"),sep=",")



# Plot PDF ----------------------------------------------------------------
  
  xmer = unique(design.dan$PDEF$GFP) 
  xmer_new = c(xmer[1:4],xmer[c(13,14,16)],xmer[7:9],xmer[c(11,12)],xmer[5:6],xmer[10])
  xmer = xmer_new
  
  scatterlistd= vector("list")
  for (x in xmer) {
    scatterlistd[[x]] = plotphase.dan(x) 
  }
  
  scatterlist.wfoci= vector("list")
  for (x in xmer) {
    scatterlist.wfoci[[x]] = plotphase.dan.wfoci(x) 
  }
  
  pdf(file= paste(folderpath,"Rplots/phase_wo_foci_calibrated.pdf",sep="/"),useDingbats = F)
  scatterlistd
  dev.off()
  
  pdf(file= paste(folderpath,"Rplots/phase_w_foci_calibrated.pdf",sep="/"),useDingbats = F)
  scatterlist.wfoci
  dev.off()
  
  pdf(file= paste(folderpath,"Rplots/overlay_plots.pdf",sep="/"),useDingbats = F)
  plotphase.dan.overlap(bottom_layer = dat.p53_no_foci_reduced, top_layer = dat.601_no_foci_reduced)
  plotphase.dan.overlap(bottom_layer = dat.p53_no_foci_reduced, top_layer = dat.602_no_foci_reduced)
  plotphase.dan.overlap(bottom_layer = dat.p53_no_foci_reduced, top_layer = dat.603_no_foci_reduced)
  dev.off()
  
  pdf(file= paste(folderpath,"Rplots/overlay_plot_v601.pdf",sep="/"),useDingbats = F)
  plotphase.overlay2(bottom_layer = dat.p53_no_foci_reduced, top_layer = dat.601_no_foci_reduced)
  dev.off()
  
  fn_plot2pdf <- function(path,plot_function,data,...){
    pdf(file=path,useDingbats = F)
    print(plotphase.nofoci(data,...))
    dev.off()
  }

  
  
# Plot Individual ---------------------------------------------------------
  #rewrite this to for loop just as everything above
  
  pdf(file= paste(folderpath,"Rplots/phase_wo_foci_p53.pdf",sep="/"),useDingbats = F)
  plotphase.nofoci(dat.p53_no_foci_reduced,y_label="[Tetramer] (nM)")
  dev.off()
  
  pdf(file= paste(folderpath,"Rplots/phase_wo_foci_v601.pdf",sep="/"),useDingbats = F)
  plotphase.nofoci(dat.601_no_foci_reduced,y_label="[Hexamer] (nM)")
  dev.off()
  
  pdf(file= paste(folderpath,"Rplots/phase_wo_foci_v602.pdf",sep="/"),useDingbats = F)
  plotphase.nofoci(dat.602_no_foci_reduced,y_label="[Hexamer] (nM)")
  dev.off()
  
  pdf(file= paste(folderpath,"Rplots/phase_wo_foci_v603.pdf",sep="/"),useDingbats = F)
  plotphase.nofoci(dat.603_no_foci_reduced,y_label="[Hexamer] (nM)")
  dev.off()
  
  pdf(file= paste(folderpath,"Rplots/phase_wo_foci_empty.pdf",sep="/"),useDingbats = F)
  plotphase.nofoci(y_label="[Hexamer or Tetramer] (nM)")
  dev.off()
  
# Plot 4mer individual ----------------------------------------------------
  
  for(tetramer in tetramer_names){
    #data_temp = subset(dat_tetramers_no_foci_reduced, xmer == tetramer)
    fn_plot2pdf(
      path = paste(folderpath, "/Rplots/", tetramer, "_no_foci_reduced.pdf", sep = ""),
      plot_function = plotphase.nofoci,
      data = list_tetramer_no_foci_reduced[[tetramer]],
      y_label = "[Tetramer] (nM)"
    )
  }
  
# Arrange Plots in 2x2 ----------------------------------------------------

plot.empty <- plotphase.nofoci()
plot.400 <-  plotphase.nofoci(dat.p53_no_foci_reduced,y_label="[Tetramer] (nM)",plot_color = "green")
plot.601 <- plotphase.nofoci(dat.601_no_foci_reduced,y_label="[Hexamer] (nM)", plot_color = "red")
plot.602 <- plotphase.nofoci(dat.602_no_foci_reduced,y_label="[Hexamer] (nM)", plot_color = "red")
plot.603 <- plotphase.nofoci(dat.603_no_foci_reduced,y_label="[Hexamer] (nM)", plot_color = "red")
plot.overlay <- plotphase.dan.overlap(bottom_layer = dat.p53_no_foci_reduced, top_layer = dat.601_no_foci_reduced)


plot_grid(plot.400,plot.empty,plot.overlay,plot.601,scale = 0.9)
plot2x2 <- plot_grid(plot.400,plot.empty,plot.overlay,plot.601,scale = 0.9)
ggsave(paste(folderpath,"Rplots/plot_2x2.pdf",sep="/"))


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



# DIAGONAL -----------


######## Cells per biol_repeat, get minimum
ncells_by_colony <- table(dat.hexamers_no_foci$xmer, dat.hexamers_no_foci$biol_replica)
ncells_by_colony
min(ncells_by_colony)
sort(ncells_by_colony)


# COMBINE POOL+BY.COLONY --------------------------------------------------

########################### FUNCTION DEFINITION ************************

binodal_calc_plot <- function(bandwidth = 0.35, diagonal_distance = 1, x_range_hexamer = c(3.5, 5.5), x_range_tetramer = c(5,7), x_range_pool = c(3.5,7)){

# CREATE LIST ************
  l_binodal <- vector(mode = "list", length = 4)
  names(l_binodal) <- unique(dat.hexamers_no_foci_reduced$xmer)
  
  for (multimer in names(l_binodal)) {
    temp <- subset(dat.hexamers_no_foci_reduced, xmer == multimer)
    
    l_binodal[[multimer]] <- list("x" = log2(temp$RFPconc), "y" = log2(temp$GFPconc))
    l_binodal[[multimer]][["close_to_diag"]]  <-  abs(l_binodal[[multimer]]$x-l_binodal[[multimer]]$y)<diagonal_distance  #TRUE / FALSE vector depending on whether a point is less than "1" from the diagonal
    l_binodal[[multimer]][["avg_dist"]]  <-  (l_binodal[[multimer]]$x+l_binodal[[multimer]]$y)/2
    l_binodal[[multimer]][["density"]]  <- density(l_binodal[[multimer]][["avg_dist"]][l_binodal[[multimer]][["close_to_diag"]]],  #logic: <- density(avg_dist[close_to_diag], from = , to = , bw = )
                                                from = .5, to = 8, bw = bandwidth)
    
    l_binodal[[multimer]][["peak_and_half"]]  <-  get_peak(l_binodal[[multimer]]$density, fromto = x_range_pool)
  }

  

  names_not_reduced <- c()
  for (multimer in unique(dat.hexamers_no_foci$xmer)) {
    for (colony in unique(dat.hexamers_no_foci$biol_replica)) {
      names_not_reduced <- append(names_not_reduced,paste(multimer,"c",colony,sep = "_"))
    }
  }
  

##  by(data = dat.no_foci, INDICES = list("xmer" = dat.no_foci$xmer, "biol_replica" = dat.no_foci$biol_replica),FUN = print)  
    
l_binodal_not_reduced <- vector(mode = "list", length = 24)
names(l_binodal_not_reduced) <- names_not_reduced
for (multimer in unique(dat.hexamers_no_foci$xmer)) {
  for (colony in unique(dat.hexamers_no_foci$biol_replica)) {

    
    entry <- paste(multimer,"c", colony, sep = "_")
    
      temp <- subset(dat.hexamers_no_foci, xmer == multimer & biol_replica == colony)
      
      l_binodal_not_reduced[[entry]] <- list("x" = log2(temp$RFPconc), "y" = log2(temp$GFPconc))
      l_binodal_not_reduced[[entry]][["close_to_diag"]]  <-  abs(l_binodal_not_reduced[[entry]]$x-l_binodal_not_reduced[[entry]]$y)<diagonal_distance  #TRUE / FALSE vector depending on whether a point is less than "1" from the diagonal
      l_binodal_not_reduced[[entry]][["avg_dist"]]  <-    (l_binodal_not_reduced[[entry]]$x+l_binodal_not_reduced[[entry]]$y)/2
      l_binodal_not_reduced[[entry]][["density"]]  <- density(l_binodal_not_reduced[[entry]][["avg_dist"]][l_binodal_not_reduced[[entry]][["close_to_diag"]]],  #logic: <- density(avg_dist[close_to_diag], from = , to = , bw = )
                                                      from = .5, to = 8, bw = bandwidth)
      
      if (substr(entry,1,3) == "p53") {
        l_binodal_not_reduced[[entry]][["peak_and_half"]]  <-  get_peak(l_binodal_not_reduced[[entry]]$density, fromto = x_range_tetramer)
      } else{
        l_binodal_not_reduced[[entry]][["peak_and_half"]]  <-  get_peak(l_binodal_not_reduced[[entry]]$density, fromto = x_range_hexamer)
      }
      
    
  }
}


bycolony_peaks <- as.data.frame(t(unlist(sapply(l_binodal_not_reduced,'[[',"peak_and_half"))))
bycolony_peaks[,"xmer"] <- as.factor(sub("_c.*","",rownames(peaks_not_reduced)))

pool_peaks <- as.data.frame(t(sapply(l_binodal,'[[',"peak_and_half"))[,c("x_peak","x_at_half")])
pool_peaks <- pool_peaks[order(rownames(pool_peaks)),]
bycolonies_mean <- aggregate(x = bycolony_peaks[,c("x_peak","x_at_half")], by = list(XMER = bycolony_peaks$xmer), FUN = mean)
bycolonies_sd <- aggregate(x = bycolony_peaks[,c("x_peak","x_at_half")], by = list(XMER = bycolony_peaks$xmer), FUN = sd)
bycolonies_se <- bycolonies_sd
bycolonies_se[,2:3] <- bycolonies_se[,2:3]/sqrt(6)

if(rownames(pool_peaks) == bycolonies_mean$XMER && colnames(pool_peaks) == colnames(bycolonies_mean)[2:3]) {
  x_diff <- abs(pool_peaks$x_peak-bycolonies_mean[,2])
  names(x_diff) <- rownames(pool_peaks)} else {print("Sort Order Error")}

result_list <- vector(mode = "list", length = 4)
names(result_list) <- c("bycolony_peaks", "pool_peaks", "bycolonies_mean", "x_diff")
result_list[["bycolony_peaks"]] <- bycolony_peaks
result_list[["pool_peaks"]] <- pool_peaks
result_list[["bycolonies_mean"]] <- bycolonies_mean
result_list[["x_diff"]] <- x_diff


# PLOT ******************
pdf(file= paste(folderpath,"Rplots/binodal_analysis.pdf",sep="/"),useDingbats = F)
par(mfrow=c(2,2))
### plot pooled *****

counter = 1
for(entry in names(l_binodal)){
  
  plot(l_binodal[[entry]]$density,main = "")
  
  x_peak <- l_binodal[[entry]]$peak_and_half["x_peak"]
  y_peak <- l_binodal[[entry]]$peak_and_half["y_peak"]
  x_peak_half <- l_binodal[[entry]]$peak_and_half["x_at_half"]
  y_peak_half <- l_binodal[[entry]]$peak_and_half["y_at_half"]
  
  points(x = x_peak, y = y_peak, col = "red")
  points(x = x_peak_half, y = y_peak_half, col = "blue")
  
  mtext(paste("x =",round(x_peak,2),"x_half =",round(x_peak_half,2)), side = 3)
  
  plot(l_binodal[[entry]]$x, 
       l_binodal[[entry]]$y, 
       main = paste(entry,"/ n =", length(l_binodal[[entry]]$x)), 
       pch = 20, xlab = "", ylab = "", xlim = c(0,14), ylim = c(0,14)
  )
  
  points(x = l_binodal[[entry]]$x[ l_binodal[[entry]]$close_to_diag], y = l_binodal[[entry]]$y[l_binodal[[entry]]$close_to_diag], pch = 20, col = 3)
  points(x= x_peak_half, y= x_peak_half, bg="blue", pch=22, cex = 2)
  points(x= x_peak, y= x_peak, bg="red", pch=22, cex = 2)
  counter = counter + 1
}


### plot by colony *****
counter = 1
for(entry in names(l_binodal_not_reduced)){
  plot(l_binodal_not_reduced[[entry]]$density,main = "")
  x_peak <- l_binodal_not_reduced[[entry]]$peak_and_half["x_peak"]
  y_peak <- l_binodal_not_reduced[[entry]]$peak_and_half["y_peak"]
  x_peak_half <- l_binodal_not_reduced[[entry]]$peak_and_half["x_at_half"]
  y_peak_half <- l_binodal_not_reduced[[entry]]$peak_and_half["y_at_half"]
  points(x = x_peak, y = y_peak, col = "red")
  points(x = x_peak_half, y = y_peak_half, col = "blue")
  mtext(paste("x =",round(x_peak,2),"x_half =",round(x_peak_half,2)), side = 3)
  
  plot(l_binodal_not_reduced[[entry]]$x, 
       l_binodal_not_reduced[[entry]]$y, 
       main = paste(entry,"/ n =", length(l_binodal_not_reduced[[entry]]$x)), 
       pch = 20, xlab = "", ylab = "", xlim = c(0,14), ylim = c(0,14)
       )
  
  points(x = l_binodal_not_reduced[[entry]]$x[ l_binodal_not_reduced[[entry]]$close_to_diag], y = l_binodal_not_reduced[[entry]]$y[l_binodal_not_reduced[[entry]]$close_to_diag], pch = 20, col = 3)
  points(x= x_peak_half, y= x_peak_half, bg="blue", pch=22, cex=2)
  points(x= x_peak, y= x_peak, bg="red", pch=22, cex=2)
  counter = counter + 1
}

dodge <-position_dodge(width = .9)
limits_x <- aes(ymax = bycolonies_mean$x_peak + bycolonies_se$x_peak,
              ymin = bycolonies_mean$x_peak - bycolonies_se$x_peak)

limits_x_half <- aes(ymax = bycolonies_mean$x_at_half + bycolonies_se$x_at_half,
                     ymin = bycolonies_mean$x_at_half - bycolonies_se$x_at_half)

print(ggplot(data = bycolonies_mean, aes(x=XMER, y=x_peak)) + 
  geom_bar(stat = "identity", position = dodge, fill = c("darkblue","darkorange","brown","orange")) + 
  geom_errorbar(limits_x, position = dodge,  width = .25))

print(ggplot(data = bycolonies_mean, aes(x=XMER, y=x_at_half)) + 
        geom_bar(stat = "identity", position = dodge, fill = c("darkblue","darkorange","brown","orange")) + 
        geom_errorbar(limits_x_half, position = dodge,  width = .25))

dev.off()
x_diff
}

diff_list <- list()
for(dist in seq(.5,1,.1)) diff_list[[as.character(dist)]] <-  binodal_calc_plot(bandwidth = .53)
dat.diff <- as.data.frame(diff_list)
min(colMeans(dat.diff))

dodge <-position_dodge(width = .9)
limits <- aes(ymax = bycolonies_mean$x_peak + bycolonies_se$x_peak,
              ymin = bycolonies_mean$x_peak - bycolonies_se$x_peak)

ggplot(data = bycolonies_mean, aes(x=XMER, y=x_peak))+geom_bar(stat = "identity", position = dodge)+geom_errorbar(limits, position = dodge,  width = .25)
