# Source Data and Functions -----------------------------------------------
library(scales) 
library(ggplot2)

#WORKING_DIR has to contain data_ready_to_be_analyzed.RData
WORKING_DIR <-  "/Users/d/211116_K42_K43"

file2source <-
  list.files(
    path = WORKING_DIR,
    pattern = "data_ready_to_analyze.RData",
    full.names = TRUE,
    recursive = TRUE
  )

source(file2source)
source("/Users/d/phase_public_final/phase_plot_functions.R")

PLOT_PATH = paste(opt$working_dir,"Rplots",sep = "/")

if( ! file.exists ( PLOT_PATH ) ){
  dir.create(PLOT_PATH, showWarnings = TRUE)
} 

data0 <- dat.ready_to_analyze



# Tidy Up Data ------------------------------------------------------------
data = vector("list")
for (w in well) {
  
  data[[w]] = data.frame (data0[[w]]$cell,
                           data0[[w]]$area,
                              data0[[w]]$pic,
                          
                              data0[[w]]$RFP_int_b5,
                              data0[[w]]$GFP_int_b0,
                          data0[[w]]$BFP_int_b5,
                          
                              data0[[w]]$f1_inRFP_toRFPmed,
                              data0[[w]]$f1_inGFP_toGFPmed,
                          data0[[w]]$f1_inBFP_toBFPmed,
                          
                              data0[[w]]$f1_inGFParea,
                          data0[[w]]$f1_inRFParea,
                          data0[[w]]$f1_inBFParea,
                          
                           data0[[w]]$inRFPnfoci,
                          data0[[w]]$inBFPnfoci,
                          data0[[w]]$inGFPnfoci,
                          
                           data0[[w]]$x,
                           data0[[w]]$y,
                          
                        rep (design.dan$PDEF$rep[design.dan$PDEF$well==w], times = length(data0[[w]]$cell)),
                        rep (design.dan$PDEF$BFP[design.dan$PDEF$well==w], times = length(data0[[w]]$cell)),
                        rep (design.dan$PDEF$RFP[design.dan$PDEF$well==w], times = length(data0[[w]]$cell)),
                        rep (design.dan$PDEF$GFP[design.dan$PDEF$well==w], times = length(data0[[w]]$cell)),
                        rep (design.dan$PDEF$well[design.dan$PDEF$well==w], times = length(data0[[w]]$cell))
                        )
  
  colnames(data[[w]]) = c(
    "cell",
    "cellsize",
    "pic",
    
    "RFP_b5",
    "GFP_b0",
    "BFP_b5",
    
    "RFP_foci",
    "GFP_foci",
    "BFP_foci",
    
    "GFP_foci_size",
    "RFP_foci_size",
    "BFP_foci_size",
    
    "nfoci_RFP",
    "nfoci_BFP",
    "nfoci_GFP",
    
    "x",
    "y",
    
    "rep",
    "tetramer",
    "dimer",
    "kinase",
    "well"
  )
}

tab = do.call(rbind, data)


# Plot and Save -----------------------------------------------------------


pdf(file= PLOT_PATH,useDingbats = F)
plot.overlap.gfp(xmer = "405", kianse_variant = "K42")
dev.off()

