#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-d","--working_dir"), type="character", default=NULL, 
              help="Absolute Path of Working Directory", metavar="character"),
  make_option(c("-n", "--nd_folder_path"), type="character", default=NULL, 
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

print("This script will automatically look for all required input files in the current dir:")
print(getwd())
print("Use options to specify custom input files/paths.")
print("To see all available options run: ./process_results.R --help")
print("It is sufficient to specify individual custom options. The other options will be determined automatically from the current dir.")

cat("Continue? [y/n] ")
ANSWER <- readLines("stdin",n=1)
if (ANSWER != "y" & ANSWER != "Y") {stop("User interrupt")}


#Check which OS ("Darwin" is MacOS)
if (Sys.info()['sysname'] == "Darwin") {
  os.foldername <- "/Volumes"
  source("/Users/d/phase/Multiwell.R")
  source("/Users/d/phase/general.R")
  source("/Users/d/phase/MicroscopeToolBox.R")
  
}else if (Sys.info()['sysname'] == "Linux") {
  os.foldername <- "/media"
  source("/home/d/phase_private_final/uscope_tools/Multiwell.R")
  source("/home/d/phase_private_final/uscope_tools/general.R")
  source("/home/d/phase_private_final/uscope_tools/MicroscopeToolBox.R")
}

#source("/home/d/phase_public_final/postprocess_result_txt.R")
#source("/home/d/phase_public_final/valvar_diploid_plot_functions.R")


if(is.null(opt$working_dir)) {
  opt$working_dir <- getwd()
}

if(is.null(opt$nd_folder_path)){
  #search for .nd file. save to opt$nd_folder_path
  opt$nd_folder_path <- list.files(path = opt$working_dir, pattern = ".nd$", full.names = TRUE, recursive = TRUE)
  if (length(opt$nd_folder_path)>1) stop("More than 1 .nd file detected. Please keep only 1 .nd file in the directory or specify it manually")
  if (length(opt$nd_folder_path)<1) stop("No .nd file detected. Please keep exactly 1 file in the directory. The name has to contain the word .nd")
  }

if (is.null(opt$pdef_path)) {
  #get all .csv files that have 'pdef' in their name
  opt$pdef_path <- list.files(path = opt$working_dir, pattern = ".*pdef.*\\.csv", full.names = TRUE, recursive = TRUE)
  if (length(opt$pdef_path)>1) stop("More than 1 pdef file detected. Please keep only 1 pdef file in the directory or specify it manually")
  if (length(opt$pdef_path)<1) stop("No pdef file detected. Please keep exactly 1 file in the directory. The name has to contain the word 'pdef' and it must be a .csv file.")
}

if (is.null(opt$results_txt_path)) {
  #get the folder that contains the result txt files
  opt$results_txt_path <- list.files(path = opt$working_dir, pattern = "images_output", full.names = TRUE, recursive = TRUE, include.dirs = TRUE, ignore.case = TRUE)
  if (length(opt$results_txt_path)>1) stop("More than 1 results folder detected. Please keep only 1 folder in the directory or specify it manually")
  if (length(opt$results_txt_path)<1) stop("No results folder detected. Please keep exactly 1 folder in the directory. The name has to be images_output (not case sensitive)")
}

cat(
  "The current working directory is: ", opt$working_dir, "\n",
  "The .nd file path is: ",opt$nd_folder_path, "\n",
  "The pdef file path is: ",opt$pdef_path, "\n",
  "The path to the folder containing the result txt files is: ",opt$results_txt_path, "\n",
  sep = ""
)

# pdef_head <- colnames(read.csv(opt$pdef_path,sep=",")) #get column names of pdef file
# CHANNELS_INPUT <- unlist(lapply(pdef_head,function(x) if(x=="GFP" | x=="RFP" | x=="BFP") x))
# 
# design = microscope.get.design(
#   F = opt$nd_folder_path,
#   D = c("comp"),
#   PDEF = opt$pdef_path,
#   FORMAT = opt$wells,
#   OUT = opt$results_txt_path,
#   CHANELS = CHANNELS_INPUT,
#   MEASURE.COL = opt$column_to_measure,
#   DIR.res = opt$out_dir
# )
# 
# postprocess_result_txt(design.file = design, min.cell.size = 1000, max.cell.size = 2500, brightfield.cutoff = 0.8)
# #do normalization now, but ask for it with if condition
# 
# 
# save(list=ls(),file = paste0(data.name,"_ready_to_analyze.RData"))
