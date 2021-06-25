##source the pdef file from the directory
#get all .csv files that have 'pdef' in their name

f.get.pdef <- function(directory){
  pdef_files <- list.files(path = directory,pattern = ".*pdef.*\\.csv")
  if (length(pdef_files)>1) stop("More than 1 pdef file detected. Please keep only 1 pdef file in the directory or specify it manually")
  if (length(pdef_files)<1) stop("No pdef file detected. Please keep exactly 1 file in the directory. The name has to contain the word 'pdef' and it must be a .csv file.")
  pdef_name <- pdef_files
}
