source("https://raw.githubusercontent.com/elevylab/RToolBox2/main/Multiwell.R")
source("https://raw.githubusercontent.com/elevylab/RToolBox2/main/general.R")
## this is a test

###
## Returns the design of an experiment 
##
## or modif 2
microscope.get.design = function(
    F=c("/media/elmicro/Joe/Duplex"),                                        # To be updated by user
    D=c("Duplex1"),                                                          # To be updated by user
    PDEF=c("/media/elusers/data/microscope/plate_defs/Duplex_384_2016.csv"), # To be updated by user
    FORMAT=384,                                                              # To be updated by user
    CHANELS=c("GFP"),                                                        # To be updated by user
    DIR.res = "/data/elevy/00_Papers/008_codon/duplexAna/figs/",             # To be updated by user
    KEY=c(),                                                                 # To be updated by user or leave empty "/media/elusers/data/microscope/plate_defs/Duplex_384_key.csv"
    MEASURE.COL = "_int_b5",                                                 # Should not change
    OUT=c("_output"),                                                         # Should not change
    MINCELLNUM=20
){
  
  if(length(F) != length(D)){
      warning("There is a problem, there should be one short description/code (given in D) for each folder provided")
  }
  
  microscope.design = list()
  
  ## 
  ## Background value
  ## 
  # microscope.design$BACKGROUND     = BACKGROUND
  
  ##
  ##
  microscope.design$CHANELS        = CHANELS
  
  ##
  ##
  microscope.design$MEASURE.COL    = MEASURE.COL
  
  ##
  ##
  microscope.design$DIR.res = DIR.res
  if( ! file.exists ( DIR.res ) ){
      dir.create(DIR.res, showWarnings = TRUE)
      print("the directory given for the results did not exist and was created") }
  
  ## Where files are located --
  ## if multiple plates it may be several folders
  ## /!\ acquisition setting must be the same for all plates)    
  microscope.design$F     = F

  ## Minimum number of cells to consider a picture
  microscope.design$MINCELLNUM     = MINCELLNUM
  
  ##
  ## String that must be appended to the original folder to locate the results folder.
  ##
  microscope.design$OUT     = OUT
  print(paste("OUT=",OUT))
  K=1
  for(each.out in microscope.design$OUT){
    res.files = list.files(path=each.out, pattern='esults.txt')
    print(paste("In plate ",K," there are ",length(res.files)," result files ",sep=""))
    K=K+1
  }
  
  ### TO BE REMOVED LATER, THIS IS ONLY BECAUSE AVITAL NAMED A DIRECTORY INCORECTLY
  # microscope.design$OUT     = gsub(design$OUT,pattern="151",rep="1151")
  
  ## Description of what's in the folders
  ## e.g., c("SD.1","SD.2","YPD.1","YPD.2")
  microscope.design$D    = D
  
  microscope.design$PDEF=NULL
  print(paste("Plate definition file (PDEF) : ",PDEF))
  microscope.design$PDEF = tryCatch(expr = { res=NA; cat("Reading csv plate definition file..."); res = read.csv(PDEF, stringsAsFactors=FALSE) },
                                    error=function(e){ cat('\rReading csv plate definition file...FAIL\n'); return(NA) },
                                    finally = { if( !is.null(dim(res)) ){ cat("\rReading csv plate definition file...OK\n"); res }
                                      else if( is.na(res) ){ stop("File not accessible for reading. Make sure it is not opened on any computer !",call.=FALSE) }
                                    })
  
  if( nrow(microscope.design$PDEF) == 0 ){ stop("the plate definition file is empty") }

  ## Checks that the column names are OK
  if( sum(grepl("well",names(microscope.design$PDEF))) == 0){
      stop("A column should be named 'well' -- perhaps you also need to check the case, e.g., Well will not work")
  }
    
    ## Adds zeros to well ids, e.g., changes A1 to A01 if needed
    microscope.design$PDEF$well = gsub(microscope.design$PDEF$well,pattern="([A-Z])([0-9]$)", replacement="\\10\\2")

    microscope.design$PDEF$Row    = gsub(microscope.design$PDEF$well, pattern="[0-9]+", rep="")
    microscope.design$PDEF$Column = as.numeric(gsub(microscope.design$PDEF$well, pattern="[A-Z]+", rep=""))    
    if( nrow(microscope.design$PDEF) > FORMAT ){
                
      if( sum ( colnames(microscope.design$PDEF)=="plate") ){
          microscope.design$PDEF.split = split(microscope.design$PDEF,microscope.design$PDEF$plate)
          empty.plate = data.frame(well=get.384.by.row(), ord=1:384, stringsAsFactors=FALSE)

          ## Now makes sure that lines in each plate match a well, even if empty or if no information is available.
          for( each.plate in names(microscope.design$PDEF.split)){
              microscope.design$PDEF.split[[each.plate]] = merge(empty.plate, microscope.design$PDEF.split[[each.plate]], by="well", all.x=TRUE)
          }
          
      } else {
          microscope.design$PDEF.split = list()
          empty.plate = data.frame(well=get.384.by.row(), ord=1:384, stringsAsFactors=FALSE)
          microscope.design$PDEF.split[[1]] = merge(empty.plate, microscope.design$PDEF, by="well", all.x=TRUE)

          warning(paste("the format specified (",FORMAT," wells) does not correspond to this plate defition file"), sep="")
      }

    } else {

        if(  nrow(microscope.design$PDEF) == FORMAT ){ 
            microscope.design$PDEF.split = list()
            for( each.plate in 1:length(microscope.design$D)){
                microscope.design$PDEF.split[[each.plate]] = microscope.design$PDEF
            }
        } else {
            stop("The plate definition file should contain at least as many lines as the plate format")
        }
      
  }
  
  microscope.design$KEY=NULL
  
  if(length(KEY)>0){
      print(paste("Plate definition KEY file (KEY) : ",KEY))
      
      microscope.design$KEY = tryCatch(
          expr = { res=NA; cat("Reading csv plate KEY file..."); res = read.csv(KEY, stringsAsFactors=FALSE) },
          error=function(e){ cat('\rReading csv plate KEY file...FAIL\n'); return(NA) },
          finally = { if( !is.null(dim(res)) ){ cat("\rReading csv plate KEY file...OK\n"); res }
          else if( is.na(res) ){ stop("File not accessible for reading. Make sure it is not opened on any computer !",call.=FALSE) }
                  })
      
      if( nrow(microscope.design$KEY) == 0 ){ stop("the plate definition file is empty") }

      for(each.ch in CHANELS){
          if( sum(grepl(each.ch,names(microscope.design$KEY))) == 0 ){
              warning(paste0("The channel ",each.ch," is given in the experiment design but does not appear in the plate KEY file"))
          }
      }
  }
      
  for(each.ch in CHANELS){
    if( sum(grepl(each.ch,names(microscope.design$PDEF))) == 0 ){
      warning(paste0("The channel ",each.ch," is given in the experiment design but does not appear in the plate definition file"))
    }
  }
  
  ## PLATE FORMAT 96/384
  microscope.design$FORMAT    =FORMAT  ##
  
  ## First, looks for the nd files in each folder
  microscope.design$ND=c()
  microscope.design$BASENAMES=c()

  
  ## Problem to resolve in the future: if there are multiple folders, we may need to use multiple plate def. files.
  K=1
  for( each.folder in F){
    
    print(paste("Going through folder",each.folder))
    if( !file.exists(each.folder) ){ stop("the directory given for the pictures does not exist") }
    
    if(substr(each.folder,nchar(each.folder),nchar(each.folder)) == "/"){
      microscope.design$F[K] = substr(each.folder,1,nchar(each.folder)-1)
    }
    
    ## "/media/elusers/data/microscope/or/06noise/151109_dip1_take1_SDfullselection/"
    # nd.file = list.files(path="/media/elusers/data/microscope/or/06noise/151109_dip1_take1_SDfullselection/", pattern=".nd")
    nd.file = list.files(path=each.folder, pattern=".nd$")
    
    print(paste("nd file is", nd.file))
    microscope.design$BASENAMES = c(microscope.design$BASENAMES,sub(nd.file, pattern=".nd$",replacement=""))
    
    #nd.file = list.files(path=each.folder, pattern=".nd")
    if(length(nd.file) > 1){
      warning(paste("Problem: there is more than one nd file in", each.folder, " (",length(nd.file),")" ))
    }
    microscope.design$ND = c(microscope.design$ND, paste(each.folder,"/",nd.file, sep=""))     
    
    nd = read.table(microscope.design$ND[K], sep="\n", stringsAsFactors=FALSE)
    nd = nd[grep("^Stage",nd[,1]),]
    
    microscope.design$S.positions[[K]] = sub("Stage([0-9]+), [A-Z][0-9]+\\.([0-9]+)", "\\1", nd,perl=TRUE)
    microscope.design$WELLS[[K]]       = sub("Stage([0-9]+), ([A-Z][0-9]+)\\.([0-9]+)", "\\2", nd,perl=TRUE)
    microscope.design$WELLS[[K]]       = gsub(microscope.design$WELLS[[K]],pattern="([A-Z])([0-9]$)", replacement="\\10\\2")
    microscope.design$PICNUM[[K]]      = sub("Stage([0-9]+), ([A-Z][0-9]+)\\.([0-9]+)", "\\3", nd,perl=TRUE)

    ## Then we loop over all well names, see if the weel exists and if it does not, we concatenate an "empty" position to these vectors.
    all.wells = c()
    if(microscope.design$FORMAT==96){
        all.wells = get.96.by.row()          
    } else {
        all.wells = get.384.by.row()          
    }
    
    for(each.well in all.wells){

        if( sum(each.well %in% microscope.design$WELLS[[K]])>0){
            # do nothing
        } else {
            microscope.design$WELLS[[K]] = c(microscope.design$WELLS[[K]], each.well)
            microscope.design$S.positions[[K]] = c(microscope.design$S.positions[[K]], 0)
            microscope.design$PICNUM[[K]] = c(microscope.design$PICNUM[[K]], 0)            
        }
    }
      
    K=K+1
    
  }
  
  ## Extracts the plate definition information to establish well to strain correspondancies 
  
  return(microscope.design)
}

###
### This function cleans up the plate def file.
### In particular, it add "empty" positions where info is missing
###
cleanup.platedef = function(plateDEF){

    ## If a column is called ord, we ignore it
    ord.col = grep("ord",colnames(plateDEF))
    if( length(ord.col)>0 ){        
        plateDEF = plateDEF[, -c(ord.col)]
    }
    
    ## First it checks how many plates there are
    plate.colnames = colnames(plateDEF)
    
    plate.col = grep(".*plate.*",plate.colnames, ignore.case=TRUE, perl=TRUE)
    
    if( length(plate.col) > 0){
        plate.colnames[plate.col[1]]="plate"

        if(length(plate.col) > 1){
            warning(paste("More than one column was identified with the string 'plate', this is likely problematic -- restart the function after proper naming of the excel file", sep=""))
        }

        colnames(plateDEF) = plate.colnames        
        NPLATES = max(plateDEF$plate)

    } else {
        NPLATES = 1
        plateDEF$plate=1
    }

    ## Second reformat the well numbers
    plateDEF$well = gsub(plateDEF$well, pattern="([A-Z])([0-9]$)", replacement="\\10\\2")

    ## Third it fills the empty positions
    plateBASE = data.frame(plate=rep(1:NPLATES,rep(384,NPLATES)), well=rep(get.384.by.row(),NPLATES), ord = 1:( length(rep(get.384.by.row(),NPLATES))) )
    
    ## Fourth it merges the plate def with the base plate.
    plateFINAL = merge(plateBASE, plateDEF, by=c("plate","well"), all.x=TRUE)

    plateFINAL = plateFINAL[order(plateFINAL$ord),]

    plateFINAL[is.na(plateFINAL)]="empty"
    return(plateFINAL)
}


###

microscope.load.data = function(design){
  
  data = list()
  K=1
  cat("\n")
  NCOL_ALL = 0
  HEADER = c()
  
  for(each.plate in design$F){
    Nb.wells.witohut.res = 0
    
    cat(paste("Processing plate ",each.plate,"\n"))
    
    data[[K]]=list()
    
    ## First extracts the screen information to establish "s" number to well correspondancies.
    S.positions = design$S.positions[[K]]
    WELLS       = design$WELLS[[K]]
    PICNUM      = design$PICNUM[[K]]
    
    cat(paste("Loading plate ",each.plate,", there are ",length(S.positions)," positions\n"));
    
    WELL.tmp = WELLS[1]
    WELL.index   = 1
    
    ### First we search for a non-empty file to create a blank template
    file.num=1
    FILE.NOT.FOUND=1
    while(FILE.NOT.FOUND & file.num < length(S.positions) ){
      
      res.file.name = paste0(design$OUT[K],"/",design$BASENAMES[K],"s",S.positions[file.num],"Results.txt")
      
      if( file.exists(res.file.name) ){                
        res.file = read.csv(res.file.name, sep="\t", row.names=c(1), stringsAsFactors=FALSE)
        HEADER = colnames(res.file)
        cat(paste("It found a non empty file (",res.file.name,") and the number of columns is ",NCOL(res.file),"\n"))
        res.file$pic=1

        if(K==1){
            NCOL_ALL=NCOL(res.file)
        }
        
        res.file = res.file[1,]
        res.file = res.file[ -c(1),]
        ## print(res.file)
        write.csv(file="/tmp/empty_file.csv", x=res.file, row.names=FALSE)
        FILE.NOT.FOUND=0
      }
      file.num = file.num+1
    }
    
    if( file.num == length(S.positions)){
        stop(paste("in plate ",each.plate,": all the result files were empty, there is no data to load"))
    }
    
    ### Then we load the data
    for( N in 1:length(S.positions) ){
      
      if(N %% 100 == 0){
          cat(paste("\rProcessing position number ",N,"    "))
      }

      res.file.name = paste(design$OUT[K],"/",design$BASENAMES[K],"s",S.positions[N],"Results.txt",sep="")
      
      res.file=matrix(2,2)
      
      # an example file to load without X/Y
      #res.file = read.csv("/media/elusers/data/microscope/or/06noise/151109_dip1_take1_SDfullselectionoutput/noise_dip1_tak11s1Results.txt",
      #    , sep="\t", row.names=c(1))
      
      ## an example file to load WITH X/Y
      # res.file = read.csv("/media/elusers/data/microscope/or/06noise/151117_YMD_selection_saturated_output/noise_dip1_tak11s1Results.txt", sep="\t", row.names=c(1))
      #    , sep="\t", row.names=c(1))
      
      
      if( file.exists(res.file.name) ){                
        res.file = read.csv(res.file.name, sep="\t", row.names=c(1), stringsAsFactors=FALSE)
        
        ## It seems that sometimes the column names get mixed up by the javascript, so lets check this
        if( sum(colnames(res.file) == HEADER) != length(HEADER)  ){
            print(paste("There is a problem with file s",N, " as the file header is different from the others and the number of cells in that file is", nrow(res.file),sep=""))
        }
        
        res.file$pic=PICNUM[N]
      } else {
        Nb.wells.witohut.res = Nb.wells.witohut.res+1
        res.file = read.csv("/tmp/empty_file.csv", stringsAsFactors=FALSE)
        #res.file$pic=PICNUM[N]
      }
      
      ## If the file was not read, or if the file read had a bad format
      ## we use a generic file instead
      if(NROW(res.file)<design$MINCELLNUM){
        Nb.wells.witohut.res = Nb.wells.witohut.res+1
        res.file = read.csv("/tmp/empty_file.csv", stringsAsFactors=FALSE)                            
      }
      
      if(WELLS[N] != WELL.tmp) {
        WELL.index = WELL.index+1
      }                      
      
      if( length(data[[K]])==0 | (WELLS[N] != WELL.tmp) ){
        data[[K]][[WELL.index]] = res.file
      } else {
        data[[K]][[WELL.index]] = rbind(data[[K]][[WELL.index]],res.file)
      }

      if( NCOL(data[[K]][[WELL.index]]) != (NCOL_ALL) ) {
          cat(paste("There is a problem with plate number",K," the number of columns recorded in the first plate (",(NCOL_ALL),") does not match this one (",NCOL(data[[K]][[WELL.index]]),")\n"))
      }
      
      WELL.tmp = WELLS[N]
      
    }

      #for(each.well in 1:design$FORMAT){

       #   if(data[[K]][[each.well]]
      #}
      
    cat("\n")
    names(data[[K]]) = unique(WELLS)
    print(paste("For plate",K," there are ",Nb.wells.witohut.res," pictures without data"))
    K=K+1
  }
  
  names(data) = design$D
  return(data)
}


uscope.process.reorder = function(data, design){
  
    for(K in 1:length(data)){
        wells.data = names(data[[K]])

        if(length(design$PDEF.split)>1){
            wells.to.match = design$PDEF.split[[K]]$well
            data[[K]] = data[[K]][match(wells.to.match, wells.data)]
        } else {
            wells.to.match = design$PDEF.split[[1]]$well
            data[[K]] = data[[K]][match(wells.to.match, wells.data)]
        }
    }    
    return(data)  
}

uscope.process.estimate.background = function(data, design){
  
    design$BACKGROUND=list()
    
    for(K in 1:length(data)){
        
        design$BACKGROUND[[K]]=list()
        
        for(each.ch in design$CHANELS){
            
            cols.fluo = intersect(grep("_int_b9",colnames(data[[K]][[1]])),grep(paste("^",each.ch,"_",sep=""),colnames(data[[K]][[1]])))
            
            tmp=sapply(data[[K]],function(x){
                if(length(x[,cols.fluo]) <= 20){
                    return(NA)
                } else {
                    return( min(x[,cols.fluo]) )
                }
            })            

            if(median(tmp, na.rm=TRUE) < 80){
                print(paste("Plate ",K," and chanel ",each.ch,": strange, median below 80 -- tmp="))
                print(tmp)
            }
            
            design$BACKGROUND[[K]][[each.ch]] = median(tmp, na.rm=TRUE)
            
    }    
    }    
    return(design)  
}


uscope.process.add.ncells = function(data){

    cur.cols = colnames(data[[1]][[1]])
    if(length(grep("ncells",cur.cols))>0){
        new.col="ncells.2"
    } else {
        new.col="ncells"
    }
    
  for(K in 1:length(data)){
    
      for(L in 1:length(data[[K]])){
          
          ncells = matrix(data=NROW(data[[K]][[L]]), ncol=1, nrow=NROW(data[[K]][[L]]))
          data[[K]][[L]] = cbind(data[[K]][[L]], ncells)
          colnames(data[[K]][[L]])[ ncol(data[[K]][[L]]) ] = new.col          
      }
  }
  return(data)
}

###

uscope.process.remove.first.pic = function(data, N=1){
  
  for(K in 1:length(data)){
    
    for(L in 1:length(data[[K]])){
        #print(paste(K,L))
        if(dim(data[[K]][[L]])[1] > 0){
 #           print(paste(K,L))
            data[[K]][[L]] = data[[K]][[L]][data[[K]][[L]]$pic>N,]
        }
    }
  }
  return(data)
}

###


uscope.process.define.empty.pos = function(data){
  
  for(K in 1:length(data)){
    
    for(L in 1:length(data[[K]])){

        if(is.null(data[[K]][[L]])){
            print(paste("filling in position",L))
            data[[K]][[L]] = data[[K]][[1]][1,]
            data[[K]][[L]] = data[[K]][[L]][-c(1),]
        }
    }
  }
  return(data)
}

###


uscope.process.remove.small = function(data, MIN.size=1000, MAX.size=2000){
     
  for(K in 1:length(data)){
    
    for(L in 1:length(data[[K]])){

        if( sum(colnames(data[[1]][[1]])=="Area") ){
            data[[K]][[L]] = data[[K]][[L]][data[[K]][[L]]$Area>MIN.size & data[[K]][[L]]$Area < MAX.size,]
        } else {
            data[[K]][[L]] = data[[K]][[L]][data[[K]][[L]]$area>MIN.size & data[[K]][[L]]$area < MAX.size,]
        }
    }
  }
  return(data)
}

uscope.process.add.meanArea = function(data){
  
    for(K in 1:length(data)){
    
        for(L in 1:length(data[[K]])){
            if(nrow(data[[K]][[L]]) > 0) {
                data[[K]][[L]]$meanArea = round(mean(data[[K]][[L]]$area))
            }
        }
    }
    return(data)
}

uscope.process.add.meanArea.per.pic = function(data){
  
    for(K in 1:length(data)){
    
        for(L in 1:length(data[[K]])){
            if(nrow(data[[K]][[L]]) > 0) {

                for( each.pic in unique(data[[K]][[L]]$pic)){
                    data[[K]][[L]]$meanArea[data[[K]][[L]]$pic==each.pic] = round(mean(data[[K]][[L]]$area[data[[K]][[L]]$pic==each.pic]))
                }
            }
        }
    }
    return(data)
}


uscope.process.remove.empty.rows = function(data){

      for(K in 1:length(data)){
    
          for(L in 1:length(data[[K]])){

              if(nrow(data[[K]][[L]]) > 0) {
                  data[[K]][[L]] = data[[K]][[L]][ data[[K]][[L]]$x > 0,]
              }
          }
      }
    return(data)    
}

###

##
## The brightfield image serve as a source of information to gate on "GOOD" cells.
## The idea is that each cell has a number of very dark / dark / medium / light / very light grey pixels.
## We assume that most cells are good, and so keep those cells which profile of pixel colors is average
##
uscope.process.BF = function(data){
  
  for(K in 1:length(data)){

      if( length(grep("InFocBin[0-9]",colnames(data[[1]][[1]])))>0){

          cols.out.of.focus = grep("OutofFocBin[0-9]",colnames(data[[K]][[1]]))
          cols.in.focus     = grep("InFocBin[0-9]",colnames(data[[K]][[1]]))
          AREA = "Area"
          
      } else {

          cols.out.of.focus = grep("BFO_lev_b[0-9]",colnames(data[[K]][[1]]))
          cols.in.focus     = grep("BFI_lev_b[0-9]",colnames(data[[K]][[1]]))
          AREA = "area"
      }
      
    for(L in 1:length(data[[K]])){
      
      PICS = as.numeric(unique(data[[K]][[L]]$pic))

      data[[K]][[L]]$in.focus.d = double( length(data[[K]][[L]]$pic) )
      data[[K]][[L]]$out.focus.d = double( length(data[[K]][[L]]$pic) )
      
      for(PIC in PICS){

        BF.info   = data[[K]][[L]][data[[K]][[L]]$pic==PIC, cols.in.focus]
        areas     = data[[K]][[L]][data[[K]][[L]]$pic==PIC, AREA]
        area.info = matrix(areas, ncol=NCOL(BF.info), nrow=NROW(BF.info), byrow=FALSE)
        BF.norm   = BF.info / area.info
        BF.mean   = colMeans(BF.norm)
        BF.dist   = as.vector(apply(BF.norm,1,function(x){ sum((x-BF.mean)^2) }))
        data[[K]][[L]]$in.focus.d[data[[K]][[L]]$pic==PIC] = (rank(BF.dist) / length(BF.dist))
        
        BF.info   = data[[K]][[L]][data[[K]][[L]]$pic==PIC, cols.out.of.focus]
        BF.norm   = BF.info / area.info
        BF.mean   = colMeans(BF.norm)
        BF.dist   = as.vector(apply(BF.norm,1,function(x){ sum((x-BF.mean)^2) }))
        data[[K]][[L]]$out.focus.d[data[[K]][[L]]$pic==PIC] = (rank(BF.dist) / length(BF.dist))
      }
      
    }
  }
  return(data)
}

##
## The fluorescence intensity across the field of view is not homogeneous
## Thus we remove the areas with extreme intensity different, that is the center and the corners
##
uscope.estimate.centroid.center = function(data){
  
  centroids = list()
  
  for(K in 1:length(data)){
    
    
    X.cent = 1900
    Y.cent = 1100
    means = c() #matrix(ncol=1,nrow=4)
    for(L in 1:length(data[[K]])){
      
      tmp = data[[K]][[L]]            
      ## 3 => top right
      ## 2 => bottom right
      ## 1 => top left
      ## 0 => bottom left
      if(NROW(tmp) > 2000){
        tmp$Q = 0 +   3*(tmp$xCentroid > X.cent & tmp$yCentroid > Y.cent) +
          2*(tmp$xCentroid > X.cent & tmp$yCentroid < Y.cent) +
          1*(tmp$yCentroid > Y.cent & tmp$xCentroid < X.cent)
        
        if(length(table(tmp$Q))==4){
          agg = aggregate(tmp$GFPavgIntBin7,by=list(Q = tmp$Q), mean)
          means = cbind(means, agg$x/max(agg$x))
        }
      }
    }
  }
  return(data)
}


##
## The fluorescence intensity across the field of view is not homogeneous
## Thus we remove the areas with extreme intensity different, that is the center and the corners
##
uscope.process.centroid = function(data, X.center=1004, Y.center=1024){
  
  for(K in 1:length(data)){
    
    for(L in 1:length(data[[K]])){

        if(dim(data[[K]][[L]])[1] > 0){
            
            data[[K]][[L]]$radial.dist =
                (
                    (data[[K]][[L]]$x-X.center)^2 +
                    (data[[K]][[L]]$y-Y.center)^2
                )^0.5
            }
    }
  }
  return(data)
}

##
## The data structure is filtered to remove all pictures where there are less than a given number of cells
##
uscope.process.remove.pics.without.cell = function(data, MIN.cells=5){
  
  for(K in 1:length(data)){
    
    for(L in 1:length(data[[K]])){
      
      n.cell.per.pic = table(data[[K]][[L]]$pic)
      to.remove = as.numeric(names(n.cell.per.pic)[which(n.cell.per.pic<MIN.cells)])
      if(length(to.remove)>0){
        data[[K]][[L]] = data[[K]][[L]][ ! data[[K]][[L]]$pic %in% to.remove,]
      }
    }
  }
  return(data)
}

###
### Filters outlier cells
###
uscope.process.remove.centroid.outliers = function(data, min=200,max=900){
  
  for(K in 1:length(data)){
    
    for(L in 1:length(data[[K]])){
      
        data[[K]][[L]] = data[[K]][[L]][which(data[[K]][[L]]$radial.dist > min & data[[K]][[L]]$radial.dist < max),]
    }
  }
  return(data)
}

###
### After the BF has been processed we discard a fraction of cells that look different from the average
###
uscope.process.remove.BF.outliers = function(data, cutoff){
  
  for(K in 1:length(data)){
    
    cols.out.of.focus = grep("BFO_lev_b[0-9]",colnames(data[[K]][[1]]))
    cols.in.focus     = grep("BFI_lev_b[0-9]",colnames(data[[K]][[1]]))
    
    for(L in 1:length(data[[K]])){
      
      PICS = as.numeric(unique(data[[K]][[L]]$pic))
      
      for(PIC in PICS){
        #print(paste("PIC number =",PIC) )
        to.remove = which( data[[K]][[L]]$pic==PIC & (data[[K]][[L]]$in.focus.d > cutoff | data[[K]][[L]]$out.focus.d > cutoff) )
        #print(paste("length=",length(to.remove)))
        if(length(to.remove)>0){
          data[[K]][[L]] = data[[K]][[L]][-to.remove,]
        }
      }
    }
  }
  return(data)
}

## 
## Removes the background value to all values
##
uscope.process.remove.background = function(data, design){
  
  for(K in 1:length(data)){
    
    for(L in 1:length(data[[K]])){
      
      for( each.ch in design$CHANELS){
        
        col.of.interest =
                    ( grepl(paste("^",each.ch,"_",sep=""),colnames(data[[K]][[L]])) |
                      grepl(paste("to",each.ch,"m",sep=""),colnames(data[[K]][[L]]))
                    ) & !grepl("_bg$",colnames(data[[K]][[L]]))
          
        col.of.interest = colnames(data[[K]][[L]])[col.of.interest]
        
        for( each.col in col.of.interest){

          data[[K]][[L]][,each.col] = data[[K]][[L]][,each.col]- design$BACKGROUND[[K]][[each.ch]]
          below.zero = which(data[[K]][[L]][,each.col] < 0)
          if(length(below.zero)>0){
            # print("Warning, background value may be set too high as values are below it - these are set to 0");
            data[[K]][[L]][below.zero,each.col]=0
          }
        }
      }
    }
  }
  return(data)
}


## 
## Removes the background value to each cell based on the background image computed (requires >April 2021 JS script from emmanuel ) 
##
uscope.process.remove.background.per.cell = function(data, design){
  
    for(K in 1:length(data)){
        
        for(L in 1:length(data[[K]])){
            
            for( each.ch in design$CHANELS){
                
                col.of.interest =
                    ( grepl(paste("^",each.ch,"_",sep=""),colnames(data[[K]][[L]])) |
                      grepl(paste("to",each.ch,"m",sep=""),colnames(data[[K]][[L]]))
                    ) & !grepl("_bg$",colnames(data[[K]][[L]]))
                
                col.of.interest = colnames(data[[K]][[L]])[col.of.interest]
                
                col.of.background = paste0(each.ch,"_int_bg")
                
                for( each.col in col.of.interest){

                    if(K==1 & L==1){
                        print(paste0("For channel ",each.ch," COL=",each.col," we are substrcting the background in ",col.of.background))
                    }
                    data[[K]][[L]][,each.col] = data[[K]][[L]][,each.col] - data[[K]][[L]][,col.of.background]
                    
                    below.zero = which(data[[K]][[L]][,each.col] < 0)
                    
                    if(length(below.zero)>0 & substr(each.col, 1, 1) != "f" & substr(each.col, nchar(each.col)-2, nchar(each.col)) != "_b8" & substr(each.col, nchar(each.col)-2, nchar(each.col)) != "_b9" & substr(each.col, nchar(each.col)-2, nchar(each.col)) != "b10"){
                        print(paste0("In plate=",K,", well=",L,", channel", each.ch, ", col=",each.col, ", some cells end-up with negtive value (this may happen although rare) - it has been set to 0"));
                        data[[K]][[L]][below.zero,each.col]=0
                    }
                }
            }
        }
    }
    return(data)
}


## 
## Removes the background value to all values
##
uscope.process.rename.undefined = function(data){
  
  for(K in 1:length(data)){
    
      for( each.well in 1:length(data[[K]])){
        for( each.col in 1:ncol(data[[K]][[each.well]])){
          
          undef.pos = grep("undefined",data[[K]][[each.well]][,each.col])
          if( length(undef.pos)>0 ){
            data[[K]][[each.well]][undef.pos,each.col] = NA
          }
          
          if(typeof( typeof(data[[K]][[each.well]][1,each.col])) == "character" ){
            data[[K]][[each.well]][,each.col] = as.numeric(as.vector(data[[K]][[each.well]][,each.col]))
          }
        }
      }
  }
  return(data)
}

uscope.process.remove.lowfluo.cells = function(data, fraction=0.8){
  
  for(K in 1:length(data)){
    
    for(L in 1:length(data[[K]])){
      
      for( each.ch in design$CHANELS){
        my.col = grep(paste(each.ch,design$MEASURE.COL,"$",sep=""), colnames(data[[K]][[1]]))
        CUTOFF = quantile( data[[K]][[L]][,my.col] , prob=c(0,fraction,1))[2]                 
        data[[K]][[L]] = data[[K]][[L]][ data[[K]][[L]][,my.col] > CUTOFF,]                
      }
    }
  }
  return(data)    
}

uscope.process.remove.highfluo.cells = function(data, fraction=0.8){
  
    for(K in 1:length(data)){
    
        for(L in 1:length(data[[K]])){
      
            for( each.ch in design$CHANELS){
                my.col = grep(paste(each.ch,design$MEASURE.COL,"$",sep=""), colnames(data[[K]][[1]]))
                CUTOFF = quantile( data[[K]][[L]][,my.col] , prob=c(0,(1-fraction),1))[2]                    
                data[[K]][[L]] = data[[K]][[L]][ data[[K]][[L]][,my.col] < CUTOFF,]                
            }
        }
    }
    return(data)    
}

## 
## Removes the extreme few percents from the data
##
uscope.process.filter.fluo = function(data, design, PERCENT){
  
    for(K in 1:length(data)){
        
        for(L in 1:length(data[[K]])){
            
            bads = rep(0, length(data[[K]][[L]][,1]))
            
            for ( each.CH in design$CHANELS){
                
                cols.of.interest = grep( paste(each.CH,"_int_b5",sep=""), colnames(data[[K]][[1]]))
                
                GFPm = data[[K]][[L]][,cols.of.interest]
                Q.GFP = quantile(GFPm,prob=seq(0,1,len=101))
                
                bads =   bads | (GFPm <= Q.GFP[PERCENT]) |
                    (GFPm >= Q.GFP[100-PERCENT])
            }
            
            data[[K]][[L]] = data[[K]][[L]][-which(bads),]
        }
    }
    return(data)
}

##
## We assume that the mean fluorescence should be the same across different pics
## If it's not, then most likely that cells were out of focus.
##
## This is not really a good way to filter things, I should keep it simple (see remove.oof.pics instead).
##
uscope.process.remove.outlier.pics = function(data, TOLERANCE=3){
  
  for(K in 1:length(data)){
    
    cols.fluo = grep("_int_b",colnames(data[[K]][[1]]))
    
    for(L in 1:length(data[[K]])){
      
      #PICS = as.numeric(unique(data[[K]][[L]]$pic))
      
      not.converged=1
      max.cycles = 3
      while(NROW(data[[K]][[L]])>0 & not.converged & max.cycles > 0){
        
        agg = aggregate(data[[K]][[L]][,cols.fluo], by=list(pic=data[[K]][[L]]$pic), mean)
        
        ## Check the differences using the chanel with strongest signal
        if(NCOL(agg)==2){
          index=2
        } else {
          index = which.max(colMeans(agg[,-c(1)]))+1
        }
        
        ## SD across pictures
        my.sd = sd(agg[,index])
        
        ## if SD across pictures significantly decrease when removing a particular pic, we remove it.
        to.remove=c()
        for(each.sample in 1:nrow(agg)){
          
          my.new.sd = sd(agg[-each.sample,index])
          if(!is.na(my.new.sd) & my.new.sd*TOLERANCE < my.sd ){
            to.remove=c(to.remove, agg$pic[each.sample])
          }                 
        }
        
        if(length(to.remove)>0){
          
          data[[K]][[L]] = data[[K]][[L]][ ! data[[K]][[L]]$pic %in% to.remove,]
          
          warning(paste("In plate ",K,"and well ",L, " the picture(s) ",paste(to.remove,",")," was/were removed"))
        } else {                    
          not.converged=0
        }
        max.cycles = max.cycles-1
      }
      
      if(NROW(data[[K]][[L]])>0){
        agg = aggregate(data[[K]][[L]][,cols.fluo], by=list(pic=data[[K]][[L]]$pic), mean)
        
        ## Check the differences using the chanel with strongest signal
        if(NCOL(agg)==2){
          index=2
        } else {
          index = which.max(colMeans(agg[,-c(1)]))+1
        }
        my.sd = sd(agg[,index])
        my.mean = mean(agg[,index])
        
        data[[K]][[L]]$pic.sd = my.sd/my.mean
      } else {
        data[[K]][[L]][1,]=NA
        data[[K]][[L]]$pic.sd = NA
        data[[K]][[L]] = data[[K]][[L]][-c(1),]
      }
    }
  }
  return(data)
}


##
## We assume that the mean fluorescence should be the same across different pics
## If it's not, then most likely that cells were out of focus.
##
## This is not really a good way to filter things, I should keep it simple (see remove.oof.pics instead).
##
uscope.process.remove.nofluo.pics = function(data){
  
  for(K in 1:length(data)){
    
    cols.fluo = grep("int_b",colnames(data[[K]][[1]]))
    
    for(L in 1:length(data[[K]])){
      
      PICS = as.numeric(unique(data[[K]][[L]]$pic))
      
      if(NROW(data[[K]][[L]])>2){
        
        agg = aggregate(data[[K]][[L]][,cols.fluo], by=list(pic=data[[K]][[L]]$pic), mean)
        
        low.fluo = which( rowSums(agg[,-c(1)] < 5)>0 )
        
        if(length(low.fluo)>0){
          pics.to.remove = agg$pic[low.fluo]
          warning(paste("In plate ",K,"and well ",L, " the picture(s) ",paste(pics.to.remove,",")," was/were removed"))
          data[[K]][[L]] = data[[K]][[L]][ ! (data[[K]][[L]]$pic %in% pics.to.remove) , ]
        }
      }
    }
  }
  return(data)
}


###
### We know that focus can be lost between weel but not within wells.
### Therefore, we check if the intensity of the first pic(s) are lower than subsequent ones
### and if so, we remove them.
###
uscope.process.remove.oof.pics = function(data){
  
  for(K in 1:length(data)){
    
    cols.fluo = grep("Brigh30Per",colnames(data[[K]][[1]]))
    
    for(L in 1:length(data[[K]])){           
      
      if(NROW(data[[K]][[L]])>50 & length(unique(data[[K]][[L]]$pic))>1 ){
        
        agg = aggregate(data[[K]][[L]][,cols.fluo], by=list(pic=data[[K]][[L]]$pic), mean)
        
        n.rows = NROW(agg)
        
        ## Check the differences using the chanel with strongest signal
        if(NCOL(agg)==2){
          index=2
        } else {
          index = which.max(colMeans(agg[,-c(1)]))+1
        }
        
        ## I need to compare the intensity to a robust measure
        ## If there are not enough cells in the last picture ... it's a problem
        ## So I take the last 100 cells of the well 
        if( sum(data[[K]][[L]]$pic==agg[n.rows,1]) > 50){
          my.mean.index = which(data[[K]][[L]]$pic==agg[n.rows,1])
        } else {
          my.mean.index = which(data[[K]][[L]]$pic==agg[n.rows,1] | data[[K]][[L]]$pic==agg[n.rows-1,1])
        }
        
        print(paste("K=",K," - L=",L, " and l=",length(my.mean.index)))
        
        my.mean.fluo = mean(data[[K]][[L]][my.mean.index, cols.fluo[index-1]], na.rm=TRUE)
        
        to.remove=c()
        for(each.sample in 1:n.rows){
          
          if( agg[each.sample, index]*1.5 < my.mean.fluo ){
            to.remove=c(to.remove, agg$pic[each.sample])
          } else {
            break;
          }                    
        }
        
        if(length(to.remove)>0){
          data[[K]][[L]] = data[[K]][[L]][ ! data[[K]][[L]]$pic %in% to.remove,]
          warning(paste("In plate ",K,"and well ",L, " the picture(s) ",paste(to.remove,",")," was/were removed"))
        }
      }
    }           
  }
  return(data)
}


## Returns the S number of pictures corresponding to a particular well and picture.
## If 0 is given for the pic number, S numbers are returned for all pictures of the well.
##
uscope.get.image.index = function(plate=1, well, design, pic=0){

    if( typeof(well) == "double"){
        well = design$PDEF$well[well]
    }

    if(pic == 0){
        s.positions = which(design$WELLS[[plate]] == well)
    } else {
        s.positions = which(design$WELLS[[plate]] == well & design$PICNUM[[plate]] == pic)
    }

    s.pos.res = rbind( as.numeric(design$PICNUM[[plate]][s.positions]), as.numeric(design$S.positions[[plate]][s.positions]))
    rownames(s.pos.res) = c("PIC no.", "S # ")
    
    return( s.pos.res )

}

uscope.process.get.means = function(data, design){
  
  n.samp = length(data)
  mat.res = matrix(nrow=NROW(design$PDEF), ncol=3*n.samp)
  colnames(mat.res)= paste(rep(c("GFP","RFP","cor"),n.samp),rep(c(1:n.samp),rep(3, n.samp)), sep=".")
  
  for(K in 1:length(data)){
    
    wells = names(data[[K]])
    
    for(L in 1:length(data[[K]])){
      
      if(length(data[[K]][[L]]$GFPavgIntBin7)>100){
        mean.GFP = median(data[[K]][[L]]$GFPavgIntBin7, na.rm=TRUE)
        mean.RFP = median(data[[K]][[L]]$RFPavgIntBin7, na.rm=TRUE)
        my.cor   = cor(data[[K]][[L]]$GFPavgIntBin7, data[[K]][[L]]$RFPavgIntBin7, met="pearson")
      } else {
        mean.GFP = NA
        mean.RFP = NA
        my.cor   = NA
      }
      K2 = 1+((K-1)*3)
      
      well.POS = which(design$PDEF$well==wells[L])
      mat.res[well.POS,K2] = round(mean.GFP)
      mat.res[well.POS,K2+1] = round(mean.RFP)
      mat.res[well.POS,K2+2] = round(my.cor,3)
      
      #design$PDEF[ which(design$PDEF$well==wells[L]),"GFP.mean"]=
      #design$PDEF[ which(design$PDEF$well==wells[L]),"RFP.mean"]=mean.RFP
    }
  }
  mat.res = data.frame(mat.res)
  return(mat.res)
}

uscope.process.get.one.stat.allp = function(design, data, col.of.interest, fun2use=mean, info){

    my.res = uscope.process.get.one.stat(design = design, data = data, plateNUM=1, col.of.interest=col.of.interest, fun2use=fun2use, info=info)

    print(dim(my.res))
    
    if(length(data)>1){

        for(K in 2:length(data)){
            my.res = rbind( my.res, uscope.process.get.one.stat(design = design, data = data, plateNUM=K, col.of.interest=col.of.interest, fun2use=fun2use, info=info) )
            print(dim(my.res))
        }
    }
    return(my.res)
}


uscope.process.get.one.stat = function(design, data, plateNUM, col.of.interest, fun2use=mean, info){
  
  
    my.col = grep(col.of.interest, colnames(data[[plateNUM]][[1]]))
  
    if( length(my.col)==0 | length(my.col) > 1){
        if(length(my.col)==0){
            stop(paste("no column was found to match the chanel name", my.col))
        } else {
            stop(paste("more than one column was found to match the chanel name", my.col))
        }
    }       
  
    res = sapply(data[[plateNUM]], function(x){ return( fun2use(x[,my.col]) ) } )
    constructs = design$PDEF.split[[plateNUM]][, info]
    
    my.RES = data.frame( mean = res, id = constructs, stringsAsFactors = FALSE, row.names = NULL)
    ##print(dim(my.RES))
    return(my.RES)
}

uscope.process.get.cellStat = function(design, data, plateNUM, col.of.interest, info){
  
  my.col = grep(col.of.interest, colnames(data[[plateNUM]][[1]]))
  
  if( length(my.col)==0 | length(my.col) > 1){
    if(length(my.col)==0){
      stop(paste("no column was found to match the chanel name", my.col))
    } else {
      stop(paste("more than one column was found to match the chanel name", my.col))
    }
  }       
  
  res = sapply(data[[plateNUM]], function(x){ return( x[,my.col] ) } )
  constructs = design$PDEF.split[[plateNUM]][, info]
  names(res) = constructs
  return(res)
}

uscope.process.get.one.mean = function(design, data, plateNUM){
  
  my.RES = list()
  
  for ( each.CH in design$CHANELS){
    
    my.col = grep(paste(each.CH,design$MEASURE.COL,"$",sep=""), colnames(data[[plateNUM]][[1]]))
    
    if( length(my.col)==0 | length(my.col) > 1){
      if(length(my.col)==0){
        stop(paste("no column was found to match the chanel name", each.CH))
      } else {
        stop(paste("more than one column was found to match the chanel name", each.CH))
      }
    }        
    my.RES[[each.CH]]$mean = sapply(data[[plateNUM]], function(x){ return( mean(x[,my.col], na.rm=TRUE) ) } )
    my.RES[[each.CH]]$sd   = sapply(data[[plateNUM]], function(x){ return( sd(x[,my.col], na.rm=TRUE) ) } )
    my.RES[[each.CH]]$n    = sapply(data[[plateNUM]], function(x){ return( sum(x[,my.col]>0, na.rm=TRUE) ) } )
  }
  #CH = paste(paste(rep(design$CHANELS, rep(length(design$CHANELS), length(design$CHANELS)))), c("",".sd"),sep="")
  #colnames(my.RES) = CH
  return(my.RES)
}

get.color.id = function(ch.name){
  
  if( length(grep("RFP",ch.name))>0){
    return("dark red")
  }
  if( length(grep("GFP",ch.name))>0){
    return("dark green")
  }
  if( length(grep("BFP",ch.name))>0){
    return("dark blue")
  }
}


### presnting correlations

cor.show = function(prom1="GPD",prom2="GPD",prot1="no",prot2="no",plate=2,design,data.means){
  
  prom.ind1 = intersect(grep(prom1,design$PDEF$GFP),grep(prom2,design$PDEF$RFP))
  prom.ind2 = intersect(grep(prom2,design$PDEF$GFP),grep(prom1,design$PDEF$RFP))
  if (prom1==prom2){
    prom.ind=prom.ind1
  }else{
    prom.ind = c(prom.ind1,prom.ind2)
  }
  
  prot.ind1 = intersect(grep(prot1,design$PDEF$GFP),grep(prot2,design$PDEF$RFP))
  prot.ind2 = intersect(grep(prot1,design$PDEF$RFP),grep(prot2,design$PDEF$GFP))
  if (prot1==prot2){
    prot.ind = prot.ind1
  }else{
    prot.ind = c(prot.ind1,prot.ind2)
  }
  all.ind = intersect(prom.ind,prot.ind)
  if (is.numeric(plate)){
    col = paste("cor",plate,sep=".")
  }else{
    col = paste("cor",which(design$D == plate),sep=".")
  }
  out = cbind(design$PDEF[all.ind,],data.means[all.ind,col])
  print(paste("correlations for plate ",plate,sep=": "))
  print(out)
  
}

### plot all reps of the same variants in different colors

plot.reps = function(YFPpr="GPD",CHpr="GPD",YFPad="no",CHad="no",plate=2,design,data){
  
  YFP.ind = intersect(grep(pattern=YFPpr,x=design$PDEF$GFP),grep(pattern=YFPad,x=design$PDEF$GFP))
  CH.ind = intersect(grep(pattern=CHpr,x=design$PDEF$RFP),grep(pattern=CHad,x=design$PDEF$RFP))
  ind = intersect(YFP.ind,CH.ind)
  
  plot(data[[plate]][[ind[1]]]$GFPavgIntBin7,data[[plate]][[ind[1]]]$RFPavgIntBin7,main=paste("plate:",plate,"construct:",design$PDEF$GFP[ind[1]],design$PDEF$RFP[ind[1]],sep=" "),cex=0.7)
  for (i in 2:length(ind)){
    if (i < (length(ind)/2+1)){
      points(x=data[[plate]][[ind[i]]]$GFPavgIntBin7,y=data[[plate]][[ind[i]]]$RFPavgIntBin7,col=i,cex=0.7)
    }else{
      points(x=data[[plate]][[ind[i]]]$GFPavgIntBin7,y=data[[plate]][[ind[i]]]$RFPavgIntBin7,col=i-length(ind)/2,pch=19,cex=0.7)
    }
  }
  
}

####
#### Plots any property on the plates to help detect plate-specific effects.
#### TODO: colorscale should be added but that will require using a different layout for the results to make space on one side
####       alternatively, colorscale could be written in a different file -- that's probably the easiest.
####
diagnostic.intensity = function(data, design, col.of.interest="GFP_int_b9", fun2use=mean, cols.to.get=c("GFP","Row","Column"), K=1 , use.quantile=TRUE){
    if(length(data)==1){
        par(mfrow=c(1,1), mai=c(1,1,1,1))
    }
    if(length(data)==2){
        par(mfrow=c(1,2), mai=c(1,1,1,1))
    }
    if(length(data) > 2 & length(data) <= 4){
        par(mfrow=c(2,2), mai=c(1,1,1,1))
    }
    if(length(data) > 4 & length(data) <= 6){
        par(mfrow=c(2,3), mai=c(1,1,1,1))
    }
    if(length(data) > 6 & length(data) <= 9){
        par(mfrow=c(3,3), mai=c(1,1,1,1))
    }
    if(length(data) > 9 & length(data) <= 12){
        par(mfrow=c(3,4), mai=c(0.2,0.5,0.5,0.15))
    }
    if(length(data) > 12 & length(data) <= 16){
        par(mfrow=c(4,4), mai=c(0.2,0.5,0.5,0.15))
    }
    if(length(data) > 16 & length(data) <= 20){
        par(mfrow=c(4,5), mai=c(0.2,0.5,0.5,0.15))
    }
    if(length(data) > 20 & length(data) <= 25){
        par(mfrow=c(5,5), mai=c(0.2,0.5,0.5,0.15))
    }

    check.range = c()
    for( K in 1:length(data)){
        check.range = c(check.range, uscope.process.get.one.stat(design = design, data = data, plateNUM=K, col.of.interest=col.of.interest, info=cols.to.get, fun2use=fun2use)$mean)
    }

    print(check.range[1:100])
    check.range = check.range[which(is.finite(check.range))]
    print(check.range[1:100])

    if(use.quantile==TRUE){
        my.range = quantile(check.range, prob=seq(0,1,len=21), na.rm=TRUE)
        print(my.range)
        my.COL = c(rgb(0.6,0.1,0.1), B2O.COL(20), rgb(1,1,1))     
    } else {
        my.range = seq(min(check.range,na.rm=TRUE), max(check.range,na.rm=TRUE),len=21)
        print(my.range)
        my.COL = c(rgb(0.6,0.1,0.1), B2O.COL(20), rgb(1,1,1))     
    }
    
    for( K in 1:length(data)){
        res.tmp = uscope.process.get.one.stat(design=design, data=data, plateNUM=K, col.of.interest=col.of.interest, info=cols.to.get, fun2use=fun2use)
        res.tmp$mean[is.na(res.tmp$mean)]=0
        to.plot = matrix(res.tmp$mean, byrow=TRUE, ncol=24, nrow=16)
        if(design$FORMAT == 96){
            image(rotate.matrix(to.plot), axes=FALSE, main=paste("Plate  ",K), zlim=c(0,10))
            axis(2, at=seq(0,1,length=8), lab=LETTERS[8:1])
            axis(3, at=seq(0,1,length=12), lab=1:12)
        }
        if(design$FORMAT == 384){
            image(rotate.matrix(to.plot), axes=FALSE, main=paste("Plate  ",K), breaks=c(-2^16, my.range, +2^16) , col=my.COL )
            axis(2, at=seq(0,1,length=16), lab=LETTERS[16:1])
            axis(3, at=seq(0,1,length=24), lab=1:24)
        }
    }
}



###
### finding the noise and noise component using formulas in the supplementarry of Elowitz et al. science 2002
###
noise.by.elowitz = function(data,design,GFP="norm.GFP",RFP="norm.RFP"){
  
  #formulas in Elowitz et al.
  #n.tot.sq = (mean(RFP^2+GFP^2) - 2*mean(RFP)*mean(GFP)) / 2*mean(RFP)*mean(GFP)
  #n.ext.sq = (mean(RFP*GFP) - mean(RFP)*mean(GFP)) / mean(RFP)*mean(GFP)
  #n.int.sq = mean((RFP-GFP)^2) / 2*mean(RFP)*mean(GFP)
  
  RFPcol = RFP
  GFPcol = GFP
  noise.mat = c()
  for (plate in design$D){
    noise.plate = matrix(NA,design$FORMAT,3)
    for (well in 1:design$FORMAT){
      AvSqAd = mean(data[[plate]][[well]][,RFPcol]^2 + data[[plate]][[well]][,GFPcol]^2)  # mean of addition of the squares
      MultAvs = mean(data[[plate]][[well]][,RFPcol]) * mean(data[[plate]][[well]][,GFPcol])  # multiplication of the means
      AvMult = mean(data[[plate]][[well]][,RFPcol] * data[[plate]][[well]][,GFPcol])   # mean of the multiplication 
      AvSubSq = mean((data[[plate]][[well]][,RFPcol] - data[[plate]][[well]][,GFPcol])^2)  # mean of the substruction
      #n.tot.sq = (AvSqAd - 2*MultAvs) / 2*MultAvs
      #n.ext.sq = (AvMult - MultAvs) / MultAvs
      #n.int.sq = AvSubSq / 2*AvMult
      n.tot = sqrt((AvSqAd - 2*MultAvs) / (2*MultAvs))
      n.ext = sqrt((AvMult - MultAvs) / MultAvs)
      n.int = sqrt(AvSubSq / (2*AvMult))
      noise.plate[well,] = c(n.tot,n.ext,n.int)
    }
    colnames(noise.plate) = paste(c("n.tot","n.ext","n.int"),plate,sep=".")
    noise.mat = cbind(noise.mat,noise.plate)
  } 
  return(cbind(design$PDEF,noise.mat))
}

###
### find noise by emmanuel's idea
###

noise.by.emmanuel = function(data,design,GFP="norm.GFP",RFP="norm.RFP"){
  
  noise.mat = c()
  for (plate in design$D){
    noise.plate = matrix(NA,design$FORMAT,8)  # the matrix has 3*2+1+1 columns tot,ext,int for each color. cells number and cor
    for (well in 1:design$FORMAT){
      cells.num = length(data[[plate]][[well]][,GFP])
      if (cells.num>50){
        n.tot.G = var(data[[plate]][[well]][,GFP]) / mean(data[[plate]][[well]][,GFP])^2 # in the case of the normalized data
        n.tot.R = var(data[[plate]][[well]][,RFP]) / mean(data[[plate]][[well]][,RFP])^2 # mean=1 but I wanted to be able to 
        # work with other data as well
        cor.GR = cor(data[[plate]][[well]][,GFP],data[[plate]][[well]][,RFP])
        reg.G = lm(data[[plate]][[well]][,GFP]~data[[plate]][[well]][,RFP])
        reg.R = lm(data[[plate]][[well]][,RFP]~data[[plate]][[well]][,GFP])
        n.ext.G = cor.GR^2 * n.tot.G
        n.ext.R = cor.GR^2 * n.tot.R
        n.int.G = var(reg.G$residuals) / mean(data[[plate]][[well]][,GFP])^2
        n.int.R = var(reg.R$residuals) / mean(data[[plate]][[well]][,RFP])^2
        
        noise.plate[well,] = c(n.tot.G,n.ext.G,n.int.G,n.tot.R,n.ext.R,n.int.R,cells.num,cor.GR)
        noise.plate[well,] = round(noise.plate[well,],3)
      }else{
        noise.plate[well,] = rep(NA,8)
      }
    }
    colnames(noise.plate) = paste(c("n.tot.G","n.ext.G","n.int.G","n.tot.R","n.ext.R","n.int.R","cells.num","correlation"),plate,sep=".")
    noise.mat = cbind(noise.mat,noise.plate)
  } 
  return(cbind(design$PDEF,noise.mat))
}

###
### make table of the precentage of intrinsic and extrinsic noise from the total noise
###

noise.contribution = function(data,design,GFP="norm.GFP",RFP="norm.RFP"){
  
  noise.mat = c()
  for (plate in design$D){
    noise.plate = matrix(NA,design$FORMAT,6)  # the matrix has 2+2+1+1 columns tot,ext,int for each color. cells number and cor
    for (well in 1:design$FORMAT){
      cells.num = length(data[[plate]][[well]][,GFP])
      if (cells.num>50){
        n.tot.G = var(data[[plate]][[well]][,GFP]) / mean(data[[plate]][[well]][,GFP])^2 # in the case of the normalized data
        n.tot.R = var(data[[plate]][[well]][,RFP]) / mean(data[[plate]][[well]][,RFP])^2 # mean=1 but I wanted to be able to 
        # work with other data as well
        cor.GR = cor(data[[plate]][[well]][,GFP],data[[plate]][[well]][,RFP])
        reg.G = lm(data[[plate]][[well]][,GFP]~data[[plate]][[well]][,RFP])
        reg.R = lm(data[[plate]][[well]][,RFP]~data[[plate]][[well]][,GFP])
        n.ext.G = cor.GR^2 * n.tot.G
        n.ext.R = cor.GR^2 * n.tot.R
        n.int.G = var(reg.G$residuals) / mean(data[[plate]][[well]][,GFP])^2
        n.int.R = var(reg.R$residuals) / mean(data[[plate]][[well]][,RFP])^2
        ext.cont.G = 100*(n.ext.G / n.tot.G)
        int.cont.G = 100*(n.int.G / n.tot.G)
        ext.cont.R = 100*(n.ext.R / n.tot.R)
        int.cont.R = 100*(n.int.R / n.tot.R)
        noise.plate[well,] = c(n.tot.G,n.tot.R,ext.cont.R,int.cont.R,cells.num,cor.GR)
        noise.plate[well,] = round(noise.plate[well,],3)
      }else{
        noise.plate[well,] = rep(NA,6)
      }
    }
    colnames(noise.plate) = paste(c("n.tot.G","n.tot.R","ext.cont","int.cont","cells.num","correlation"),plate,sep=".")
    noise.mat = cbind(noise.mat,noise.plate)
  } 
  return(cbind(design$PDEF,noise.mat))
}

###
### making boxplots to present the data
###

make.boxplot = function(data,design,plate,proms = c("TEF","GPD"),prots = c("no","D16","DH"),marks = c("optY","optCh"), PDF){
  
  print(paste("Now doing plate ",plate,"..."))
  pdf(PDF, width=8, height=8)
  par(mfcol=c(1,1), mai=c(3,0.5,0.5,0.5))
  list.res = list()
  
  ## same prom, same prot
  PROMS1 = rep(proms, c(length(prots),length(prots)))
  PROMS2 = PROMS1
  PROTS1 = rep(prots, length(proms))
  PROTS2 = PROTS1
  
  print("Same promoter, same protein");    print(paste("Promoters1 = ",paste(PROMS1,collapse="; ")));    print(paste("Promoters2 = ",paste(PROMS2,collapse="; ")));    print(paste("Proteins 1 = ",paste(PROTS1,collapse="; ")));    print(paste("Proteins 2 = ",paste(PROTS2,collapse="; ")))
  
  K=1
  col.border=c()
  col.box   =c()
  
  ## ADD data to list.res ...
  for(N in 1:length(PROMS1)){
    construct1 = paste(PROMS1[N],"-optY-",PROTS1[N],sep="")
    construct2 = paste(PROMS2[N],"-optCh-",PROTS2[N],sep="")
    
    indexes = which(grepl(construct1, design$PDEF$GFP) & grepl(construct2,design$PDEF$RFP))
    
    list.res[[K]] = cor(data[[plate]][[ indexes[1] ]]$GFPavgIntBin7,data[[plate]][[ indexes[1] ]]$RFPavgIntBin7)
    for(each.index in indexes[-c(1)]){
      list.res[[K]] = c(list.res[[K]], cor(data[[plate]][[each.index]]$GFPavgIntBin7,data[[plate]][[each.index]]$RFPavgIntBin7))
    }
    names(list.res)[[K]] = paste(construct1,"  vs ",construct2,sep="")
    K=K+1
    col.border=c(col.border,1)
    col.box   = c(col.box,"dark red")
  }
  
  ## same prom, diff prot
  ## makes all prots combinations
  if(length(prots)==2){
    tmp = expand.grid(0:1, 0:1)
  }
  if(length(prots)==3){
    tmp = expand.grid(0:1, 0:1, 0:1)
  }
  if(length(prots)==4){
    tmp = expand.grid(0:1, 0:1, 0:1, 0:1)
  }
  tmp = tmp[rowSums(tmp)==2,]
  PROTS1 = c()
  PROTS2 = c()
  for(i in 1:nrow(tmp)){
    PROTS1 = c(PROTS1, prots[which(tmp[i,]==1)[1]])
    PROTS2 = c(PROTS2, prots[which(tmp[i,]==1)[2]])
  }
  PROMS1 = c(rep(proms[1],length(PROTS1)),rep(proms[2],length(PROTS1)))
  PROMS2 = PROMS1
  PROTS1 = rep(PROTS1,2)
  PROTS2 = rep(PROTS2,2)
  
  print("Same promoter, DIFF protein");    print(paste("Promoters1 = ",paste(PROMS1,collapse="; ")));    print(paste("Promoters2 = ",paste(PROMS2,collapse="; ")));    print(paste("Proteins 1 = ",paste(PROTS1,collapse="; ")));    print(paste("Proteins 2 = ",paste(PROTS2,collapse="; ")))
  ## ADD data to list.res ...
  for(N in 1:length(PROMS1)){
    construct1.1 = paste(PROMS1[N],"-optY-",PROTS1[N],sep="")
    construct1.2 = paste(PROMS2[N],"-optCh-",PROTS2[N],sep="")            
    indexes1 = which(grepl(construct1.1, design$PDEF$GFP) & grepl(construct1.2,design$PDEF$RFP))
    
    construct2.1 = paste(PROMS1[N],"-optCh-",PROTS1[N],sep="")
    construct2.2 = paste(PROMS2[N],"-optY-",PROTS2[N],sep="")            
    indexes2 = which(grepl(construct2.1, design$PDEF$RFP) & grepl(construct2.2,design$PDEF$GFP))
    
    list.res[[K]]   = cor(data[[plate]][[ indexes1[1] ]]$GFPavgIntBin7,data[[plate]][[ indexes1[1] ]]$RFPavgIntBin7)
    list.res[[K+1]] = cor(data[[plate]][[ indexes2[1] ]]$GFPavgIntBin7,data[[plate]][[ indexes2[1] ]]$RFPavgIntBin7)
    for( each.i in 2:length(indexes1) ){
      list.res[[K]] = c(list.res[[K]], cor(data[[plate]][[ indexes1[each.i] ]]$GFPavgIntBin7,data[[plate]][[ indexes1[each.i] ]]$RFPavgIntBin7))
      list.res[[K+1]] = c(list.res[[K+1]], cor(data[[plate]][[ indexes2[each.i] ]]$GFPavgIntBin7,data[[plate]][[ indexes2[each.i]]]$RFPavgIntBin7))            
    }
    names(list.res)[[K]] = paste(construct1.1,"  vs ",construct1.2,sep="")
    names(list.res)[[K+1]] = paste(construct2.2,"  vs ",construct2.1,sep="")
    col.border=c(col.border,rep(col.border[length(col.border)]+1,2))
    col.box   = c(col.box, "dark blue", "dark blue")
    K=K+2
  }  
  
  ## diff prom, same prot
  PROTS1 = prots
  PROTS2 = prots
  PROMS1 = rep(proms[1], length(prots))
  PROMS2 = rep(proms[2], length(prots))
  
  print("DIFF promoter, same protein");    print(paste("Promoters1 = ",paste(PROMS1,collapse="; ")));    print(paste("Promoters2 = ",paste(PROMS2,collapse="; ")));    print(paste("Proteins 1 = ",paste(PROTS1,collapse="; ")));    print(paste("Proteins 2 = ",paste(PROTS2,collapse="; ")))
  ## ADD data to list.res ...
  for(N in 1:length(PROMS1)){
    construct1.1 = paste(PROMS1[N],"-optY-",PROTS1[N],sep="")
    construct1.2 = paste(PROMS2[N],"-optCh-",PROTS2[N],sep="")            
    indexes1 = which(grepl(construct1.1, design$PDEF$GFP) & grepl(construct1.2,design$PDEF$RFP))
    
    construct2.1 = paste(PROMS1[N],"-optCh-",PROTS1[N],sep="")
    construct2.2 = paste(PROMS2[N],"-optY-",PROTS2[N],sep="")            
    indexes2 = which(grepl(construct2.1, design$PDEF$RFP) & grepl(construct2.2,design$PDEF$GFP))
    
    list.res[[K]]   = cor(data[[plate]][[ indexes1[1] ]]$GFPavgIntBin7,data[[plate]][[ indexes1[1] ]]$RFPavgIntBin7)
    list.res[[K+1]] = cor(data[[plate]][[ indexes2[1] ]]$GFPavgIntBin7,data[[plate]][[ indexes2[1] ]]$RFPavgIntBin7)
    for( each.i in 2:length(indexes1) ){
      list.res[[K]] = c(list.res[[K]], cor(data[[plate]][[ indexes1[each.i] ]]$GFPavgIntBin7,data[[plate]][[ indexes1[each.i] ]]$RFPavgIntBin7))
      list.res[[K+1]] = c(list.res[[K+1]], cor(data[[plate]][[ indexes2[each.i] ]]$GFPavgIntBin7,data[[plate]][[ indexes2[each.i]]]$RFPavgIntBin7))            
    }
    names(list.res)[[K]] = paste(construct1.1,"  vs ",construct1.2,sep="")
    names(list.res)[[K+1]] = paste(construct2.2,"  vs ",construct2.1,sep="")
    col.border=c(col.border,rep(col.border[length(col.border)]+1,2))
    col.box   = c(col.box, "orange", "orange")
    K=K+2
  }  
  
  ## diff prom, diff prot
  PROTS1 = c()
  PROTS2 = c()
  for(i in 1:nrow(tmp)){
    PROTS1 = c(PROTS1, prots[which(tmp[i,]==1)[1]])
    PROTS2 = c(PROTS2, prots[which(tmp[i,]==1)[2]])
  }
  for(i in 1:nrow(tmp)){
    PROTS2 = c(PROTS2, prots[which(tmp[i,]==1)[1]])
    PROTS1 = c(PROTS1, prots[which(tmp[i,]==1)[2]])
  }
  PROMS1 = rep(proms[1], length(PROTS1))
  PROMS2 = rep(proms[2], length(PROTS2))    
  
  print("DIFF promoter, DIFF protein");    print(paste("Promoters1 = ",paste(PROMS1,collapse="; ")));    print(paste("Promoters2 = ",paste(PROMS2,collapse="; ")));    print(paste("Proteins 1 = ",paste(PROTS1,collapse="; ")));    print(paste("Proteins 2 = ",paste(PROTS2,collapse="; ")))
  ## ADD data to list.res ...
  for(N in 1:length(PROMS1)){
    construct1.1 = paste(PROMS1[N],"-optY-",PROTS1[N],sep="")
    construct1.2 = paste(PROMS2[N],"-optCh-",PROTS2[N],sep="")            
    indexes1 = which(grepl(construct1.1, design$PDEF$GFP) & grepl(construct1.2,design$PDEF$RFP))
    
    construct2.1 = paste(PROMS1[N],"-optCh-",PROTS1[N],sep="")
    construct2.2 = paste(PROMS2[N],"-optY-",PROTS2[N],sep="")            
    indexes2 = which(grepl(construct2.1, design$PDEF$RFP) & grepl(construct2.2,design$PDEF$GFP))
    
    list.res[[K]]   = cor(data[[plate]][[ indexes1[1] ]]$GFPavgIntBin7,data[[plate]][[ indexes1[1] ]]$RFPavgIntBin7)
    list.res[[K+1]] = cor(data[[plate]][[ indexes2[1] ]]$GFPavgIntBin7,data[[plate]][[ indexes2[1] ]]$RFPavgIntBin7)
    for( each.i in 2:length(indexes1) ){
      list.res[[K]] = c(list.res[[K]], cor(data[[plate]][[ indexes1[each.i] ]]$GFPavgIntBin7,data[[plate]][[ indexes1[each.i] ]]$RFPavgIntBin7))
      list.res[[K+1]] = c(list.res[[K+1]], cor(data[[plate]][[ indexes2[each.i] ]]$GFPavgIntBin7,data[[plate]][[ indexes2[each.i]]]$RFPavgIntBin7))            
    }
    names(list.res)[[K]] = paste(construct1.1,"  vs ",construct1.2,sep="")
    names(list.res)[[K+1]] = paste(construct2.2,"  vs ",construct2.1,sep="")
    col.border=c(col.border,rep(col.border[length(col.border)]+1,2))
    col.box   = c(col.box, "dark green", "dark green")
    K=K+2
  }  
  
  boxplot(list.res, las=2, border=col.border, col=col.box, ylim=c(-0.4,1))
  abline(v=1:length(list.res), lty=1,col="light grey")
  abline(h=pretty(c(0,1)), lty=1,col="grey")
  dev.off()    
}

#
###
### assisting the function boxplot
###

one.boxplot = function(PROMS1,PROMS2,PROTS1,PROTS2,data){ # data is noise table 
  
  ROWnames1 = paste(PROMS1,marks[1],PROTS1,sep="-")
  ROWnames2 = paste(PROMS2,marks[2],PROTS2,sep="-")
  ROWnames.total = (paste(ROWnames1,ROWnames2,sep="_"))
  
  indices = as.vector(sapply(X=ROWnames.noHead,FUN = function(x) {which(x==data$GFP.RFP)}))
  part.data = data[indices,]
  part.data = cbind(ROWnames.total,part.data)
  par(mar=c(11,3,3,3))
  boxplot(correlation.YMD~ROWnames.total,data=part.data,las=2,cex.axis=0.7,main="correlation boxplot")
  
}

find.ind = function(PROM1,PROM2,PROT1,PROT2,data){
  ROWname1 = c(paste(PROM1,"optY",PROT1,sep="-"))
  ROWname2 = c(paste(PROM2,"optCh",PROT2,sep="-"))
  ROWname.total = (paste(ROWname1,ROWname2,sep="_"))
  indices = which(ROWname.total==data$GFP.RFP)
  return(indices)
}

###
### Make a boxplot with the mean values for each protein in all contexts
###
boxplot.means = function(data, design){
  
  data.means = uscope.process.get.means(data, design)    
  
  pdf("/media/elusers/users/emmanuel/GFP2-2.pdf", width=3*length(design$D), height=5)
  par(mfcol=c(1,length(design$D)), mai=c(1.2,0.1,0.3,0.1))
  K=1
  for(each.design in design$D){
    data.GFP = split(data.means[,paste("GFP.",K, sep="")], design$PDEF$GFP)
    boxplot(data.GFP, las=2, col="green", main=each.design)
    K=K+1
  }
  dev.off()
  
  pdf("/media/elusers/users/emmanuel/RFP2-2.pdf", width=3*length(design$D), height=5)
  par(mfcol=c(1,length(design$D)), mai=c(1.2,0.1,0.3,0.1))
  K=1
  for(each.design in design$D){
    data.GFP = split(data.means[,paste("RFP.",K, sep="")], design$PDEF$RFP)
    boxplot(data.GFP, las=2, col="red",main=each.design)
    K=K+1
  }
  dev.off()
}


###
### make boxplot for noise. segment by GFP or RFP only
###

make.box.simple = function(proms = c("TEF","GPD"),prots = c("no","cODC","D16","DH"),
                           plate="YMD",param="correlation",flProt="optY",data=noiseVal.table.emmanuel,PDF){
  list.res = list()
  col.arg = c()
  prot.num = length(prots)
  prom.num = length(proms)
  for (PROM in 1:prom.num){
    for (PROT in 1:prot.num){
      ind = grep(x=data$GFP.RFP,pattern=paste(proms[PROM],flProt,prots[PROT],sep="-"))
      list.res[[(PROM-1)*prot.num+PROT]] = data[ind,paste(param,plate,sep=".")]
      names(list.res)[(PROM-1)*prot.num+PROT] = paste(proms[PROM],flProt,prots[PROT],sep="-")
      col.arg[(PROM-1)*prot.num+PROT] = (PROM-1)*prot.num+PROT
    }
  }
  pdf(PDF)
  boxplot(list.res,names=names(list.res),las=2,cex.axis=0.7,main=paste(param,flProt,sep="."),col=col.arg)
  dev.off()
}


get.well = function(PROM1, PROT1, PROM2, PROT2, design, plate=1){
  construct1 = paste(PROM1,"-optY-", PROT1,sep="")
  construct2 = paste(PROM2,"-optCh-",PROT2,sep="")
  indexes = which(grepl(construct1, design$PDEF$GFP) & grepl(construct2,design$PDEF$RFP))
  wells = design$PDEF$well[indexes]
  s.positions = design$S.positions[[plate]][which(design$WELLS[[plate]] %in% wells)]
  return(list( DEF=design$PDEF[indexes,], pos=s.positions) )
}


get.sNUM = function(design, plate=1, well){
    positions = design$WELLS[[plate]] %in% c(well)
    s.positions = design$S.positions[[plate]][positions]
    return(as.numeric(s.positions))
}


###
### normalize the data
###

uscope.process.norm = function(data,colG,colR){
  
  norm.data = vector("list",3)
  for (plate in 1:length(data)){
    norm.data[[plate]] = vector("list",length(data[[plate]]))
    for (well in 1:length(data[[plate]])){
      norm.GFP = data[[plate]][[well]][,colG] / mean(data[[plate]][[well]][,colG])
      norm.RFP = data[[plate]][[well]][,colR] / mean(data[[plate]][[well]][,colR])
      norm.data[[plate]][[well]] = cbind(data[[plate]][[well]],norm.GFP,norm.RFP)
    }
  }
  return(norm.data)
  
  
}

## 
## Takes:
## - a directory path where a screen is continued ()
## - start position (e.g., C01 or E01) for the new set of images
## - the screen name (the part of the image names that's before the underscore)
## - PPW, the number of picture per well for that screen
## - DOIT = 1 --> must be set to 1 to actually do the renaming, if not it will just print how things will be renamed.
##
## --> it will rename all the file in the current directory with the correct sNUM
## 
## AFTER RENAMING, the files can be moved to the main directory.
##

#uscope.merge.screen.rename.sNum(FOLDER="/Volumes/elmicro/meta/phase/Conditions_phase_diagram/KO_library_with_P5O6/Ghil/full_screen/plate_9/ImagesM1", START="M1", SCREEN_NAME = "image1", PPW=12, DOIT=1)
#uscope.merge.screen.rename.sNum(FOLDER="/Volumes/elmicro/meta/phase/Conditions_phase_diagram/KO_library_with_P5O6/Ghil/full_screen/plate_12/ImagesG1", START="G1", SCREEN_NAME = "images1", PPW=12, DOIT=0)


#uscope.merge.screen.rename.sNum(FOLDER="/Volumes/elmicro/meta/phase/Conditions_phase_diagram/KO_library_with_P5O6/Ghil/full_screen/plate_12/ImagesE1", START="E1", SCREEN_NAME = "images1", PPW=12, DOIT=1)

#uscope.merge.screen.rename.sNum(FOLDER="/Volumes/elmicro/meta/phase/Conditions_phase_diagram/KO_library_with_P5O6/Ghil/full_screen/plate_14/ImagesG1", START="G1", SCREEN_NAME = "image1", PPW=12, DOIT=0)

uscope.merge.screen.rename.sNum = function(FOLDER="/Volumes/elmicro/meta/phase/Conditions_phase_diagram/KO_library_with_P5O6/Ghil/full_screen/plate_13/ImagesE1", 
                                           START="E1",
                                           SCREEN_NAME = "image1",
                                           PPW=12,
                                           DOIT=0) {
  print(paste0("Please make sure to run this function with 'DOIT=0' at first, and check the output by hand in the file ~/merge_tmp1.txt"))
  print(paste0(""))
  print(paste0("If you are happy with the renaming described in that file, proceed to run this function with DOIT=1"))
  print(paste0(""))
  if(grep(substring(START,1,1), c("C","E","G","I","K","M","O"))){
    
      all.files = list.files(path=FOLDER)

      print(paste0("Now processing folder",FOLDER, " and ", length(all.files), " files were found"))
      
      if( length(grep(paste0("^",SCREEN_NAME,"[_s]"), all.files[1]))==0 ) {
        print(paste0("Looks like you gave a wrong screen name: ",SCREEN_NAME, " doesn't match ",all.files[1] ))
        return("Exiting")
      }
      
      s.number.order = c()
      
      for( each.file in all.files){
        s.number = as.numeric(sub(pattern=paste0("^(",SCREEN_NAME,".*_s)([0-9]+)([\\.Ri\\-].*)$"), replacement="\\2", x=each.file, perl=TRUE))
        if(is.na(s.number)){
          s.number.order = c(s.number.order,0)
        } else {
          s.number.order = c(s.number.order,s.number)
        }
      }
      
      all.files = all.files[order(s.number.order, decreasing=TRUE)]
      
      WELLS = get.plate.wells(PPW)
      S.first = grep(paste0(START,".1$"), WELLS)
      print(paste("This folder contains images starting at well ",START))
      print(paste("The s number corresponding to this position should have been",S.first))
      print(paste("The s number of all the files will be changed and ",(S.first-1)," will be added"))
      msg = readline(prompt="Are you OK with this? Y/N ")
      if(msg == "Y"){
        
        sink(file="~/merge_tmp1.txt")
        
        for( each.file in all.files){
          s.number = as.numeric(sub(pattern=paste0("^(",SCREEN_NAME,".*_s)([0-9]+)([\\.Ri\\-].*)$"), replacement="\\2", x=each.file, perl=TRUE))
          if(!is.na(s.number)){
            new.num = s.number + S.first-1
            new.file = sub(pattern=paste0("^(",SCREEN_NAME,".*_s)([0-9]+)([\\.Ri\\-].*)$"), replacement=paste0("\\1",new.num,"\\3"), x=each.file, perl=TRUE)
            
            print(paste0("RENAMING:",each.file, " ==> ", new.file ))
            
            if(DOIT){
              file.rename( paste0(FOLDER,"/",each.file), to= paste0(FOLDER, "/", new.file))
            }
          }
        }
        sink()
      }
  } else {
    print("A new screen must start at odd-numbered rows otherwise it cannot be merged")
  }
}

## 
## Takes:
## - A directory where all images will be copied (typically that's the "Images" or "Images_output" folder)
## - A directory where a screen is continued -- YOU MUST HAVE APPLIED uscope.merge.screen.rename.sNum on that directory
## - The screen name (the part of the image names that's before the underscore)
## - DOIT = 1 --> must be set to 1 to actually do the renaming, if not it will just print how things will be renamed.
##
## --> it will MOVE all the files to the master Images directory. If images already exist with a particular sNUM, they will be moved to a BACKUP directory
## 

#uscope.merge.screen.merge.two.dir(FOLDER_DEST="/Volumes/elmicro/meta/phase/Conditions_phase_diagram/KO_library_with_P5O6/Ghil/full_screen/plate_9/Images/",
#                                  FOLDER_SOURCE="/Volumes/elmicro/meta/phase/Conditions_phase_diagram/KO_library_with_P5O6/Ghil/full_screen/plate_9/ImagesM1/",
#                                  SCREEN_NAME = "image1",
#                                  DOIT=1) 

#uscope.merge.screen.merge.two.dir(FOLDER_DEST="/Volumes/elmicro/meta/phase/Conditions_phase_diagram/KO_library_with_P5O6/Ghil/full_screen/plate_14/Images/",
#                                  FOLDER_SOURCE="/Volumes/elmicro/meta/phase/Conditions_phase_diagram/KO_library_with_P5O6/Ghil/full_screen/plate_14/ImagesG1/",
#                                SCREEN_NAME = "image1",
#                                DOIT=0) 

uscope.merge.screen.merge.two.dir = function(FOLDER_DEST="/Volumes/elmicro/meta/phase/Conditions_phase_diagram/KO_library_with_P5O6/Ghil/full_screen/plate_9/Images/",
                                             FOLDER_SOURCE="/Volumes/elmicro/meta/phase/Conditions_phase_diagram/KO_library_with_P5O6/Ghil/full_screen/plate_9/ImagesM1/",
                                           SCREEN_NAME = "p21",
                                           DOIT=0) {
  
  print(paste0("This function must be run after 'uscope.merge.screen.rename.sNum'"))
  print(paste0(""))
  print(paste0("Please make sure to run this function with 'DOIT=0' at first, and check the output by hand in the file ~/merge_tmp2.txt"))
  print(paste0(""))
  print(paste0("If you are happy with the file moves described in that file, proceed and run this function again, with the option DOIT=1"))
  print(paste0(""))
  
  ## If folder ends with "/", we trim it
  if(substr(FOLDER_DEST,nchar(FOLDER_DEST),nchar(FOLDER_DEST)) == "/"){
    FOLDER_DEST = substr(FOLDER_DEST,1,nchar(FOLDER_DEST)-1)
  }
  
  ## If folder ends with "/", we trim it
  if(substr(FOLDER_SOURCE,nchar(FOLDER_SOURCE),nchar(FOLDER_SOURCE)) == "/"){
    FOLDER_SOURCE = substr(FOLDER_SOURCE,1,nchar(FOLDER_SOURCE)-1)
  }
  
  print(paste("The images in folder",FOLDER_SOURCE," will be moved to the folder",FOLDER_DEST))
  print(paste("Images in",FOLDER_DEST," that already exist will be moved to ", paste0(FOLDER_SOURCE,"_BAKMERGE") ))
  msg = readline(prompt="Are you OK with this? Y/N ")
  
  if(msg == "Y"){
    
    ### We will move all doubles to a new directory 
    if( ! file.exists ( paste0(FOLDER_SOURCE,"_BAKMERGE") ) ){
      print(paste0("Creating directory: ", FOLDER_SOURCE,"_BAKMERGE"))
      dir.create(paste0(FOLDER_SOURCE,"_BAKMERGE"), showWarnings = TRUE)
    } else {
      print(paste0("Directory: ",FOLDER_SOURCE,"_BAKMERGE already exists"))
    }
    
    sink(file="~/merge_tmp2.txt")
    
    all.files.source = list.files(path=FOLDER_SOURCE)
    
    for( each.file in all.files.source){
      
      if(length(grep('.nd$', each.file))==0){
      
        if( file.exists(paste0(FOLDER_DEST,"/",each.file)) ){
          
          print(paste0("Moving ",FOLDER_DEST,"/",each.file, " ===> ", FOLDER_SOURCE, "_BAKMERGE/", each.file))
          print(paste0("Moving ",FOLDER_SOURCE,"/",each.file, " ===> ", FOLDER_DEST, "/", each.file))
          if(DOIT){
            file.rename( paste0(FOLDER_DEST,"/",each.file), to= paste0(FOLDER_SOURCE, "_BAKMERGE/", each.file))
            file.rename( paste0(FOLDER_SOURCE,"/",each.file), to= paste0(FOLDER_DEST, "/", each.file))
          }
        } else {
          print(paste0("Just moving ",FOLDER_SOURCE,"/",each.file, " ===> ", FOLDER_DEST, "/", each.file))
          if(DOIT){
            file.rename( paste0(FOLDER_SOURCE,"/",each.file), to= paste0(FOLDER_DEST, "/", each.file))
          }
        }
      }
    }
    sink()
  }
}


get.position.for.merge = function(well="P1", ppw=12){
  
  
}


merge.folders.old = function(folders=
                           c("/media/elusers/data/microscope/or/06noise/151213-YPD_005_Sat-1",
                             "/media/elusers/data/microscope/or/06noise/151213-YPD_005_Sat-1.2"),
                         positions.stop=c(769)){
  

    ## Make sure the position is correct with
  
  ## the first element of position.stop will be the first number taken by s1 from the next directory.  
  print(paste("Going through folder",folders[1]))
  
  for( each.F in 1:length(folders)){
    if(substr(folders[each.F],nchar(folders[each.F]),nchar(folders[each.F])) == "/"){
      folders[each.F] = substr(folders[each.F],1,nchar(folders[each.F])-1)
    }
  }
  
  ### First I need to move all doubles to a new directory (i.e., all files from first screen equal or above the number given
  if( ! file.exists ( paste(folders[1],"_BAKMERGE",sep="") ) ){
    print(paste("Creating directory", paste(folders[1],"_BAKMERGE",sep="")))
    dir.create(paste(folders[1],"_BAKMERGE",sep=""), showWarnings = TRUE)
  } else {
    print(paste("Directory ",paste(folders[1],"_BAKMERGE",sep="")," already exists"))
  }
  
  stk.files = list.files(path=folders[1])
  stk.files = stk.files[ grepl(pattern=".*s([0-9]+)\\.tif$",x=stk.files) | grepl(pattern=".*s([0-9]+)\\.stk$",x=stk.files)  ]
  new.index = positions.stop[1]
  
  for( each.file in stk.files){
    s.number = as.numeric(sub(pattern="^.*_s([0-9]+)\\..*$", replacement="\\1", x=each.file, perl=TRUE))
    print(paste(each.file, "S=",s.number))
    if(s.number >= new.index){            
      #print(paste("FROM: ",paste(folders[1],"/",each.file,sep="") , " TO=", paste(folders[1],"_BAKMERGE/",each.file,sep="")))
      file.rename(paste(folders[1],"/",each.file,sep=""), to=paste(folders[1],"_BAKMERGE/",each.file,sep=""))
    }
  }
  
  ### Then I need to move all files from second directory to first directory
  stk.files = list.files(path=folders[2])
  stk.files = stk.files[ grepl(pattern=".*s([0-9]+)\\.tif$",x=stk.files) | grepl(pattern=".*s([0-9]+)\\.stk$",x=stk.files)  ]
  for( each.file in stk.files){
    #print(each.file)
    s.number = as.numeric(sub(pattern="^.*_s([0-9]+)\\..*$", replacement="\\1", x=each.file, perl=TRUE))
    
    ## Either this replacement matches 
    new.file = sub("_s[0-9]+\\.",paste("_s",(s.number+new.index-1),"\\.",sep=""), each.file)
    ## Or that one does
    ## new.file = sub("_s[0-9]+-",paste("_s",(s.number+new.index-1),"-",sep=""), each.file)        
    
    #print(paste("FROM: ",paste(folders[2],"/",each.file,sep="") , " TO=", paste(folders[1],"/",new.file,sep="")))
    
    if( ! file.exists(paste(folders[1],"/",new.file,sep="")) ){
      file.copy(from=paste(folders[2],"/",each.file,sep=""), to=paste(folders[1],"/",new.file,sep=""))
    } else {
      print(paste("Problem: file",new.file," already exists in DIR1"))
    }
  }
}


## 
## After two folders have been merged, files should be renamed so that they are consistent.
## 
## --> THIS DOESN"T MAKE SENSE
##
rename.basename.folder = function(folder, basename = "MayaGFP_plate91"){
  
  stk.files = list.files(path=folder, pattern=".[stk|tif]$")
  for (each.file in stk.files) {
    file.second.part = sub(pattern="^(.*)(_w[0-9].*)$", replacement="\\2", x=each.file, perl=TRUE)
    new.file = paste(basename, file.second.part, sep="")
    file.rename(from=paste(folder,"/",each.file,sep=""), to=paste(folder,"/",new.file,sep=""))
    #print(paste("FROM: ",paste(folder,"/",each.file,sep="") , " TO=", paste(folder,"/",new.file,sep="")))        
  }    
}

toto = function(){
  
  for( each.F in 1:length(folders)){
    stk.files = list.files(path=folders[each.F], pattern=".stk")
  }
  
  print(paste("nd file is", nd.file))
  microscope.design$BASENAMES = c(microscope.design$BASENAMES,sub(nd.file, pattern=".nd",replacement=""))
  
  #nd.file = list.files(path=each.folder, pattern=".nd")
  if(length(nd.file) > 1){
    warning(paste("Problem: there is more than one nd file in", each.folder, " (",length(nd.file),")" ))
  }
  microscope.design$ND = c(microscope.design$ND, paste(each.folder,"/",nd.file, sep=""))     
  
  nd = read.table(microscope.design$ND[K], sep="\n", stringsAsFactors=FALSE)
  nd = nd[grep("^Stage",nd[,1]),]
  
  microscope.design$S.positions[[K]] = sub("Stage([0-9]+), [A-Z][0-9]+\\.([0-9]+)", "\\1", nd,perl=TRUE)
  microscope.design$WELLS[[K]]       = sub("Stage([0-9]+), ([A-Z][0-9]+)\\.([0-9]+)", "\\2", nd,perl=TRUE)
  microscope.design$PICNUM[[K]]      = sub("Stage([0-9]+), ([A-Z][0-9]+)\\.([0-9]+)", "\\3", nd,perl=TRUE)
  
  K=K+1
  
}    

get.plate.wells = function(ppw=10){
    
  WELLS.NUM = 1:24
  WELLS.LET = 1:16
  WELLS.NUM.rev = WELLS.NUM[ length(WELLS.NUM):1]  
  WELLS.times = ppw
  
  #print(paste("WELLS N=",WELLS.NUM))
  #print(paste("WELLS L=",WELLS.LET))
  
  WELL.IDs = c()
  PIC.NUM=0
  EVEN.PASS = TRUE ## To ZIG-ZAG among wells
  
  for( well.letter in WELLS.LET ){
    
    if(EVEN.PASS){
      WELLS.NUM.GO = WELLS.NUM
      EVEN.PASS = FALSE
    } else {
      WELLS.NUM.GO = WELLS.NUM.rev
      EVEN.PASS = TRUE
    }
    
    for (well.num in WELLS.NUM.GO){
        
        for ( well.pic in 1:WELLS.times){
          #WELL.IDs = c(WELL.IDs, paste(LETTERS[well.letter],well.num,'.',well.pic , sep="") )
          
          WELL.ID = paste(LETTERS[well.letter],well.num,'.',well.pic , sep="")
          #print(WELL.ID)
          WELL.IDs = c(WELL.IDs, WELL.ID)
        }
      }
  }
  return(WELL.IDs)
}

microscope.position.MM = function(file.path="./384_Joe_matrical.stg", PLATE = 384, NUM.from = 1, NUM.to=8, LET.from=1, LET.to=4, exclude=c() , TIMES=8, OFFSET=0){
  
  WELLS.NUM = NUM.from:NUM.to
  WELLS.LET = LET.from:LET.to
  WELLS.NUM.rev = WELLS.NUM[ length(WELLS.NUM):1]  
  WELLS.times = TIMES
  
  print(paste("WELLS N=",WELLS.NUM))
  print(paste("WELLS L=",WELLS.LET))
  
  file = c()
  file = c(file, '"Stage Memory List", Version 5.0')
  file = c(file, '0, 0, 0, 0, 0, 0, 0, "microns", "microns"' )
  file = c(file, '0' )
  file = c(file, length(WELLS.NUM)*length(WELLS.LET)*WELLS.times )
  
  X.init = 0
  Y.init = 0  
  
  if(PLATE == 96){
    WELL.sep = 9000
  } else {
    WELL.sep = 4500
  }
  
  PIC.size  = 300
  
  ##    
  #x.pos.in.well = rep( c(-4.5*PIC.size, -1.5*PIC.size, 1.5*PIC.size, 4.5*PIC.size,  4.5*PIC.size, 1.5*PIC.size, -1.5*PIC.size, -4.5*PIC.size) ,3)
  #y.pos.in.well = rep( c( 4.5*PIC.size,  1.5*PIC.size,-1.5*PIC.size,-4.5*PIC.size), rep(4,4) )
  
  ##  Goes like this within each well:
  ##      ------- |
  ##     | -----| |
  ##     | | >--| |
  ##     | |------|
  ##     |---------    
  
  x.pos.in.well = c(0 , 1.5*PIC.size, 1.5*PIC.size, 0, -1.5*PIC.size, -1.5*PIC.size, -1.5*PIC.size, 0, 1.5*PIC.size, 3*PIC.size, 3*PIC.size, 3*PIC.size, 3*PIC.size, 1.5*PIC.size, 0 , -1.5*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size, -1.5*PIC.size, 0*PIC.size, 1.5*PIC.size, 3*PIC.size)
  y.pos.in.well = c(0, 0, 1.5*PIC.size, 1.5*PIC.size, 1.5*PIC.size,0, -1.5*PIC.size, -1.5*PIC.size, -1.5*PIC.size, -1.5*PIC.size,0, 1.5*PIC.size, 3*PIC.size, 3*PIC.size, 3*PIC.size, 3*PIC.size, 3*PIC.size, 1.5*PIC.size, 0*PIC.size, -1.5*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size)
  
  
  WELL.IDs = c()
  
  PIC.NUM=0
  
  EVEN.PASS = TRUE ## To ZIG-ZAG among wells
  
  for( well.letter in WELLS.LET ){
    
    if(EVEN.PASS){
      WELLS.NUM.GO = WELLS.NUM
      EVEN.PASS = FALSE
    } else {
      WELLS.NUM.GO = WELLS.NUM.rev
      EVEN.PASS = TRUE
    }
    
    for (well.num in WELLS.NUM.GO){
      
      Y = -(  Y.init + ( (-1+well.letter) * WELL.sep))
      X = -(  X.init + ( (-1+well.num   ) * WELL.sep))
      
      if( ! paste(LETTERS[well.letter],well.num,sep="") %in% exclude){
        
        for ( well.pic in 1:WELLS.times){
          
          X.tmp.center = X
          Y.tmp.center = Y
          
          X.tmp = X.tmp.center + x.pos.in.well[well.pic]
          Y.tmp = Y.tmp.center + y.pos.in.well[well.pic]
          
          WELL.IDs = c(WELL.IDs, paste(LETTERS[well.letter],well.num,'.',well.pic , sep="") )
          
          WELL.ID = paste(LETTERS[well.letter],well.num,'.',well.pic , sep="")
          
          file = c(file, paste('"',WELL.ID,'"',", ",X.tmp,", ", Y.tmp,", 0, 0, 0, FALSE, -9999, TRUE, TRUE, 0",sep=""))
          
        }
      }
    }
  }
  #write.table(file, file="/home/elevy/.gvfs/elusers on mata.weizmann.ac.il/data/microscope/R/384-matrical.stg", quote=FALSE, eol="\r", row.names=FALSE, col.names=FALSE)
  write.table(file, file=file.path, quote=FALSE, eol="\r", row.names=FALSE, col.names=FALSE)
}


## Create Excel reports
##
## It will show avg per well, per strain, and per construct
##
uscope.get.report.mean = function(data, design, res.name=""){
  
  ##
  ## Gets the mean
  ##
  res.all           = list()
  res.all.construct = list()
    
  for( each.ch in design$CHANELS){
    res.all[[each.ch]] = data.frame(construct=design$PDEF[[each.ch]])
    #res.all.construct[[each.ch]] = data.frame(constructs = )
  }
  
  N.samples = length(data[[1]])
  
  for(K in 1:length(design$D)) {
      
    all.means = uscope.process.get.one.mean(design, data, K)
    
    ##
    ## FIRST WE print things per WELL
    ##
    ## For each channel --> creates a barplot showing avg intensity/well
    ##
    for( each.ch in design$CHANELS){          
      
      if(N.samples > 90){
          pdf(file=paste(design$DIR.res, "/MEAN.per.well-",design$D[K],"-",res.name,"-",each.ch,".pdf",sep=""), width=12*(N.samples/90), height=9)
      } else {
          
          pdf(file=paste(design$DIR.res, "/MEAN.per.well-",design$D[K],"-",res.name,"-",each.ch,".pdf",sep=""), width=12, height=9)
      }
        
      par(mai=c(3,1,1,1))
      
      col.name = get.color.id(each.ch)
      tabbedMeans = all.means[[each.ch]]$mean
      tabbedSE    = all.means[[each.ch]]$sd/2
      
      names(tabbedMeans) = design$PDEF[[each.ch]]
      
      barCenters = barplot(tabbedMeans, col=col.name, ylim=c(0, max(tabbedMeans+tabbedSE, na.rm=TRUE) ), las=2, cex.lab=0.5 )
      
      arrows(barCenters, tabbedMeans - tabbedSE, barCenters,
             tabbedMeans + tabbedSE, lwd = 1.5, angle = 90,
             code = 3, length = 0.05)
      text(barCenters, tabbedMeans+(mean(tabbedMeans, na.rm=TRUE)/10), all.means[[each.ch]]$n)
      
      res.all[[each.ch]] = cbind(res.all[[each.ch]], round(all.means[[each.ch]]$mean), round(all.means[[each.ch]]$sd), all.means[[each.ch]]$n)
      
      colnames(res.all[[each.ch]])[ (NCOL(res.all[[each.ch]])-2):NCOL(res.all[[each.ch]])] = paste(paste(each.ch,".",design$D[K],sep=""), c("",".sd",".n"),sep="")
      dev.off()
    }
    
    ##
    ## THEN WE SEARCH FOR IDENTICAL CONSTRUCTS
    ##
    for( each.ch in design$CHANELS){
      
      pdf(file=paste(design$DIR.res, "/MEAN.per.construct-",design$D[K],"-",res.name,"-",each.ch,".pdf",sep=""), width=12, height=8)
      
      par(mai=c(3,1,0.5,0.5))
      
      col.name = get.color.id(each.ch)
      
      
      tabbedMeans2      = aggregate(x=data.frame( x=all.means[[each.ch]]$mean), by=list(id=design$PDEF[[each.ch]]), mean, na.rm=TRUE)
      tabbedMeans.order = aggregate(x=data.frame( x=design$PDEF$well), by=list(id=design$PDEF[[each.ch]]), function(x){return(x[1]) })
      tabbedMeans2$ord = tabbedMeans.order$x
      tabbedMeans2 = tabbedMeans2[ order(tabbedMeans2$ord),]
      
      tabbedMeans = as.matrix(tabbedMeans2[,2])           
      rownames(tabbedMeans) = tabbedMeans2$id
      
      tabbedSE    = aggregate(x=all.means[[each.ch]]$mean, by=list(id=design$PDEF[[each.ch]]), sd, na.rm=TRUE)[,2]
      
      ## We cannot calculate a SD if the number of samples is 1 ..
      tabbedSE[ is.na(tabbedSE) | is.nan(tabbedSE)]=0
      
      barCenters = barplot(t(tabbedMeans), col=col.name, ylim=c(0, max(tabbedMeans[,1]+tabbedSE, na.rm=TRUE)) ,las=2)
      
      arrows(barCenters, tabbedMeans[,1] - tabbedSE, barCenters,
             tabbedMeans[,1] + tabbedSE, lwd = 1.5, angle = 90,
             code = 3, length = 0.05)
      #text(barCenters, tabbedMeans+(mean(tabbedMeans, na.rm=TRUE)/10), all.means[[each.ch]]$n)
      dev.off()
    }              
  }
  
  for( each.ch in design$CHANELS){
    write.csv(file=paste(design$DIR.res, "/tab-report-",res.name,"-",each.ch,".csv",sep=""), x=res.all[[each.ch]])
    #res.all.construct[[each.ch]] = data.frame(constructs = )

  }
  
}


## Create Excel reports
##
## It will show avg per well, per strain, and per construct
##
uscope.get.report.custom = function(data, design, res.name=""){ ##### TO BE COMPLETED
  
  ##
  ## Gets the mean
  ##
  res.all           = list()
  res.all.construct = list()
    
  for( each.ch in design$CHANELS){
    res.all[[each.ch]] = data.frame(construct=design$PDEF[[each.ch]])
    #res.all.construct[[each.ch]] = data.frame(constructs = )
  }
  
  N.samples = length(data[[1]])
  
  for(K in 1:length(design$D)) {
      
      all.means = uscope.process.get.one.stat(design, data, K, myCOL, myFUN)
    
    ##
    ## FIRST WE print things per WELL
    ##
    ## For each channel --> creates a barplot showing avg intensity/well
    ##
    for( each.ch in design$CHANELS){          
      
      if(N.samples > 90){
          pdf(file=paste(design$DIR.res, "/MEAN.per.well-",design$D[K],"-",res.name,"-",each.ch,".pdf",sep=""), width=12*(N.samples/90), height=9)
      } else {
          pdf(file=paste(design$DIR.res, "/MEAN.per.well-",design$D[K],"-",res.name,"-",each.ch,".pdf",sep=""), width=12, height=9)
      }
        
      par(mai=c(3,1,1,1))
      
      col.name = get.color.id(each.ch)
      tabbedMeans = all.means[[each.ch]]$mean
      tabbedSE    = all.means[[each.ch]]$sd/2
      
      names(tabbedMeans) = design$PDEF[[each.ch]]
      
      barCenters = barplot(tabbedMeans, col=col.name, ylim=c(0, max(tabbedMeans+tabbedSE, na.rm=TRUE) ), las=2, cex.lab=0.5 )
      
      arrows(barCenters, tabbedMeans - tabbedSE, barCenters,
             tabbedMeans + tabbedSE, lwd = 1.5, angle = 90,
             code = 3, length = 0.05)
      text(barCenters, tabbedMeans+(mean(tabbedMeans, na.rm=TRUE)/10), all.means[[each.ch]]$n)
      
      res.all[[each.ch]] = cbind(res.all[[each.ch]], round(all.means[[each.ch]]$mean), round(all.means[[each.ch]]$sd), all.means[[each.ch]]$n)
      
      colnames(res.all[[each.ch]])[ (NCOL(res.all[[each.ch]])-2):NCOL(res.all[[each.ch]])] = paste(paste(each.ch,".",design$D[K],sep=""), c("",".sd",".n"),sep="")
      dev.off()
    }
    
    ##
    ## THEN WE SEARCH FOR IDENTICAL CONSTRUCTS
    ##
    for( each.ch in design$CHANELS){
      
      pdf(file=paste(design$DIR.res, "/MEAN.per.construct-",design$D[K],"-",res.name,"-",each.ch,".pdf",sep=""), width=12, height=8)
      
      par(mai=c(3,1,0.5,0.5))
      
      col.name = get.color.id(each.ch)
      
      
      tabbedMeans2      = aggregate(x=data.frame( x=all.means[[each.ch]]$mean), by=list(id=design$PDEF[[each.ch]]), mean, na.rm=TRUE)
      tabbedMeans.order = aggregate(x=data.frame( x=design$PDEF$well), by=list(id=design$PDEF[[each.ch]]), function(x){return(x[1]) })
      tabbedMeans2$ord = tabbedMeans.order$x
      tabbedMeans2 = tabbedMeans2[ order(tabbedMeans2$ord),]
      
      tabbedMeans = as.matrix(tabbedMeans2[,2])           
      rownames(tabbedMeans) = tabbedMeans2$id
      
      tabbedSE    = aggregate(x=all.means[[each.ch]]$mean, by=list(id=design$PDEF[[each.ch]]), sd, na.rm=TRUE)[,2]
      
      ## We cannot calculate a SD if the number of samples is 1 ..
      tabbedSE[ is.na(tabbedSE) | is.nan(tabbedSE)]=0
      
      barCenters = barplot(t(tabbedMeans), col=col.name, ylim=c(0, max(tabbedMeans[,1]+tabbedSE, na.rm=TRUE)) ,las=2)
      
      arrows(barCenters, tabbedMeans[,1] - tabbedSE, barCenters,
             tabbedMeans[,1] + tabbedSE, lwd = 1.5, angle = 90,
             code = 3, length = 0.05)
      #text(barCenters, tabbedMeans+(mean(tabbedMeans, na.rm=TRUE)/10), all.means[[each.ch]]$n)
      dev.off()
    }              
  }
  
  for( each.ch in design$CHANELS){
    write.csv(file=paste(design$DIR.res, "/tab-report-",res.name,"-",each.ch,".csv",sep=""), x=res.all[[each.ch]])
    #res.all.construct[[each.ch]] = data.frame(constructs = )

  }
  
}





## General directory to save images and reports
init.result.folders = function(design){
  
  if( ! file.exists ( design$DIR.res ) ){ stop("the directory given for the results does not exist") }
  
  for( each.plate in design$D){
    if( ! file.exists ( paste(design$DIR.res,"/", each.plate,sep="") ) ){
      dir.create(paste(design$DIR.res,"/", each.plate,sep=""), showWarnings = TRUE)
    }
  }
}


## Counts the total number of cells in the given data structure
uscope.count.cells = function(data){

    for(K in 1:length(data)){

        cols.to.sum = 0
        
        for(L in 1:length(data[[K]])){

            cols.to.sum = cols.to.sum + NROW(data[[K]][[L]])
        }

        print(paste("For plate",K, "there is a total of ",cols.to.sum," cells"));
    }
}


.BDtestOK = function(){
  
  # GOOD EXAMPLE FROM OR's EXPERIMENTS
  design.or = microscope.get.design(
    F=c("/media/elusers/data/microscope/or/02STAB/160120-STAB_all_002"), ## folder of PICS
    D=c("STAB"),
    PDEF=c("/media/elusers/data/microscope/plate_defs/160120-STAB_New_collection.csv"),
    FORMAT=384,
    OUT="_output",
    CHANELS=c("GFP"),
    MEASURE.COL = "_int_b5",
    DIR.res = "/data/benjamin/"
  )
  
  data.or = microscope.load.data(design.or)
  
  design.or = uscope.process.estimate.background(data.or,design.or)
  
  data.or    = uscope.process.reorder(data.or)
  data.or    = uscope.process.remove.first.pic(data.or)
  data.or1   = uscope.process.remove.background(data.or, design.or)
  data.or2   = uscope.process.remove.small(data.or1, MIN.size=800,MAX.size=2000)
  data.or2.0 = uscope.process.BF(data.or2)
  data.or2.1 = uscope.process.remove.BF.outliers(data.or2.0, cutoff=0.8)
  
  stat = uscope.process.get.one.stat(design = design.or, data = data.or2.1, col.of.interest = 'GFP_int_b5', plateNUM = "STAB" ) 
  
  M = aggregate(x = data.frame( mean = stat$mean ), by = list(id = stat$id), mean, na.rm=TRUE )
  S = aggregate(x = data.frame( sd = stat$mean ), by = list(id = stat$id), sd,  na.rm=TRUE )
  bp = merge( M, S, by='id', all.x=TRUE, all.y=TRUE, sort=FALSE)
  bp = bp[order(bp$mean),]
  
  # BARS WILL BE ORDERED IN ASCENDING VALUE ON Y-AXIS
  bp$order = factor(x = bp$id, levels = bp$id)
  
  require(ggplot2)
  ggplot( bp, aes(x=order,y=mean)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
}


.BDtestERROR1 = function(){
  
  design.correct = design = microscope.get.design(
    F=c("/media/elusers/data/microscope/or/02STAB/160120-STAB_all_002"), ## folder of PICS
    D=c("STAB"),
    PDEF=c("/media/elusers/data/microscope/plate_defs/160120-STAB_New_collection.csv"),
    FORMAT=384,
    OUT="_output",
    CHANELS=c("GFP"),
    MEASURE.COL = "_int_b5",
    DIR.res = "/data/benjamin/"
  )
  
  ##### BAD EXAMPLES OF DESIGNS #####
  # - fake directory for PICS
  design$F = "media/elusers/data/microscope/or/02STAB/160120-STAB_all_002"
  
  # - No decription corresponding to the picture folders
  design = design.correct  
  design$D = c()
  
  # - PLATE DEF WITH NO DATA BUT WITH HEADER
  design = design.correct  
  design$PDEF = "/media/elusers/data/microscope/plate_defs/Empty.csv"
  
  # - PLATE DEF DOES NOT CONTAIN ANYTHING
  design = design.correct  
  design$PDEF = "/media/elusers/data/microscope/plate_defs/None.csv"
  
  # - WRONG FORMAT SPECIFIED
  design = design.correct  
  design$FORMAT = 96
  
  # - CHANELS NOT USED IN PLATE DEFINITION FILE
  design = design.correct  
  design$CHANELS = c("blabla")
  
  # - FAKE RESULTS DIRECTORY
  design = design.correct  
  design$DIR.res = "data/benjamin/"
}  


###
### Function to generate positions for the microscope
###
get.microscope.positionList = function(file.path="./384_Joe_matrical.stg", PLATE = 384, NUM.from = 1, NUM.to=8, LET.from=1, LET.to=4, exclude=c() , TIMES=8, OFFSET=0){

    WELLS.NUM = NUM.from:NUM.to
    WELLS.LET = LET.from:LET.to
    WELLS.NUM.rev = WELLS.NUM[ length(WELLS.NUM):1]  
    WELLS.times = TIMES
  
    print(paste("WELLS N=",WELLS.NUM))
    print(paste("WELLS L=",WELLS.LET))
    
    file = c()
    file = c(file, '"Stage Memory List", Version 5.0')
    file = c(file, '0, 0, 0, 0, 0, 0, 0, "microns", "microns"' )
    file = c(file, '0' )
    file = c(file, length(WELLS.NUM)*length(WELLS.LET)*WELLS.times )
       
    X.init = 0
    Y.init = 0  

    if(PLATE == 96){
        WELL.sep = 9000
    } else {
        WELL.sep = 4500
    }

    PIC.size  = 300
    
    ##    
    #x.pos.in.well = rep( c(-4.5*PIC.size, -1.5*PIC.size, 1.5*PIC.size, 4.5*PIC.size,  4.5*PIC.size, 1.5*PIC.size, -1.5*PIC.size, -4.5*PIC.size) ,3)
    #y.pos.in.well = rep( c( 4.5*PIC.size,  1.5*PIC.size,-1.5*PIC.size,-4.5*PIC.size), rep(4,4) )

    ##  Goes like this within each well:
    ##      ------- |
    ##     | -----| |
    ##     | | >--| |
    ##     | |------|
    ##     |---------    

    x.pos.in.well = c(0 , 1.5*PIC.size, 1.5*PIC.size, 0, -1.5*PIC.size, -1.5*PIC.size, -1.5*PIC.size, 0, 1.5*PIC.size, 3*PIC.size, 3*PIC.size, 3*PIC.size, 3*PIC.size, 1.5*PIC.size, 0 , -1.5*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size, -1.5*PIC.size, 0*PIC.size, 1.5*PIC.size, 3*PIC.size)
    y.pos.in.well = c(0, 0, 1.5*PIC.size, 1.5*PIC.size, 1.5*PIC.size,0, -1.5*PIC.size, -1.5*PIC.size, -1.5*PIC.size, -1.5*PIC.size,0, 1.5*PIC.size, 3*PIC.size, 3*PIC.size, 3*PIC.size, 3*PIC.size, 3*PIC.size, 1.5*PIC.size, 0*PIC.size, -1.5*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size)

    
    WELL.IDs = c()
    
    PIC.NUM=0
    
    EVEN.PASS = TRUE ## To ZIG-ZAG among wells
  
    for( well.letter in WELLS.LET ){
    
        if(EVEN.PASS){
            WELLS.NUM.GO = WELLS.NUM
            EVEN.PASS = FALSE
        } else {
            WELLS.NUM.GO = WELLS.NUM.rev
            EVEN.PASS = TRUE
        }
    
        for (well.num in WELLS.NUM.GO){
            
            Y = -(  Y.init + ( (-1+well.letter) * WELL.sep))
            X = -(  X.init + ( (-1+well.num   ) * WELL.sep))

            if( ! paste(LETTERS[well.letter],well.num,sep="") %in% exclude){
      
                for ( well.pic in 1:WELLS.times){

                    X.tmp.center = X
                    Y.tmp.center = Y
        
                    X.tmp = X.tmp.center + x.pos.in.well[well.pic]
                    Y.tmp = Y.tmp.center + y.pos.in.well[well.pic]
                    
                    WELL.IDs = c(WELL.IDs, paste(LETTERS[well.letter],well.num,'.',well.pic , sep="") )

                    WELL.ID = paste(LETTERS[well.letter],well.num,'.',well.pic , sep="")
                    
                    file = c(file, paste('"',WELL.ID,'"',", ",X.tmp,", ", Y.tmp,", 0, 0, 0, FALSE, -9999, TRUE, TRUE, 0",sep=""))
                                        
                }
            }
        }
    }
    #write.table(file, file="/home/elevy/.gvfs/elusers on mata.weizmann.ac.il/data/microscope/R/384-matrical.stg", quote=FALSE, eol="\r", row.names=FALSE, col.names=FALSE)
    write.table(file, file=file.path, quote=FALSE, eol="\r", row.names=FALSE, col.names=FALSE)
}




###
### Function to generate positions following the order of the microscope
###
get.microscope.positions = function(file.path="~/384_screen_order.csv", PLATE = 384, NUM.from = 1, NUM.to=20, LET.from=1, LET.to=5, exclude=c() , TIMES=1, OFFSET=0){

    WELLS.NUM = NUM.from:NUM.to
    WELLS.LET = LET.from:LET.to
    WELLS.NUM.rev = WELLS.NUM[ length(WELLS.NUM):1]  
    WELLS.times = TIMES
  
    print(paste("WELLS N=",WELLS.NUM))
    print(paste("WELLS L=",WELLS.LET))
    
    file = c()
    
    X.init = 0
    Y.init = 0  

    if(PLATE == 96){
        WELL.sep = 9000
    } else {
        WELL.sep = 4500
    }

    PIC.size  = 300
    
    ##    
    #x.pos.in.well = rep( c(-4.5*PIC.size, -1.5*PIC.size, 1.5*PIC.size, 4.5*PIC.size,  4.5*PIC.size, 1.5*PIC.size, -1.5*PIC.size, -4.5*PIC.size) ,3)
    #y.pos.in.well = rep( c( 4.5*PIC.size,  1.5*PIC.size,-1.5*PIC.size,-4.5*PIC.size), rep(4,4) )

    ##  Goes like this within each well:
    ##      ------- |
    ##     | -----| |
    ##     | | >--| |
    ##     | |------|
    ##     |---------    

    x.pos.in.well = c(0 , 1.5*PIC.size, 1.5*PIC.size, 0, -1.5*PIC.size, -1.5*PIC.size, -1.5*PIC.size, 0, 1.5*PIC.size, 3*PIC.size, 3*PIC.size, 3*PIC.size, 3*PIC.size, 1.5*PIC.size, 0 , -1.5*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size, -1.5*PIC.size, 0*PIC.size, 1.5*PIC.size, 3*PIC.size)
    y.pos.in.well = c(0, 0, 1.5*PIC.size, 1.5*PIC.size, 1.5*PIC.size,0, -1.5*PIC.size, -1.5*PIC.size, -1.5*PIC.size, -1.5*PIC.size,0, 1.5*PIC.size, 3*PIC.size, 3*PIC.size, 3*PIC.size, 3*PIC.size, 3*PIC.size, 1.5*PIC.size, 0*PIC.size, -1.5*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size, -3*PIC.size)

    
    WELL.IDs = c()
    
    PIC.NUM=0
    
    EVEN.PASS = TRUE ## To ZIG-ZAG among wells
  
    for( well.letter in WELLS.LET ){
    
        if(EVEN.PASS){
            WELLS.NUM.GO = WELLS.NUM
            EVEN.PASS = FALSE
        } else {
            WELLS.NUM.GO = WELLS.NUM.rev
            EVEN.PASS = TRUE
        }
    
        for (well.num in WELLS.NUM.GO){
            
            Y = -(  Y.init + ( (-1+well.letter) * WELL.sep))
            X = -(  X.init + ( (-1+well.num   ) * WELL.sep))

            if( ! paste(LETTERS[well.letter],well.num,sep="") %in% exclude){
      
                for ( well.pic in 1:WELLS.times){

                    X.tmp.center = X
                    Y.tmp.center = Y
        
                    X.tmp = X.tmp.center + x.pos.in.well[well.pic]
                    Y.tmp = Y.tmp.center + y.pos.in.well[well.pic]
                    
                    WELL.IDs = c(WELL.IDs, paste(LETTERS[well.letter],well.num,'.',well.pic , sep="") )                    
                                        
                    WELL.ID = paste(LETTERS[well.letter],well.num , sep="")

                    WELL.ID = gsub(WELL.ID, pattern="([A-Z])([0-9]$)", replacement="\\10\\2")

                    file = c(file, WELL.ID)
                                        
                }
            }
        }
    }
    #write.table(file, file="/home/elevy/.gvfs/elusers on mata.weizmann.ac.il/data/microscope/R/384-matrical.stg", quote=FALSE, eol="\r", row.names=FALSE, col.names=FALSE)
    write.csv(cbind(1:length(file),file), file=file.path, quote=FALSE, row.names=FALSE, col.names=FALSE)
}
