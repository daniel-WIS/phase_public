#' This does the initial postprocessing of the txt files that are usually in images_output after 
#' running uscope_analyze_FLUO
#' 
#'
#' @param design.file 
#' @param min.cell.size 
#' @param max.cell.size 
#' @param brightfield.cutoff 
#'
#' @return
#' @export
#'
#' @examples
postprocess_result_txt <-
  function(design.file,
           min.cell.size = 1000,
           max.cell.size = 2500,
           brightfield.cutoff = 0.8) {
    data.raw = microscope.load.data(design.file)
    data.raw = uscope.process.add.ncells(data = data.raw)
    lsos() ## check how much memory the object takes
    
    design.file = uscope.process.estimate.background(data.raw, design.file)
    data.1    = uscope.process.reorder(data.raw, design = design.file)
    data.1    = uscope.process.remove.first.pic(data.1)
    data.1    = uscope.process.remove.background(data.1, design.file)
    
    n_cells_before = uscope.count.cells(data.1)
    
    data.1    = uscope.process.remove.small(data.1, MIN.size = min.cell.size, MAX.size =
                                              max.cell.size)
    data.1   = uscope.process.BF(data.1)
    data.12   = uscope.process.remove.BF.outliers(data.1, cutoff = brightfield.cutoff)
    
    data.12dan   = uscope.process.add.ncells(data = data.12)
    
    n_cells_after = uscope.count.cells(data.12)
    
    print(
      n_cells_before,
      "cells were loaded and",
      n_cells_after,
      "cells remained after initial postprocessing."
    )
    
  }