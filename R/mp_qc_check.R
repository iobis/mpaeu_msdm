mp_qc_check <- function(dataset,
                        qc_steps = c(
                          "on_sea",
                          "on_depth",
                          "outlier_spatial",
                          "outlier_distcoast"
                        )
                        sea_layer = NULL) {
  
  if ("on_sea" %in% qc_steps) {
    if (!is.null(sea_layer)) {
      
    } else {
      obistools
    }
  }
  
  if ("dist_coast") {
    
  }
  
}