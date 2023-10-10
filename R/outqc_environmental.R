#' Detect geographical outliers based on distance between points
#'
#' @param pts the points for which detect outliers
#' @param variables the name of the variables that should be considered. One or more of
#' "bathymetry", "sstemperature", "sssalinity" or "shoredistance". Alternatively,
#' you can supply a SpatRaster with environmental layers. In that case, the function
#' will extract the information from those points.
#' @param methods which methods to use to detect outliers. One of "mdist_iqr",
#' "mdist_mad" or "isoforest" (see details)
#' @param iqr_mltp the multiplier of the interquantile range to use
#' @param mad_mltp the multiplier of the MAD to use
#' @param isoforest_th the threshold for the isoforest detection
#' @param iso_single if \code{TRUE}, then the Isolation Forest algorithm will run in
#' single variable mode, across each variable. If \code{FALSE} (default), it will
#' consider all variables together in a single model.
#' @param limit_rem a number between 0 and 1 limiting the percentage of points that can
#' be considered outliers. The default value (0.01) says that no more than 1% of points
#' can be outliers
#'
#' @return a `matrix` with columns equal the number of methods * variables used. Each line is 
#' one observation; values of 0 indicates normal values, while 1 indicates outliers.
#' @export
#' 
#' @details
#' Methods available are:
#' - iqr: see if the observation is within \code{iqr_mltp} IQRs from the first and third quartile
#' - mad: see if the observation is within \code{mad_mltp} MADs from the distance to the median value
#' - isoforest: uses an extended isolation forest model (or single model, if \code{iso_single = TRUE}) 
#' to detect outliers. The higher the value of the output of such model, the more
#' likely it is to be an outlier. What value to use as a
#' threshold is difficulty to compute, as we don't have any control information
#' of what is indeed an outlier (otherwise we could have removed it!). After
#' tests with simulated and real data, the 0.8 value seems to be a reasonable
#' threshold.
#' 
#' If variables names are supplied, the [obistools::lookup_xy()] function is used
#' to get the environmental information.
#' 
#' \code{NA} values are not considered in the calculations, and are marked as NA
#' in the final outputs.
#' 
#'
#' @examples
#' \dontrun{
#' 
#' }
outqc_env <- function(pts,
                      variables = c("bathymetry", "sstemperature", "sssalinity", "shoredistance"),
                      methods = c("iqr", "mad", "isoforest"),
                      iqr_mltp = 3,
                      mad_mltp = 6,
                      isoforest_th = 0.8,
                      iso_single = FALSE,
                      limit_rem = 0.01) {
  
  
  if (class(variables)[1] == "SpatRaster") {
    data_info <- terra::extract(variables, pts, ID = F)
  } else {
    if (any(!variables %in% c("bathymetry", "sstemperature", "sssalinity", "shoredistance"))) {
      stop('Variables must be one of "bathymetry", "sstemperature", "sssalinity" or "shoredistance".')
    }
    data_info <- obistools::lookup_xy(pts)
    data_info <- data_info[,variables]
    if (length(variables) == 1) {
      eval(parse(text = paste0("data_info <- data.frame(", variables, " = data_info)")))
    }
  }
  
  if (any(grepl("iqr", methods))) {
    
    data_iqr_res <- apply(data_info, 2, function(x) {
      
      to <- order(x, decreasing = T)
      
      iqr <- .iqr_out(x, threshold = iqr_mltp, id = -999999)
      
      iqr_res <- ifelse(iqr == -999999, 1, 0)
      iqr_res[-to[1:ceiling(length(iqr_res)*limit_rem)]] <- 0
      
      iqr_res
      
    })
    
    colnames(data_iqr_res) <- paste0(colnames(data_iqr_res), "_iqr")
  }
  
  if (any(grepl("mad", methods))) {
    
    data_mad_res <- apply(data_info, 2, function(x) {
      
      to <- order(x, decreasing = T)
      
      mad <- .mad_out(x, threshold = mad_mltp, id = -999999)
      
      mad_res <- ifelse(mad == -999999, 1, 0)
      mad_res[-to[1:ceiling(length(mad_res)*limit_rem)]] <- 0
      
      mad_res
      
    })
    
    colnames(data_mad_res) <- paste0(colnames(data_mad_res), "_mad")
  }
  
  if (any(grepl("isoforest", methods))) {
    
    if (iso_single) {
      
      data_isoforest_res <- apply(data_info, 2, function(x) {
        
        isof_m <- isotree::isolation.forest(
          data.frame(var = x),
          ndim = 1,
          ntrees = 500,
          nthreads = 1,
          penalize_range = FALSE,
          prob_pick_pooled_gain = 0,
          prob_pick_avg_gain = 0)
        
        isof_vals <- predict(isof_m, data.frame(var = x))
        
        to <- order(isof_vals, decreasing = T)
        
        isof_vals <- ifelse(isof_vals >= isoforest_th, 1, 0)
        isof_vals[-to[1:ceiling(length(isof_vals)*limit_rem)]] <- 0
        
        isof_vals
        
      })
      
      colnames(data_isoforest_res) <- paste0(colnames(data_isoforest_res), "_iso")
      
    } else {
      
      isof_m <- isotree::isolation.forest(
        data_info,
        ndim = 2,
        ntrees = 500,
        nthreads = 1,
        penalize_range = FALSE,
        prob_pick_pooled_gain = 0,
        prob_pick_avg_gain = 0)
      
      isof_vals <- predict(isof_m, data_info)
      
      to <- order(isof_vals, decreasing = T)
      
      isof_vals <- ifelse(isof_vals >= isoforest_th, 1, 0)
      isof_vals[-to[1:ceiling(length(isof_vals)*limit_rem)]] <- 0
      
      data_isoforest_res <- data.frame(isoforest = isof_vals)
      
    }
  }
  
  for (i in 1:length(methods)) {
    if (i==1) {
      results <- eval(parse(text = paste0("data_", methods[i], "_res")))
    } else {
      results <- cbind(results,
                       eval(parse(text = paste0("data_", methods[i], "_res"))))
    }
  }
  
  return(results)
}