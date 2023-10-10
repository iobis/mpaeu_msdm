#' Detect geographical outliers based on distance between points
#'
#' @param pts the points for which detect outliers
#' @param dist_folder the folder where the distance layers produced with [outqc_get_distances()]
#' are stored
#' @param methods which methods to use to detect outliers. One of "mdist_iqr",
#' "mdist_mad" or "isoforest" (see details)
#' @param k the number of nearest points to consider when computing the metrics.
#' If \code{NULL}, then nrow(pts) - 1 is used.
#' @param iqr_mltp the multiplier of the interquantile range to use
#' @param mad_mltp the multiplier of the MAD to use
#' @param isoforest_th the threshold for the isoforest detection
#' @param mdist_method which method to use to get the mean distance. 
#' Usually mean, but can be "median"
#' @param limit_rem a number between 0 and 1 limiting the percentage of points that can
#' be considered outliers. The default value (0.01) says that no more than 1% of points
#' can be outliers
#' @param parallel enable parallel computation on the extraction of values from
#' the distance layers (see [outqc_query_distances()])
#' @param mc_cores number of cores to be used by parallel extraction
#'
#' @return a `matrix` with columns equal the number of methods used. Each line is 
#' one observation; values of 0 indicates normal values, while 1 indicates outliers.
#' @export
#' 
#' @details
#' Methods available are:
#' - mdist_iqr: see if the observation is within \code{iqr_mltp} IQRs from the first and third quartile
#' - mdist_mad: see if the observation is within \code{mad_mltp} MADs from the distance to the median value
#' - isoforest: uses an extended isolation forest model to detect outliers. The higher the value of the output
#' of such model, the more likely it is to be an outlier. What value to use as a threshold is difficulty
#' to compute, as we don't have any control information of what is indeed an outlier (otherwise we could have removed it!).
#' After tests with simulated and real data, the 0.8 value seems to be a reasonable threshold.
#' 
#'
#' @examples
#' \dontrun{
#' 
#' }
outqc_geo <- function(pts,
                      dist_folder = "distances",
                      methods = c("mdist_iqr", "mdist_mad", "isoforest"),
                      k = 15,
                      iqr_mltp = 3,
                      mad_mltp = 6,
                      isoforest_th = 0.8,
                      mdist_method = "mean",
                      limit_rem = 0.01,
                      parallel = TRUE,
                      mc_cores = NULL) {
  
  if (is.null(k)) {
    mdk <- nrow(pts)-1
  } else {
    mdk <- k
  }
  
  dists <- outqc_query_distances(pts,
                                 kdist = mdk,
                                 distfolder = dist_folder,
                                 parallel = parallel,
                                 mc_cores = mc_cores,
                                 returnid = T)
  
  if (any(grepl("mdist", methods))) {
    
    m_dists <- apply(dists$dist, 1, mdist_method)
    
    to <- order(m_dists, decreasing = T)
    
    if (any(grepl("mdist_iqr", methods))) {
      
      iqr <- .iqr_out(m_dists, threshold = iqr_mltp, id = -999999)
      
      mdist_iqr_res <- ifelse(iqr == -999999, 1, 0)
      mdist_iqr_res[-to[1:ceiling(length(mdist_iqr_res)*limit_rem)]] <- 0
    }
    
    if (any(grepl("mdist_mad", methods))) {
      
      mad <- .mad_out(m_dists, threshold = mad_mltp, id = -999999)
      
      mdist_mad_res <- ifelse(mad == -999999, 1, 0)
      mdist_mad_res[-to[1:ceiling(length(mdist_mad_res)*limit_rem)]] <- 0
    }
    
  }
  
  if (any(grepl("isoforest", methods))) {
    
    isof_m <- isotree::isolation.forest(
      dists$dist,
      ndim = 2,
      ntrees = 500,
      nthreads = 1,
      penalize_range = FALSE,
      prob_pick_pooled_gain = 0,
      prob_pick_avg_gain = 0)
    
    isof_vals <- predict(isof_m, dists$dist)
    
    to <- order(isof_vals, decreasing = T)
    
    isof_vals <- ifelse(isof_vals >= isoforest_th, 1, 0)
    isof_vals[-to[1:ceiling(length(isof_vals)*limit_rem)]] <- 0
    
    isoforest_res <- isof_vals
  }
  
  results <- matrix(nrow = nrow(pts), ncol = length(methods))
  colnames(results) <- methods
  
  for (i in 1:length(methods)) {
    results[,i] <- eval(parse(text = paste0(methods[i], "_res")))
  }
  
  return(results)
}

#' @export
.iqr_out <- function(values, threshold, id = NA) {
  
  q1 <- quantile(values, 0.25, na.rm = T)
  q3 <- quantile(values, 0.75, na.rm = T)
  iqr <- IQR(values, na.rm = T)
  
  lower_bound <- q1 - (threshold * iqr)
  upper_bound <- q3 + (threshold * iqr)
  
  values[values > upper_bound] <- id
  values[values < lower_bound] <- id
  
  return(values)
  
}

#' @export
.mad_out <- function(values, threshold, id = NA) {
  
  mad <- mad(values, na.rm = T)
  med_val <- median(values, na.rm = T)
  
  lower_bound <- med_val - (threshold * mad)
  upper_bound <- med_val + (threshold * mad)
  
  values[values > upper_bound] <- id
  values[values < lower_bound] <- id
  
  return(values)
  
}
