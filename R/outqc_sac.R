#' Calculates spatial autocorrelation between points
#' 
#' This function calculates the spatial autocorrelation between a set of points and,
#' optionally, prune the records based on the minimum non-significant distance.
#'
#' @param pts a \code{data.frame} with the coordinates to the points. Should
#'   have at least two columns indicating the longitude and latitude (in that
#'   order)
#' @param env_layers environmental layers in SpatRaster format
#' @param autocor_classdist the steps (in km) to calculate the autocorrelation
#' @param autocor_maxdist the maximum distance (in km) to which calculate the
#'   autocorrelation
#' @param autocor_signif the significance level
#' @param return_pruned if \code{TRUE}, it will return the pruned records (if
#'   the first non-correlated distance is higher or equal than \code{prune_threshold},
#'   otherwise the original dataset is returned)
#' @param prune_threshold the minimum distance (in km) to use for prunning the
#'   records
#' @param plot_result plot a graphical representation of the results
#' @param verbose print control messages
#'
#' @return the first non-correlated distance (in km) or prunned records
#' @export
#'
#' @details For each distance class, a linear model tests the effect of
#' correlation with geographic distance. This procedure finds the minimum
#' non-significant autocorrelated distance which is optionally used to prune the
#' occurrence records
#' 
#' Note that if the minimum distance between points is higher than the maximum
#' distance being considered (\code{autocor_maxdist}), the function will abort
#' as it's impossible to assess the SAC. In that case, it makes no sense in prunning
#' the points.
#' 
#' If \code{return_pruned=TRUE}, the function [spThin::thin()] is used.
#' 
#' @references
#'   [https://jorgemfa.medium.com/reducing-spatial-autocorrelation-in-species-distribution-models-fe84d4269cee](https://jorgemfa.medium.com/reducing-spatial-autocorrelation-in-species-distribution-models-fe84d4269cee)
#' @references Boavida, J., Assis, J., Silva, I. et al. (2016) Overlooked
#'   habitat of a vulnerable gorgonian revealed in the Mediterranean and Eastern
#'   Atlantic by ecological niche modelling. Scientific Reports. 6, 36460.
#' 
#' @author Original version by Jorge Assis (jmassis@ualg.pt), adapted to use terra functions.
#'
#' @examples
#' \dontrun{
#' new_pts <- outqc_sac(pts, env)
#' }
outqc_sac <- function(pts,
                      env_layers,
                      autocor_classdist = 5,
                      autocor_maxdist = 200,
                      autocor_signif = 0.05,
                      return_pruned = TRUE,
                      prune_threshold = 10,
                      plot_result = TRUE,
                      verbose = TRUE) {
  
  pts_original <- pts
  
  if(nrow(pts) > 1000) { 
    if (verbose) cli::cli_alert_info("More than 1000 records. Using a maximum of 1000 random records.")
    pts <- pts[sample(1:nrow(pts), 1000, replace = FALSE), ]
  }
  
  pts <- pts[!duplicated(pts),]
  
  presences_env <- data.frame(terra::extract(env_layers, pts[,1:2], ID = F))
  to_keep <- apply(presences_env, 1, sum)
  presences_env <- presences_env[!is.na(to_keep),]
  pts <- pts[!is.na(to_keep),]
  
  pts_vec <- terra::vect(as.matrix(pts[,1:2]), crs = terra::crs(env_layers))
  
  space <- terra::distance(pts_vec, pts_vec)
  space <- space/1000
  teste <- sp::spDists(as.matrix(pts[,1:2]),as.matrix(pts[,1:2]),longlat=TRUE)
  
  if (min(space[space > 0]) > autocor_maxdist) {
    cli::cli_abort("Impossible to comput the SAC, {.var autocor_maxdist} is lower
                   then the minimum distance between points.
                   You probably don't need to prune the points.")
  }
  
  env_space <- ecodist::distance(presences_env, method = "euclidean")
  env_space <- as.matrix(env_space)
  
  n_class <- round(autocor_maxdist / autocor_classdist)
  
  results_mat <- data.frame(classdist_from = seq(0, autocor_maxdist - autocor_classdist, autocor_classdist),
                            classdist_to = seq(autocor_classdist, autocor_maxdist, autocor_classdist),
                            R = NA, pval = NA)
  
  for (i in 1:nrow(results_mat)) {
    
    d1 <- results_mat[i, 1]
    d2 <- results_mat[i, 2]
    
    env_space_d <- as.vector(env_space)
    space_d <- as.vector(space)
    
    remove <- which(space_d < d1 | space_d > d2)
    env_space_d <- env_space_d[-remove]
    space_d <- space_d[-remove]
    
    lm_mod <- try(lm(space_d ~ env_space_d), silent = T)
    
    if (class(lm_mod)[1] != "try-error") {
      f <- summary(lm_mod)$fstatistic
    }
    p <- 0
    
    tryCatch(p <- pf(f[1], f[2], f[3], lower.tail = FALSE), error = function(e) {Error <<- TRUE })
    
    if (class(lm_mod)[1] != "try-error") {
      results_mat[i, 3] <- summary(lm_mod)$adj.r.squared
    }
    results_mat[i, 4] <- p
  }
  
  vec_sign <- as.numeric(results_mat[,4] > autocor_signif)
  if(vec_sign[1] == 1) { vec_sign[1] <- 0}
  
  distance <- round(results_mat[which(results_mat[,4] >= autocor_signif), 2][1])
  if(is.na(distance)) {distance <- autocor_maxdist}
  
  if (plot_result) {
    
    vec_sign_pch <- c(19, 1)
    
    par(mar = c(4.5, 5.5, 4.5, 4.5) , bg = "#FFFFFF")
    plot(results_mat[,2], results_mat[,3], axes = FALSE, xlab = "Distance (Km)", ylab = "")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#F6F6F6")
    lines(results_mat[,2], results_mat[,3], lty = 2, col="#000000", type = "l", xlab = "Distance (Km)", ylab = "")
    
    points(results_mat[, 2], results_mat[, 3], pch = vec_sign_pch[vec_sign + 1],
           col = "#5B5B5B")
    axis(2, las = 2, col = "White", col.ticks = "Black")
    axis(1, las = 0, col = "White", col.ticks = "Black")
    box()
    title(ylab = "Correlation (R)", mgp = c(4, 1, 0))
    abline(h = 0, lty = 1)
    abline(v = distance, lty = 3, col = "grey70")
    text(x = (distance + 0.2), y = round(max(results_mat[,3], na.rm = T) - 0.01, 3),
         labels = paste("\U2192 First non-correlated distance (in km):", distance),
         pos = 4, col = "grey70")
    legend("topright", legend = c("Significant", "Non significant"), pch = c(19, 1), col = "#5B5B5B")
    
  }

  if (verbose) cli::cli_alert_success("First non-correlated distance: {distance} km")
  
  if (return_pruned) {
    
    if (distance >= prune_threshold) {
      if (ncol(pts_original) < 3) {
        pts_original$species <- "species1"
      }
      
      pruned <- spThin::thin(pts_original,
                             lat.col = colnames(pts_original)[2],
                             long.col = colnames(pts_original)[1],
                             spec.col = colnames(pts_original)[ncol(pts_original)],
                             thin.par = distance,
                             locs.thinned.list.return = TRUE,
                             reps = 1,
                             write.files = F,
                             write.log.file = F,
                             verbose = verbose)
      pruned <- pruned[[1]]
      colnames(pruned) <- c("decimalLongitude", "decimalLatitude")
      
      if (verbose) cli::cli_alert_success("Input records: {nrow(pts_original)} | Final records: {nrow(pruned)}")
    } else {
      if (verbose) cli::cli_alert_danger("Minimum distance is lower than {.var prune_threshold} ({.val {prune_threshold}km}). Returning original dataset.")
      pruned <- pts_original
    }
    
    return(pruned)
  } else {
    return(distance)
  }
}
