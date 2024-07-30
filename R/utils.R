############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# July of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################### Utilities ######################################

#' Retrieve the most recent file
#' 
#' @param path path to search for the file
#' @param name name of the file to be searched (uses [grepl()])
#' @param inverse if `TRUE` return the oldest one
#' @param by_date if `TRUE` uses the last modified date instead
#' 
#' @details 
#' In this project, many files are named with the date on the YMD format, like:
#' 'all_splist_20240723.csv'. If `by_date` is `FALSE`, the function will simply
#' order the files (because is alphabetic, the most recent date will come first)
#' If the target file is not in that format, then it is better to use `by_date`
#' 
#' @returns file path
#' @export
#' 
#' @examples
#' \dontrun{
#' recent_file("data", "all_splist")
#' }
recent_file <- function(path, name, inverse = FALSE, by_date = FALSE) {
  f <- list.files(path, full.names = T)
  f_sel <- f[grepl(name, f)]
  
  if(length(f_sel) < 1) {
    cli::cli_alert_warning("No file found.")
    fr <- NULL
  } else if (length(f_sel) == 1) {
    fr <- f_sel
  } else {
    if(by_date) {
      f_sel_info <- fs::file_info(f_sel)
      f_sel_info <- f_sel_info[order(f_sel_info$modification_time,
                                     decreasing = !inverse),]
      fr <- f_sel_info$path[1]
    } else {
      fr <- f_sel[order(f_sel, decreasing = !inverse)][1]
    }
  }
  return(fr)
}