############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
################################## SDM modules #################################

#' Split single file parquet dataset into multiple files
#'
#' @param local_file path to the local parquet dataset file/folder.
#' @param database_name the name of the original database. Will be used as 
#' a grouping key in the final folder format.
#' @param grouping_key \code{character} indicating the key used for spliting 
#' (e.g. AphiaID, for OBIS datasets)
#' @param sel_keys an optional \code{vector} containing a set of keys for which 
#' you want the split. This is specially relevant if you have a full export and
#' want the split for just a set of species. Note that if keys are substituted, 
#' \code{sel_keys} should reflect the new key, not the older one.
#' @param change_key an optional \code{data.frame} containing two columns (key, new_key) 
#' for changing the final grouping key (see details).
#' @param sel_columns an optional character vector of column names which should be
#' selected for the final files. Should include at least the grouping key.
#' @param run_in_batches if \code{sel_keys} are supplied and this option is \code{TRUE},
#' then the dataset is split in several batches according to the \code{batch_size}. The function
#' loop through the batches to save all groups. This option is strongly recommended if
#' the dataset is huge, as the split by Arrow function can fail in those cases.
#' @param batch_size an optional size for the batch. If it's too big, then the 
#' function may still fail. If too small, then it may take too much time to run. If \code{NULL}, 
#' the default of 100 is used.
#' 
#' @details
#' If you have a single parquet file or a parquet dataset that have a grouping key
#' and want to split in the standard format used in the project, then it's more
#' recommended to use this function then the [obissdm::mp_get_local()]. This function
#' use the base Arrow capability to split the file, so it's really fast.
#' 
#' ## How Files are saved?
#' 
#' Files will be saved on "data/species" (if the folder does not exist on the
#' working directory it will be created) following a hive structure:
#' 
#' - key='grouping key of the species'
#' - date='date being generated'
#' - ftype='type of file, i.e., which database generated'
#' 
#' So, for a species with AphiaID 12345, from OBIS,
#' downloaded on 2023-01-05, a folder would be structured as that:
#' 
#' data/species/key=12345/date=20230105
#' 
#' With the following folder inside:
#' ftype=obis/
#' 
#' With a file called spdata0.parquet.
#' 
#' This structure enables easy indexing and querying with multiple species, specially
#' when working with parquet datasets. Each of those keys will become a filtering feature.
#' For more details see [arrow::open_dataset()]
#' 
#' ## Changing the key
#' 
#' The dataset will be split by the grouping key. However, sometimes you want to consistently
#' save all the files using a single identifier. For example, GBIF files will be saved
#' by specieskey (or other grouping), but you may want to identify all by the AphiaID. 
#' In that case you can supplly the optional \code{change_key} argument which should be
#' a \code{data.frame} with two columns: the original grouping key (name "key") and the
#' equivalent new key (named "new_key"). Files will be renamed according to the new key.
#' In this case, if you want to also select just part of the dataset, the \code{sel_keys}
#' should be based on the "new_key" column, as the filter is performed after the key
#' substitution.
#' 
#'
#' @return saved files
#' @export
#' 
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' \dontrun{
#' split_dataset("local_dataset.parquet", "obis", "AphiaID")
#' }
split_dataset <- function(local_file,
                          database_name,
                          grouping_key,
                          sel_keys = NULL,
                          change_key = NULL,
                          sel_columns = NULL,
                          run_in_batches = FALSE,
                          batch_size = NULL) {
  
  if (!is.null(change_key)) {
    if (!all(c("key", "new_key") %in% colnames(change_key)) | length(colnames(change_key)) != 2) {
      stop("change_key have two columns named key and new_key")
    }
    if (!is.null(sel_keys)) {
      if (!all(sel_keys == change_key$new_key)) {
        stop("change_key[['new_key']] have to be equal (in the same order) thand sel_keys")
      }
    }
  }
  
  if (is.null(sel_columns)) {
    add_sel <- FALSE
  } else {
    sel_columns <- paste(sel_columns, collapse = ",")
    add_sel <- TRUE
  }
  
  if (tools::file_ext(local_file) == "parquet") {
    database <- open_dataset(local_file) %>%
      {if(add_sel) select(., eval(strsplit(sel_columns, ",")[[1]])) else .}
  } else {
    database <- open_csv_dataset(local_file) %>%
      {if(add_sel) select(., eval(strsplit(sel_columns, ",")[[1]])) else .}
  }
  
  if (!is.null(sel_keys)) {
    
    if (run_in_batches) {
      cli::cli_alert_info("Running in batches...")
      if (is.null(batch_size)) {
        batch_size <- 100
      }
      cuts <- seq(batch_size, length(sel_keys), by = 100)
      if (max(cuts) != length(sel_keys)) {
        cuts <- c(cuts, length(sel_keys))
      }
    } else {
      cuts <- length(sel_keys)
    }
    
    st <- 1
    
    size_init <- length(list.files("data/species", recursive = T))
    
    for (z in 1:length(cuts)) {
      
      size_init_l <- length(list.files("data/species", recursive = T))
      
      sel_keys_k <- sel_keys[st:cuts[z]]
      
      if (!is.null(change_key)) {
        
        change_key_k <- change_key[st:cuts[z],]
        
        eval(parse(text = glue::glue(
          'database %>%
    mutate(key = [grouping_key]) %>%
    mutate(date = format(Sys.Date(), "%Y%m%d")) %>%
    mutate(ftype = "[database_name]") %>%
    left_join(change_key_k, by = "key") %>%
    mutate(key = new_key) %>%
    select(-new_key) %>%
    filter(key %in% sel_keys_k) %>%
    group_by(key, date, ftype) %>%
    write_dataset("data/species/",
                  format = "parquet",
                  basename_template = "spdata{i}.parquet")',
          .open = "[", .close = "]"
        )))
      } else {
        eval(parse(text = glue::glue(
          'database %>%
    mutate(key = [grouping_key]) %>%
    filter(key %in% sel_keys_k) %>%
    mutate(date = format(Sys.Date(), "%Y%m%d")) %>%
    mutate(ftype = "[database_name]") %>%
    group_by(key, date, ftype) %>%
    write_dataset("data/species/",
                  format = "parquet",
                  basename_template = "spdata{i}.parquet")',
          .open = "[", .close = "]"
        )))
      }
      
      collect <- length(list.files("data/species", recursive = T))
      batch_info <- paste("Batch", z, "of", paste0(length(cuts), "."))
      cat(batch_info, collect-size_init_l, "new files. Total =",  collect-size_init, "\n")
      st <- (cuts[z] + 1)
    }
    
  } else {
    if (!is.null(change_key)) {
      eval(parse(text = glue::glue(
        'database %>%
    mutate(key = [grouping_key]) %>%
    mutate(date = format(Sys.Date(), "%Y%m%d")) %>%
    mutate(ftype = "[database_name]") %>%
    left_join(change_key, by = "key") %>%
    mutate(key = new_key) %>%
    select(-new_key) %>%
    group_by(key, date, ftype) %>%
    write_dataset("data/species/",
                  format = "parquet",
                  basename_template = "spdata{i}.parquet")',
        .open = "[", .close = "]"
      )))
    } else {
      eval(parse(text = glue::glue(
        'database %>%
    mutate(key = [grouping_key]) %>%
    mutate(date = format(Sys.Date(), "%Y%m%d")) %>%
    mutate(ftype = "[database_name]") %>%
    group_by(key, date, ftype) %>%
    write_dataset("data/species/",
                  format = "parquet",
                  basename_template = "spdata{i}.parquet")',
        .open = "[", .close = "]"
      )))
    }
  }
  
  return(invisible(NULL))
}
