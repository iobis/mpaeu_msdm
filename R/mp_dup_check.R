#' Remove duplicate datasets
#'
#' @description
#' Removes duplicate datasets based on geohash, species, and year, 
#' and calculating cosine similarities between datasets. This function is
#' a wrap around the code developed by Pieter Provoost.
#' 
#' Details on the method are available at https://iobis.github.io/notebook-duplicates/
#' 
#' @author Provoost, Pieter
#'
#' @param data_a A named `list` of datasets which are to be compared or a `data.frame`.
#' Should have at least the collumns ScientificName, decimalLongitude,
#' decimalLatitude and year.
#' @param data_b (data.frame) If `data_a` is a `data.frame`, then `data_b` should be supplied.
#' @param exclude Should the function automatically exclude the duplicate records? If `FALSE`
#' a list with the 2 datasets and 1 additional `data.frame` with duplicate records will be returned.
#' @param as_single If `TRUE`, the function will merge the final datasets. Only works when exclude is `TRUE`.
#' @param gh_precision precision parameter for the geohash.
#' (see details at [geohashTools::gh_encode()]).
#' @param sim_threshold similarity threshold (0 to 1) at which values are suspect and should be removed or tagged.
#'
#' @return A `data.frame` with cleaned occurrence records or a 
#' `list` with cleaned datasets or a `list` with duplicate records.
#' @export
#'
#' @examples
#' 
#' 
mp_dup_check <- function(data_a,
                         data_b = NULL,
                         exclude = T,
                         as_single = T,
                         gh_precision = 2,
                         sim_threshold = 0.85) {
  
  if (is.null(data_b)) {
    if (is.list(data_a)) {
      data_a <- lapply(1:length(data_a), function(x){
        data_a[[x]]$dataset_id <- names(data_a)[x]
        data_a[[x]] %>%
          select(scientificName, decimalLongitude, decimalLatitude, year, dataset_id) %>%
          mutate(unique_id = paste(dataset_id, 1:nrow(.), sep = "_")) %>%
          filter(!is.na(year)) %>%
          mutate(year = as.integer(year))
      })
      
      data_a <- bind_rows(data_a)
    } else {
      stop("Second dataset was not supplied. data_a should be a list of datasets or data_b should be supplied.")
    }
  } else {
    data_a %>%
      mutate(dataset_id = "dataset_a") %>%
      select(scientificName, decimalLongitude, decimalLatitude, year, dataset_id) %>%
      mutate(unique_id = paste(dataset_id, 1:nrow(.), sep = ""))
    data_b %>%
      mutate(dataset_id = "dataset_a") %>%
      select(scientificName, decimalLongitude, decimalLatitude, year, dataset_id) %>%
      mutate(unique_id = paste(dataset_id, 1:nrow(.), sep = ""))
    data_a <- bind_rows(data_a, data_b)
  }
  
  # Filter and get geohash
  stats <- data_a %>%
    filter(!is.na(year) & decimalLatitude < 90 & decimalLongitude < 180) %>%
    mutate(geohash = gh_encode(decimalLatitude, decimalLongitude, gh_precision)) %>%
    mutate(cell = factor(paste(geohash,
                               tolower(gsub(" ", "_", scientificName)),
                               year, sep = "_"))) %>%
    group_by(dataset_id, cell) %>%
    summarize(records = n())
  
  n_cells <- length(levels(stats$cell))
  dataset_ids <- unique(data_a$unique_id)
  
  vectors <- list()
  
  for (id in dataset_ids) {
    message(id)
    vector <- rep(0, n_cells)
    dataset <- stats %>%
      filter(dataset_id == id)
    for (i in 1:nrow(dataset)) {
      vector[as.numeric(dataset$cell[i])] <- dataset$records[i]
    }
    vectors[[id]] <- as(vector, "sparseVector")
  }
  
  write("x y similarity", file = "similarity.txt", append = FALSE)
  
  parallel::mclapply(1:length(dataset_ids), function(i) {
    dataset_x <- dataset_ids[i]
    x <- as.vector(vectors[[dataset_x]])
    for (j in (i + 1):length(dataset_ids)) {
      dataset_y <- dataset_ids[j]
      y <- as.vector(vectors[[dataset_y]])
      similarity <- coop::cosine(x, y)
      line <- paste(dataset_x, dataset_y, format(similarity, scientific = FALSE))
      write(line, file = "similarity.txt", append = TRUE)
    }
  }, mc.cores = 6)
  
  similarity <- fread("similarity.txt", sep = " ", header = TRUE)
  
  datasets <- robis::dataset() %>%
    tidyr::unnest(statistics) %>%
    select(id, url, title, records = Occurrence)
  
  suspect <- similarity %>%
    filter(similarity > 0.85) %>%
    left_join(datasets, by = c("x" = "id")) %>%
    left_join(datasets, by = c("y" = "id"), suffix = c("_x", "_y")) %>%
    arrange(desc(similarity), desc(records_x + records_y)) %>%
    as_tibble()
  
}