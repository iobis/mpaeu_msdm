library(robis)
library(sf)
# Settings
sf_use_s2(FALSE)

st.area <- st_read("data/shapefiles/mpa_europe_starea.shp")

st.area.simp <- st_simplify(st_buffer(st.area, dist = 3), dTolerance = 2)

st_covered_by(st.area, st.area.simp)

sp.onarea <- checklist(geometry = st_as_text(st_geometry(st.area.simp[1,])))

head(sp.onarea)


#### TESTE DOWNLOAD E PRATICA
mp_get_occ <- function(sp,
                       database = c("obis", "gbif"),
                       local_file = NULL,
                       study_area = NULL,
                       range_date = NULL,
                       save_results = T,
                       save_folder = "data/species/",
                       save_format = "parquet",
                       print.summary = F,
                       verbose = F
){
  
  if (verbose) cli_inform(c("i" = "Downloading data for {sp}.",
                            "{symbol$pointer} Data will be downloaded from {database}."))
  
  # Verify if range date is correct
  if (!is.null(range_date)) {
    if (length(range_date) < 2) {
      range_date <- c(range_date, Sys.Date())
    }
    range_date <- as.Date(range_date)
  }
  
  # Get data from each provider
  
  # Start a list to hold results
  occs <- list()
  
  # Assign databse type for local sources
  dbtype <- database
  dbtype[grepl("local", dbtype)] <- "local"
  
  for (i in 1:length(dbtype)) {
    
    if (verbose) cli_progress_step("Downloading from {database[i]}")
    
    # Obtain occurrence data from each provider
    occs[[i]] <- switch (dbtype[i],
                        obis = robis::occurrence(
                          scientificname = sp,
                          startdate = range_date[1],
                          enddate = range_date[2],
                          geometry = study_area
                        ),
                        gbif = rgbif::occ_data(
                          scientificName = sp,
                          geometry = study_area,
                          eventDate = range_date
                        )$data,
                        local = .mp_get_occ_local(
                          local_file[database[i]],
                          scientificName = sp,
                          startdate = range_date[1],
                          enddate = range_date[2],
                          geometry = study_area
                        )
    )
    
    if (verbose) cli_progress_done()
    
  }
  
  # Assign names
  names(occs) <- database
  
  # Check for duplicates
  occs_nodups <- mp_dup_check(occs, exclude = T, as_single = T)
  
  # Perform quality control
  occs_qc <- mp_qc_check(occs_nodups)
  
  if (print.summary) {
    .mp_print_occ_get(sp,
                      dbs = database,
                      full = occs,
                      dups = occs_nodups,
                      qc = occs_qc)
  }
  
  # Save
  if (save_results) {
    id <- worrms::wm_name2id(sp)
    
    if (any(save_format %in% c("parquet", "csv")) & length(save_format) == 1) {
      switch (save_format,
        parquet = write_parquet(occs_qc, glue("{save_folder}spdat_{id}.{save_format}")),
        csv = write_csv_arrow(occs_qc, glue("{save_folder}spdat_{id}.{save_format}"))
      )
    }
    
    return(invisible(NULL))
  } else {
    return(occs_qc)
  }
}

match(
  database,
  c("obis", "gbif", "local*")
)


switch("teste",
       teste = print("FOI"),
       casa = print("MINHA CASA"),
       grepl("local", ))
