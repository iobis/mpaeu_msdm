.mp_print_occ_get <- function(sp,
                              dbs,
                              full,
                              dups,
                              qc) {
  
  nr_full <- unlist(lapply(full, nrow))
  
  nr_dups <- nrow(dups)
  nr_qc <- nrow(qc)
  
  cli_h1("Data download for {sp}")
  cli_inform("")
  cli_inform("Sources: {dbs}")
  cli_inform("")
  
  fmt <- ansi_columns(
    glue('{names(nr_full)} = {nr_full}'),
    width = 50,
    fill = "rows",
    max_cols=3,
    align = "center",
    sep = "   "
  )
  
  print(boxx(fmt, padding = c(1,0,1,0), header = col_cyan("Occurrence found in each source"),
       footer = col_yellow(glue("Total number of records: {sum(nr_full)}"))))
  
  cli_inform("")
  cli_inform("Number of points retained after duplicate checking: {.strong {col_cyan(nr_dups)}}")
  cli_inform("Number of points retained after quality control: {.strong {col_cyan(nr_qc)}}")
 
  return(invisible(NULL)) 
}