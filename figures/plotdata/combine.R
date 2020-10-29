library( tidyverse )

## Load individual .csv files
c("Fig1C", "Fig2B", "Fig3",
  "Fig4B", "Fig4C", "Fig5A", "Fig5B",
  str_c("Suppl", 1:6)) %>%
    set_names( str_c(., ".csv"), . ) %>%
    map( read_csv, col_types=cols() ) %>%
    openxlsx::write.xlsx( file="plotdata.xlsx" )
