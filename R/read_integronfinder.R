# Read one IntegronFinder 2.x `.integrons` file into a tidy data frame.
# Columns are mapped BY NAME, never by position.
suppressPackageStartupMessages({library(readr); library(dplyr)})

read_integronfinder <- function(path) {
  raw <- readLines(path, warn = FALSE)
  raw <- raw[nchar(trimws(raw)) > 0 & !startsWith(trimws(raw), "#")]
  if (length(raw) < 2) stop("No data rows in ", path)
  df <- readr::read_tsv(paste(raw, collapse = "\n"),
                        show_col_types = FALSE, progress = FALSE)
  expected <- c("ID_integron","ID_replicon","pos_beg","pos_end",
                "strand","type_elt","annotation","model")
  missing <- setdiff(expected, names(df))
  if (length(missing))
    warning("columns not found (continuing): ", paste(missing, collapse = ", "))
  df |>
    mutate(pos_beg  = suppressWarnings(as.numeric(pos_beg)),
           pos_end  = suppressWarnings(as.numeric(pos_end)),
           strand   = suppressWarnings(as.integer(strand)),
           type_elt = tolower(as.character(type_elt)))
}
