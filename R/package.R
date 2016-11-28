.packageName <- "MassInSpectoR"

.packageEnv <- new.env()
.packageEnv$db = data.frame()

# convenience functions:

#' Load MassBank DB table
#'
#' @param path 
#'
#' @return
#' @export
#'
#' @examples
loadMassBankDB <- function(path)
{
  db.fi <- list.files("C:/Daten/Support/Jen - uncleaned spectra/tables", "_pos.csv", full.names = TRUE)
  db.tot <- lapply(db.fi, read.csv)
  db <- do.call(rbind, db.tot)
  rm(db.tot)
  .packageEnv$db <- db
  invisible(db)
}

loadMISettings <- function(file)
{
  o <- yaml.load_file(file)
  l <- list()
  l[[.packageName]] <- o
  do.call(options, l)
}