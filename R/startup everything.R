source('C:/Daten/R scripts/KetcheR/spectrumMerging.R')
source('C:/Daten/R scripts/KetcheR/similarityFunctions.R')
source('C:/Daten/R scripts/KetcheR/shiftmatch.R')
source('C:/Daten/R scripts/KetcheR/sameCEdata_func.R')
#source('C:/Daten/R scripts/KetcheR/R/KetcheR.R')
source('C:/Daten/R scripts/KetcheR/extractFragments.R')
source('C:/Daten/R scripts/KetcheR/cdBinding.R')

setwd("C:/Daten/R scripts/KetcheR/")
#devtools::install_github("zachcp/chemdoodle")
#Sys.unsetenv("JAVA_HOME")
library(chemdoodle)
library(RMassBank)

# convert smiles to sdf
getSdf <- function(smiles)
{
  mol <- getMolecule(smiles)
  tf <- tempfile(fileext=".sdf")
  write.molecules(mol, tf)
  sdf <- readLines(tf)
  paste(sdf, collapse="\r\n")
}
# convert smiles to sdf

# or this?
getSdf <- function(smiles)
{
  #mol <- getMolecule(smiles)
  sdf <- getCactus(smiles, "sdf")
  
  tf <- tempfile(fileext=".sdf")
  writeLines(sdf, tf)
  mol <- load.molecules(tf)[[1]]
  mol <- remove.hydrogens(mol)
  write.molecules(mol, tf)
  sdf <- readLines(tf)
  paste(sdf, collapse="\r\n")
}

# Given a vector or list, drop all the NULL items in it
dropNulls <- function(x) {
  x[!vapply(x, is.null, FUN.VALUE=logical(1))]
}

# Send an update message to a URL input on the client.
# This update message can change the value and/or label.
updateCDInput <- function(session, inputId,
                               label = NULL, value = NULL) {
  message <- dropNulls(list(label = label, value = value))
  session$sendInputMessage(inputId, message)
}

gfpath <- "C:/Software/GenForm/GenForm.exe" 
