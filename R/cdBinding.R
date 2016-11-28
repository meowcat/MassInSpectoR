chemdoodle_sketcherInput <- function(inputId, label = "Get SMILES") {
  tagList(
    singleton(includeScript(
      system.file("js/cdBinding.js", package=.packageName))),
    span(class="chemdoodle_sketcher",
         chemdoodle_sketcher(),
         
         tags$button(id = inputId, class="chemdoodle_sketcher btn btn-default", type="button", label)
    )
  )
}


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
