ketcherInput <- function(inputId, label = "Get SMILES") {
  tagList(
    singleton(includeScript("inst/js/ketcherBinding.js")),
    span(class="ketcher",
         tag("iframe", list(id="ketcherFrame", name="ketcherFrame", 
                            style="min-width:410px;min-height:510px;width:100%;border-style:none",
                            src="ketcher/ketcher.html", scrolling="no",seamless="seamless"
         )
         ),
         tags$button(id = inputId, class="ketcher btn btn-default", type="button", label)
    )
  )
}


# Given a vector or list, drop all the NULL items in it
dropNulls <- function(x) {
  x[!vapply(x, is.null, FUN.VALUE=logical(1))]
}

# Send an update message to a URL input on the client.
# This update message can change the value and/or label.
updateKetcherInput <- function(session, inputId,
                               label = NULL, value = NULL) {
  message <- dropNulls(list(label = label, value = value))
  session$sendInputMessage(inputId, message)
}


# convert sdf to smiles
getSmiles <- function(sdf)
{
  tf <- tempfile(fileext=".sdf")
  writeLines(sdf, tf)
  mol <- load.molecules(tf)[[1]]
  get.smiles(mol)
}


# convert smiles to sdf
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


getMolecule.sh <- function (smiles) 
{
  capture.output(mol <- parse.smiles(smiles)[[1]])
  do.aromaticity(mol)
  #convert.implicit.to.explicit(mol)
  do.aromaticity(mol)
  do.typing(mol)
  do.isotopes(mol)
  return(mol)
}