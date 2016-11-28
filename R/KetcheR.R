ketcherInput <- function(inputId, label = "get smiles") {
  tagList(
    singleton(includeScript("inst/js/ketcherBinding.js")),
    span(class="ketcher",
      tag("iframe", list(id="ketcherFrame", name="ketcherFrame", 
                         style="min-width:310px;min-height:410px;width:100%;border-style:none",
                         src="ketcher/ketcher.html", scrolling="no"
                         )
      ),
      tags$button(id = inputId, class="ketcher btn btn-default", type="button", label)
    )
  )
}

# Send an update message to a URL input on the client.
# This update message can change the value and/or label.
updateKetcherInput <- function(session, inputId,
                           label = NULL, value = NULL) {
  session$sendInputMessage(inputId, message)
}


# convert smiles to sdf
getSdf <- function(smiles)
{
  mol <- getMolecule(smiles)
  tf <- tempfile(fileext=".sdf")
  write.molecules(mol, tf)
  sdf <- readLines(tf)
  paste(sdf, collapse="\r\n")
}