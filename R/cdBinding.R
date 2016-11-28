chemdoodle_sketcherInput <- function(inputId, label = "Get SMILES") {
  tagList(
    singleton(includeScript("cdBinding.js")),
    span(class="chemdoodle_sketcher",
         chemdoodle_sketcher(),
         
         tags$button(id = inputId, class="chemdoodle_sketcher btn btn-default", type="button", label)
    )
  )
}