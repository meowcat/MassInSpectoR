library(shiny)
library(shinyjs)
library(MSnbase)
library(RMassBank)
library(rcdk)
library(splashR)

ketcherDir <- "C:/Daten/R scripts/KetcheR/inst/html"
jsDir <- "C:/Daten/R scripts/KetcheR/inst/js"

addResourcePath("ketcher", ketcherDir)
addResourcePath("js", jsDir)
addResourcePath("icons", "C:/Daten/R scripts/KetcheR/inst/html/icons")

cfmPath <- "C:/Software/cfm-id-2.0_win32/cfm-predict.exe"
cfmSettings <- "0.001 C:/Software/cfm-id-2.0_win32/param_output.log C:/Software/cfm-id-2.0_win32/param_config.txt 1"

server <- function(input, output, session)
{
  #output$smiles <- reactive({input$smi})
  container <- reactiveValues(
    w = newMsmsWorkspace()
  )
  
  observe({
    shiny::validate(
      need(!is.null(input$smi), "No molecule set"),
      need(input$smi != "", "No molecule set")
    )
    
    #updateTextInput(session, "smiles", value=getSmiles(input$smi))
    updateTextInput(session, "smiles", value=input$smi)
  })
  
  smiles <- reactive({
    input$smiles
  })
  
  observe({
    if(input$setSmiles == 0)
      return()
    isolate({
        updateKetcherInput(session, "smi", value=getSdf(input$smiles))        
      })
    
  })
  
  mz <- reactive({
    if(input$process == 0)
      return()
    isolate(
      findMass(input$smiles) + findMz.formula("", mode="pH")$mzCenter
    )
  })
  
  output$mz <- renderText(mz())
  
  observe({
    if(input$shift == 0)
      return()
    isolate({
      updateNumericInput(session, "shiftTo", value=mz())
    })
    
    
  })
  
  # The three compound selection fields: keep up to date when new compounds are added
  updateCompounds <- function(id = c("1", "2", "V"))
  {
    isolate({
      names(container$w@spectra) -> nms
      print(container$w@spectra)
      if("V" %in% id)
        updateSelectInput(session, "compoundV", choices=nms, selected = input$compoundV)
      if("1" %in% id)
        updateSelectInput(session, "compound1", choices=nms, selected = input$compound1)
      if("2" %in% id)
        updateSelectInput(session, "compound2", choices=nms, selected = input$compound2)
    })
  }
  
  
  # The three input index fields: keep up to date when new compounds are selected
  observe({
    print(input$compoundV)
    rmbCpd <- container$w@spectra[[input$compoundV]]
    shiny::validate(
      need(!is.null(rmbCpd), label="Compound")
    )
    updateSelectInput(session, "inputIndexV", choices=as.character(seq_len(length(rmbCpd@children))))
  })
  observe({
    rmbCpd <- container$w@spectra[[input$compound1]]
    shiny::validate(
      need(!is.null(rmbCpd), label="Compound")
    )
    updateSelectInput(session, "inputIndex1", choices=as.character(seq_len(length(rmbCpd@children))))
  })
  observe({
    rmbCpd <- container$w@spectra[[input$compound2]]
    shiny::validate(
      need(!is.null(rmbCpd), label="Compound")
    )
    updateSelectInput(session, "inputIndex2", choices=as.character(seq_len(length(rmbCpd@children))))
  })
  
  # Simulate with CFM-ID
  observe(
    {
      if(input$process == 0)
        return()
      isolate({
        progress <- shiny::Progress$new(session, min=1, max=15)
        on.exit(progress$close())
        progress$set(message="Processing...")
        tmpOut <- tempfile(fileext=".mgf")
        smi <- input$smiles
        cl <- paste(
          cfmPath,
          smi,
          cfmSettings,
          tmpOut,
          sep=" "
        )
        system(cl)
        specAll <- readMgfData(tmpOut)
        spectra <- lapply(rownames(specAll@featureData), function(name) specAll@assayData[[name]])
        rmbCpd <- new("RmbSpectraSet")
        rmbCpd@children <- as(lapply(
          spectra, function(sp) new("RmbSpectrum2", sp, versions=c("RmbSpectrum2" = "0.1.1"))
        ), "SimpleList")
        name <- input$nameSim
        rmbCpd@mz <- mz()
        if(name=="")
          name <- smi
        rmbCpd@name <-name
        container$w@spectra[[name]] <- filterMerge(rmbCpd)
        updateCompounds()
       })
    }
  )
  
  # spectum upload
  observe({
    fi <- input$refspec
    
    isolate({
      shiny::validate(need(!is.null(fi), message="File not available"))
      path <- fi$datapath
      
      specAll <- readMgfData(path)
      spectra <- lapply(rownames(specAll@featureData), function(name) specAll@assayData[[name]])
      spectra <- spectra[
        unlist(lapply(spectra, function(sp) sp@precursorMz > 0))
      ]
      rmbCpd <- new("RmbSpectraSet")
      rmbCpd@children <- as(lapply(
        spectra, function(sp) new("RmbSpectrum2", sp, versions=c("RmbSpectrum2" = "0.1.1"))
      ), "SimpleList")
      
      rmbCpd@mz <- rmbCpd@children[[1]]@precursorMz
      
      name <- input$nameUpload
      if(name=="")
        name <- fi$name
      rmbCpd@name <-name
      container$w@spectra[[name]] <- filterMerge(rmbCpd)
      print(container$w@spectra)
      updateCompounds()
    })
    
  })
  
  # spectrum DB import
  observe({
    if(input$import == 0) return()
    
    isolate({
      id <- as.numeric(input$cpdID)
      cpd <- getCompound(db, id)
      container$w@spectra[[as.character(id)]] <- filterMerge(cpd)
      updateCompounds()
    })
  })
  
  
  cpdV <- reactive(
    {
      input$compoundV
      cpd <- isolate({
        shiny::validate(
          need(!is.null( container$w@spectra[[input$compoundV]]), label="Compound")
        )
        container$w@spectra[[input$compoundV]]
      })
      
      shiftTo <- ifelse(input$shiftToV, cpd@mz, 0)
      #if(input$shift1 > 0)
      cpd@children <- cpd@children + (input$shiftV - shiftTo)
      container$newMzV <- cpd@mz +(input$shiftV - shiftTo)

      cpd
    }
  )
      
  cpd1 <- reactive(
    {
      input$compound1
      rp <- input$removeParent
      cpd <- isolate({
        shiny::validate(
          need(!is.null( container$w@spectra[[input$compound1]]), label="Compound")
        )
        container$w@spectra[[input$compound1]]
      })
      # Apply modifications: Remove parent, or shift
      if(rp)
      {
        
        m1 <- cpd@mz - input$tolerance
        m2 <- cpd@mz + input$tolerance
#         cpd@children <- as(lapply(cpd@children, function(chi)
#           {
#           property(chi, "m1",  addNew=TRUE, "numeric") <- m1
#           property(chi, "m2",  addNew=TRUE, "numeric") <- m2
#           chi
#         }), "SimpleList")
        cpd@children <- selectPeaks(cpd@children, (mz < m1) | (mz > m2))
      }
      shiftTo <- ifelse(input$shiftTo1, cpd@mz, 0)
      #if(input$shift1 > 0)
        cpd@children <- cpd@children + (input$shift1 - shiftTo)
      container$newMz1 <- cpd@mz +(input$shift1 - shiftTo)
        
      cpd
    }
  )
  cpd2 <- reactive(
    {
      input$compound2
      rp <- input$removeParent
      cpd <- isolate({
        shiny::validate(
          need(!is.null( container$w@spectra[[input$compound2]]), label="Compound")
        )
        container$w@spectra[[input$compound2]]
      })
      # Apply modifications: Remove parent, or shift
      if(rp)
      {
        m1 <- cpd@mz - input$tolerance
        m2 <- cpd@mz + input$tolerance
#         cpd@children <- as(lapply(cpd@children, function(chi)
#         {
#           property(chi, "m1",  addNew=TRUE, "numeric") <- m1
#           property(chi, "m2",  addNew=TRUE, "numeric") <- m2
#           chi
#         }), "SimpleList")
        cpd@children <- selectPeaks(cpd@children, (mz < m1) | (mz > m2))
      }
      shiftTo <- ifelse(input$shiftTo2, cpd@mz, 0)
      #if(input$shift2 > 0)
      cpd@children <- cpd@children + (input$shift2 - shiftTo)
      container$newMz2 <- cpd@mz + (input$shift2 - shiftTo)
      cpd
    }
  )
  
  specV <- reactive(
    {
      input$inputIndexV
      cpdV()@children[[as.numeric(input$inputIndexV)]]
      
    }
  )
  
  spec1 <- reactive(
    {
      input$inputIndex1
      cpd1()@children[[as.numeric(input$inputIndex1)]]
      
    }
  )
  spec2 <- reactive(
    {
      input$inputIndex2
      cpd2()@children[[as.numeric(input$inputIndex2)]]
      
    }
  )
  
  output$mz1 <- renderText(
    {
      paste(cpd1()@mz, "->", container$newMz1)
    }
  )
  output$mz2 <- renderText(
    {
      paste(cpd2()@mz, "->", container$newMz2)
    }
  )
  output$mzV <- renderText(
    {
      paste(cpdV()@mz, "->", container$newMzV)
    }
  )

  # display the tables
  specTable <- reactive({
    spec <- specV()
    shiny::validate(need(!is.null(spec), "Data is not yet available"))
    isolate(
      {
        getData(spec)
      }
    )
  })
  
  output$vTab <- renderDataTable(specTable())
  output$vSpec <- renderPlot({
    #t1 <- specTable()
    #t2 <- inputSpecTable
    plot(intensity ~ mz, data=specTable(), type="h")
  })
  

  simMatrix <- reactive({
    if(input$calcMatrix == 0) return()
    
    isolate({
      progress <- shiny::Progress$new(session, min=1, max=15)
      on.exit(progress$close())
      progress$set(message="Processing...")
      
      s1 <- cpd1()@children
      s2 <- cpd2()@children
      tabs1 <- lapply(s1, getData)
      tabs2 <- lapply(s2, getData)
      as.data.frame(
        outer(seq_len(length(tabs1)), seq_len(length(tabs2)), Vectorize(function(i1, i2)
        {
          t1 <- tabs1[[i1]]
          t2 <- tabs2[[i2]]
          round(OrgMSsim(t1, t2, t = input$tolerance, b = input$cutoff)[[2]],2)
        }
        , c("i1", "i2"))))
    })
  })
  
  observe({simMatrix()})
  
  output$simMatrix <- renderDataTable({simMatrix()})
  
  # Similarity
  similarity <- reactive({
    shiny::validate(
      need(cpd1(), "no spec 1"),
      need(cpd2(), "no spec 2")
    )
    outFile <- tempfile(fileext=".png")
    spec1 <- getData(spec1())
    spec2 <- getData(spec2())

    
    # calculate similarity and capture the output in a PNG file
    png(outFile, width=600, height=600)
    par(mar=c(2,1,1,1)+0.1)
    res <- OrgMSsim(spec1, spec2, t = input$tolerance, b = input$cutoff, top.label = cpd1()@name, bottom.label = cpd2()@name,
                    xlim=range(spec1$mz, spec2$mz))
    dev.off()
    res[[3]] <- outFile
    res
  })
  
  shifts <- reactive({
    shiny::validate(
      need(cpd1(), "no spec 1"),
      need(cpd2(), "no spec 2")
    )
    shiftMatch <- stepwiseShiftMatch(spec1(), spec2(), cutoff=0.05)
    shiftMatch$ratio <- shiftMatch$addedIntensity / shiftMatch$intensity
    shiftMatch[order(shiftMatch$totalIntensity, decreasing = TRUE),,drop=FALSE]
  })
  
  
  searchRes <- reactive({
    if(input$searchV == 0)
      return(data.frame())
    isolate({
      specV <- normalize(specV(), slot="relint", scale=1, precision=3)
      spec <- getData(specV)
      hits <- lapply(seq_len(nrow(spec)), function(n)
      {
        mz <- spec[n, "mz"]
        hits.mz <- which(abs(db$FragmentMZ - mz) < 0.001)
        ret <- db[hits.mz, c("CompoundID", "CollisionEnergy", "Resolution", "FragmentRelativeIntensity")]
        hit <- cbind(mz = rep(mz, nrow(ret)), QueryRI = rep(spec[n, "relint"], nrow(ret)), ret)
        hit
      })
      hits.all <- do.call(rbind, hits)
      hits.all$hitSpectrum <- interaction(hits.all$CompoundID , hits.all$CollisionEnergy)
      hits.all$CRI <- hits.all$FragmentRelativeIntensity * hits.all$QueryRI
      hits.sum <- unclass(as.matrix(xtabs(CRI ~ hitSpectrum, data=hits.all)))
      #hits.sum <- cbind(total = rowSums(hits.sum), hits.sum)
      hits.sum <- round(hits.sum, 2)
      hits.sum <- cbind(hitSpectrum.c = row.names(hits.sum), as.data.frame(hits.sum))
      hits.sum$cpd <- sapply(strsplit(as.character(hits.sum[,"hitSpectrum.c"]), "\\."), "[", 1)
      hits.sum <- merge(hits.sum, compoundTable, all.x=TRUE)
      #hits.sum <- hits.sum[order(hits.sum[,"total"], decreasing=TRUE), ]
    })
  })
  
  searchSummary <- reactive({
    if(is.null(searchRes()))
      return(data.frame())
    if(nrow(searchRes() == 0))
      return(data.frame())
    return(data.frame())
    res <- searchRes()
    summ <- unclass(as.matrix(xtabs( ~ CompoundID + CollisionEnergy, data=res)))
    summ <- cbind(as.data.frame(summ), total = rowSums(summ))
    summ <- summ[order(summ[,"total"], decreasing=TRUE),,drop=FALSE]
    summ
  })
  
  output$simPlot <- renderImage({list(src=similarity()[[3]])})
  output$similarity <- renderText({similarity()[[2]]})
  output$simTab <- renderDataTable({similarity()[[1]]})
  
  output$shiftTab <- renderDataTable({shifts()})
  
  output$searchRes <- renderDataTable({searchRes()})
  output$searchSummary <- renderDataTable({searchSummary()})
  
}

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      #fileInput("refspec", "Comparison")
      
      
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Upload",
                 fileInput("refspec", "Comparison"),
                 textInput("nameUpload", "Name", "")
                 ),
        tabPanel("MassBankDB",
                 selectInput("cpdID", "compound ID", as.character(compounds)),
                 actionButton("import", "Import"),
                 #textInput("nameUpload", "Name", "")
                 tabsetPanel(
                   tabPanel("Search results", 
                            dataTableOutput("searchRes")),
                   tabPanel("Search summary",
                            dataTableOutput("searchSummary"))
                 )
        ),
        tabPanel("Simulate",
                 includeScript("inst/js/ketcherBinding.js"),
                 #numericInput('n', 'Number of obs', 100),
                 ketcherInput("smi"),
                 textInput("smiles", "SMILES", ""),
                 actionButton("setSmiles", "Set SMILES"),
                 actionButton("process", "Process"),
                 textOutput("mz"),
                 textInput("nameSim", "Name", ""),
                 actionButton("shift", "Shift input spectra to mz")
        ), # tabPanel input
        tabPanel("spectra",
                 selectInput("compoundV", "Compound",c()),
                 selectInput("inputIndexV", "Index",c()),
                 "m/z:", textOutput("mzV"),
                 checkboxInput("shiftToV", "Shift to absolute mz", FALSE),
                 numericInput("shiftV", "mz shift", 0),
                 actionButton("searchV", "Search in MassBank DB"),
                 #textOutput("specOut"),
                 plotOutput("vSpec"),
                 dataTableOutput("vTab")
        ), # tabPanel CFM-ID
#         tabPanel("Input spectra",
#                  selectInput("compound2", "Compound",c()),
#                  #textOutput("specOut"),
#                  plotOutput("inputSpec"),
#                  dataTableOutput("inputTab")
#         ), # tabPanel Input Spectra
        tabPanel("Comparison",
                 selectInput("compound1", "Compound 1",c()),
                 selectInput("inputIndex1", "Index", c()),
                 "m/z:", textOutput("mz1"),
                 checkboxInput("shiftTo1", "Shift to absolute mz", FALSE),
                 numericInput("shift1", "mz shift", 0),
                 
                 selectInput("compound2", "Compound 2",c()),
                 selectInput("inputIndex2", "Index", c()),
                 "m/z:", textOutput("mz2"),
                 checkboxInput("shiftTo2", "Shift to absolute mz", FALSE),
                 numericInput("shift2", "mz shift", 0),
                 
                 checkboxInput("removeParent", "Remove parent peak from MS2 for comparison", FALSE),
                 numericInput("tolerance", "Tolerance (mz)", 0.001),
                 numericInput("cutoff", "Relative cutoff %", 5),
                 actionButton("calcMatrix", "calculate matrix")
                 
        ),
        tabPanel("Similarity matrix",
                 dataTableOutput("simMatrix")
        ),
        tabPanel("Similarity",
                 textOutput("similarity"),
                 imageOutput("simPlot", width=600, height=600),
                 tabsetPanel(
                   tabPanel("Aligned spectra", 
                            dataTableOutput("simTab")),
                   tabPanel("Shift summary", 
                            dataTableOutput("shiftTab"))
                 )
                 )
        
      ) # tabsetPanel
    ) # mainPanel
  ) # sidebarLayout
) # ui
                

shinyApp(ui=ui, server=server)

