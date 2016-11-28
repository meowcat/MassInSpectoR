library(shiny)
library(shinyjs)
library(MSnbase)
library(RMassBank)
library(rcdk)
library(splashR)
library(zoo)

# ketcherDir <- "C:/Daten/R scripts/KetcheR/inst/html"
# jsDir <- "C:/Daten/R scripts/KetcheR/inst/js"
# 
# addResourcePath("ketcher", ketcherDir)
# addResourcePath("js", jsDir)
# addResourcePath("icons", "C:/Daten/R scripts/KetcheR/inst/html/icons")


.Workbench.server <- function(input, output, session)
{
  #output$smiles <- reactive({input$smi})
  container <- reactiveValues(
    w = newMsmsWorkspace()
  )
  
  # for GenForm:
  pool <- reactiveValues()
  pool$cli <- ""
  
  # from ketcher: this observed input$smi
  
  # observe({
  #   shiny::validate(
  #     need(!is.null(input$smi), "No molecule set"),
  #     need(input$smi != "", "No molecule set")
  #   )
  #   
  #   #updateTextInput(session, "smiles", value=getSmiles(input$smi))
  #   updateTextInput(session, "smiles", value=input$smi)
  # })
  
  smiles <- reactive({
    input$smiles
  })
  
  # from ketcher: this was for getting smiles into ketcher
  
  # observe({
  #   if(input$setSmiles == 0)
  #     return()
  #   isolate({
  #       updateKetcherInput(session, "smi", value=getSdf(input$smiles))        
  #     })
  #   
  # })
  
  # Process the chemdoodle sketcher click:
  # read out the smiles
  observe({
                 input$getSmiles
                 if(is.null(input$getSmiles))
                   return()
                 
                 moljson <- input$getSmiles
                 mol <- processChemDoodleJson(moljson)
                 smi <- toSmiles(mol)
                 updateTextInput(session, "smiles", value=smi)
                 
                 })
  
  # Set chemdoodle structure from SMILES:
  observe({
    input$setSmiles
    if(input$setSmiles == 0)
      return()
    isolate({
      smi <- input$smiles
      sdf <- getSdf(smi)
      updateCDInput(session, "getSmiles", value=sdf)
      
    })
    #var mol = ChemDoodle.readMOL(fileContent);
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
      
      # This is an ugly hack. MSP+ will solve this.
      titles <- strsplit(as.character(specAll@featureData$TITLE), ";")
      msLevel <- unlist(lapply(titles, function(title)
        {
        if(title[[1]] == "msLevel 1")
          return(1)
        else
          return(2)
      }))
      
      
      # this generates names X1..X10 for 10 spectra, so they are in the right order...
      # Alternatively the scan no could be used, but since it doesn't correctly
      # extract MS1, this is still not a good solution.
      
      specNames <- paste0("X", seq_len(
        nrow(specAll@featureData)
      ))
      
      spectra <- lapply(specNames, function(name) specAll@assayData[[name]])
      # spectra <- spectra[
      #   unlist(lapply(spectra, function(sp) sp@precursorMz > 0))
      # ]
      rmbCpd <- new("RmbSpectraSet")
      if(!is.na(match(1, msLevel)))
        rmbCpd@parent <- new("Spectrum1", mz=spectra[[match(1, msLevel)]]@mz, intensity=spectra[[match(1, msLevel)]]@intensity)
      else
        rmbCpd@parent <- new("Spectrum1")
       
      rmbCpd@children <- as(lapply(
        spectra[which(msLevel == 2)], function(sp) new("RmbSpectrum2", sp, versions=c("RmbSpectrum2" = "0.1.1"))
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
  
  parentV <- reactive(
    {
      if(is.null(cpdV()))
        return(data.frame(mz=numeric(0), i=numeric(0)))
      if(is.null(cpdV()@parent))
        return(data.frame(mz=numeric(0), i=numeric(0)))
      as.data.frame(cpdV()@parent)
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
  
  ## GenForm processing
  
  # build CLI from user input
  cli <- reactive({
    input[["cli.build"]]
    isolate({
      #"ion=+H ff=C0-30H0-60N0-10O0-6Cl0-5 ppm=5 acc=3 rej=10"
      sprintf("ion=%s ff=%s ppm=%s acc=%s rej=%s %s",
              input[["cli.ion"]],
              input[["cli.ff"]],
              input[["cli.ppm"]],
              input[["cli.acc"]],
              input[["cli.rej"]],
              input[["cli.args"]])
    })
  })
  
  observe({
    d <- cli()
    updateTextInput(session, "cli", value=d)
  })
  r.ms1 <- reactive({
    parentV()
    if(input$m == "") return()
    
    ms1.t <- parentV()[,,drop=FALSE]
    m <- as.numeric(input$m)
    ms1.t <- ms1.t[(as.numeric(ms1.t[,1]) > (m -1)) & (as.numeric(ms1.t[,1]) < (m+ 10)),c(1,2),drop=FALSE]
  })
  
  
  gfresult <- reactive({
    input$processGF
    if(input$processGF > 0)
    {
      

        isolate({
          #try({
          # format the output values
          progress <- shiny::Progress$new(session, min=1, max=15)
          on.exit(progress$close())
          progress$set(message="Processing...")
          cli <- input$cli
          m <- input$m
          ms1.t <- r.ms1()
          ms1.temp <- tempfile("ms1-")
          write.table(ms1.t[,,drop=FALSE], file=ms1.temp, sep="\t", row.names=FALSE, col.names=FALSE)
          
          ms2.t <- getData(specV())
          ms2.temp <- tempfile("ms2-")
          write.table(ms2.t[,,drop=FALSE], file=ms2.temp, sep="\t", row.names=FALSE, col.names=FALSE)
          # call GenForm
          # is MS1 given?
          if(nrow(ms1.t) > 0)
          {
            cl <- sprintf(
              "%s ms=%s msms=%s %s out analyze",
              gfpath, ms1.temp, ms2.temp, cli)
            r1 <- system(cl, intern=TRUE)
            cl <- sprintf(
              "%s ms=%s msms=%s %s out analyze loss",
              gfpath, ms1.temp, ms2.temp, cli)
            r2 <- system(cl, intern=TRUE)
          }
          else
          {
            cl <- sprintf(
              "%s m=%s msms=%s %s out analyze",
              gfpath, m, ms2.temp, cli)
            r1 <- system(cl, intern=TRUE)
            cl <- sprintf(
              "%s m=%s msms=%s %s out analyze loss",
              gfpath, m, ms2.temp, cli)
            r2 <- system(cl, intern=TRUE)
          }
          pool$cli <- cl
          return(list(r1=r1, r2=r2))
          #})
        })
      }
      else(return(NA))
  })
  
  
  gfprocessed <- reactive({
    d <- gfresult()
    
    if(is.na(d))
      return(NA)
    #save(d, file="d.RData")
    d.processed <- lapply(d, function(r)
    {
      # cut the first 4 lines and the last line out
      r<-r[5:(length(r)-1)]
      # all lines not starting with whitespace are a new formula
      lines.nf <- which( substring(r,1,1) != " ")
      # assign all other lines to a formula
      assigned <- unlist(lapply(1:length(r), function(n) max(c(0,lines.nf)[c(0,lines.nf) < n], na.rm=FALSE)))
      formulae <- r[lines.nf]
      msms <- r[-lines.nf]
      msms.split <- split(msms, assigned[-lines.nf])
      formulae.t <- read.table(text=formulae, header=FALSE, sep="\t")
      if(ncol(formulae.t) < 5)
        formulae.t <- data.frame(character(), numeric(), numeric(), numeric(), numeric())
      colnames(formulae.t) <- c("formula", "dppm", "msmv", "msmsmv", "cmv")
      rownames(formulae.t) <- formulae.t$formula
      formulae.t$lineindex <- as.character(lines.nf)
      msms.split.t <- lapply(msms.split, function(msms.formula)
      {
        msms.f.t <- read.table(text=msms.formula, header=FALSE, sep="\t")
        if(ncol(msms.f.t) < 3)
          msms.f.t <- data.frame(numeric(),numeric(),numeric())
        colnames(msms.f.t) <- c("mz", "formula", "dppm")
        msms.f.t[,"mz"] <- na.locf(msms.f.t[,"mz"])
        msms.f.t
      })
      names(msms.split.t) <- names(msms.split)
      return(list(formulae = formulae.t, msms = msms.split.t))
    })
    d.formulae <- d.processed[[1]]$formulae
    d.formulae$sign <- ifelse(sign(d.formulae$dppm)==1, "+", "-")
    d.formulae$dppm <- abs(d.formulae$dppm)
    d.formulae <- d.formulae[,c("formula", "sign", "dppm", "msmv", "msmsmv", "cmv", "lineindex")]
    d.msms <- lapply(1:length(d.processed[[1]]$msms), function(n)
    {
      r <- cbind(d.processed[[1]]$msms[[n]], d.processed[[2]]$msms[[n]])
      if(!all(r[,1] == r[,4])) stop("Error in merging MSMS tables")
      if(!all(r[,3] == r[,6])) stop("Error in merging MSMS tables")
      r <- r[,c(1,2,5,3)]
      colnames(r) <- c("mz", "formula", "loss", "dppm")
      r[,"sign"] <- ifelse(sign(r$dppm)==1, "+", "-")
      r$dppm <- abs(r$dppm)
      r <- r[,c("mz", "formula", "loss", "sign", "dppm")]
      r
    })
    names(d.msms) <- names(d.processed[[1]]$msms)
    list(formulae = d.formulae, msms=d.msms)
  })
  
  
  observe({
    input$m.find
    if(input$m.find == 0) return()
    
    isolate({
      # find from precursor instead!
      #cpdVxx <<- cpdV()
      m <- cpdV()@mz
      
      # ms1 <- parentV()
      # print(ms1)
      # ms1xx <<- ms1
      # m <- ms1[which.max(ms1[,2,drop=FALSE]),1]
      updateTextInput(session,"m", value=m[[1]])
    })
  })
  
  output$formulae <- renderDataTable(gfprocessed()$formulae[,c("formula", "sign", "dppm", "msmv", "msmsmv", "cmv")])
  
  observe({
    d <- gfprocessed()
    if(!is.na(d))
    {
      l <- names(d$msms)
      names(l) <- d$formulae[match(l,d$formulae$lineindex),"formula"]
      updateSelectInput(session,"formula",choices=l, selected=l[[1]])
    }
  })
  
  msmsTable <- reactive({
    input$formula
    
    isolate({
      d <- gfprocessed()
      d$msms[[input$formula]]
    })
  })
  
  output$msms <- renderDataTable(msmsTable())
  
  output$cli <- renderText(pool$cli)
  #output$gfresults <- renderText(gfnice())
  
  output$msmsPlot <- renderPlot({
    ms2 <- r.ms2()
    msms <- msmsTable()
    
    plot(ms2[,2] ~ ms2[,1], data=ms2, type='h', col="black")
    abline(v=msms$mz, col="green")
    lines(ms2[,2] ~ ms2[,1], data=ms2, type='h', col="black")
  })
  
  
}



mPanel <-   tabsetPanel(
  tabPanel(
    "input",
    textInput("cli", "CLI parameters", ""),
    textInput("m","Parent ion",0),
    actionButton("m.find", "Find parent ion"),
    # p("MS1 (2 or 3 row with tab or space separator"),
    # checkboxInput("ms1.header", "with header row", TRUE),
    # tags$textarea(id="ms1", rows=10, cols=50),
    # p("MS2 (2 or 3 row with tab or space separator"), 
    # checkboxInput("ms2.header", "with header row", TRUE),
    # tags$textarea(id="ms2", rows=10, cols=50),
    actionButton("processGF", "Process")
  ),
  
  tabPanel(
    "output",
    textOutput("cli"),
    dataTableOutput("formulae"),
    selectInput("formula", "Formula",c()),
    dataTableOutput("msms")
  ),
  tabPanel(
    "msms",
    plotOutput("msmsPlot"))
)




.Workbench.ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      #fileInput("refspec", "Comparison")
      
      conditionalPanel(
        "input.mainTab == 'GenForm'",
          selectInput("cli.ion", "Ion", c("+H", "-H", "+e", "-e", "+Na"), "+H"),
          textInput("cli.ppm","MS1 ppm", 5),
          textInput("cli.acc","MS2 accept ppm", 3),
          textInput("cli.rej","MS2 reject ppm", 10),
          textInput("cli.ff","Fuzzy formula", "C0-30H0-60N0-10O0-6Cl0-5"),
          textInput("cli.args","arguments", "oei exist"),
          actionButton("cli.build", "Make CLI")
        )
    ),
    
    mainPanel(
      tabsetPanel(id="mainTab",
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
                 #includeScript("cdBinding.js"),
                chemdoodle_sketcherInput("getSmiles"),
                 textInput("smiles", "SMILES", ""),
                 #actionButton("getSmiles", "SMILES from viewer"),
                 actionButton("setSmiles", "SMILES to viewer (WIP)"),
                 actionButton("process", "Process"),
                 textOutput("mz"),
                 textInput("nameSim", "Name", ""),
                 actionButton("shift", "Shift input spectra to mz")
                

        ), # tabPanel input
        
        tabPanel("formula",
                 mPanel,
                 value="GenForm"),
        
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
                 tabsetPanel(
                 tabPanel("Settings",
                          
                          
                          fluidRow(
                            column(width=5,
                                   selectInput("compound1", "Compound 1",c()),
                                   "m/z:", textOutput("mz1"),
                                   checkboxInput("shiftTo1", "Shift to absolute mz", FALSE),
                                   numericInput("shift1", "mz shift", 0)
                            ),
                            column(width=5,
                                   selectInput("compound2", "Compound 2",c()),
                                   "m/z:", textOutput("mz2"),
                                   checkboxInput("shiftTo2", "Shift to absolute mz", FALSE),
                                   numericInput("shift2", "mz shift", 0)
                            )
                          ),
                          
                          checkboxInput("removeParent", "Remove parent peak from MS2 for comparison", FALSE),
                          numericInput("tolerance", "Tolerance (mz)", 0.001),
                          numericInput("cutoff", "Relative cutoff %", 5),
                          actionButton("calcMatrix", "calculate matrix")),
                 
                 
                 tabPanel("Similarity matrix",
                          dataTableOutput("simMatrix")
                 ),
                 tabPanel("Spectra similarity",
                          textOutput("similarity"),
                          fluidRow(
                            column(5, selectInput("inputIndex1", "cpd 1 spectrum", c())),
                            column(5, selectInput("inputIndex2", "cpd 2 spectrum", c()))
                          )
                          ,
                          imageOutput("simPlot", width=600, height=600),
                          tabsetPanel(
                            tabPanel("Aligned spectra", 
                                     dataTableOutput("simTab")),
                            tabPanel("Shift summary", 
                                     dataTableOutput("shiftTab"))
                          )
                 )
                 ) # tabsetPanel
        ) # tabPanel Comparison
      ) # tabsetPanel
    ) # mainPanel
  ) # sidebarLayout
) # ui
                


runBench <- function(options = getOption(.packageName))
{
  
  if(is.null(options))
  {
    error(paste0("Set options with option(", .packageName, ") or pass optionList to function alternatively. ?", .packageName,
                 "-options"))
  }
  
  gfpath <- options$path.GenForm
  cfmPath <- "C:/Software/cfm-id-2.0_win32/cfm-predict.exe"
  cfmSettings <- "0.001 C:/Software/cfm-id-2.0_win32/param_output.log C:/Software/cfm-id-2.0_win32/param_config.txt 1"
  
  
  shinyErr <- options("shiny.error")
  options("shiny.error" = NULL) 
  
  if(is.null("gfpath")){
    message("-----------------------------------------------------")
    message("Path to GenForm.exe is not set!")
    message("Before running this shiny app, set the GenForm path like this:")
    message("gfpath <- \"C:/Software/GenForm/GenForm.exe\" ")
    #gfpath <- "C:/Software/GenForm/GenForm.exe"
  }
  message("-----------------------------------------------------")
  message("Ignore all warnings and errors here. They are normal.")
  message("-----------------------------------------------------")
  
  
  
  shinyApp(ui = .Workbench.ui, server = .Workbench.server)
  options("shiny.error" = shinyErr)
}

