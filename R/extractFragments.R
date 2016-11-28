
getCompound <- function(db, id)
{
  subDb <- db[db$CompoundID == id,,drop=FALSE]
  ces <- unique(subDb$CollisionEnergy)
  cpd <- new("RmbSpectraSet")
  spectra.db <-   lapply(ces, function(ce) extractFrags(db, id, ce))
  cpd@children <- as(
    lapply(spectra.db, function(sp) 
      {
        s <- new("RmbSpectrum2",
                 mz=sp$spectrum[,1],
                 intensity=sp$spectrum[,2],
                 precursorMz=as.numeric(sp$info$mz),
                 collisionEnergy=as.numeric(sp$info$CE)
                 )
    }),
    "SimpleList")
  cpd@name <- as.character(id)
  cpd@mz <- spectra.db[[1]]$info$MZ
  cpd
}


extractFrags <- function(all.frag, id, opCE)
{
  info <- list(
    found=c(),
    MZ=c(),
    Ion=c(),
    Origin=c(),
    Res=c(),
    CE=opCE
  )
  #   sameCE[[j]]$result <- c()
  #   sameCE[[j]]$parMZ <- c()
  #   sameCE[[j]]$parIon <- c()
  #   sameCE[[j]]$parOrigin <- c()
  #   sameCE[[j]]$parRes <- c()
  #   sameCE[[j]]$parCE <- c()
  #   sameCE[[j]]$TPMZ <- c()
  #   sameCE[[j]]$TPIon <- c()
  #   sameCE[[j]]$TPOrigin <- c()
  #   sameCE[[j]]$TPRes <- c()
  #   sameCE[[j]]$TPCE <- c()
  #   sameCE[[j]]$pairID <- pair.ID
  
  
  
  
  
  wheret <- which(all.frag$CompoundID == id)
  if(length(wheret) == 0){
    y <- as.data.frame(c(0,0))
    info$found <- FALSE
    return(list(info=info))
  }
  
  tpmass <- as.numeric(paste(all.frag[wheret[1], "PrecursorMZ"]))
  info$MZ <- tpmass
  
  tpion <- paste(all.frag[wheret[1], "IonMode"])
  info$Ion <- tpion
  #   ption <- paste(all.frag[wherep[1], "IonMode"])
  #   sameCE[[j]][i,"TPIon"] <- ption
  
  #   if(ption != tpion){
  #     sameCE[[j]][i,"result"] <- paste("Different ionization modes!!")
  #     y[[i]] <- as.data.frame(c(0,0))
  #     z[[i]] <- as.data.frame(c(0,0))
  #     next
  #   }
  
  
  ### Extraction of TP fragments
  
  tpdat <- subset(all.frag, all.frag$CompoundID == id)
  tpdat <- subset(tpdat, tpdat$CollisionEnergy == opCE)
  
  if(nrow(tpdat) == 0){
    y <- as.data.frame(c(0,0))
    info$found <- FALSE
    return(list(info=info))
  }
  
  res <- unique(tpdat$Resolution)
  
  if(length(res) > 1){
    
    take <- max(res)
    tpdat <- subset(tpdat, tpdat$Resolution == take)
    info$Res <- take
  }else{
    info$Res <- res
    
  }
  
  if(mean(table(tpdat$FragmentRank)) != 1 ){
    
    take <- unique(tpdat$origin)[1]
    tpdat <- subset(tpdat, tpdat$origin == take)
    info$Origin <- take
    # break
    
  }else{
    info$Origin <- unique(tpdat$origin)
  }
  
  tp <- as.data.frame(cbind(tpdat$FragmentMZ, tpdat$FragmentIntensity))
  y <- as.data.frame(tp[order(tp$V1), ])
  y <- unique(y)
  info$found <- TRUE
  return(list(info=info, spectrum=y))
}
