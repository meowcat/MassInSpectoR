
### Producing the list for comparing MS/MS spectra from same CEs for the pairs

sameCEdata <- function(all.frag, pairs){
  
  
  TP.ID <- as.integer(pairs$TP_UchemID)
  par.ID <- as.integer(pairs$par_UchemID)
  pair.ID <- as.integer(pairs$pair_ID)
  y <- list()
  z <- list()
  nowp <- list(0)
  nowt <- list(0)
  
  nCE <- dim(table(all.frag$CollisionEnergy))
  sameCE <- as.list(rep("", nCE))
  names(sameCE) <- unique(all.frag$CollisionEnergy)
  
  
  for(j in 1:nCE){
    
    opCE <- names(table(all.frag$CollisionEnergy))[j]
    sameCE[[j]] <- data.frame(matrix(NA, ncol = 1, nrow = nrow(pairs)))
    
    #colnames(sameCE) <- names(table(all.frag$CollisionEnergy))
    rownames(sameCE[[j]]) <- rownames(pairs)
    
    sameCE[[j]]$result <- c()
    sameCE[[j]]$parMZ <- c()
    sameCE[[j]]$parIon <- c()
    sameCE[[j]]$parOrigin <- c()
    sameCE[[j]]$parRes <- c()
    sameCE[[j]]$parCE <- c()
    sameCE[[j]]$TPMZ <- c()
    sameCE[[j]]$TPIon <- c()
    sameCE[[j]]$TPOrigin <- c()
    sameCE[[j]]$TPRes <- c()
    sameCE[[j]]$TPCE <- c()
    sameCE[[j]]$pairID <- pair.ID
    
    for(i in 1:length(TP.ID)){
      
      
      tpid <- as.integer(paste(TP.ID[i]))
      ptid <- as.integer(paste(par.ID[i]))
      tp <- extractFragments(all.frag, tpid, opCE)
      pt <- extractFragments(all.frag, ptid, opCE)
      
      tpInfo <- tp$info
      ptInfo <- pt$info
      sameCE[[j]][i,c("parMZ", "parIon", "parOrigin", "parRes", "parCE")] <- ptInfo[c("MZ", "Ion", "Origin", "Res", "CE")]
      sameCE[[j]][i,c("TPMZ", "TPIon", "TPOrigin", "TPRes", "TPCE")] <- tpInfo[c("MZ", "Ion", "Origin", "Res", "CE")]
      if(tpInfo$found & ptInfo$found)
      {
        sameCE[[j]][i,"result"] <- paste("both measured")
        y[[i]] <- tp$spectrum
        z[[i]] <- pt$spectrum
        if(tpInfo$Ion != ptInfo$Ion){
          sameCE[[j]][i,"result"] <- paste("Different ionization modes!!")
          y[[i]] <- as.data.frame(c(0,0))
          z[[i]] <- as.data.frame(c(0,0))
        }
      }
      else if(tpInfo$found & !ptInfo$found)
      {
        y[[i]] <- as.data.frame(c(0,0))
        z[[i]] <- as.data.frame(c(0,0))
        sameCE[[j]][i,"result"] <- paste("No par")
      }
      else if(!tpInfo$found & ptInfo$found)
      {
        y[[i]] <- as.data.frame(c(0,0))
        z[[i]] <- as.data.frame(c(0,0))
        sameCE[[j]][i,"result"] <- paste("No TP")
      }
      
      
      
    }
    
    nowp[[j]] <- z
    nowt[[j]] <- y
    
  }
  
  #rm(pt, ptdat, tp, tpdat, opCE, ptid, ption, ptmass, res, take, tpid, tpion, tpmass, wherep, wheret)
  
  return(list(sameCE, nowt, nowp))
}