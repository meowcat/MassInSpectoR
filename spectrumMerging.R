

filterMerge <- function(cpd, dmz = 0.001)
{
  peaks <- lapply(cpd@children, getData)
  peaks <- do.call(rbind, peaks)
  merged <- .filterMerge(peaks, dmz)
  spec <- new("RmbSpectrum2")
  #, mz=merged$mz, intensity=merged$intensity)
  spec <- setData(spec, merged)
  property(spec, "count", addNew=TRUE, "numeric") <- merged$count
  property(spec, "addedIntensity", addNew=TRUE, "numeric") <- merged$addedIntensity
  property(spec, "totalIntensity", addNew=TRUE, "numeric") <- merged$totalIntensity
  cpd@children <- as(c(spec, unlist(cpd@children)), "SimpleList")
  cpd
}


.filterMerge <- function (peaks, dmz, int = 1e20) 
{
  cutoff_int_limit <- int
  cutoff_mz_limit <- dmz
  peaks_o <- peaks[order(peaks$intensity, decreasing = TRUE), ]
  n <- 1
  peaks_o$count <- rep(0, nrow(peaks_o))
  while (n <= nrow(peaks_o)) {
    
    windowpeaks <- which( (peaks_o$mz > peaks_o[n, "mz"] - cutoff_mz_limit) & (peaks_o$mz < peaks_o[n, "mz"] + cutoff_mz_limit))
    if(length(windowpeaks) > 1) {
      peaks_o[n, "count"] <- length(windowpeaks)
      peaks_o[n, "addedIntensity"] <- sum(peaks_o[windowpeaks[-1], "intensity"])
      peaks_o <- peaks_o[-windowpeaks[-1],]
    }
    n <- n + 1
  }
  if(sum(peaks_o$count == 0)>0)
  {
    peaks_o[peaks_o$count == 0,c("count", "addedIntensity")] <- cbind(rep(1, sum(peaks_o$count == 0)),
                                                                      rep(0, sum(peaks_o$count == 0)))
    peaks_o[,"totalIntensity"] <- rowSums(peaks_o[,c("intensity", "addedIntensity")])
  }
  return(peaks_o[order(peaks_o$mz), ])
}
