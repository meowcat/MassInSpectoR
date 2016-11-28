
OrgMSsim <- function(spec.top, spec.bottom, y = 1, z = 1, t = 0.005, b = 0.5, plotHtT = TRUE,
                     top.label = NULL, bottom.label = NULL, xlim = c(50,1000)){
  
  SimOutput <- list()
  
  top_tmp <- data.frame(mz = spec.top[, 1], intensity = spec.top[,2] ) 
  top_tmp$normalized <- round((top_tmp$intensity/max(top_tmp$intensity)) * 
                                100)
  top_tmp <- subset(top_tmp, top_tmp$mz >= xlim[1] & top_tmp$mz <= 
                      xlim[2])
  
  top_plot <- data.frame(mz = top_tmp$mz, intensity = top_tmp$normalized)
  top <- subset(top_plot, top_plot$intensity >= b)
  bottom_tmp <- data.frame(mz = spec.bottom[, 1], intensity = spec.bottom[,2] )
  
  bottom_tmp$normalized <- round((bottom_tmp$intensity/max(bottom_tmp$intensity)) * 
                                   100);
  bottom_tmp <- subset(bottom_tmp, bottom_tmp$mz >= xlim[1] & 
                         bottom_tmp$mz <= xlim[2])
  bottom_plot <- data.frame(mz = bottom_tmp$mz, intensity = bottom_tmp$normalized)
  bottom <- subset(bottom_plot, bottom_plot$intensity >= b)
  
  if((nrow(top) == 0) | (nrow(bottom) == 0))
  {
    plot.new()
    return(list(data.frame(), 0))
  }
  
  for (i in 1:nrow(bottom)) top[, 1][bottom[, 1][i] >= top[, 
                                                           1] - t & bottom[, 1][i] <= top[, 1] + t] <- bottom[, 
                                                                                                              1][i]
  alignment <- merge(top, bottom, by = 1, all = TRUE)
  if (length(unique(alignment[, 1])) != length(alignment[, 
                                                         1])) 
    warning("the m/z tolerance is set too high")
  alignment[, c(2, 3)][is.na(alignment[, c(2, 3)])] <- 0
  
  alignment[,1] <- sapply(alignment[,1], function(x) { x^y } )
  alignment[,2] <- sapply(alignment[,2], function(x) { x^z } )
  alignment[,3] <- sapply(alignment[,3], function(x) { x^z } )
  
  names(alignment) <- c("mz", "intensity.top", "intensity.bottom")
  match <-  as.vector(rep("NA", nrow(alignment)))
  for(i in 1:nrow(alignment)){
    if(alignment[i, "intensity.top"] != 0 && alignment[i, "intensity.bottom"] != 0){
      match[i] <- as.integer(1)
    }else{
      match[i] <- 0
    }
  }
  
  alignment$match <- as.numeric(match)
  SimOutput[[1]] <- as.data.frame(alignment)
  

 # print(alignment)
  u <- alignment[, 2]
  v <- alignment[, 3]
  similarity_score <- as.vector((u %*% v)/(sqrt(sum(u^2)) * 
                                             sqrt(sum(v^2))))
  SimOutput[[2]] <- as.numeric(similarity_score)
  
  if(plotHtT == TRUE){
  
  #pdf(file=pdftitle, width=11, height=6)
    
  plot.new()
  plot.window(xlim = xlim, ylim = c(-125, 125))
  ticks <- c(-100, -50, 0, 50, 100)
  for (i in 1:length(top_plot$mz)) lines(rep(top_plot$mz[i], 
                                             2), c(0, top_plot$intensity[i]), col = "blue")
  for (i in 1:length(bottom_plot$mz)) lines(rep(bottom_plot$mz[i], 
                                                2), c(0, -bottom_plot$intensity[i]), col = "red")
  axis(2, at = ticks, labels = abs(ticks), pos = xlim[1], ylab = "intensity")
  axis(1, pos = -125)
  lines(xlim, c(0, 0))
  rect(xlim[1], -125, xlim[2], 125)
  mtext("m/z", side = 1, line = 2)
  mtext("intensity (%)", side = 2, line = 2)
  plot.window(xlim = c(0, 20), ylim = c(-10, 10))
  text(10, 9, top.label)
  text(10, -9, bottom.label)
  #dev.off()
  }

  return(SimOutput)
}


