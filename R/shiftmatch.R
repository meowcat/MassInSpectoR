stepwiseShiftMatch <- function(s1, s2, order=c("count", "totalIntensity", "intensity"), maxHits = 20, cutoff = 0.05)
{
  s1 <- normalize(s1, slot="relint", scale=1, precision=3)
  s2 <- normalize(s2, slot="relint", scale=1, precision=3)
  s1 <- selectPeaks(s1, relint > cutoff)
  s2 <- selectPeaks(s2, relint > cutoff)
  df1 <- getData(s1)
  df2 <- getData(s2)
  shift <- as.vector(outer(df2$mz, df1$mz, "-"))
  int <- as.vector(outer(df1$i, df2$i, "*"))
  df <- data.frame(mz=shift, intensity=int)
  df.merged <- .filterMerge(df, 0.001)
}


