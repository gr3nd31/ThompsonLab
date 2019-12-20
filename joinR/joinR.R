# This script will take two files: a 'tagger' file (which will how cells are numbered) and a 'fraction' file (which has the fraction that needs identifying)
# It then will save the new, tagged fraction file.
# Of note, this works by a 'nearest neighbor' principle using the centerpoints of the tagger cells compared to the adjusted centerpoints of the fractions
# This means that it is ENTIRELY possible that the wrong cell will be called, but this is unlikely and shouldn't happen frequently
# Also, depending on unforseen possibilities, this script takes only the first possible "shortest" distance. Something to keep in mind

joinR <- function(y, x = "dna.csv"){
  tagger <- read.csv(x)
  fraction <- read.csv(y)
  if (!"Number" %in% names(tagger)){
    names(tagger)[1] <- "Number"
    write.csv(tagger, file = x, row.names = FALSE)
  }
  if (!"Number" %in% names(fraction)){
    names(fraction)[1] <- "Number"
    write.csv(fraction, file = x, row.names = FALSE)
  }
  
  fraction$cell <- "unknown"
  tagger$distance <- 0
  for (i in 1:nrow(fraction)){
    tagger$distance <- sqrt((tagger$X-fraction$XM[i])^2+(tagger$Y-fraction$YM[i])^2)
    closest <- min(tagger$distance)
    fraction$cell[i] <- subset(tagger, distance == closest)$Number[1]
    if (length(closest) > 1){
      print("PING!")
    }
  }
  write.csv(fraction, file = y, row.names = FALSE)
}