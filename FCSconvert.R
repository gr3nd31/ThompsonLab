FCSconvert <- function(dirFCS = "./",
                      limit = F,
                      concat = T){
  if(!require('flowCore')) {install.packages('flowCore')}
  setwd(dirFCS)
  fiList <- list.files(pattern = ".fcs")
  for(a in fiList){
    interimFCS <- flowCore::exprs(flowCore::read.FCS(a, transformation = F))
    interim <- as.data.frame(interimFCS)
    filname <- gsub("fcs", "csv", a)
    write.csv(interim, filname, row.names = F)
    if (concat==T){
      if(limit != F){
        interim <- interim[sample(nrow(interim), limit),]
      }
      interim$file <- a
      if(!exists("cells")){
        cells <- interim
      } else{
        cells <- rbind(cells, interim)
      }
      write.csv(cells, "all_cells.csv", row.names = F)
    }
  }
}