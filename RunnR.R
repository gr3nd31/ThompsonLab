#Welcome to RunnR version 0.5 - Runnin' Down a Dream

#This began as a theoretical package to analyze DNA fiber assays and, while it still may used for that in the future, it has been developed for cell line analyses instead
# Use:
#   1) Bulk analyses can be done by using the oneAndDone() function within a parent directory containing subdirectories which themselves hold the line txt files
#   2) Individual analysis can be done by using the lData() and lineAnalysis() functions

# Changes from versions 0.4:
#   -Added a RunnR oneAndDone() function
#   -Pearson's correlation and significance check is performed following graph generation
#   -Added legend to graphs (double check this, as things didn't quite line up during beta testing)

# Planned projects to be updated:
#   - Fix peak calling to automate.... everything?


#This loads the appropraite packages
suppressPackageStartupMessages(library(ggplot2))

dafault <<- FALSE


#This function loads the .csv file into the appropriate dataframe and adds the 'pix' variable for graphing
lData <- function(x, modelT=T){
  pointz <<- read.table(file=x, header=TRUE, fileEncoding = "latin1", sep = ",")
  if(!"pix" %in% colnames(pointz)){
    pointz$pix <<- c(1:nrow(pointz))
    write.table(pointz, file = x, sep = ",")
  }
  if (dafault == FALSE){
    pointz$redNam <<- readline(prompt = "What is the name of the red color: ")
    pointz$greenNam <<- readline(prompt = "What is the name of the green color: ")
    pointz$blueNam <<- readline(prompt = "What is the name of the blue color: ")
  }
  lineAnalysis()
  if (modelT==T){
    modelIt()
  }
}

#This function is used to save the figure as a high-res tiff. You can alter the dimensions if you wish.
savetif <- function(){
  name <- readline(prompt = "What would you like the name of the file to be: ")
  name <- paste0(name, ".tif")
  tiff(name, res = 300, units = "in", width = 12, height = 4)
  print(linera+theme_general)
  dev.off()
}

#----------------------------------------------------------------------------------------------------------------------
#This function graphs a line analysis for all three colors, then performs a multiple linear regression of R~B and R~G.
lineAnalysis <- function(){
  pointz$G <<- pointz$G/(mean(pointz$G)+0.1)
  pointz$R <<- pointz$R/(mean(pointz$R)+0.1)
  pointz$B <<- pointz$B/(mean(pointz$B)+0.1)
  print("pointz adjusted")
  linera <<- ggplot(data = pointz, aes(x = pix))+
    geom_step(data = pointz, aes(x=pix, y = G, colour = pointz$greenNam), size = 1.5)+
    geom_step(data = pointz, aes(x=pix, y = R, colour = pointz$redNam), size = 1.5)+
    geom_step(data = pointz, aes(x=pix, y = B, colour = pointz$blueNam), size = 1.5)+
    scale_colour_manual("", 
                       breaks = c(pointz$redNam, pointz$greenNam, pointz$blueNam),
                       values = c("blue", "red", "darkgreen"))+
    xlab("Relative pixel position")+
    ylab("Relative pixel intentsity")+
    theme(axis.line = element_line(colour = "black", size = 2))+
    theme(axis.text = element_text(face = "bold", color = "black", size = 12))+
    theme(axis.title = element_text(face = "bold", color = "black", size = 20))+
    theme_classic()+
    theme_general
  print(linera)
}

modelIt <- function(df=pointz){
  #Change the next lines to alter the regression analysis
  model <<- lm(df$G~df$R, data = df)
  coorz <<- cor.test(x = df$G, y = df$R, method = "pearson")
  
  # Significance is tested and, if found, the linear regression model is displayed
  if (coorz[3] < 0.05){
    print(paste0("Pearson's correlation finds significance (p = ", as.numeric(coorz[3]), ")."))
    print(summary(model))
  } else{
    print(paste0("Pearson's correlation does not find significance (p = ", as.numeric(coorz[3]), ")."))
  }
  print("Call 'summary(model)' to see the regression results again")
}

#----------------------------------------------------------------------------------------------------------------------
#This function graphs a line analysis for 'G', then puts a line at the midpoint for peak calling.
grapho <- function(){
  linera<<-ggplot(data=pointz, aes(x=pix, y=G))+
    geom_step(size = 1, colour = "darkgreen")+
    theme_classic()+
    xlab("Relative pixel position")+
    ylab("Relative pixel intensity")+
    theme(axis.line = element_line(colour = "black", size = 2))+
    theme(axis.ticks = element_line(colour = "black", size = 2))+
    geom_hline(yintercept = 127, linetype = "dashed", size = 1)+
    theme(axis.text = element_text(face = "bold", color = "black", size = 12))+
    theme(axis.title = element_text(face = "bold", color = "black", size = 20))+
    theme_general
  print(linera)
}

#This function is for peak calling the green channel. It returns the number of peaks and the average length by finding infection points across the midline and storing them in the 'runz' variable.
runner <- function(){
  cutoff<<-127
  runz <<- data.frame(pix=c(1), run_num = c(0), PoT = c(1), length=c(1))
  lenz<-0
  if(!"pix" %in% colnames(pointz)){
    pointz$pix <<- c(1:nrow(pointz))
  }
  c<-1
  for(i in 2:(nrow(pointz)-1)) {
    j<-abs(pointz$G[i]-pointz$G[i-1])
    if (pointz$G[i]>=cutoff & pointz$G[i-1]<cutoff){
      lenz<-lenz+1
      x<-c(pointz$pix[i], lenz, 2, 1)
      runz<<-rbind(runz, x)
    } else if(pointz$G[i]<=cutoff & pointz$G[i-1]>cutoff){
      lenz<-lenz+1
      x<-c(pointz$pix[i], lenz, 1, 1)
      runz<<-rbind(runz, x)
    }
    c<-c+1
  }
  runz<<-rbind(runz, c(nrow(pointz), lenz+1, if(pointz$G[nrow(pointz)] >127){1} else{2}, 1))
  for (i in 1:nrow(runz)){
    if (runz$PoT[i]=="1"){
      runz$PoT[i]<<-"peak"
    } else{
      runz$PoT[i]<<-"trough"
    }
  }
  for (i in 2:nrow(runz)){
    runz$length[i]<<-runz$pix[i]-runz$pix[i-1]
  }
  print("The number of peaks and troughs are:")
  print(table(runz$PoT))
  cat("\n")
  print(paste("The average peak is", as.character(round(mean(runz$length[runz$PoT=="peak"]), digits=2)), "pixals."))
  print(paste("The peak standard deviation is", as.character(round(sd(runz$length[runz$PoT=="peak"]), digits=2)), "pixals."))
  cat("\n")
  print(paste("The average trough is", as.character(round(mean(runz$length[runz$PoT=="trough"]), digits=2)), "pixals."))
  print(paste("The trough standard deviation is", as.character(round(sd(runz$length[runz$PoT=="trough"]), digits=2)), "pixals."))
  grapho()
}


#Here is the general theme used by Jason to darken things and things
theme_general <- theme(legend.position = "top", 
                       legend.text = element_text(size = 12), 
                       legend.title = element_text(size = 16), 
                       axis.line = element_line(color = "black", size = 1.5), 
                       axis.ticks = element_line(color = "black", size = 1.5), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       axis.title = element_text(size = 21, family = "sans", color = "black"), 
                       axis.text = element_text(size = 16, family = "sans", color = "black"),
                       legend.key.size = unit(1.5, "cm")
)

#----------------------------------------------------------------------------------------------------------------------
# This is the RunnR version of oneAndDone, which, starting from a parent directory, iterates through subdirectories containing the txt files of the line analyses and generates a graph of each
# Instead of asking for names for each graph, however, it will use the same names for everything

oneAndDone <- function(){
  dafault <<- TRUE
  redNam <<- readline(prompt = "What is the name of the red color for all graphs: ")
  greenNam <<- readline(prompt = "What is the name of the green color for all graphs: ")
  blueNam <<- readline(prompt = "What is the name of the blue color for all graphs: ")
  folds <<- list.dirs()
  for (i in folds[2:length(folds)]){
    setwd(i)
    filz <- list.files()
    for (j in filz){
      lData(j)
      pointz$redNam <<- redNam
      pointz$greenNam <<- greenNam
      pointz$blueNam <<- blueNam
      namen <- paste0(substring(j, 1, nchar(j)-4), "_line.tif")
      tiff(namen, res = 300, units = "in", width = 12, height = 4)
      print(linera+theme_general)
      dev.off()
    }
    setwd("../")
  }
}

#----------------------------------------------------------------------------------------------------------------------

print("Welcome to RunnR, the script to perform line analyses and fiber assays. To start, call 'lData(x)' to load a text file and perform a line analysis with multiple linear regression if signifance is found (default is Pearson's correlation), or use oneAndDone() for bulk graph generation. If you would like to quantify a fiber assay, call 'runner()' instead. Note - runner is currently only set up to analyze single color.")