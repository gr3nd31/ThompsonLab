#SirMixaPlot version 7.0 - Electric Plot-a-loo

#----------------------------------------------------------------------------------------------------------

#This program is designed to take a csv (comma separated values) file and generate scatterplots based on columns.

#Changes from SirMixaPlot version 6.2:
#   -Takes incoming data from imaGen()
#   -Got rid of the savetif() since ggsave() is prebuilt into ggplot2
#   -color_density() now works with everything automatically based on the current qq
#   -The gate -> negPop -> normalizer cascade now works without being prompted for variables and only if EdU has been corrected
#   -Removed show_color_names() as colors no longer matter and targets are labeled in the headers


#Things still that need to be coded:
#   -I dunno. Think of some things

#-----------------------------------------------------------------------------------------------------------

#Program opens required packages

#If these packages aren't install, uncomment the following lines to install them.

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("MASS")
#install.packages("viridis")
#install.packages("stringr")

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(stringr))

#------------------------------------------------------------------------------------------------------------

#Functions for your pleasure

#--------------------------------------------------------------------------------------

#The imaGen() script takes any number of csv files that are within your working directory and combines them. Take caution, however, as this script is designed for a particular imagej csv outpt
#It is important that the jetData uses the Set Measurements with the following settings active:
#   Area    Standard Deviation    Min & Max gray value    Center of Mass    Mean gray value   Perimenter    Display label
# Additionally, for the best results, decimal places should be set to 9
imaGen <- function(directory="./"){
  setwd(directory)
  filez <- list.files(pattern = ".csv")
  cells <- read.table(filez[1], sep = ",", header = TRUE)
  cat("REMEMBER: One of these csv's must be labeled 'dna'")
  colo <- readline(prompt= paste0("What color is in ", filez[1], ": "))
  names(cells) <- paste0(names(cells), "_", colo)
  names(cells)[1]<-"Number"
  names(cells)[3]<-"Area"
  names(cells)[10]<-"Perimeter"
  names(cells)[9]<-"Y_location"
  names(cells)[8]<-"X_location"
  cells<-cells[,-2]
  cells <- cbind(cells[,1:2], cells[,7:9], cells[,3:6])
  cells$Area <- cells$Area*100
  for (i in filez[2:length(filez)]){
    if (!grepl("_all", i)){
      interim <- read.table(i, sep = ",", header = TRUE)
      colo <- readline(prompt = paste0("What color is in ", i, ": "))
      names(interim) <- paste0(names(interim), "_", colo)
      cells <- cbind(cells, interim[,4:6])
    }
  }
  filnam <- readline(prompt = "What should the file be named: ")
  write.csv(cells, file = paste0(filnam,"_all.csv"), row.names = FALSE)
}

#--------------------------------------------------------------------------------------
#joiner() connects the sirmixaplot() with grapho() and can be used to graph new plots without running through sirmixaplot(). Default dataset is cells but any dataset can be used
joiner <- function(df=cells){
  #newAxis is reset for the changeAxis() function
  newAxis <<- FALSE
  
  #cats the dataset columns
  cat("\n")
  cat("Data categories detected:")
  cat("\n")
  print(sort(names(df)))
  cat("\n")
  
  #This sets the x and y variables
  px<<-readline(prompt = "px - What parameter is the x-axis: ")
  xname<<-readline(prompt = "What is the name of the x axis: ")
  cat("\n")
  py<<-readline(prompt = "py - What parameter is the y-axis: ")
  yname<<-readline(prompt = "What is the name of the y axis: ")
  cat("Working...")
  cat("\n")
  
  #calls the graphing function
  grapho(df)
}

#grapho() actually makes the scatter plots
grapho <- function(df=cells){
  
  #backup dataframe is generated
  interim <<- cells

  cat("Data frame generated")
  cat("\n")
  
  #Initial scatter plot is generated
  qq<<-ggplot(data=df, aes_string(px, py))+
    theme_linedraw()+
    theme(plot.title = element_text(size=40, face="bold"))+
    theme(axis.title = element_text(size=30, face="bold"))+
    geom_point(size=2)+
    scale_x_continuous(name=xname)+
    scale_y_continuous(name=yname)+
    theme(axis.text = element_text(face='bold', size=16), axis.ticks = element_blank())
  
  if(newAxis == TRUE){
    changeAxis(store[1], store[2], store[3], store[4])
  }
  
  #Program finishes up, generates the scrollable HTML file
  cat("Plot generated")
  cat("\n")
  print(qq)

  cat("Call 'modthequad(df)' to apply quadrant analysis")
  cat("\n")
  cat("Call 'cull()' to remove certain populations from original dataset (DANGER, but less than before")
  cat("\n")
  cat("Call 'eGod(df)' to apply ergodic analysis of S phase cells")
  cat("\n")
  cat("Call 'color_denisty()' to add density gradient to qq plot")
  cat("\n")
  cat("Call 'changeAxis(min X, max, X, min Y, max Y)' to change the axes of qq plots")
  cat("\n")
  cat("Call 'parser(df, df$variable)' to create a gated population")
  cat("\n")
}

#----------------------------------------------------------------

sirmixaplot <- function(filo){
  #The script, m'lord
  
  #This is the beginning of the readline prompts for the program. User inputs desired output conditionals. Also, certain columns are added for safety...
  thingee <<- filo
  cat(paste("Opening file: ", thingee, sep=""))
  cat("\n")
  
  #This checks if the file has the '_cells', which signals it has been previously processed. Otherwise it extracts the appropriate rows and creates a cells file
  if (grepl("_cells.csv", thingee, fixed = TRUE)){
    cells <<- read.table(file=thingee,header=TRUE, fileEncoding = "latin1", sep = ",")
  } else {
    cells <<- read.table(file=thingee,header=TRUE, fileEncoding = "latin1", sep = ",")
    cat("Creating a new file, sir or madam.")
    cat("\n")
    thingee <<- paste0(substr(thingee, 1, nchar(thingee)-4), "_cells.csv")
    cat(paste(thingee, "created."))
    for (i in 1:ncol(cells)){
      if (grepl("Mean_", names(cells)[i])){
        cells[paste0("I", names(cells)[i])] <<- cells[,i]*cells$Area
      }
    }
    cells$log2_dna <<- log(cells$IMean_dna, 2)
    for (i in 1:ncol(cells)){
      if (grepl("Mean_", names(cells)[i])){
        cells[paste0("L", names(cells)[i])] <<- log(cells[,i], 10)
      }
    }
  }
  
  #writes the file
  write.csv(cells, file = thingee, row.names = FALSE)
  
  #Checks that some size correctiong has occured and begins making the graph
  if(TRUE %in% str_detect(names(cells), "ALIMean_.")){
    joiner(cells)
  } else{
    joiner(cells)
    gate()
  }
}

#----------------------------------------------------------------------------------
#This function creates quadrants and calculates the percentage in each quadrant
#   Geom_line function can be altered to change the asthetics of the quadrant lines

#quadrify() assesses if a point is above or below the assigned x- and y-intercepts as True/False
quadrify <<- function(df){
  for (i in df[px]){
    lr <<- (i>as.numeric(xi))
    xmin <<- min(df[px])
    xmax <<- max(df[px])
  }
  for (i in df[py]){
    ud <<- (i>as.numeric(yi))
    ymin <<- min(df[py])
    ymax <<- max(df[py])
  }
}

#quadricate() takes the True/False pairs for each point and assigns it a quadrant
quadricate <<- function(df){
  for (i in 1:nrow(df)){
    if(df[i, "lr"]==FALSE & df[i, "ud"]==FALSE){
      quadrant <<- append(quadrant, "Q3")
    }else if(df[i, "lr"]==TRUE & df[i, "ud"]==FALSE){
      quadrant <<- append(quadrant, "Q4")
    }else if(df[i, "lr"]==FALSE & df[i, "ud"]==TRUE){
      quadrant <<- append(quadrant, "Q1")
    }else {
      quadrant <<- append(quadrant, "Q2")
    }
  }
}

#CountR() gets the percentage of each quadrant
countR <<- function(vecky, df){
  quad <<- c("Q1", "Q2", "Q3", "Q4")
  todos <<- nrow(df)
  qt1 <<- sum(vecky == "Q1")
  qt2 <<- sum(vecky == "Q2")
  qt3 <<- sum(vecky == "Q3")
  qt4 <<- sum(vecky == "Q4")
  PercentTotal <<- c(format(round(qt1/todos*100, 2), nsmall = 2), format(round(qt2/todos*100, 2), nsmall = 2), format(round(qt3/todos*100, 2), nsmall = 2), format(round(qt4/todos*100, 2), nsmall = 2))
  quads <<- data.frame(quad, PercentTotal)
}

modthequad <<- function(df=cells){
  quadrant <<- c()
  xi<<-readline(prompt = "What is the x-intercept: ")
  yi<<-readline(prompt = "What is the y-intercept: ")
  quadrify(df)
  fu <<- data.frame(lr, ud)
  quadricate(fu)
  countR(quadrant, df)
  qqmod <<- qq
  qqmod <<- qqmod+
    geom_hline(yintercept = as.numeric(yi), color = "red", linetype = "dashed", size = 2)+
    geom_vline(xintercept = as.numeric(xi), color = "red", linetype = "dashed", size = 2)+
    annotate("text", x = xmin, y = ymax, label = paste0("Q1: ", as.character(format(round((qt1/todos*100), 2), nsmall = 2)), "%"))+
    annotate("text", x = xmax, y = ymax, label = paste0("Q2: ", as.character(format(round((qt2/todos*100), 2), nsmall = 2)), "%"))+
    annotate("text", x = xmin, y = ymin, label = paste0("Q3: ", as.character(format(round((qt3/todos*100), 2), nsmall = 2)), "%"))+
    annotate("text", x = xmax, y = ymin, label = paste0("Q4: ", as.character(format(round((qt4/todos*100), 2), nsmall = 2)), "%"))
  print(quads)
  print(qqmod+geom_density_2d())
  cat("\n")
  cat("call 'qqmod' to see plot with quadrants.")
}

#---------------------------------------------------------------------------
#This function will remove certain points from the cells dataframe.

cull <- function(){
  print(colnames(cells))
  para <<- readline(prompt = "Which parameter should be culled: ")
  bORs <<- readline(prompt = "Data with values (b)igger, (s)maller, or (=) to the target: ")
  vale <<- readline(prompt = "What is the target value: ")
  remover<<-c()
  if (bORs=="b"){
    for (i in 1:nrow(cells)){
      checkr <<- cells[i, para]
      if (checkr > as.numeric(vale)){
        remover <<- c(remover, i)
      }
    }
    cells<<-cells[-remover,]
  } else if (bORs=="s"){
    for (i in 1:nrow(cells)){
      checkr <<- cells[i, para]
      if (checkr < as.numeric(vale)){
        remover <<- c(remover, i)
      }
    }
    cells<<-cells[-remover,]
  } else if (bORs=="="){
    for (i in 1:nrow(cells)){
      checkr <<- cells[i, para]
      if (checkr == vale){
        remover <<- c(remover, i)
      }
    }
    cells<<-cells[-remover,]
  } else{
    grapho(cells)
  }
  grapho(cells)
  cat("\n")
  cat("Previous 'cells' has been saved as 'interim'")
}

#----------------------------------------------------------------------------
#This function performs ergodic analysis on the graph, but can currently only do it along the x-axis, with a y-axis cut-off value.
#eGod will take your scatterplot and perform ergodic analysis on an x-axis

eGod <- function(df=cells){
  rate <<- as.integer(readline(prompt = "What is the rate value: "))
  bind <<- as.integer(readline(prompt = "How many bins: "))
  bn <<- as.integer(readline(prompt = "Remove number of final bins: "))
  hig <<- max(df[px])
  low <<- min(df[px])
  divine <<- (hig-low)/(bind-1)
  yt <<- as.integer(readline(prompt = "What is your S phase cutoff: "))
  binner <<- matrix(nrow = nrow(df), ncol = 1)
  
  sPhase <<- 0
  binner <<- c()
  binSize <<- c()
  
  for (i in 1:(bind-bn)){
    binSize[i] <<- low+(i*divine)
  }
  for (i in 1:nrow(df)){
    if (df[i,py]>yt){
      binner[i] <<- floor((df[i,px]-low)/divine)+1
      sPhase <<- sPhase+1
    }else{
      binner[i] <<- 0
    }
  }
  tot <- nrow(df)
  
  CpB <<- c()
  for (i in 1:(bind-bn)){
    ghetto <- sum(binner == i)
    CpB[i] <<- ghetto
  }
  
  a <- log(2)/rate
  Ft <- sPhase/tot
  bRate <<- c()
  for (i in 1:length(CpB)){
    bRate[i] <<- a*((2-Ft)/(CpB[i]/sPhase))
  }
  binNum <- c()
  for (i in 1:length(CpB)){
    binNum[i] <- i
  }
  justice <<- data.frame(binNum, binSize, CpB, bRate)
  hh <<- ggplot(justice, aes(binSize))+
    geom_bar(aes(weight=bRate))+
    xlab("Size of binning variable")+
    ylab("Rate")
  print(justice)
  print(hh)
}

#----------------------------------------------------------------
#This function creates new dataframes of gated populations.

gate <- function(df=cells, gateName="negPop"){
  print(qq+geom_density_2d())
  cat(paste0("Gating out the ", gateName, " population:"))
  cat("\n")
  x1 <<- readline(prompt = "Which value of pX should the gating begin: ")
  x2 <<- readline(prompt = "Which value of pX should the gating end: ")
  y1 <<- readline(prompt = "Which value of py should the gating begin: ")
  y2 <<- readline(prompt = "Which value of py should the gating end: ")
  xSet <- c(x1, x2)
  ySet <- c(y1, y2)
  
  #gates<<-append(gates, gateName)
  
  workaround<<- subset(df, as.numeric(get(px)) > as.numeric(xSet[1]) & as.numeric(get(px)) < as.numeric(xSet[2]))
  workaround <<- subset(workaround, as.numeric(get(py)) > as.numeric(ySet[1]) & as.numeric(get(py)) < as.numeric(ySet[2]))
  assign(gateName, workaround, envir = .GlobalEnv)
  
  cells$gatePop<<-cells$Number %in% workaround$Number
  
  #grapho(workaround)
  if(gateName == "negPop"){
    negPop_correct()
  } else{
    colnames(cells)[match("gatePop", names(cells))]<<-gateName
  }
}

#This function uses a true negative population to adjust the background
negPop_correct <- function(){
  x <- strsplit(px, "_")[[1]][2]
  y <- strsplit(py, "_")[[1]][2]
  px_background <<- mean(negPop[names(negPop)==paste0("Mean_", x)][,1])
  py_background <<- mean(negPop[names(negPop)==paste0("Mean_", y)][,1])
  if (x != "dna"){
    cells[paste0("AIMean_", x)] <<- cells[paste0("IMean_", x)]-(cells$Area*px_background)
    cells[paste0("AIMean_", x)] <<- cells[paste0("AIMean_", x)] + abs(min(cells[paste0("AIMean_", x)]))+1
    cells[paste0("ALIMean_", x)] <<- log(cells[paste0("AIMean_", x)], 10)
  }
  if (y != "dna"){
    cells[paste0("AIMean_", y)] <<- cells[paste0("IMean_", y)]-(cells$Area*py_background)
    cells[paste0("AIMean_", y)] <<- cells[paste0("AIMean_", y)] + abs(min(cells[paste0("AIMean_", y)]))+1
    cells[paste0("ALIMean_", y)] <<- log(cells[paste0("AIMean_", y)], 10)
  }
  #if(grepl("CMean_", px) | grepl("CMean_", py)){
   # if (grepl("CMean_", px)){
  #    cells[paste0("AICMean_", x)] <<- cells[paste0("ICMean_", x)]-(cells$Area*px_background)
  #    cells[paste0("AICMean_", x)] <<- cells[paste0("AICMean_", x)] + abs(min(cells[paste0("AIMean_", x)]))+1
  #    cells[paste0("ALICMean_", x)] <<- log(cells[paste0("AICMean_", x)], 10)
  #  } else if(grepl("CMean_", py)){
  #    cells[paste0("AICMean_", y)] <<- cells[paste0("ICMean_", y)]-(cells$Area*py_background)
  #    cells[paste0("AICMean_", y)] <<- cells[paste0("AICMean_", y)] + abs(min(cells[paste0("AIMean_", y)]))+1
  #    cells[paste0("ALICMean_", y)] <<- log(cells[paste0("AICMean_", y)], 10)
  #  }
  #}
  write.csv(cells, file = thingee, row.names = FALSE)
  
  if (y == "edu" | x == "edu"){
    normalizer()
  }
}

#----------------------------------------------------------------
#This function will automatically pseudo-color density of a cell cycle profile and provides the basic functions for denisty pseudo-coloring
get_density <- function(px, py, ...){
  theme_set(theme_bw(base_size = 16))
  dens <<- MASS::kde2d(px, py, ...)
  ix <<- findInterval(px, dens$x)
  iy <<- findInterval(py, dens$y)
  ii <<- cbind(ix, iy)
  return(dens$z[ii])
}

color_density<-function(){
  cells$density <<- get_density(cells[px][,1], cells[py][,1], n=100)
  qq <<- ggplot(data = cells, aes_string(px, py)) + geom_point(aes(color = density)) + scale_color_viridis()+xlab(xname)+ylab(yname)
  
  if(newAxis == TRUE){
    changeAxis(store[1], store[2], store[3], store[4])
  }
  print(qq)
}

#----------------------------------------------------------------
#This function changes the axis of qq files via changeAxis(X0, X1, Y0, Y1)

changeAxis <- function(w,x,y,z){
  newAxis <<- TRUE
  store <<- c(w,x,y,z)
  qq<<-qq+coord_cartesian(xlim = c(w,x), ylim = c(y, z))
  print(qq+geom_density_2d())
}
#----------------------------------------------------------------

#----------------------------------------------------------------
#This function asks you to define a G1 peak and re-assigns that value to  '1'
normalizer <- function(){
  
  #create a column holder (just in case...)
  cells$edu<<-"holder"
  cells$ploidy<<-"holder"
  cat("Determining the G1 population now.")
  cat("\n")
  
  #This part calls an EdU negative popuation and creates the normalized values
  bigBin <- which.max(density(cells$ALIMean_edu)$y)
  edu_neg <- density(cells$ALIMean_edu)$x[bigBin]+1
  edu_hist <- ggplot(cells, aes(ALIMean_edu))+geom_density()+geom_vline(xintercept = edu_neg)+xlab("Log EdU")
  print(edu_hist)
  g1_good <- readline(prompt = paste0("Is ", format(round(edu_neg, 2), nsmall = 2), " representative of the EdU cutoff? (y/n) "))
  if (g1_good == "n"){
    edu_neg = as.numeric(readline(prompt = "What value should EdU be cutoff as negative? "))
    edu_hist <- ggplot(cells, aes(ALIMean_edu))+geom_density()+geom_vline(xintercept = edu_neg)+xlab("Log EdU")
    print(edu_hist)
  }
  cat(paste0("Using ", format(round(edu_neg, 2), nsmall = 2), " as the EdU cutoff."))
  cat("\n")
  for (i in 1:nrow(cells)){
    if (cells$ALIMean_edu[i] > edu_neg){
      cells$edu[i] <<- "Positive"
    } else {
      cells$edu[i] <<- "Negative"
    }
  }
  eduNegCells <-subset(cells, ALIMean_edu < edu_neg)
  cells$edu_norm <<- (cells$ALIMean_edu+1)-mean(eduNegCells$ALIMean_edu)
  
  #This part calls the 2N peak that is EdU-negative   
  bigBin <- which.max(density(eduNegCells$log2_dna)$y)
  diploid <- density(eduNegCells$log2_dna)$x[bigBin]
  dna_hist <<- ggplot(eduNegCells, aes(log2_dna))+geom_density()+geom_vline(xintercept = diploid)+xlab("DNA content")
  print(dna_hist)
  g1_good <- readline(prompt = paste0("Is ", format(round(diploid, 2), nsmall = 2), " representative of the 2N peak? (y/n) "))
  if (g1_good == "n"){
    diploid = as.numeric(readline(prompt = "What is the value of the 2N peak? "))
    dna_hist <<- ggplot(eduNegCells, aes(log2_dna))+geom_density()+geom_vline(xintercept = diploid)+xlab("DNA content")
    print(dna_hist)
  }
  cat(paste0("Using ", format(round(diploid, 2), nsmall = 2), " as the 2N peak value."))
  cat("\n")
  
  
  #This part bins each point's ploidy        
  cells$dna_norm <<- (cells$log2_dna+1)-diploid
  
  for (i in 1:nrow(cells)){
    if (cells$dna_norm[i] < 1.5) {
      cells$ploidy[i] <<- "2N"
    } else if (cells$dna_norm[i] > 2.5) {
      cells$ploidy[i] <<- ">4N"
    } else{
      cells$ploidy[i] <<- "4N"
    }
  }
  
  ckF <- readline(prompt = paste0("Is the filename ", thingee, "? (y/n) "))
  if (ckF == "y"){
    write.csv(cells, file = thingee, row.names = FALSE)
  } else{
    thingee <<-readline(prompt = 'What is the name of the file: ')
    write.csv(cells, file = thingee, row.names = FALSE)
  }
  px<<-"dna_norm"
  xname <<- "DNA content (Log 2)"
  py <<- "edu_norm"
  yname <<- "EdU content (Log 10)"
  changeAxis(0.5,3.5,0.5,2.5)
  grapho(cells)
}
#----------------------------------------------------------------
#this function just makes a histogram for you (yay!)

makehisto<-function(df=cells, va="log2_dna"){
  hh<<-ggplot(df, aes(x = va))+geom_density()+
    theme_classic()+
    theme(axis.line = element_line(color = "black", size = 1.5), 
          axis.ticks = element_line(color = "black", size = 1.5), 
          axis.text = element_text(size = 24, family = "sans", color = "black"),
          axis.title = element_text(size = 36, family = "sans", color = "black"))+
    xlab("Variable")
  print(hh)
}

#-----------------------------------------------------------------
parser <- function(df=cells, va){
  catNam <<- readline(prompt = "What is the name of this variable: ")
  catNum <<- as.integer(readline(prompt = "How many categories: "))
  cells$holder <<- "NA"
  catz <<- c()
  for (i in 1:catNum){
    catz[i]<-readline(prompt = paste0("What is the name of category ", i, " of ", catNam, ": "))
  }
  cat("\n")
  for (i in catz){
    owl <<- as.numeric(readline(prompt = paste0("What is the low value of ", i, ": ")))
    hige <<- as.numeric(readline(prompt = paste0("What is the high value of ", i, ": ")))
    for (j in 1:nrow(cells)){
      if(va[j] >= owl & va[j] <= hige){
        cells$holder[j]<<-i
      }
    }
    cat("\n")
  }
  if (catNam %in% colnames(cells)){
    cells[match(catNam, names(cells))] <<- NULL
    colnames(cells)[match("holder", names(cells))] <<- catNam
  } else{
    colnames(cells)[match("holder", names(cells))] <<- catNam
  }
  write.csv(cells, file = thingee, row.names = FALSE)
}

#------------------------------------------------------------------

#This function is designed to iterate through the "_cells" tagged csvs in a working directory and concatenate them together into single, csv
concater <- function(){
  ruSure <- readline(prompt = "Preparing to concanenate all csv's in the current working directory. Are you sure (y/n)? ")
  if (ruSure == "y"){
    setwd("./finished")
    path <- getwd()
    lister_cells <- dir(path, pattern = "_cells.csv")
    whole <- read.csv(lister_cells[1])
    whole$file <- lister_cells[1]
    for (i in lister_cells[2:length(lister_cells)]){
      parrt <- read.csv(i)
      parrt$file <- i
      whole <- rbind(parrt, whole)
    } 
    cat("Done!")
    cat("\n")
    fileName <- readline(prompt = "What would you like the filename to be: ")
    fileName <- paste0(fileName, ".csv")
    write.csv(whole, file = fileName, row.names = FALSE)
    whole <<- read.csv(fileName)
    cat(paste("Saved as", fileName))
  } else {
    cat("Okay, well try again later?")
  }
}

oneAndDone <- function(){
  path <- getwd()
  lister <- dir(path, pattern = ".csv")
  dir.create(path = paste0(path, "/finished"))
  for (i in 1:length(lister)){
    sirmixaplot(lister[i])
    print(lister[i])
    cells$infect <<- readline(prompt = "Infection: ")
    cells$cellType <<- readline(prompt = "Cell type: ")
    write.csv(cells, file = thingee, row.names = FALSE)
    file.rename(from = paste0(substr(thingee, 1, nchar(thingee)-10), ".csv"), to = paste0(path, "/finished/", paste0(substr(thingee, 1, nchar(thingee)-10), ".csv")))
    file.rename(from = paste0(substr(thingee, 1, nchar(thingee)-10), "_cells.csv"), to = paste0(path, "/finished/", paste0(substr(thingee, 1, nchar(thingee)-10), "_cells.csv")))
    cat("\n")
  }
  concater()
}

#-----------------------------------------------------------------
# So this is the compensating function in case (god forbid) you have spectral overlap. Make sure to run this having graphed the offending colors against eachother
# IMPORTANT - The skewed color must be graphed on the y axis

compensate <- function(df = cells){
  # First it generates a datframe to hold point. This is kept open to eventaully morph this into a shape fitting gate
  gridIron <- data.frame(X=c(1), Y=c(1))
  cat("Thank you for choosing coompensate(). To begin, lets assign a first point. Be sure the intended line follows the overlap line:")
  # Then you pick a point on the line of the overlapping population
  gridIron$X[1] <- as.numeric(readline(prompt = "What is the x-value of the first point? "))
  gridIron$Y[1] <- as.numeric(readline(prompt = "What is the x-value of the first point? "))
  cat("\n")
  print(qq+geom_point(data = gridIron, aes(x=gridIron$X, y=gridIron$Y, color = "red", size = 16)))
  cat(paste0("Great, I have added that point. Now lets move on to the second point"))
  gridIron <- rbind(gridIron, c(1, 1))
  # Then you pick a second point
  gridIron$X[2] <- as.numeric(readline(prompt = "What is the x-value of the second point? "))
  gridIron$Y[2] <- as.numeric(readline(prompt = "What is the x-value of the second point? "))
  print(qq+geom_point(data = gridIron, aes(x=gridIron$X, y=gridIron$Y, color = "red", size = 16))+
          geom_segment(aes(x = gridIron[1,1], y = gridIron[1,2], xend = gridIron[2,1], yend = gridIron[2,2], color = "red", size=16)))
  # Asks if these two points lie on the 
  #print(gridIron)
  yORn <<- readline(prompt = "Does this look right (y/n)? ")
  if (yORn == "n"){
    cat("Sorry, lets try that again.")
    cat("\n")
    compensate()
  } else {
    cat("Excellent, generating the linear equation now:")
    cat("\n")
    rise <- gridIron[2,2]-gridIron[1,2]
    run <- gridIron[2,1]-gridIron[1,1]
    slop <- rise/run
    yint <- gridIron[1,2]-(gridIron[1,1]*slop)
    print(qq+geom_point(data = gridIron, aes(x=gridIron$X, y=gridIron$Y, color = "red", size = 16))+
            geom_abline(slope = slop, intercept = yint, size = 2, color = "blue"))
    cat(paste0("Slope has been determined to be: ", slop))
    cat("\n")
    cat(paste0("Y-intercept has been determined to be: ", yint))
    yORn <- readline(prompt = "Does this look right (y/n)? ")
    if (yORn == "n"){
      cat("Sorry, lets try that again.")
      cat("\n")
      compensate()
    } else {
      cat("Applying compensation.")
      check <- df
      check[py] <- (check[py]-yint)/(slop*check[px])
      cat("Regraphing:")
      joiner(check)
      yORn <- readline(prompt = "Does this look right (y/n)? ")
      if (yORn == "n"){
        cat("Sorry, lets try that again.")
        cat("\n")
        compensate()
      } else {
        cells <<- check
        cat(paste0("Cells saved with compensated ", py))
        yORn <- readline(prompt = "Would you like to recalculate the means (y/n)? ")
        if (yORn == "y"){
          for (i in 1:ncol(cells)){
            if (grepl(py, names(cells)[i])){
              cells[paste0("I", names(cells)[i])] <<- cells[,i]*cells$Area
            }
          }
          for (i in 1:ncol(cells)){
            if (grepl(py, names(cells)[i])){
              cells[paste0("L", names(cells)[i])] <<- log(cells[,i]+1, 10)
            }
          }
          cat("Means recalculated. Thank you for your patience and patronage.")
          write.csv(cells, thingee, row.names = FALSE)
        }
      }
    }
  }
}

#-----------------------------------------------------------------
#Color palettes:

Color_blind<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb_black <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#The cellCycle_colors palette are the hexdecimal codes for: 
#                        G1,         S,       G2/M,    mitosis, & premature mitosis respecitvely, as used in Justice et al, 2019
cellCycle_colors <- c("#D4D4D4", "#98C84C", "#23B8CC", "#F16B1A", "#E5001C")


theme_general <- theme(legend.position = "right", 
                       legend.text = element_text(size = 15), 
                       legend.title = element_text(size = 18), 
                       axis.line = element_line(color = "black", size = 1.5), 
                       axis.ticks = element_line(color = "black", size = 1.5), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       axis.title = element_text(size = 36, family = "sans", color = "black"), 
                       axis.text = element_text(size = 24, family = "sans", color = "black"),
                       legend.key.size = unit(1.5, "cm")
)

cat("Welcome, SirMixaPlot! type `sirmixaplot(filename)' to get started, or call imaGen(directory) to generate an '_all.csv' file.")
