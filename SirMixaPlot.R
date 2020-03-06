#SirMixaPlot version 8.2 - Danger! High Voltage

#----------------------------------------------------------------------------------------------------------

#This program is designed to take a csv (comma separated values) file and generate scatterplots based on columns.
# Changes from version 8.1:
#   -Changed the defualt EdU neg setting on normalizer()
#   -Added a work-around for concater() should rbinding fail due to missing columns
#   -imaGen() now also creates a concatanated csv of all the ROIs in the WholeCell folder with the _ROI_ tag
#   -Added runPCA() function to take a dataframe and generate a PCA biplot; still in beta testing

#Changes from SirMixaPlot version 8.0:
#   -Added some comments and typos
#   -Added the reMap() function for re-creating a facsimile of the original image

#Changes from SirMixaPlot version 7.0:
#   -Takes data in from the YggData.imj macro
#   -imaGen() is now much bigger and generates combined CSVs for nuclear vs whole cell comparisons

#Things still that need to be coded:
#   -Exploratory data analyses function called explore() made

#-----------------------------------------------------------------------------------------------------------

#Program opens required packages

#If these packages aren't install, uncomment the following lines to install them.

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("MASS")
#install.packages("viridis")
#install.packages("stringr")
#devtools::install_github("vqv/ggbiplot")

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(stringr))
#suppressPackageStartupMessages(library(ggbioplot))

#------------------------------------------------------------------------------------------------------------

# Functions for your pleasure

#--------------------------------------------------------------------------------------

#The imaGen() script takes any number of csv files that are within your working directory and combines them. Take caution, however, as this script is designed for a particular imagej csv outpt
#It is important that the jetData uses the Set Measurements with the following settings active:
#   Area    Standard Deviation    Min & Max gray value    Center of Mass    Mean gray value   Perimenter    Display label
# Additionally, for the best results, decimal places should be set to 9

# If you have used the YggData macro, this will all be automatically set for you

#This is the joinR function used to assign nuclear number to each ROI that may or may not be localized to the nucleus
joinR <- function(y, x = "dna.csv"){
  # Files are opened
  tagger <- read.csv(x)
  fraction <- read.csv(y)
  # The number column is named, if not present
  if (!"Number" %in% names(tagger)){
    names(tagger)[1] <- "Number"
    write.csv(tagger, file = x, row.names = FALSE)
  }
  if (!"Number" %in% names(fraction)){
    names(fraction)[1] <- "Number"
    write.csv(fraction, file = y, row.names = FALSE)
  }
  #Default cellid is given
  fraction$cell <- "unknown"
  #Default distance is given
  tagger$distance <- 0
  # The closest nucleus is calculated for each ROI
  for (i in 1:nrow(fraction)){
    tagger$distance <- sqrt((tagger$X-fraction$XM[i])^2+(tagger$Y-fraction$YM[i])^2)
    closest <- min(tagger$distance)
    fraction$cell[i] <- subset(tagger, distance == closest)$Number[1]
    #If two or more nuclei are picked, it lets you know
    if (length(closest) > 1){
      print("PING! More than one nucleus at minimum distance. Something is wrong (probably)")
    }
  }
  write.csv(fraction, file = y, row.names = FALSE)
}

imaGen <- function(directory="./", colorz = T){
  # The directory is set and a list of CSV's is generated
  directoryN <- paste0(directory, "/Nuclear/")
  directoryC <- "WholeCell/"
  setwd(directoryN)
  filez <- list.files(pattern = ".csv")
  # The first file is opened to serve as a template
  cells <- read.table(filez[1], sep = ",", header = TRUE)
  # If the csv name is not simply the target name (colorz == F), then you will be asked to name a target for each image
  # REMEMBER: One of the targets must be 'dna'
  if(colorz == F){
    cat("REMEMBER: One of these color's must be labeled 'dna'")
    cat("\n")
    cat(paste0("This file is ", filez[1]))
    cat("\n")
    colo <- readline(prompt= paste0("What color is in ", filez[1], ": "))
  } else{
    # Otherwise, the target is taken from the csv name
    colo <- substr(filez[1], 1, nchar(filez[1])-4)
  }
  cells$Area <- cells$Area*100
  names(cells) <- paste0(names(cells), "_NUC_", colo)
  # The first column with designated as the cell number which is dictated by the nucleus and is consistent across all downstream applications
  names(cells)[1] <- "Number"
  # The label is removed as this isn't required after 'colo' labeling
  cells<-cells[,-2]
  # Default area is quite small, so its boosted to make the numbers more "real"
  # The rest of the csv's in the Nuclear directory are opened and cbinded to the first set
  for (i in filez[2:length(filez)]){
    # Already bound csv's (_all tagged) are ignored
    if (!grepl("_all", i)){
      interim <- read.table(i, sep = ",", header = TRUE)
      if(colorz == F){
        cat("REMEMBER: One of these color's must be labeled 'dna'")
        cat("/n")
        cat(paste0("This file is ", i))
        cat("/n")
        colo <- readline(prompt= paste0("What color is in ", i, ": "))
      } else{
        colo <- substr(i, 1, nchar(i)-4)
      }
      names(interim) <- paste0(names(interim), "_NUC_", colo)
      cells <- cbind(cells, interim[,4:8], interim[,11:12], interim[,15:16])
    }
  }
  filnam <- readline(prompt = "What should the base file be named: ")
  write.csv(cells, file = paste0(filnam,"_NUC_all.csv"), row.names = FALSE)
  xoo <- list.files(pattern = "dna")
  setwd("../")
  
  #Now that the nuclear data is collected, the whole data will be combined
  
  # First, you navigate to the wholeCell folder
  setwd(directoryC)
  # Then you get the list of files
  filez <- list.files()
  # Then you designate which tagger file to use. If colorz == T, the script assumes your file is called dna.csv and uses the same csv as found in the Nuclear folder
  if(colorz==F){
    pathy <- getwd()
    cat(paste0("Current directory is ", pathy))
    cat("\n")
    tagger <- readline(prompt = "What is the path to the dna file: ")
  } else{
    tagger <- paste0("../Nuclear/", xoo)
  }
  # Every ROI is designated a nucleus number
  for (i in filez){
    joinR(i, xoo)
  }
  # The tagger file is opened
  cells <- read.csv(tagger)
  # Default area is quite small, so its boosted to make the numbers more "real"
  cells$Area <- cells$Area*100
  # The dna values are labeled as such
  names(cells) <- paste0(names(cells), "_WC_dna")
  # The first column with designated as the cell number which is dictated by the nucleus and is consistent across all downstream applications
  names(cells)[1] <- "Number"
  # The label is removed as this isn't required after 'colo' labeling
  cells<-cells[,-2]
  #Each ROI file is summarized and added to the nucleus data
  for (i in filez){
    if(colorz == F){
      cat("REMEMBER: One of these color's must be labeled 'dna'")
      cat("\n")
      cat(paste0("This file is ", i))
      cat("\n")
      colo <- readline(prompt= paste0("What color is in ", i, ": "))
    } else{
      colo <- substr(i, 1, nchar(i)-4)
    }
    # Temporary dataset is made for easy column labeling
    interim <- cells
    # The ROI files is opened
    fraction <- read.csv(i)
    for (j in 1:nrow(fraction)){
      cellP <- subset(interim, Number == fraction$cell[j])
      fraction$NucDist[j] <- sqrt((fraction$XM[j] - mean(cellP$X_WC_dna))^2 + (fraction$YM[j] - mean(cellP$Y_WC_dna))^2)
    }
    write.csv(fraction, file = i, row.names = F)
    #Each nucleus number is analyzed for ROI's designated to it
    for (j in 1:nrow(interim)){
      interim$ROI_Num[j] <- 0
      interim$ROI_Area[j] <- 0
      interim$ROI_Mean[j] <- 0
      interim$ROI_Stdev[j] <- 0
      interim$ROI_Mode[j] <- 0
      interim$ROI_Perimeter[j] <- 0
      interim$ROI_IntDens[j] <- 0
      interim$ROI_IntTotal[j] <- 0
      interim$ROI_NucDist[j] <- 0
      if(interim$Number[j] %in% unique(fraction$cell)){
        # The ROIs that are assigned to a particular nucleus number are subsetted
        hitz <- subset(fraction, cell == interim$Number[j])
        
        # Individual ROI nuclear distances are calculated
        #The number of ROIS assigned to that cell
        interim$ROI_Num[j] <- nrow(hitz)
        
        #The average average of all the ROIs
        interim$ROI_Area[j] <- mean(hitz$Area)*100
        
        #The average Mean intensity of all the ROIs
        interim$ROI_Mean[j] <- mean(hitz$Mean)
        
        #The average standard deviation of all the ROIs
        interim$ROI_Stdev[j] <- mean(hitz$StdDev)
        
        #The average mode of pixel intensities of all the ROIs
        interim$ROI_Mode[j] <- mean(hitz$Mode)
        
        #The average perimeter of all the ROIs
        interim$ROI_Perimeter[j] <- mean(hitz$Perim.)
        
        #The average integrated mean/density of the ROIs
        interim$ROI_IntDens[j] <- mean(hitz$IntDen)*100
        
        #The total of the integrated means/densities of all the ROIs
        interim$ROI_IntTotal[j] <- sum(hitz$IntDen)*100
        
        #The average distance from the nucleus of all the ROIs
        interim$ROI_NucDist[j] <- mean(hitz$NucDist)
      } 
    }
    # A target tag is added for later identification
    names(interim) <- paste0(names(interim), "_", colo)
    # The new variables are added to their cell ID
    cells <- cbind(cells, interim[,(ncol(cells)+1):ncol(interim)])
  }
  #filnam <- readline(prompt = "What should the whole cell file be named: ")
  interim <- read.csv(filez[1])
  colo <- substr(filez[1], 1, nchar(i)-4)
  interim$roi <- colo
  for (i in filez[2:length(filez)]){
    x <- read.csv(i)
    colo <- substr(i, 1, nchar(i)-4)
    for (j in names(interim)){
      if (!j %in% names(x)){
        x[j] <- "NA"
      }
    }
    for (j in names(x)){
      if (!j %in% names(interim)){
        interim[j] <- "NA"
      }
    }
    x$roi <- colo
    interim <- rbind(interim, x)
  }
  write.csv(interim, file = paste0(filnam, "_ROI_all.csv"), row.names = F)
  write.csv(cells, file = paste0(filnam, "_WC_all.csv"), row.names = F)
  setwd("../")
  
  #Now that we have a nuclear dataset and a WC dataset, we might as well put it all together...
  nucka <- paste0("Nuclear/", filnam, "_NUC_all.csv")
  nuke <- read.csv(nucka)
  wucka <- paste0("WholeCell/", filnam, "_WC_all.csv")
  wuke <- read.csv(wucka)
  cells <- cbind(nuke, wuke[,22:ncol(wuke)])
  write.csv(cells, file = paste0(filnam, "_WN_all.csv"), row.names = F)
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
grapho <- function(df=cells, X=px, Y=py, Xn = xname,  Yn=yname){
  
  #backup dataframe is generated
  interim <<- cells

  cat("Data frame generated")
  cat("\n")
  
  #Initial scatter plot is generated
  qq<<-ggplot(data=df, aes_string(X, Y))+
    theme_linedraw()+
    theme(plot.title = element_text(size=40, face="bold"))+
    theme(axis.title = element_text(size=30, face="bold"))+
    geom_point(size=2)+
    scale_x_continuous(name=Xn)+
    scale_y_continuous(name=Yn)+
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
  cat("Call 'cull()' to remove certain populations from original dataset (DANGER, but less than before)")
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
    cells$log2_dna <<- log(cells$IntDen_NUC_dna, 2)
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
  print(sort(names(cells)))
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
  
  workaround<<- subset(df, as.numeric(get(px)) > as.numeric(xSet[1]) & as.numeric(get(px)) < as.numeric(xSet[2]))
  workaround <<- subset(workaround, as.numeric(get(py)) > as.numeric(ySet[1]) & as.numeric(get(py)) < as.numeric(ySet[2]))
  assign(gateName, workaround, envir = .GlobalEnv)
  
  cells$gatePop<<-cells$Number %in% workaround$Number
  
  if(gateName == "negPop"){
    negPop_correct()
  } else{
    colnames(cells)[match("gatePop", names(cells))]<<-gateName
  }
}

#This function uses a true negative population to adjust the background
negPop_correct <- function(){
  x <<- strsplit(px, "_")
  x <- paste0(x[[1]][length(x[[1]])-1], "_", x[[1]][length(x[[1]])])
  y <- strsplit(py, "_")
  y <- paste0(y[[1]][length(y[[1]])-1], "_", y[[1]][length(y[[1]])])
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
  write.csv(cells, file = thingee, row.names = FALSE)
  if (y == "NUC_edu" | x == "NUC_edu"){
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
  bigBin <- which.max(density(cells$ALIMean_NUC_edu)$y)
  edu_neg <- density(cells$ALIMean_NUC_edu)$x[bigBin]+0.5
  edu_hist <- ggplot(cells, aes(ALIMean_NUC_edu))+geom_density()+geom_vline(xintercept = edu_neg)+xlab("Log EdU")
  print(edu_hist)
  g1_good <- readline(prompt = paste0("Is ", format(round(edu_neg, 2), nsmall = 2), " representative of the EdU cutoff? (y/n) "))
  if (g1_good == "n"){
    edu_neg = as.numeric(readline(prompt = "What value should EdU be cutoff as negative? "))
    edu_hist <- ggplot(cells, aes(ALIMean_NUC_edu))+geom_density()+geom_vline(xintercept = edu_neg)+xlab("Log EdU")
    print(edu_hist)
  }
  cat(paste0("Using ", format(round(edu_neg, 2), nsmall = 2), " as the EdU cutoff."))
  cat("\n")
  for (i in 1:nrow(cells)){
    if (cells$ALIMean_NUC_edu[i] > edu_neg){
      cells$edu[i] <<- "Positive"
    } else {
      cells$edu[i] <<- "Negative"
    }
  }
  eduNegCells <-subset(cells, ALIMean_NUC_edu < edu_neg)
  cells$edu_norm <<- (cells$ALIMean_NUC_edu+1)-mean(eduNegCells$ALIMean_NUC_edu)
  
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
  hh<<-ggplot(df, aes(x = df[va]))+geom_density()+
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
    lister_cells <- dir(path, pattern = "all.csv")
    whole <- read.csv(lister_cells[1])
    whole$file <- lister_cells[1]
    for (i in lister_cells[2:length(lister_cells)]){
      parrt <- read.csv(i)
      parrt$file <- i
      for (j in names(whole)){
        if (!j %in% names(parrt)){
          parrt[j] <- "NA"
        }
      }
      for (j in names(parrt)){
        if (!j %in% names(whole)){
          whole[j] <- "NA"
        }
      }
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

# oneAndDone() can be performed in a parent direcroty containing subdirectories that each contain 1 or more  all.csv tagged files
oneAndDone <- function(directory = "./"){
  setwd(directory)
  dirs <- list.files()
  dirs <- dirs[1:length(dirs)]
  dir.create(path = "./finished")
  for (j in dirs){
    setwd(j)
    path <- getwd()
    lister <- dir(path, pattern = "all.csv")
    sirmixaplot(lister[1])
    print(lister[1])
    
    #Alter/comment these lines to add identifier information for the parent directory
    cells$infect <<- readline(prompt = "What is the infection type of this file: ")
    cells$knockdown <<- readline(prompt = "What is the knockdown type of this file: ")
    cells$treatment <<- readline(prompt = "What is the treatment type of this file: ")
    cells$bacteria <<- readline(prompt = "What is the bacteria type of this file: ")
    cells$time <<- readline(prompt = "What is the timepoint of this file: ")
    
    setwd("../finished")
    write.csv(cells, file = thingee, row.names = FALSE)
    setwd("../")
    #file.rename(from = paste0(substr(thingee, 1, nchar(thingee)-10), ".csv"), to = paste0("finished/", paste0(substr(thingee, 1, nchar(thingee)-10), ".csv")))
    #file.rename(from = paste0(substr(thingee, 1, nchar(thingee)-10), "_cells.csv"), to = paste0(path, "/", paste0(substr(thingee, 1, nchar(thingee)-10), "_cells.csv")))
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
  gridIron$Y[1] <- as.numeric(readline(prompt = "What is the y-value of the first point? "))
  cat("\n")
  print(qq+geom_point(data = gridIron, aes(x=gridIron$X, y=gridIron$Y, color = "red", size = 16)))
  cat(paste0("Great, I have added that point. Now lets move on to the second point"))
  gridIron <- rbind(gridIron, c(1, 1))
  # Then you pick a second point
  gridIron$X[2] <- as.numeric(readline(prompt = "What is the x-value of the second point? "))
  gridIron$Y[2] <- as.numeric(readline(prompt = "What is the y-value of the second point? "))
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
      check[py] <- check[py]-(slop*check[px]+yint)
      cat("Regraphing:")
      grapho(check)
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

#------------------------------------------------------

# This function take the nuclear locations and sizes and re-maps the original image
#   df is the dataframe, defaulted to 'cells'
reMap <- function(df = cells, 
                  #   highlight is a vector containing the Number(s) to be highlighted in remapped image
                  highlight = F, 
                  #   X and Y refer to the beginning of the X and Y position names in the dataframe, defaulted to the nuclear geometric position
                  X = "X_NUC_", 
                  Y = "Y_NUC_", 
                  #   S is the size of the ROI, defaulted to nucelar size
                  S = "Area_NUC_", 
                  #   Xn and Yn are the axis names and don't really need to be changed
                  Xn = "X position", 
                  Yn = "Y position",
                  #   colz is a vector of colors to use for the highlights variable. Be sure to have enough colors for the number of populations in the highlights variable
                  colz = cb_black) {
  #These set the variable names using regex to anchor tot he front of the variable name
  X <- names(df)[str_detect(names(df), paste0("^",X))]
  Y <- names(df)[str_detect(names(df), paste0("^",Y))]
  S <- names(df)[str_detect(names(df), paste0("^",S))]
  #This creates a ggplot object for graphing
  test <- ggplot(data = df, aes(y = unlist(df[X]), x = unlist(df[Y]), size = unlist(df[S])))
  #If the highlight variable is a vector and is not F, its graphed here
  if (is.vector(highlight) & highlight != F) {
    #print(1)
    check <- df
    check$highlight <- check$Number %in% highlight
    qq <<- test+geom_point(aes(color = check$highlight))+
      ylab(Yn)+
      xlab(Xn)+
      scale_color_manual(values = colz)+
      theme_classic()+
      theme(axis.text = element_text(face='bold', size=16), axis.ticks = element_blank())+
      theme(plot.title = element_text(size=40, face="bold"))+
      theme(axis.title = element_text(size=30, face="bold"))
  } 
  # If highlight is a not a vector (such as variable called) and isn't F, then its graphed
  else if (length(highlight) != 1 & highlight != F){
    #print(2)
    qq <<- test+geom_point(aes(color = highlight))+
      ylab(Yn)+
      xlab(Xn)+
      theme_classic()+
      theme(axis.text = element_text(face='bold', size=16), axis.ticks = element_blank())+
      theme(plot.title = element_text(size=40, face="bold"))+
      theme(axis.title = element_text(size=30, face="bold"))+
      scale_color_manual(values = colz)
  }
  # Finally, if highlight IS FALSE, then its graphed without a color option
  else {
    #print(3)
    qq <<- test+geom_point()+
      ylab(Yn)+
      xlab(Xn)+
      theme_classic()+
      theme(axis.text = element_text(face='bold', size=16), axis.ticks = element_blank())+
      theme(plot.title = element_text(size=40, face="bold"))+
      theme(axis.title = element_text(size=30, face="bold"))
  }
  #The plot is graphed
  print(qq)
}

#-----------------------------------------------------------------
# Beta run of PCA function

runPCA <- function(df=cells, 
                   choicez = c(1,2), 
                   subz = F, 
                   groupz = F, 
                   circz=F,
                   axe=F,
                   elly=F){  
  
  # get rid of non-numerics
  num <- unlist(lapply(df, is.numeric))
  cellist <- df[ , num]
  
  # Get rid of columns with no variance
  x <- foo(cellist)
  cellist <- cellist[,-x]
  
  # Get rid of NAs and Infs
  cellist[is.na(cellist)] <- 0
  cellist[cellist == Inf] <- 0
  cellist[cellist == -Inf] <- 0 
  
  # Get only the nuclear columns & remove NA
  if (subz != F){
    cellist <- cellist[str_detect(names(cellist), subz) | str_detect(names(cellist), "Number")]
  }
  
  # Make the pca object
  cell.pca <<- prcomp(cellist, center = TRUE, scale = TRUE)

  #Report how much of the data is shown
  howManyLeft <- round((nrow(cellist)*ncol(cellist))/(nrow(df)*ncol(df))*100, digits = 1)
  cat(paste0("Around ", howManyLeft, "% of the data is represented after trimming."))
  cat("\n")
  
  # Make the biplot
  if (groupz ==  T){
    elly <- T
    print(names(cellist))
    grope <- readline(prompt = "What should the PCA be grouped by from the 'cellist' dataframe: ")
    n <- which(colnames(cellist)==grope)
    #if (is.character(cellist[grope])){
    #  cat("Converting group to numeric")
    #  cellist[grope] <- as.list(cellist[grope])
    #}
    pp <<- ggbiplot(cell.pca,
                    ellipse = elly, 
                    choices = choicez, 
                    var.axes = axe, 
                    circle = circz, 
                    groups = as.factor(cellist[,n]))+
      #geom_density_2d(aes(color = cellist[grope]))+
      theme_classic2()
  } else {
    pp <<- ggbiplot(cell.pca, 
                    ellipse = elly, 
                    choices = choicez, 
                    var.axes = axe, 
                    circle = circz)+
    geom_density_2d()+
    theme_classic2()
  }
  if (elly == F){
    print(pp)
  } else {
    print(pp)
  }
}

foo <- function(dat) {
  out <- lapply(dat, function(x) length(unique(x)))
  want <- which(!out > 1)
  unlist(want)
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
