library(ggpubr)

glycoRun <- function(
  # Whole protein gel profile ladder file
  ladderWP = "ladder.csv",
  # Glycosylation protein gel profile ladder file
  ladderGP = "ladder.csv",
  # Whole protein profile ladder values
  ladderValuesWP = c(250, 150, 100, 75, 50, 37, 25, 20, 15, 10),
  # Glycolsyation ladder values
  ladderValuesGP = c(75, 25),
  # Whole protein gel profile file
  wp = "wp.csv",
  # Glycosylation gel profile file
  gp = "gp.csv",
  # Output file name
  outputCSV = "allData.csv",
  # Correct for background variation along the gel
  corrections = T,
  # Location of the peak data
  peakCSV = "csv/allData.csv",
  # Numerical threshold for peak calling
  givenThreshold = F,
  # The number of iterations attempted to get the ladder right before it gives up
  ladderAttempts = 10,
  # If the ladder attempts max out, should the peaks called be assumed as the upper or lower portion of the given ladder sizes?
  upperOrLowerLadder = "lower") {
  
  #Lists all the files in the main folder
  currentSetup <- list.files()
  # If 'figures', 'bands', or 'csv' aren't present, these directories get made
  if(!"figures" %in% currentSetup){
    dir.create("figures")
    pngList <- list.files(pattern = ".png")
    for (a in pngList){
      file.rename(a, paste0("figures/", a))
    }
  }
  if(!"bands" %in% currentSetup){
    dir.create("bands")
  }
  if(!"csv" %in% currentSetup){
    dir.create("csv")
  }
  
  # If the corrections option is TRUE and the ladder file is present, the data gets corrected
  correctionsCheck <- list.files()
  # First the ladder gets corrected
  if (corrections == T & ladderWP %in% correctionsCheck){
    # The csv's are listed
    cList <- list.files(pattern = ".csv")
    # The csvs are read into a dataframe
    pxWP <- read.csv(ladderWP)
    # The dataframe is given a name and type on ladder
    # For non-ladder samples, these are taken from the file name: 'name'_'type'.csv
    pxWP$name <- "Ladder"
    pxWP$type <- "Ladder"
    # Dataframe names are changed
    names(pxWP)[1] <- "position"
    names(pxWP)[2] <- "value"
    
    # An initial linear relationship is determined and the slope and intercepts are gathered
    sampleLine <- lm(value~position, data = pxWP)
    intercept <- sampleLine[1]$coefficients[1]
    slope <- sampleLine[1]$coefficients[2]
    # The value at each position is adjusted by the initial linear regression
    pxWP$adjValue <- pxWP$value-(slope*pxWP$position+intercept)
    
    # The adjusted values are used to generate a new regression
    # So long as the slope is above above 0.001, this process is iterated
    while (abs(slope) > 0.001){
      sampleLine <- lm(adjValue~position, data = pxWP)
      intercept <- sampleLine[1]$coefficients[1]
      slope <- sampleLine[1]$coefficients[2]
      pxWP$adjValue <- pxWP$adjValue-(slope*pxWP$position+intercept)
    } 
    
    if (ladderGP != ladderWP){
      pxGP <- read.csv(ladderGP)
      # The dataframe is given a name and type on ladder
      # For non-ladder samples, these are taken from the file name: 'name'_'type'.csv
      pxGP$name <- "Ladder"
      pxGP$type <- "Ladder"
      # Dataframe names are changed
      names(pxGP)[1] <- "position"
      names(pxGP)[2] <- "value"
      
      # An initial linear relationship is determined and the slope and intercepts are gathered
      sampleLine <- lm(value~position, data = pxGP)
      intercept <- sampleLine[1]$coefficients[1]
      slope <- sampleLine[1]$coefficients[2]
      # The value at each position is adjusted by the initial linear regression
      pxGP$adjValue <- pxGP$value-(slope*pxGP$position+intercept)
      
      # The adjusted values are used to generate a new regression
      # So long as the slope is above above 0.001, this process is iterated
      while (abs(slope) > 0.001){
        sampleLine <- lm(adjValue~position, data = pxGP)
        intercept <- sampleLine[1]$coefficients[1]
        slope <- sampleLine[1]$coefficients[2]
        pxGP$adjValue <- pxGP$adjValue-(slope*pxGP$position+intercept)
      }
    } else {
      pxGP <- pxWP
    }

    # The same process is done with the other csv files
    # First the WP gel
    mName <- strsplit(wp, "_")[[1]][1]
    interim <- read.csv(wp)
    interim$name <- mName
    interim$type <- "Whole Protein"
    names(interim)[1] <- "position"
    names(interim)[2] <- "value"
    sampleLine <- lm(value~position, data = interim)
    intercept <- sampleLine[1]$coefficients[1]
    slope <- sampleLine[1]$coefficients[2]
    interim$adjValue <- interim$value-(slope*interim$position+intercept)
    while (abs(slope) > 0.001){
      sampleLine <- lm(adjValue~position, data = interim)
      intercept <- sampleLine[1]$coefficients[1]
      slope <- sampleLine[1]$coefficients[2]
      interim$adjValue <- interim$adjValue-(slope*interim$position+intercept)
    }
    pxWP <- rbind(pxWP, interim)
    
    # Then the GP gel
    mName <- strsplit(gp, "_")[[1]][1]
    interim <- read.csv(gp)
    interim$name <- mName
    interim$type <- "Gylcosylated Protein"
    names(interim)[1] <- "position"
    names(interim)[2] <- "value"
    sampleLine <- lm(value~position, data = interim)
    intercept <- sampleLine[1]$coefficients[1]
    slope <- sampleLine[1]$coefficients[2]
    interim$adjValue <- interim$value-(slope*interim$position+intercept)
    while (abs(slope) > 0.001){
      sampleLine <- lm(adjValue~position, data = interim)
      intercept <- sampleLine[1]$coefficients[1]
      slope <- sampleLine[1]$coefficients[2]
      interim$adjValue <- interim$adjValue-(slope*interim$position+intercept)
    }
    if (ladderWP == ladderGP){
      pxGP <- rbind(pxWP, interim)
      tick <- "Good"
    } else {
      pxGP <- rbind(pxGP, interim)
      tick <- "Bad"
    }

    # The dataframe is saved as the output file
    write.csv(pxWP, paste0("csv/WP_", outputCSV), row.names = F)
    write.csv(pxGP, paste0("csv/GP_", outputCSV), row.names = F)
    # A seperate 'tag' ID is generated to make iterating through name/type subsets easier
    # This is not saved into any final csv
    px$tag <- paste0(px$name, "_", px$type)
    
    if (tick == "bad"){
      
    }
    
    # Now we make some initial graphs
    # First, we iterate through the names, which will have both 'types' of data in them
    for(a in unique(px$name)){
      # This makes sure we don't pull ladder data
      if (a != "Ladder"){
        # Only that name data is selected as an interim file
        interim <- subset(px, name == a | name == "Ladder")
        # A ggplot object is created using the interim data
        # This will be a line overlay between the types of data coom and glyco
        draft <- ggplot(data = interim, aes(
          x = position,
          y = adjValue,
          color = tag
        ))+
          geom_line()+
          theme_classic()+
          theme(legend.position = "top")+
          facet_wrap(~name, nrow = 2)
        print(draft)
        # The graph is saved
        ggsave(paste0("figures/",a, ".pdf"))
        
        # Now a correlation scatter plot is generated and saved
        interim <- subset(px, name == a)
        coom <- subset(interim, type == "coom")
        glyco <- subset(interim, type != "coom")
        interim <- data.frame("PositionFromTop" = coom$position)
        interim$coom <- coom$adjValue
        interim$glyco <- glyco$adjValue
        draft <- ggplot(data = interim, aes(
          x = coom,
          y = glyco,
        ))+
          geom_point(size = 4, aes(color = PositionFromTop))+
          geom_density_2d()+
          geom_abline(slope = 1, intercept = 0, size = 2, linetype = 2)+
          theme_classic()+
          xlab("Coomasie signal")+
          ylab("Glycosylation signal")
        print(draft)
        ggsave(paste0("figures/",a, "_signalCorr.pdf"))
        
        # Now we generate a glycosylation signal score by subtracting the glyco score from the coomassie score
        interim$coom <- interim$coom+abs(min(interim$coom))
        interim$glyco <- interim$glyco+abs(min(interim$glyco))
        interim$relGylcoSignal <- interim$glyco-interim$coom
        interim$relGylcoSignal <- interim$relGylcoSignal+abs(min(interim$relGylcoSignal))
        # And now we'll normalize it by lower quartile normalization because why not?
        interim$relGylcoSignal <- interim$relGylcoSignal/quantile(interim$relGylcoSignal)[2]
        
        #Now we make the super cool heat map histograms using the newly calculated diagrams
        ladder <- subset(px, name == "Ladder")
        interim$ladder <- ladder$adjValue+abs(min(ladder$adjValue))
        draft <- ggplot(data = interim, aes(
          x = PositionFromTop,
          y = coom,
          color = relGylcoSignal,
          fill = relGylcoSignal
        ))+
          scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = mean(interim$relGylcoSignal))+
          scale_fill_gradient2(low = "blue", mid = "green", high = "red", midpoint = mean(interim$relGylcoSignal))+
          geom_bar(stat = "identity")+
          theme_classic2()+
          theme(legend.position = "top")+
          xlab("Relative position")+
          ylab("Coomassie gel intensity")
        draft
        ggsave(paste0("figures/", a, "_relGlyco.pdf"))
        draft+geom_line(aes(y = ladder), color = "#5c5c5cff")
        ggsave(paste0("figures/", a, "_relGlyco_ladder.pdf"))
        write.csv(interim, paste0("csv/", a, ".csv"), row.names = F)
      }
    }
  } else {
    # If the correction parameter FALSE, then the data is assumed corrected and opened from its save location for peak calling
    px <- read.csv(peakCSV)
    # The temp tag IDs are given for peak calling iteration
    px$tag <- paste0(px$name, "_", px$type)
  }

  #Peaks are called in the various tag ID sets
  for (a in unique(px$tag)){
    # interim file based on tag ID is generated
    tagSet <- subset(px, tag == a)
    print(paste0("Calling peaks on ", a))
    # If the interim set is from the ladder, a special peak calling is done
    if (a == "Ladder_Ladder"){
      peakCaller(df = tagSet, dfName = "ladder", 
                 threshold = F,
                 expectedPeaks = T,
                 expectedNumber = length(ladderValues),
                 ladder = ladderValues,
                 attemptsForPeaks = ladderAttempts, 
                 ladderUpOrDown = upperOrLowerLadder)
      
      # Uses the bands on the ladder to calculate kD based on position but this doesn't really work... great
      # In the future, linear regression should be done between points, but this is... okay
      ladderBands <- read.csv("bands/ladder.csv")
      ladderEq <- lm(size~position, data = ladderBands)
      intercept <- ladderEq$coefficients[1]
      slope <- ladderEq$coefficients[2]
      
      # Draws a QC plot of the regression model
      draft <- ggplot(data = ladderBands, aes(
        x = position,
        y = size
      ))+geom_point(size = 4)+
        geom_smooth(method = "lm",)+
        theme_classic()
      print(draft)
      # Saves the QC plot
      ggsave("figures/ladderQC.pdf")
      
      # Saves the size data of the whole dataset
      px$size <- px$position*slope+intercept
      write.csv(px, peakCSV, row.names = F)
      
    } else {
      # Peak calling is done on the data based on the entered parameters
      peakCaller(df = tagSet, 
                 dfName = a, 
                 threshold = givenThreshold,
                 expectedPeaks = F,
                 expectedNumber = length(ladderValues),
                 ladder = ladderValues)
      # peaks are read in from the peak calling program
      peakBands <- read.csv(paste0("bands/", a, ".csv"))
      # peak sizes are calculated
      peakBands$size <- peakBands$position*slope+intercept
      # New file with developed
      write.csv(peakBands, paste0("bands/", a, ".csv"), row.names = F)
    }
  }
  
  # Now peaks are assinged between the coomassie and glycosylation data
  print("Combining called peaks by name...")
  banderSnatch()
  
  for(a in unique(px$tag)){
    tagSet <- subset(px, tag == a)
    
    # Redraws the profile plot by size rather than position
    draft <- ggplot(data = tagSet,
                    aes(x = size,
                        y = adjValue))+
      geom_line()+
      theme_classic()
    print(draft)
    ggsave(paste0("figures/", a, "_sized.pdf"))
  }
}

peakCaller <- function(df = px,
                       dfName,
                       threshold = 0,
                       expectedPeaks = F,
                       expectedNumber = 9,
                       ladder,
                       attemptsForPeaks=10,
                       ladderUpOrDown = "lower"){
  peaks <- df
  if(threshold != F){
    print(paste0("Threshold given of ", threshold, "."))
    peaks <- subset(peaks, adjValue > threshold)
    
    bandNumber <- 1
    peaks$band <- bandNumber
    for(a in 2:nrow(peaks)){
      if (peaks[a,]$position != peaks[a-1,]$position+1){
        bandNumber <- bandNumber+1
      }
      peaks[a,]$band <- bandNumber
    }
    print(paste0("Total peaks found: ", bandNumber, "."))
    
  } else {
    if(expectedPeaks == T){
      print(paste0("Peaks expected: ", expectedNumber))
      holder <- peaks
      threshold <- mean(peaks$adjValue)+sd(peaks$adjValue)
      peaks <- subset(peaks, adjValue > threshold)
      
      bandNumber <- 1
      peaks$band <- bandNumber
      for(a in 2:nrow(peaks)){
        if (peaks[a,]$position != peaks[a-1,]$position+1){
          bandNumber <- bandNumber+1
        }
        peaks[a,]$band <- bandNumber
      }
      
      if (max(peaks$band) != expectedNumber){
        counter <- 1
        print("Incorrect peaks found. Exploring other thresholds...")
        while (max(peaks$band) != expectedNumber){
          if (max(peaks$band) > expectedNumber){
            threshold <- threshold+1
          } else {
            threshold <- threshold-1
          }
          print(paste0("Setting threshold at ", round(threshold, 2)))
          peaks <- holder
          peaks <- subset(peaks, adjValue > threshold)
          bandNumber <- 1
          peaks$band <- bandNumber
          for(a in 2:nrow(peaks)){
            if (peaks[a,]$position != peaks[a-1,]$position+1){
              bandNumber <- bandNumber+1
            }
            peaks[a,]$band <- bandNumber
          }
          print(paste0("Found ", bandNumber, " bands"))
          print(paste0("Expected: ", expectedNumber))
          if(bandNumber==expectedNumber){
            break
          }
          counter <- counter+1
          if (counter > attemptsForPeaks){
            print("Max number of attempts reached...")
            print(paste0("Using ", bandNumber, " bands set for the ", ladderUpOrDown, " set of bands."))
            break
          }
        }
      }
    } else {
      print("Threshold not given.")
      threshold <- mean(peaks$adjValue)+sd(peaks$adjValue)
      print(paste0("Using mean+sd threshold of: ",threshold,"."))
      peaks <- subset(peaks, adjValue > threshold)
      
      bandNumber <- 1
      peaks$band <- bandNumber
      for(a in 2:nrow(peaks)){
        if (peaks[a,]$position != peaks[a-1,]$position+1){
          bandNumber <- bandNumber+1
        }
        peaks[a,]$band <- bandNumber
      }
      print(paste0("Total peaks found: ", bandNumber))
    }
  }
  
  draft <- ggplot(data = df, aes(
    x = position,
    y = adjValue,
  ))+geom_line()+
    geom_hline(yintercept = threshold, size = 1, color = "red", linetype = 2)+
    theme_classic()
  
  print(draft)
  ggsave(paste0("figures/", dfName, "_threshold.pdf"))
  
  for (peakID in unique(peaks$band)){
    interim <- subset(peaks, band == peakID)
    if(!exists("bands")){
      bands <- data.frame(
        "name" = unique(interim$name),
        "type" = unique(interim$type),
        "bandID" = peakID,
        "position" = median(interim$position),
        "width" = nrow(interim),
        "intensity" = sum(interim$adjValue),
        "size" = 0
      )
    } else {
      turkey <- data.frame(
        "name" = unique(interim$name),
        "type" = unique(interim$type),
        "bandID" = peakID,
        "position" = median(interim$position),
        "width" = nrow(interim),
        "intensity" = sum(interim$adjValue),
        "size" = 0
      )
      bands <- rbind(bands, turkey)
    }
    
  }
  if (dfName == "ladder"){
    
    bands$size <- 0
    if (bandNumber != length(ladder)){
      if (bandNumber > length(ladder)){
        bands <- bands[1:length(ladder),]
      }
      if (ladderUpOrDown == "lower"){
        ladder <- ladder[(1+length(ladder)-nrow(bands)):length(ladder)]
      } else {
        ladder <- ladder[1:length(bands)]
      }
    }
    for (a in 1:length(ladder)){
      bands[bands$bandID == a,]$size <- ladder[a]
    }
  }
  write.csv(bands, paste0("bands/", dfName, ".csv"), row.names = F)
}

banderSnatch <- function(){
  bList <- list.files(path = "bands", pattern = "csv")
  for (a in bList){
    if (!grepl("bands", a) & !grepl("ladder", a)){
      a <- paste0("bands/", a)
      if(!exists("allBands")){
        allBands <- read.csv(a)
      } else {
        someBands <- read.csv(a)
        allBands <- rbind(allBands, someBands)
      }
      for (gName in unique(allBands$name)){
        interim <- subset(allBands, name == gName)
        write.csv(interim, paste0("bands/",gName, "_bands.csv"), row.names = F)
        draft <- ggplot(data = interim, aes(
          x = size,
          y = intensity,
          color = type
        ))+
          geom_point(size = 4)+
          theme_classic()
        print(draft)
        ggsave(paste0("figures/",gName, "_bands.pdf"))
      }
    }
  }
}
