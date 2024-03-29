library(ggpubr)

glycoRun <- function(
  # Whole protein gel profile ladder file
  ladderWP = "ladder_wp.csv",
  # Glycosylation protein gel profile ladder file
  ladderGP = "ladder_wp.csv",
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
  
  print("Running GlycoRunner...")
  cat("\n")
  #Lists all the files in the main folder
  currentSetup <- list.files()
  # If 'figures', 'bands', or 'csv' aren't present, these directories get made
  if(!"figures" %in% currentSetup){
    print("Setting up directories...")
    dir.create("figures")
    pngList <- list.files(pattern = ".png")
    for (a in pngList){
      file.rename(a, paste0("figures/", a))
    }
  }
  # if(!"bands" %in% currentSetup){
  #   dir.create("bands")
  # }
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
    # Dataframe names are changed
    names(pxWP)[1] <- "position"
    names(pxWP)[2] <- "Ladder_value"

    # An initial linear relationship is determined and the slope and intercepts are gathered
    sampleLine <- lm(Ladder_value~position, data = pxWP)
    intercept <- sampleLine[1]$coefficients[1]
    slope <- sampleLine[1]$coefficients[2]
    # The value at each position is adjusted by the initial linear regression
    pxWP$adjValue <- pxWP$Ladder_value-(slope*pxWP$position+intercept)
    
    # The adjusted values are used to generate a new regression
    # So long as the slope is above above 0.001, this process is iterated
    while (abs(slope) > 0.001){
      sampleLine <- lm(adjValue~position, data = pxWP)
      intercept <- sampleLine[1]$coefficients[1]
      slope <- sampleLine[1]$coefficients[2]
      pxWP$adjValue <- pxWP$adjValue-(slope*pxWP$position+intercept)
    }
    pxWP$Ladder_value <- pxWP$adjValue
    pxWP <- pxWP[,1:2]

    # The same process is done with the other csv files
    # First the WP gel
    interim <- read.csv(wp)
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
    write.csv(interim, "csv/wp_corrected.csv", row.names = F)
    # The corrected WP points are then bound to the WP ladder dataframe
    if (nrow(interim) != nrow(pxWP)){
      print("The WP ladder is not the same length as the WP sample file. Script will now self-destruct...")
    }
    pxWP$WP_value <- interim$adjValue

    # Now we do the same for the GP ladder and sample
    # If the ladders are different files, then the ladders and samples will be corrected
    if (ladderGP != ladderWP){
      print("Your ladder files are different and this version isn't fully vetted. Prepare for possible failure.")
      cat("\n")
      # The GP ladder is read into a CSV
      pxGP <- read.csv(ladderGP)
      # Dataframe names are changed
      names(pxGP)[1] <- "position"
      names(pxGP)[2] <- "Ladder_value"
      
      # An initial linear relationship is determined and the slope and intercepts are gathered
      sampleLine <- lm(Ladder_value~position, data = pxGP)
      intercept <- sampleLine[1]$coefficients[1]
      slope <- sampleLine[1]$coefficients[2]
      # The value at each position is adjusted by the initial linear regression
      pxGP$adjValue <- pxGP$Ladder_value-(slope*pxGP$position+intercept)
      
      # The adjusted values are used to generate a new regression
      # So long as the slope is above above 0.001, this process is iterated
      while (abs(slope) > 0.001){
        sampleLine <- lm(adjValue~position, data = pxGP)
        intercept <- sampleLine[1]$coefficients[1]
        slope <- sampleLine[1]$coefficients[2]
        pxGP$adjValue <- pxGP$adjValue-(slope*pxGP$position+intercept)
      }
      pxGP$Ladder_value <- pxGP$adjValue
      pxGP <- pxGP[,1:2]

      # Then the GP gel is read, corrected, and bound to the GP ladder
      interim <- read.csv(gp)
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
      # The corrected GP values are saved
      write.csv(interim, "csv/gp_corrected.csv", row.names = F)
      if (nrow(interim) != nrow(pxGP)){
        print("The GP ladder is not the same length as the GP sample file. Script will now self-destruct...")
      }
      # The corrected GP values are saved to the GP ladder
      pxGP$GP_value <- interim$adjValue

      #Now we have to correct everything
      # First, lets assign ladder peaks to their corresponding values
      print("Calling peaks on the whole protein ladder...")
      pxWP <- peakCaller(df = pxWP,
                 value = "Ladder_value",
                 threshold = F,
                 expectedPeaks = T,
                 expectedNumber = length(ladderValuesWP),
                 ladder = ladderValuesWP,
                 attemptsForPeaks = ladderAttempts, 
                 ladderUpOrDown = upperOrLowerLadder)
      # Now that ladder peaks are assigned, we'll fill out the rest of the positions by linear regression between peaks
      # First the pxWP values
      ladderValuesWP <- unique(pxWP[pxWP$size != 0,]$size)
      for (peakz in 1:length(ladderValuesWP)){
        x1 <- pxWP[pxWP$size == ladderValuesWP[peakz],]$position[1]
        y1 <- ladderValuesWP[peakz]
        
        if(peakz == 1){
          x2 <- pxWP[pxWP$size == ladderValuesWP[peakz+1],]$position[1]
          y2 <- ladderValuesWP[peakz+1]
          peakSlope <- (y2-y1)/(x2-x1)
          peakIntercept <- y1-x1*peakSlope
          pxWP[pxWP$size == 0 & pxWP$position < x2,]$size <- (pxWP[pxWP$size == 0 & pxWP$position < x2,]$position*peakSlope)+peakIntercept
        } else if (peakz == length(ladderValuesWP)){
          x2 <- max(pxWP$position)
          y2 <- 0.1
          peakSlope <- (y2-y1)/(x2-x1)
          peakIntercept <- y1-x1*peakSlope
          pxWP[pxWP$size == 0 & pxWP$position > x1,]$size <- (pxWP[pxWP$size == 0 & pxWP$position > x1,]$position*peakSlope)+peakIntercept
        } else {
          x2 <- pxWP[pxWP$size == ladderValuesWP[peakz+1],]$position[1]
          y2 <- ladderValuesWP[peakz+1]
          peakSlope <- (y2-y1)/(x2-x1)
          peakIntercept <- y1-x1*peakSlope
          pxWP[pxWP$size == 0 & pxWP$position > x1 & pxWP$position < x2,]$size <- (pxWP[pxWP$size == 0 & pxWP$position > x1 & pxWP$position < x2,]$position*peakSlope)+peakIntercept
        }
      }
      
      # Now to assign the ladder peaks to the pxGP dataset
      cat("\n")
      print("Calling peaks on the glycoprotein ladder...")
      pxGP <- peakCaller(df = pxGP,
                 value = "Ladder_value",
                 threshold = F,
                 expectedPeaks = T,
                 expectedNumber = length(ladderValuesGP),
                 ladder = ladderValuesGP,
                 attemptsForPeaks = ladderAttempts, 
                 ladderUpOrDown = upperOrLowerLadder)

      # Now we adjust the positions so that the GP values and WP values match based on the ladder peaks
      # First, we'll align the GP position based on the ladder peaks 
      # Since not all ladder marks may have been found, we'll restrict our ladder bands to those that were assigned
      ladderValuesGP <- unique(pxGP[pxGP$size != 0,]$size)
      # Now we assign a placeholder value newPosition
      pxGP$newPosition <- 0.1
      # Then we give the max GP position a size equal to the smallest WP pixel size
      pxGP[pxGP$position == max(pxGP$position),]$size <- min(pxWP$size)
      # And we'll give the minimal GP position the maximal WP pixel size
      pxGP[pxGP$position == min(pxGP$position),]$size <- max(pxWP$size)
      # These will act as our upper and lower bounds or 'anchors' (see below)
      # While we could just use these anchors to bin or expand the GP data to fit WP positions, its better to do this on a stack-by-stack basis where
      # stacks are the pixels between known boundaries. We've got two boundaries, but the GP ladder peaks can provide us with other boundaries
      
      # For each band in the GP ladder we'll assign is a newPosition based on the position of the pixel in the pxWP with the closest size
      for (gpPeakNumber in length(ladderValuesGP)){
        gpPeak <- ladderValuesGP[gpPeakNumber]
        # Then we find the position of the WP pixel that is closest in size to the GP ladder peak and assign the GP ladder peak that position as a newPosition
        pxGP[pxGP$size == gpPeak,][1,]$newPosition <- pxWP[abs(pxWP$size-gpPeak) == min(abs(pxWP$size-gpPeak)),]$position
      }
      # Now that we have some newposition anchors, we'll use math to calculate the newPosition of the interpeak pixels
      # first, we'll make an array of the anchor points
      anchors <- append(c(max(pxWP$size)), ladderValuesGP)
      anchors <- append(anchors, c(min(pxWP$size)))
      
      pxWP$GP_value <- 0
      
      # Now we iterate through the stacks
      for (anchorNum in 2:length(anchors)){
        # First we generate wp and gp lists of the pixels that fit within the stack bounded by our known sizes
        # The WP ipList is easy, since sizes have already been determined for every position
        ipList_wp <- pxWP[pxWP$size <= anchors[anchorNum-1] & pxWP$size >= anchors[anchorNum],]
        # But since only certain positions have assigned sizes to the GP values, we take pixels with positions bound by the known points
        # So first we get the maximal position bound by the bigger size defining the stack
        upperGP <- pxGP[pxGP$size == anchors[anchorNum-1],]$position[1]
        # Then we get the minimal position bound by the smaller size defining the stack
        lowerGP <- pxGP[pxGP$size == anchors[anchorNum],]$position[1]
        # Then we just pull the GP pixels that have original position within the stack
        ipList_gp <- pxGP[pxGP$position >= upperGP & pxGP$position <= lowerGP,]
        # Now, there are three possibilities:
        # 1) The ipList_wp and the ipList_gp are the same length, wherein we can ignore increasing or decreasing pixels and just assign the GP_values directly
        if (nrow(ipList_gp) == nrow(ipList_wp)){
          pxWP[pxWP$position %in% ipList_wp$position,]$GP_value <- ipList_gp$value
          
          
        } else if (nrow(ipList_gp) < nrow(ipList_wp)){
        # 2) The ipList_wp could be larger than the correspond ipList_gp.
          # Here, we'll have to bin the wp pixels into a number of bins defined by the gp pixels
          gpBin <- split(ipList_gp, cut(ipList_gp$position, nrow(ipList_wp)))
          # Then we assign the average GP value of each bin to the corresponding wp pixel
          for (ip_wp in 1:nrow(ipList_wp)){
            pxWP[pxWP$position == ipList_wp[ip_wp,]$position,]$GP_value <- mean(gpBin[[ip_wp]]$GP_value)
          }
          
          # Finally, we'll make a linear equation of the GP_Values to fill in the gaps
          emptyNest <- pxWP[is.na(pxWP$GP_value),]
          for(nest in 1:nrow(emptyNest)){
            nestPosition <- emptyNest[nest,]$position
            upperValue <- pxWP[pxWP$position == nestPosition-1,]$GP_value
            lowerValue <- pxWP[pxWP$position == nestPosition+1,]$GP_value
            pxWP[pxWP$position == emptyNest[nest,]$position,]$GP_value <- mean(upperValue, lowerValue)
          }

        } else if (nrow(ipList_gp) > nrow(ipList_wp)){
        # 3) The ipList_wp is smaller than the ipList_gp.
          # In this final case, we'll bin the gp pixels in to a number of bins defined by the wp pixels in the size stack
          gpBin <- split(ipList_gp, cut(ipList_gp$position, nrow(ipList_wp)))
          # Then we assign the average GP value of each bin to the corresponding wp pixel
          for (ip_wp in 1:nrow(ipList_wp)){
            pxWP[pxWP$position == ipList_wp[ip_wp,]$position,]$GP_value <- mean(gpBin[[ip_wp]]$GP_value)
          }
        }
      }
      # Finally, we can bind the assign everything to px
      px <- pxWP
      
    } else {
      # if the same ladder is used for both GP and WP, then everything is simply bound together
      # Then the GP gel is read, corrected, and bound to the WP ladder
      interim <- read.csv(gp)
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
      # The corrected GP values are saved
      write.csv(interim, "csv/gp_corrected.csv", row.names = F)
      if (nrow(interim) != nrow(pxWP)){
        print("The WP ladder and WP sample file is not the same length as the GP sample file. Script may self-destruct...")
      }
      pxWP$GP_value <- interim$adjValue
      #Assign the variable to px
      px <- pxWP

      #Run the peak-assigning script
      px <- peakCaller(df = px, 
                value = "Ladder_value",
                threshold = F,
                expectedPeaks = T,
                expectedNumber = length(ladderValuesWP),
                ladder = ladderValuesWP,
                attemptsForPeaks = ladderAttempts, 
                ladderUpOrDown = upperOrLowerLadder)
      
      ladderValuesWP <- unique(px[px$size != 0,]$size)
      for (peakz in 1:length(ladderValuesWP)){
        x1 <- px[px$size == ladderValuesWP[peakz],]$position[1]
        y1 <- ladderValuesWP[peakz]
        
       if(peakz == 1){
         x2 <- px[px$size == ladderValuesWP[peakz+1],]$position[1]
         y2 <- ladderValuesWP[peakz+1]
         peakSlope <- (y2-y1)/(x2-x1)
         peakIntercept <- y1-x1*peakSlope
         px[px$size == 0 & px$position < x2,]$size <- (px[px$size == 0 & px$position < x2,]$position*peakSlope)+peakIntercept
       } else if (peakz == length(ladderValuesWP)){
         x2 <- max(px$position)
         y2 <- 0.1
         peakSlope <- (y2-y1)/(x2-x1)
         peakIntercept <- y1-x1*peakSlope
         px[px$size == 0 & px$position > x1,]$size <- (px[px$size == 0 & px$position > x1,]$position*peakSlope)+peakIntercept
       } else {
         x2 <- px[px$size == ladderValuesWP[peakz+1],]$position[1]
         y2 <- ladderValuesWP[peakz+1]
         peakSlope <- (y2-y1)/(x2-x1)
         peakIntercept <- y1-x1*peakSlope
         px[px$size == 0 & px$position > x1 & px$position < x2,]$size <- (px[px$size == 0 & px$position > x1 & px$position < x2,]$position*peakSlope)+peakIntercept
       }
      }
    }
    
    #Now that the pixels are all aligned and such, lets calculate the relative glycosylation score
    # First, lets get rid of negative values
    write.csv(px, "csv/full_values.csv", row.names = F)
    px$GP_value <- px$GP_value+abs(min(px$GP_value))
    px$WP_value <- px$WP_value+abs(min(px$WP_value))
    px$Ladder_value <- px$Ladder_value+abs(min(px$Ladder_value))
    # Now to subtract the GP_value from the WP_values and get rid of negatives
    px$relGylcoScore <- px$GP_value-px$WP_value
    px$relGylcoScore <- px$relGylcoScore +abs(min(px$relGylcoScore))
    #Now we'll lower quartile normalize because... Z scores don't look as good?
    px$relGylcoScore <- px$relGylcoScore/quantile(px$relGylcoScore)[2]
    
    # The dataframe is saved as the output file
    write.csv(px, "csv/full_values.csv", row.names = F)
    
    cat("\n")
    print("CSVs generated, creating graphs...")
    
    # First, a correlation scatter plot is generated and saved
    draft <- ggplot(data = px, aes(
      x = WP_value,
      y = GP_value,
    ))+
      geom_point(size = 4, aes(color = position))+
      geom_density_2d()+
      geom_abline(slope = 1, intercept = quantile(px$GP_value)[2], size = 2, linetype = 2)+
      theme_classic()+
      xlab("Whole protein signal")+
      ylab("Glycosylation signal")
    draft
    ggsave(paste0("figures/signalCorr_position.pdf"), width = 7.5, height = 7.5, units = "in")
    
    draft <- ggplot(data = px, aes(
      x = WP_value,
      y = GP_value,
    ))+
      geom_point(size = 4, aes(color = relGylcoScore))+
      geom_density_2d()+
      geom_abline(slope = 1, intercept = quantile(px$GP_value)[2], size = 2, linetype = 2)+
      theme_classic()+
      scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = median(px$relGylcoScore))+
      xlab("Whole protein signal")+
      ylab("Glycosylation signal")
    draft
    ggsave(paste0("figures/signalCorr_relGlycoScore.pdf"), width = 7.5, height = 7.5, units = "in")
    
    #Now we make the super cool heat map histograms using the newly calculated diagrams
    draft <- ggplot(data = px, aes(
      x = position,
      y = WP_value,
      color = relGylcoScore,
      fill = relGylcoScore
    ))+
      scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = mean(px$relGylcoScore))+
      scale_fill_gradient2(low = "blue", mid = "green", high = "red", midpoint = mean(px$relGylcoScore))+
      geom_bar(stat = "identity")+
      theme_classic2()+
      theme(legend.position = "top")+
      xlab("Relative position")+
      ylab("Whole protein gel intensity")
    draft
    ggsave("figures/relGlyco.pdf", width = 10, height = 4, units = "in")
    draft+geom_line(aes(y = Ladder_value), color = "#5c5c5cff")
    ggsave(paste0("figures/relGlyco_ladder.pdf"), width = 10, height = 4, units = "in")
    
    draft <- ggplot(data = px, aes(
      x = position,
      y = WP_value,
      color = WP_value,
      fill = WP_value
    ))+
      scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = mean(px$WP_value))+
      scale_fill_gradient2(low = "blue", mid = "green", high = "red", midpoint = mean(px$WP_value))+
      geom_bar(stat = "identity")+
      theme_classic2()+
      theme(legend.position = "top")+
      xlab("Relative position")+
      ylab("Whole protein gel intensity")
    draft
    ggsave("figures/WP.pdf", width = 10, height = 4, units = "in")
    draft+geom_line(aes(y = Ladder_value), color = "#5c5c5cff")
    ggsave(paste0("figures/WP_ladder.pdf"), width = 10, height = 4, units = "in")
    
    draft <- ggplot(data = px, aes(
      x = position,
      y = GP_value,
      color = GP_value,
      fill = GP_value
    ))+
      scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = mean(px$GP_value))+
      scale_fill_gradient2(low = "blue", mid = "green", high = "red", midpoint = mean(px$GP_value))+
      geom_bar(stat = "identity")+
      theme_classic2()+
      theme(legend.position = "top")+
      xlab("Relative position")+
      ylab("Glyco protein gel intensity")
    draft
    ggsave("figures/GP.pdf", width = 10, height = 4, units = "in")
    draft+geom_line(aes(y = Ladder_value), color = "#5c5c5cff")
    ggsave(paste0("figures/GP_ladder.pdf"), width = 10, height = 4, units = "in")
    
    cat("\n")
    print("Glycorunner complete.")
  } else{
    print("Ladder not present. Try again; maybe better?")
  }
}

peakCaller <- function(df = px,
                       value = "Ladder_value",
                       threshold = 0,
                       expectedPeaks = F,
                       expectedNumber = 9,
                       ladder,
                       attemptsForPeaks=10,
                       ladderUpOrDown = "lower"){
  peaks <- df
  boomer <- df
  if(threshold != F){
    print(paste0("Threshold given of ", threshold, "."))
    peaks <- peaks[peaks[value] > threshold,]
    
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
      threshold <- mean(unlist(peaks[value]))+sd(unlist(peaks[value]))
      peaks <- peaks[peaks[value] > threshold,]
      
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
          peaks <- peaks[peaks[value] > threshold,]
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
      print(paste0("Peaks found: ", bandNumber))
    } else {
      print("Threshold not given.")
      threshold <- mean(unlist(peaks[value]))+sd(unlist(peaks[value]))
      print(paste0("Using mean+sd threshold of: ",threshold,"."))
      peaks <- peaks[peaks[value] > threshold,]
      
      bandNumber <- 1
      peaks$band <- bandNumber
      for(a in 2:nrow(peaks)){
        if (peaks[a,]$position != peaks[a-1,]$position+1){
          bandNumber <- bandNumber+1
        }
        peaks[a,]$band <- bandNumber
      }
      print(paste0("Peaks found: ", bandNumber))
    }
  }
  boomer$size <- 0
  
  if (value == "Ladder_value"){
    for (ladderBand in unique(peaks$band)){
      bandSize <- ladder[ladderBand]
      midPoint <- median(unlist(peaks[peaks$band == ladderBand,]["position"]))
      closestPoint <- min(abs(peaks[peaks$band == ladderBand,]["position"]-midPoint))
      peakCoord <- peaks[abs(peaks["position"]-midPoint) == closestPoint & peaks$band == ladderBand,]$position
      boomer[boomer$position %in% peakCoord,]$size <- bandSize
    }
  }
  
  # for (peakID in unique(peaks$band)){
  #   interim <- subset(peaks, band == peakID)
  #   if(!exists("bands")){
  #     bands <- data.frame(
  #       "name" = unique(interim$name),
  #       "type" = unique(interim$type),
  #       "bandID" = peakID,
  #       "position" = median(interim$position),
  #       "width" = nrow(interim),
  #       "intensity" = sum(interim$adjValue),
  #       "size" = 0
  #     )
  #   } else {
  #     turkey <- data.frame(
  #       "name" = unique(interim$name),
  #       "type" = unique(interim$type),
  #       "bandID" = peakID,
  #       "position" = median(interim$position),
  #       "width" = nrow(interim),
  #       "intensity" = sum(interim$adjValue),
  #       "size" = 0
  #     )
  #     bands <- rbind(bands, turkey)
  #   }
  # }
  # 
  # if (value == "Ladder_value"){
  #   bands$size <- 0
  #   if (bandNumber != length(ladder)){
  #     if (bandNumber > length(ladder)){
  #       bands <- bands[1:length(ladder),]
  #     }
  #     if (ladderUpOrDown == "lower"){
  #       ladder <- ladder[(1+length(ladder)-nrow(bands)):length(ladder)]
  #     } else {
  #       ladder <- ladder[1:length(bands)]
  #     }
  #   }
  #   for (a in 1:length(ladder)){
  #     bands[bands$bandID == a,]$size <- ladder[a]
  #   }
  # }
  # write.csv(bands, paste0("bands/", dfName, ".csv"), row.names = F)
  return(boomer)
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

dirsList <- list.dirs()
for (deer in dirsList){
  if(deer != "."){
    setwd(deer)
    print(paste0("Running sample: ", deer))
    glycoRun(ladderWP = "ladder_wp.csv", ladderGP = "ladder_gp.csv")
    setwd("../")
  }
}
