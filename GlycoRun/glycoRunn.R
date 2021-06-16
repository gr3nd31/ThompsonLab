library(ggpubr)

glycoRun <- function(ladderValues = c(250, 150, 100, 70, 50, 37, 25, 20, 15),
                     outputCSV = "allData.csv",
                     corrections = T,
                     peakCSV = "csv/allData.csv",
                     givenThreshold = F) {
  currentSetup <- list.files()
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
  
  correctionsCheck <- list.files()
  if (corrections == T & "ladder.csv" %in% correctionsCheck){
    cList <- list.files(pattern = ".csv")
    px <- read.csv("ladder.csv")
    #file.rename("ladder.csv", "csv/ladder.csv")
    px$name <- "Ladder"
    px$type <- "Ladder"
    names(px)[1] <- "position"
    names(px)[2] <- "value"
    
    sampleLine <- lm(value~position, data = px)
    intercept <- sampleLine[1]$coefficients[1]
    slope <- sampleLine[1]$coefficients[2]
    px$adjValue <- px$value-(slope*px$position+intercept)
    
    while (abs(slope) > 0.001){
      #dSample <- px[sample(nrow(px), nrow(px)/10),]
      sampleLine <- lm(adjValue~position, data = px)
      intercept <- sampleLine[1]$coefficients[1]
      slope <- sampleLine[1]$coefficients[2]
      px$adjValue <- px$adjValue-(slope*px$position+intercept)
    }
    
    for (a in cList){
      if(a != "ladder.csv" & a != outputCSV){
        mName <- strsplit(a, "_")[[1]][1]
        mType <- strsplit(strsplit(a, "_")[[1]][2], ".csv")[[1]][1]
        
        interim <- read.csv(a)
        interim$name <- mName
        interim$type <- mType
        names(interim)[1] <- "position"
        names(interim)[2] <- "value"
        
        #file.rename(a, paste0("csv/", a))
        
        #dSample <- px[sample(nrow(px), nrow(px)/10),]
        sampleLine <- lm(value~position, data = interim)
        intercept <- sampleLine[1]$coefficients[1]
        slope <- sampleLine[1]$coefficients[2]
        interim$adjValue <- interim$value-(slope*interim$position+intercept)
        
        while (abs(slope) > 0.001){
          #dSample <- px[sample(nrow(px), nrow(px)/10),]
          sampleLine <- lm(adjValue~position, data = interim)
          intercept <- sampleLine[1]$coefficients[1]
          slope <- sampleLine[1]$coefficients[2]
          interim$adjValue <- interim$adjValue-(slope*interim$position+intercept)
        }
        px <- rbind(px, interim)
      }
    }
    write.csv(px, paste0("csv/", outputCSV), row.names = F)
    px$tag <- paste0(px$name, "_", px$type)
    
    for(a in unique(px$name)){
      if (a != "Ladder"){
        interim <- subset(px, name == a | name == "Ladder")
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
        ggsave(paste0("figures/",a, ".pdf"))
        
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
        
        interim$coom <- interim$coom+abs(min(interim$coom))
        interim$glyco <- interim$glyco+abs(min(interim$glyco))
        interim$relGylcoSignal <- interim$glyco-interim$coom
       # interim$relGylcoSignal <- log((interim$relGylcoSignal+abs(min(interim$relGylcoSignal))+0.1),2)
        
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
    px <- read.csv(peakCSV)
    px$tag <- paste0(px$name, "_", px$type)
  }

  for (a in unique(px$tag)){
    tagSet <- subset(px, tag == a)
    print(paste0("Calling peaks on ", a))
    if (a == "Ladder_Ladder"){
      peakCaller(df = tagSet, dfName = "ladder", 
                 threshold = F,
                 expectedPeaks = T,
                 expectedNumber = length(ladderValues),
                 ladder = ladderValues)
      
      # Uses the bands on the ladder to calculate kD based on position
      ladderBands <- read.csv("bands/ladder.csv")
      ladderEq <- lm(size~position, data = ladderBands)
      intercept <- ladderEq$coefficients[1]
      slope <- ladderEq$coefficients[2]
      
      # Draws a QC plot of the regession model
      draft <- ggplot(data = ladderBands, aes(
        x = position,
        y = size
      ))+geom_point(size = 4)+
        geom_smooth(method = "lm",)+
        theme_classic()
      print(draft)
      ggsave("figures/ladderQC.pdf")
      
      # Saves the size data of the whole dataset
      px$size <- px$position*slope+intercept
      write.csv(px, peakCSV, row.names = F)
      
    } else {
      peakCaller(df = tagSet, 
                 dfName = a, 
                 threshold = givenThreshold,
                 expectedPeaks = F,
                 expectedNumber = 10,
                 ladder = ladderValues)
      peakBands <- read.csv(paste0("bands/", a, ".csv"))
      peakBands$size <- peakBands$position*slope+intercept
      write.csv(peakBands, paste0("bands/", a, ".csv"), row.names = F)
    }
  }
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
                       ladder){
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
        print("Incorrect peaks found. Exploring other thresholds...")
        while (max(peaks$band != expectedNumber)){
          if (max(peaks$band) > expectedNumber){
            threshold <- threshold+0.5
          } else {
            threshold <- threshold-0.5
          }
          peaks <- holder
          peaks <- subset(peaks, adjValue > threshold)
          
          bandNumber <- 1
          for(a in 2:nrow(peaks)){
            if (peaks[a,]$position != peaks[a-1,]$position+1){
              bandNumber <- bandNumber+1
            }
            peaks[a,]$band <- bandNumber
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
