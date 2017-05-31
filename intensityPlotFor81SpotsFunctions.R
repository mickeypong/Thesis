# Read the color spot colors for all 81 spots
readCorrectSpotColors = function(correctSpotColorFilename) {
    correctSpotColors <- read.csv("correctSpotColors.csv")
    originalValues <- c("G", "O", "R", "Y")
    firstIndex <- transformVectorToNum(correctSpotColors$color1, originalValues, c(1,2,3,4))
    secondIndex <- transformVectorToNum(correctSpotColors$color2, originalValues, c(5,6,7,8))
    thirdIndex <- transformVectorToNum(correctSpotColors$color3, originalValues, c(9,10,11,12))
    fourthIndex <- transformVectorToNum(correctSpotColors$color4, originalValues, c(13,14,15,16))
    correctSpotColorMatrix <- cbind(firstIndex, secondIndex, thirdIndex, fourthIndex)
    correctSpotColorMatrix
}


# Transform a character vector to a number vector
# When vector == originalValues[i], that value is numValues[i]
transformVectorToNum = function (aVector, originalValues, numValues) {
    newVector <- vector("numeric", length(aVector))
    for(i in 1:length(originalValues)) {
        newVector <- newVector + (aVector == originalValues[i]) * numValues[i]
    }
    newVector
}


# Read color spots from a file
# names of spots include: "index", "xcoord", "ycoord", "avgIntensity", "bgIntensity",  
# "pureIntensity" = avg - bg
readColorSpots = function (spotFilename, refSpots = NULL, header = TRUE) {
    colorSpots <- read.csv(spotFilename, header)
    colorSpots$pureIntensity <- colorSpots$avgIntensity - colorSpots$bgIntensity
    colorSpots
}


# Add the relativeToRefIntensity column as pure - refPure
addRelativeToRefIntensity = function (listOfColorSpots, refColorSpots) {
    for (i in 1 : length(listOfColorSpots)) {
        listOfColorSpots[[i]]$relativeToRefIntensity <- 
            listOfColorSpots[[i]]$pureIntensity / refColorSpots$pureIntensity
    }
    listOfColorSpots
}

# Get normalization factor to calibrate for each color
# There are 16 of time, and the 4 colors at certain time should have 
# a fixed ratio assuming the fluorophore DNA binding efficiency is the same,
# and none of the fluor are bleached.
getNormalizationFactorsForColor = function(listOfColorSpots,
                                           minNormalizationFactor = 0.1,
                                           nameOfIntensitySource = "relativeToRefIntensity",
                                           lowerPercentForAvg = 5, 
                                           upperPercentForAvg = 10) {
    
    # normalziation factor is a matrix of row as time1234 and colume as colorGORY
    numOfcolors <- 4
    numOfTimes <- length(listOfColorSpots) / numOfcolors
    nf <- matrix(0, numOfTimes, numOfcolors)
    
    # 1. Get the orignal normalization factor as the average top 5 to 10 perscent
    # intensity, note R1, O2, Y3, G4 have values less than 0.1, indicating there
    # is no such fluors used in the experiment, so there factors need to be
    # derived from other factors
    for (i in 1 : length(listOfColorSpots)) {
        normalizationFactor <- getAvgTopLowerToUpperPercent(
            listOfColorSpots[[i]][[nameOfIntensitySource]], 
            lowerPercentForAvg, upperPercentForAvg)
        nf[ceiling(i / numOfcolors), (i - 1) %% numOfcolors + 1] <-
            normalizationFactor
    }
    
    # 2. Get the new normalization factor with better values for R1, O2, Y3, G4
    # Use the other fluors at the same time point, and the ratio of R2 to 
    # other fluors at time 2, etc to estimate the nf for R1.
    newNF <- nf
    times <- c(1 : numOfTimes)
    for (time in 1 : nrow(nf)) {
        noFluorColorI <- which(nf[time,] < minNormalizationFactor)
        fluorColorI <- which(nf[time,] > minNormalizationFactor)
        otherTimes <- times[!times %in% time]
        
        # Get ratios of nf of other fluor colors to the no-fluor color
        otherColorRatios <- vector("numeric", numOfcolors) 
        for (anotherColor in fluorColorI) {      
            
            # Get ratio of nf of another color to nf of the no-fluor color based 
            # on the ratio from oter time points
            anotherColorRatios <- NULL
            for (otherTime in otherTimes) {
                if (nf[otherTime, anotherColor] > minNormalizationFactor &
                    nf[otherTime, noFluorColorI] > minNormalizationFactor) {
                    anotherColorRatios <- 
                        c(anotherColorRatios, nf[otherTime, noFluorColorI] / 
                              nf[otherTime, anotherColor])
                }
            }
            if (!is.null(anotherColorRatios)) {
                anotherColorRatio <- mean(anotherColorRatios)
                otherColorRatios[anotherColor] <- anotherColorRatio
            }  
        }
        
        numOfUseableOtherColor <- sum(otherColorRatios != 0)
        newNF[time, noFluorColorI] <- 1 /
            numOfUseableOtherColor * sum(otherColorRatios * nf[time,])
    }
    newNF
}


# Get maximum top lower to upper percent relative intensities for each color in a wash
# usually set to top 5 to 10 percent
getAvgTopLowerToUpperPercent = function(intensities, lower = 5, upper = 10) {
    numOfSpots <- length(intensities)
    sortedPureIntensity <- sort(intensities, decreasing = TRUE)
    topLowerToUpperPercent <- sortedPureIntensity[
        c(ceiling(numOfSpots * lower / 100) : ceiling(numOfSpots * upper / 100))]
    avgTopLowerToUpperPercent <- mean(topLowerToUpperPercent)
    avgTopLowerToUpperPercent
}


# Add the feature of colorCalibratedRelativeToRefIntensity, use the
# colorNormalizationFactors as denominator
# input:listOfColorSpots and
# colorNormalizationFactors is a matrix of time1234 x colorGORY
# output: listOfColorSpots with colorCalibratedRelativeToRefIntensity
addColorCalibratedRelativeToRefIntensity = function (
    listOfColorSpots, colorNormalizationFactors, nameOfIntensitySource = "relativeToRefIntensity") {    
    numOfColors <- ncol(colorNormalizationFactors)
    for (i in 1 : length(listOfColorSpots)) {
        listOfColorSpots[[i]]$colorCalibratedRelativeToRefIntensity <- 
            listOfColorSpots[[i]][[nameOfIntensitySource]] /
            colorNormalizationFactors[ceiling(i / numOfColors), 
                                     (i - 1) %% numOfColors + 1]
    }
    listOfColorSpots
}


# Get a 81 x 16 matrix of intensities
# nrow = 81, which is number of spots
# ncol = 16, which is number of colors
getMatrixOfSpotToColorIntensities = function (listOfColorSpots, nameOfIntensity) {
    allColors <- length(listOfColorSpots)
    numOfSpots <- nrow(listOfColorSpots[[1]])
    spotToColorIntensities <- matrix(0, numOfSpots, allColors)
    for (i in 1 : allColors) {
        colorSpots <- listOfColorSpots[[i]]
        spotToColorIntensities[,i] <- colorSpots[[nameOfIntensity]]
    }
    spotToColorIntensities
}


# Add the feature of spotNormalizedIntensities
# Normalize all 81 spots, so the 16 colors in a spot is normalized to 1.
# input: listOfColorSpots
# output: listOfColorSpots with spotNormalizedIntensity
addSpotNormalizedIntensities = function(listOfColorSpots) {    
    spotToColorIntensities <- getMatrixOfSpotToColorIntensities(
        listOfColorSpots, "colorCalibratedRelativeToRefIntensity")
    
    # get normalization factor for each spot
    squaredSpotToColorIntensities <- spotToColorIntensities * spotToColorIntensities
    sumEachRowOfSquaredSpotToColorIntensities <- rowSums(squaredSpotToColorIntensities)
    spotNormalizationFactors <- sqrt(sumEachRowOfSquaredSpotToColorIntensities)
    
    # normalize each spot using normalization factor
    for (i in 1 : length(listOfColorSpots)) {
        listOfColorSpots[[i]]$spotNormalizedIntensity <- 
            listOfColorSpots[[i]]$colorCalibratedRelativeToRefIntensity / 
            spotNormalizationFactors
    }
    listOfColorSpots
}


# plot the 81 spots with givine intensities from G, Y, O, R color spots
plotGORY2 = function(listOfColorSpots, nameOfIntensity, washPoint = 1, plotnameSuffix = "") {
    par(mfrow = c(9,1), 
        oma = c(3,4,2,0) + 0.1,
        mar = c(0,0,1,1) + 0.5);
    
    barPlotData <- rbind(listOfColorSpots[[1 + (washPoint-1)*4]][[nameOfIntensity]], 
                         listOfColorSpots[[2 + (washPoint-1)*4]][[nameOfIntensity]],  
                         listOfColorSpots[[3 + (washPoint-1)*4]][[nameOfIntensity]], 
                         listOfColorSpots[[4 + (washPoint-1)*4]][[nameOfIntensity]]);
    
    minYLimit <- min(barPlotData)
    maxYLimit <- max(barPlotData)
    firstPortion <- 1:9;
    barplot(barPlotData[,firstPortion], col = c("green", "darkorange", "red", "yellow"), 
            beside = TRUE, names.arg = firstPortion, ylim = c(minYLimit, maxYLimit),
            main = paste(nameOfIntensity, "of 81 spots at wash", washPoint, plotnameSuffix));
    drawSubBarplot(barPlotData, maxYLimit, 10:18);
    drawSubBarplot(barPlotData, maxYLimit, 19:27);
    drawSubBarplot(barPlotData, maxYLimit, 28:36);
    drawSubBarplot(barPlotData, maxYLimit, 37:45);
    drawSubBarplot(barPlotData, maxYLimit, 46:54);
    drawSubBarplot(barPlotData, maxYLimit, 55:63);
    drawSubBarplot(barPlotData, maxYLimit, 64:72);
    drawSubBarplot(barPlotData, maxYLimit, 73:81);
}

# Plot all 16 colors of 81 spots in a plot
plotGORYForAll = function(listOfColorSpots, nameOfIntensity) {
    par(mfrow = c(9,1), 
        oma = c(3,4,2,0) + 0.1,
        mar = c(0,0,1,1) + 0.5);
    
    barPlotData <- NULL
    for (i in 1 : length(listOfColorSpots)) {
        barPlotData <- rbind(barPlotData, listOfColorSpots[[i]][[nameOfIntensity]]);
    }
  
    minYLimit <- min(barPlotData)
    maxYLimit <- max(barPlotData)
    firstPortion <- 1:9;
    barplot(barPlotData[,firstPortion], 
            col = c("green", "darkorange", "red", "yellow", "green", "darkorange", "red", "yellow",
                    "green", "darkorange", "red", "yellow", "green", "darkorange", "red", "yellow"), 
            beside = TRUE, names.arg = firstPortion, ylim = c(minYLimit, maxYLimit),
            main = paste("Relative color intensity of 81 spots"));
    drawSubBarplot(barPlotData, maxYLimit, 10:18);
    drawSubBarplot(barPlotData, maxYLimit, 19:27);
    drawSubBarplot(barPlotData, maxYLimit, 28:36);
    drawSubBarplot(barPlotData, maxYLimit, 37:45);
    drawSubBarplot(barPlotData, maxYLimit, 46:54);
    drawSubBarplot(barPlotData, maxYLimit, 55:63);
    drawSubBarplot(barPlotData, maxYLimit, 64:72);
    drawSubBarplot(barPlotData, maxYLimit, 73:81);
}


# plot hte 80-20 mix of spots
# Input: list of colorSpots at time 1, list of colorSpots at time 2
plotGORY_mixer = function(time1ColorSpots, time2ColorSpots, plotnameSuffix, percentage = 0.8) {
    listOfNewColorSpots <- list()
    for (i in 1:length(time1ColorSpots)) {
        listOfNewColorSpots[[i]] <- data.frame(
            avgIntensity = time1ColorSpots[[i]]$avgIntensity * percentage + time2ColorSpots[[i]]$avgIntensity * (1-percentage),
            bgIntensity = time1ColorSpots[[i]]$bgIntensity * percentage + time2ColorSpots[[i]]$bgIntensity * (1-percentage),
            pureIntensity = time1ColorSpots[[i]]$pureIntensity * percentage + time2ColorSpots[[i]]$pureIntensity * (1-percentage))
    }
    
    listOfNewColorSpots <- addRelativeToRefIntensity(listOfNewColorSpots, G25)
    
    maxRelativeIntensities1 <- getMaxRelativeIntensities(time1ColorSpots[[1]],
                                                         time1ColorSpots[[2]],
                                                         time1ColorSpots[[3]],
                                                         time1ColorSpots[[4]])
    maxRelativeIntensities2 <- getMaxRelativeIntensities(time2ColorSpots[[1]],
                                                         time2ColorSpots[[2]],
                                                         time2ColorSpots[[3]],
                                                         time2ColorSpots[[4]])
    print(maxRelativeIntensities1)
    print(maxRelativeIntensities2)
    maxRelativeIntensities <- apply(data.frame(maxRelativeIntensities1, maxRelativeIntensities2), 1, max)
    colorNormalizationFactors <- matrix(maxRelativeIntensities, 1, 4)
    listOfNewColorSpots <- addColorCalibratedRelativeToRefIntensity(listOfNewColorSpots, 
                                                                    colorNormalizationFactors)
    
    plotGORY2(listOfNewColorSpots, "colorCalibratedRelativeToRefIntensity", washPoint = 1, 
              paste(",", percentage, "-", (1-percentage), "mix of", plotnameSuffix))
}

# get spots that have extra colors compared to the correct colors
# return a data frame of columns as 1) spotIndex, 2) numOfExtraColor
getSpotWithExtraColor = function(spotToColorIntensities, correctSpotColor, 
                                 intensityThreshold = 0.15) {
    spotWithExtraColor <- data.frame(spotIndex = numeric(), 
                                     numOfExtraColor = numeric(),
                                     extraSpotIndex = character())
    for(i in 1: nrow(spotToColorIntensities)) {
        choosenSpotColor <- which(spotToColorIntensities[i,] > intensityThreshold)
        correctSpotColor <- correctSpotColors[i,]
        extraColors <- setdiff(choosenSpotColor, correctSpotColor)
        if (length(extraColors) != 0) {
            spotWithExtraColor <- rbind(spotWithExtraColor, 
                                        data.frame(spotIndex = i, numOfExtraColor = length(extraColors)))
        }
    }
    spotWithExtraColor
}


# get spots that have less colors compared to the correct colors
# return a data frame of columns as 1) spotIndex, 2) numOfDeficitColor
getSpotWithDeficitColor = function(spotToColorIntensities, correctSpotColor, 
                                intensityThreshold = 0.15) {
    spotWithDeficitColor <- data.frame(spotIndex = numeric(), numOfDeficitColor = numeric())
    for(i in 1: nrow(spotToColorIntensities)) {
        choosenSpotColor <- which(spotToColorIntensities[i,] > intensityThreshold)
        correctSpotColor <- correctSpotColors[i,]
        insufficientColors <- setdiff(correctSpotColor, choosenSpotColor)
        if (length(insufficientColors) != 0) {
            spotWithDeficitColor <- rbind(spotWithDeficitColor, 
                                          data.frame(spotIndex = i, 
                                                     numOfDeficitColor = length(insufficientColors)))
        }
    }
    spotWithDeficitColor
}


# get the standard deviation devided by mean of the true positive intensities in a color
# totally 16 colors
getSigmasForPositiveColors = function(listOfColorSpots, nameOfIntensity, correctSpotColors) {
    sigmas <- NULL
    for(i in 1:length(listOfColorSpots)) {
        colorSpots <- listOfColorSpots[[i]]
        numOfColors <- 4
        positiveIndex <- which(correctSpotColors[,ceiling(i / numOfColors)] == i)
        positiveIntensities <- colorSpots[[nameOfIntensity]][positiveIndex]
        sigmas <- c(sigmas, sd(positiveIntensities) / mean(positiveIntensities))
    }
    sigmas
}






# Debug
# Get maximum percent relative intensities for each color in a wash
getMaxRelativeIntensities = function(gData, oData, rData, yData) {
    c(max(gData$relativeToRefIntensity), max(oData$relativeToRefIntensity), 
      max(rData$relativeToRefIntensity), max(yData$relativeToRefIntensity))
    #c(getAvgTopLowerToUpperPercent(gData$relativeToRefIntensity), 
    #  getAvgTopLowerToUpperPercent(oData$relativeToRefIntensity), 
    #  getAvgTopLowerToUpperPercent(rData$relativeToRefIntensity), 
    #  getAvgTopLowerToUpperPercent(yData$relativeToRefIntensity))
}


plotGORY = function(gData, oData, rData, yData, pioneerData, isToNoramlWithSameColor = FALSE, 
                    limitToNormal = 0.1, maxRelativeIntensities = NA, washPoint = 1) {
    
#    gData$relativeToRefIntensity <- gData$pureIntensity / pioneerData$pureIntensity
#    oData$relativeToRefIntensity <- oData$pureIntensity / pioneerData$pureIntensity
#    rData$relativeToRefIntensity <- rData$pureIntensity / pioneerData$pureIntensity
#    yData$relativeToRefIntensity <- yData$pureIntensity / pioneerData$pureIntensity
    
    if (is.na(maxRelativeIntensities)) {
        maxRelativeIntensities <- getMaxRelativeIntensities(gData, oData, rData, yData)
    }
    
    maxRelativeG <- ifelse(maxRelativeIntensities[1] > limitToNormal, maxRelativeIntensities[1], 1);
    maxRelativeO <- ifelse(maxRelativeIntensities[2] > limitToNormal, maxRelativeIntensities[2], 1);
    maxRelativeR <- ifelse(maxRelativeIntensities[3] > limitToNormal, maxRelativeIntensities[3], 1);
    maxRelativeY <- ifelse(maxRelativeIntensities[4] > limitToNormal, maxRelativeIntensities[4], 1);
    
    maxYLimit <- max(maxRelativeG, maxRelativeO, maxRelativeR, maxRelativeY);
    
    if (isToNoramlWithSameColor) {
        maxYLimit <- 1;
        print(maxRelativeG);
        print(maxRelativeO);
        print(maxRelativeR);
        print(maxRelativeY);
        
        gData$relativeToRefIntensity <- gData$relativeToRefIntensity / maxRelativeG;
        oData$relativeToRefIntensity <- oData$relativeToRefIntensity / maxRelativeO;
        rData$relativeToRefIntensity <- rData$relativeToRefIntensity / maxRelativeR;
        yData$relativeToRefIntensity <- yData$relativeToRefIntensity / maxRelativeY;
    }
    
    #print(gData$relativeIntensity)
    #print(oData$relativeIntensity)
    #print(rData$relativeIntensity)
    #print(yData$relativeIntensity)
    
    par(mfrow = c(9,1), 
        oma = c(3,4,2,0) + 0.1,
        mar = c(0,0,1,1) + 0.5);
    
    barPlotData <- rbind(gData$relativeToRefIntensity, oData$relativeToRefIntensity, 
                         rData$relativeToRefIntensity, yData$relativeToRefIntensity);
    
    firstPortion <- 1:9;
    barplot(barPlotData[,firstPortion], col = c("green", "darkorange", "red", "yellow"), 
            beside = TRUE, names.arg = firstPortion, ylim = c(-0.2, maxYLimit),
            main = paste("Relative color intensity of 81 spots at wash", washPoint));
    drawSubBarplot(barPlotData, maxYLimit, 10:18);
    drawSubBarplot(barPlotData, maxYLimit, 19:27);
    drawSubBarplot(barPlotData, maxYLimit, 28:36);
    drawSubBarplot(barPlotData, maxYLimit, 37:45);
    drawSubBarplot(barPlotData, maxYLimit, 46:54);
    drawSubBarplot(barPlotData, maxYLimit, 55:63);
    drawSubBarplot(barPlotData, maxYLimit, 64:72);
    drawSubBarplot(barPlotData, maxYLimit, 73:81);
}

drawSubBarplot = function(barPlotData, maxYLimit, portion) {
    barplot(barPlotData[,portion], col = c("green", "darkorange", "red", "yellow"), 
            beside = TRUE, names.arg = portion, ylim = c(-0.2, maxYLimit));
}



