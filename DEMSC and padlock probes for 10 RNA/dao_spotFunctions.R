source("/Users/jackie/Downloads/SGTC microscope/2016-08-11-10mRNAs vs 4 exon mRNAs-37C-45C Amp 3rd time/4 exon mRNAs/spotPlotFunctions.R")
source("/Users/jackie/Downloads/SGTC microscope/2016-08-11-10mRNAs vs 4 exon mRNAs-37C-45C Amp 3rd time/4 exon mRNAs/spotRegisterFunctions.R")
source("/Users/jackie/Downloads/SGTC microscope/2016-08-11-10mRNAs vs 4 exon mRNAs-37C-45C Amp 3rd time/4 exon mRNAs/spotStackFunctions.R")

# Read spots from a daostorm spotFile, and use tiffFile to add intensity info
# Note that the spots Y axis is from top to down, and x axis is from left to right
# The corresponding pixel of Y is y_um / pix_XY + 1?
# The corresponding pixel of X is x_um / pix_XY + 1?
# Actually the XCENTER and YCENTER from daostorm spots match well with the tiff
# Too bad that I have to modify it for match the Pos_X and Pos_Y of FQ spots that has a unit of micron
read_daostorm_spots = function (spotsFilename, color, tifFilename) {
    # Load spots of each graph with X-cord, Y-cord, Amp, and Color
    if (file.exists(spotsFilename)) {
        all_content = readLines(spotsFilename);
        spots <- read.csv(textConnection(all_content), sep = "\t");
        spots$Color = color
        spots$isUsed <- FALSE
        names(spots)[names(spots) == "BRIGHTNESS"] <- "AMP"
    } else {
        spots <-  data.frame(YCENTER = numeric(), 
                             XCENTER = numeric(),
                             AMP = numeric(),
                             Color = character(),
                             isUsed <- logical())
    }
    
    spots <- spots %>% mutate(Pos_Y = (YCENTER - 1) * pix_XY)
    spots <- spots %>% mutate(Pos_X = (XCENTER - 1) * pix_XY)
    
    # use the tiff plot to get exact pixel brightness
    tiff = readTIFF(tifFilename, as.is = TRUE)
    
    # Daostorm's algorithm is so good that it can even find spots off the tiff image!
    # Remove the spots that has the index to be outside of index of tiff matrix
    spots <- spots %>% filter((XCENTER >= 1) & (XCENTER <= size(tiff)[2]) & (YCENTER >= 1) & (YCENTER <= size(tiff)[1]))
    
    matrixIndex <- cbind(spots$YCENTER, spots$XCENTER)
    spots$tiffAmp <- tiff[matrixIndex] 
    spots$bgTiffAmp <- getMinNearbyAmp(tiff, matrixIndex, 5)
    spots$minusBgTiffAmp <- spots$tiffAmp - spots$bgTiffAmp
    maxMinusBgTiffAmp <- max(spots$minusBgTiffAmp)
    spots$normalTiffAmp <-  spots$minusBgTiffAmp / maxMinusBgTiffAmp
    
    spots
}


# Given a vector of values, give the value that has is at the lowest X Percentage
# e.g. when X = 40 and values = c(1, 3 ,2, 5, 4) the result is 2
dao_lowestXPercentValue = function (values, lowestXPercent) {
    numOfValues = length(values)
    if (numOfValues != 0) {
        lowestXPercentIndex <- order(values)[ceil(numOfValues * lowestXPercent/100)]
        lowestXPercentValue <- values[lowestXPercentIndex]
    } else {
        lowestXPercentValue <- 0
    }
    lowestXPercentValue
}


# Translate spots in x, y direction to fit to refSpots and add the new positions as Pos_X and Pos_Y while
# the old Pos_X and Pos_Y as Old_Pos_X and Old_Pos_Y, add the Cell_Index
# Translate spots in a cell at a time. 
# maxPreCospotDist is the max distance allowed between a initial spot and a refSpot, 
# if not within the distance, it's likely that the spot is dirt
# minSpotNumPerCell is the min number of spots in a cell, if not over the number,
# it's likely there are very few spots to have a good translation, so we might as well not translate.
# and also return a data.frame that has row number as polygon number and deltaX and deltaY column.
dao_translateAllCells = function (spots, refSpots, polygons_um, maxPreCospotDist, 
                                  minusBgTiffAmpThreshold, minSpotNumPerCell = 3) {
    
    spotsToTranslate <- spots %>% mutate(Old_Pos_X = Pos_X, Old_Pos_Y = Pos_Y, Cell_Index = 0)
    for (i in 1 : length(polygons_um)) {
        print(i)
        polygon_um <- polygons_um[[i]]
        encSpotsIndice <- enclosedSpotsIndice(spots, polygon_um)
        encSpots <- enclosedSpots(spots, polygon_um)
        encRefSpots <- enclosedSpots(refSpots, polygon_um)
        if (length(encSpotsIndice) >= minSpotNumPerCell) {
            translatedEncSpots <- dao_translateACell(polygon_um, encSpots, encRefSpots, 
                                        maxPreCospotDist, minusBgTiffAmpThreshold)
            spotsToTranslate[encSpotsIndice,]$Pos_X <- translatedEncSpots$Pos_X
            spotsToTranslate[encSpotsIndice,]$Pos_Y <- translatedEncSpots$Pos_Y
            spotsToTranslate[encSpotsIndice,]$Cell_Index <- (numeric(length(encSpotsIndice)) + 1)*i     
        } 
    }
    spotsToTranslate
}


# Translate spots in a cell to refSpots.
# maxPreCospotDist is the max distance allowed between a initial good spot and a refSpot, 
# if not within the distance, it's likely that the spot is dirt
# minusBgTiffAmpThreshold is the minimum minusBgTiffAmp of the encSpots, which are bright enough 
# to be used as ancher for translations, so we don't accidentally
# include overwhelming dim fake spots, which often appear near the side of nuclei.
dao_translateACell = function (encSpots, encRefSpots, maxPreCospotDist, minusBgTiffAmpThreshold) {
    brightEncSpots <- encSpots %>% filter(minusBgTiffAmp > minusBgTiffAmpThreshold)
    
    # Get the coSpots on the refSpots
    closestRefSpots <- closestSpot(brightEncSpots, encRefSpots)
    
    # For ankering, choose the top 90 birghtEncSpots that is closest to their coSpots as candidate
    coSpotDistances <- distance(brightEncSpots, closestRefSpots)
    shortestDistanceOrder <- order(coSpotDistances)
    numOfEncSpots <- nrow(brightEncSpots)
    numOfTopShortestSpots <- 90  # used to be 10
    topEncSpots <- brightEncSpots[shortestDistanceOrder[1:min(numOfEncSpots, numOfTopShortestSpots)],]
    
    #brightEncSpotsInOrder <- brightEncSpots[shortestDistanceOrder,]
    #closestRefSpotsInOrder<- closestRefSpots[shortestDistanceOrder,]
    
    # For loop to get the best ankering spot
    currentShorestAvgDist <- avgDistance(brightEncSpots, closestRefSpots)
    currentMovedSpots <- encSpots
    for (i in seq_along(topEncSpots[,1])) {
        topEncSpot <- topEncSpots[i,]
        regEncSpots <- dao_regSpotsUsingASpot(topEncSpot, brightEncSpots, encRefSpots, maxPreCospotDist)
        shortestAvgDistBright <- avgDistanceToClosestRefSpot(regEncSpots, encRefSpots)
        
        if (shortestAvgDistBright < currentShorestAvgDist) {
            currentShorestAvgDist <- shortestAvgDistBright
            regSpotsAndAvgDist <- regEncSpots
            currentMovedSpots <- dao_regSpotsUsingASpot(topEncSpot, encSpots, encRefSpots, maxPreCospotDist)
        }
    }
    currentMovedSpots
}


# get the best coSpot of a spot within maxPreCospotDist, so that 
# after translating all the spots, the average distance between all translatedSpots
# and their closest refSpots is shortest.
# Return the translated spots
dao_regSpotsUsingASpot = function (spot, spots, refSpots, maxPreCospotDist) {
    # Get the candidate refCoSpots to spot
    neighborRefSpots <- neighborhood(spot, refSpots, maxPreCospotDist)
    
    if (nrow(neighborRefSpots) == 0) {
        movedSpots <- spots
    } else {
        averageDistance <- numeric(nrow(neighborRefSpots))
        for (i in seq_along(nrow(neighborRefSpots))) {
            neighborRefSpot <- neighborRefSpots[i,]
            movedSpots <- moveXYFromSpot1ToSpot2(spots, spot, neighborRefSpot)
            averageDistance[i] <- avgDistanceToClosestRefSpot(movedSpots, refSpots)
        }
        bestNeighborRefSpot <- neighborRefSpots[which.min(averageDistance),]
        movedSpots <- moveXYFromSpot1ToSpot2(spots, spot, bestNeighborRefSpot)
    } 
    movedSpots
}


dao_translateAllNucleis = function(spots, refSpots, nucleiPolygons_um) {
    
}



# Given refSpots and the time point, use transF1CospotIndice (if time is 1) to get 
# the closest transF1Spots for each refspot and choose the truely closest one.
# return dataFrame of transF1Spots that might truely overlap with refSpots.
dao_getClosestTransCospots = function(refSpots, transMergedSpots, mergedName) {
    
    coSpotIndiceList <- refSpots[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]]
    transMergedSpots <- transMergedSpots[[mergedName]]
    
    closestTransCospots <- spotDataFrame(transMergedSpots)
    for (r in seq_along(refSpots[,1])) {
        refSpot <- refSpots[r, ]
        coSpotIndice <- coSpotIndiceList[[r]]
        if (length(coSpotIndice) >= 1) {
            closestTransCospot <- closestSpot(refSpot, transMergedSpots[coSpotIndice,])
            closestTransCospots <- rbind(closestTransCospots, closestTransCospot)
        }
    }
    closestTransCospots
}


# Given the interested pixX and pixY location, plots the tiffMatrix within plot_area_dist
# and the spots
dao_plotNearbyTiffRaster = function (pixX_int, pixY_int, spots, normTiffMatrix, pix_XY, 
                                     plot_area_dist = 200, colorPalette = colorRampPalette(c("black", "green")) (1000),
                                     plotname = "", neighborSpotSize = 2) {
    resX <- dim(normTiffMatrix)[2]
    resY <- dim(normTiffMatrix)[1]
    pixPlotArea <- round(plot_area_dist / pix_XY)
    minPixX <- max(1, pixX_int - pixPlotArea)
    maxPixX <- min(pixX_int + pixPlotArea, resX)
    minPixY <- max(1, pixY_int - pixPlotArea)
    maxPixY <- min(pixY_int + pixPlotArea, resY)
    
    nearbyMatrix <- normTiffMatrix[minPixY : maxPixY, minPixX : maxPixX]
    rasterImage <- raster(nearbyMatrix, ymn = minPixY-0.5, ymx = maxPixY+0.5, xmn = minPixX-0.5, xmx = maxPixX+0.5) 
    plot(rasterImage, col = colorPalette, xlim=c(minPixX-0.5, maxPixX+0.5),ylim=c(minPixY-0.5, maxPixY+0.5), main = plotname)
    
    # plot the nearby spots on the tiff plot
    fakeSpot <- spots[1,]
    fakeSpot$Pos_X <- (pixX_int - 1) * pix_XY
    fakeSpot$Pos_Y <- (pixY_int - 1) * pix_XY
    
    neighborSpots <- neighbors_exclude_spot(fakeSpot, spots, plot_area_dist)
    neightborIndice2D <- spotsToIndice2D_not_integer(neighborSpots, pix_XY)
    neighborPixXOnNewMatrix <- getXIndiceFromIndice2D(neightborIndice2D)
    neighborPixYOnNewMatrix <- maxPixY - (getYIndiceFromIndice2D(neightborIndice2D) - minPixY)
    points(neighborPixXOnNewMatrix, neighborPixYOnNewMatrix, col = "white", cex = neighborSpotSize)
}

# Update the IsUsed colume and the CoRefSpotIndex for transMergedSpots that has a refCoSpot
dao_updateIsUsedForTransMergedSpots = function (transMergedSpots, refSpots) {
    for (mergedName in mergedNames) {
        transMergedSpots[[mergedName]]$IsUsed <- logical(nrow(transMergedSpots[[mergedName]]))
        #oneOnOneRefIndice <- which(refSpots[paste("transMergedNumCoSpots", mergedName, sep = "")] == 1)
        #oneOnOneFNIndice <- unlist(refSpots[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]][oneOnOneRefIndice])
        oneOnMultipleFNIndice <- unlist(refSpots[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]])
        transMergedSpots[[mergedName]]$IsUsed[oneOnMultipleFNIndice] <- TRUE
        #transMergedSpots[[mergedName]] <- transMergedSpots[[mergedName]] %>% mutate(
        #    CoRefSpotIndex = numeric(nrow(transMergedSpots[[mergedName]])))
        #transMergedSpots[[mergedName]]$CoRefSpotIndex[oneOnOneFNIndice] <- oneOnOneRefIndice
        print(sum(transMergedSpots[[mergedName]]$IsUsed))
    }
    transMergedSpots
}


# Update the IsUsed colume and the CoRefSpotIndex for transMergedSpots that has one-on-one refCoSpot
dao_updateIsUsedForOneOnOneTransMergedSpots = function (transMergedSpots, refSpots) {
    for (mergedName in mergedNames) {
        transMergedSpots[[mergedName]]$IsUsed <- logical(nrow(transMergedSpots[[mergedName]]))
        oneOnOneRefIndice <- which(refSpots[paste("transMergedNumCoSpots", mergedName, sep = "")] == 1)
        oneOnOneFNIndice <- unlist(refSpots[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]][oneOnOneRefIndice])
        transMergedSpots[[mergedName]]$IsUsed[oneOnOneRefIndice] <- TRUE
        #transMergedSpots[[mergedName]] <- transMergedSpots[[mergedName]] %>% mutate(
        #    CoRefSpotIndex = numeric(nrow(transMergedSpots[[mergedName]])))
        #transMergedSpots[[mergedName]]$CoRefSpotIndex[oneOnOneFNIndice] <- oneOnOneRefIndice
        print(sum(transMergedSpots[[mergedName]]$IsUsed))
    }
    transMergedSpots
}


# Debug
# Fill in the empty coSpot if there is a close unused coSpot within 1 um distance.
# Fill in by changing the column of transFNNumCoSpots and transFNCoSpotsIndice
dao_fillInEmptyCospotsWithUnusedSpots = function (refSpots, secondaryCospotDist = 1) {
    for (mergedName in mergedNames) {
        indexOfZeroSpots <- which(refSpots[paste("transMergedNumCoSpots", mergedName, sep = "")] == 0)
        interestedSpots <- refSpots[indexOfZeroSpots,]
        closestSpotsInFN <- closestSpot(interestedSpots, transMergedSpots[[mergedName]], 
                                        countOnlyUnusedSpots = TRUE)
        closestSpotIndexInFN <- closestSpotIndex(interestedSpots, transMergedSpots[[mergedName]], 
                                                 countOnlyUnusedSpots = TRUE)        
        # Test if the closest spot to closestSpotsInFN is interestedSpot, if not, disgard
        closestSpotInRefSpots <- closestSpot(closestSpotsInFN, refSpots)
        
        # Test if the distance between the interestedSpot and closestSpotsInFN is < secondaryCospotDist, if not, disgar
        iDistanceSmall <- which(distance(interestedSpots, closestSpotsInFN) < secondaryCospotDist &
                                    closestSpotInRefSpots$Pos_X == interestedSpots$Pos_X &
                                    closestSpotInRefSpots$Pos_Y == interestedSpots$Pos_Y) 
        
        # Update the columns of transMergedNumCoSpots and transMergedCoSpotsIndice of refSpots
        interestedSpots[iDistanceSmall, paste("transMergedNumCoSpots", mergedName, sep = "")] <- 1
        interestedSpots[iDistanceSmall, paste("transMergedCoSpotsIndice", mergedName, sep = "")] <- 
            closestSpotIndexInFN[iDistanceSmall]
        refSpots[indexOfZeroSpots, ] <- interestedSpots
    }
    refSpots
    
#     # debug to remove
#     iColumn_transFNNumCoSpots <- which(grepl("NumCoSpots", names(spots)))
#     iColumn_transFNCoSpotsIndice <- which(grepl("CoSpotsIndice", names(spots)))
#     transFNSpots <- list(transF1Spots, transF2Spots, transF3Spots, transF4Spots)
#     for (i in 1:length(iColumn_transFNNumCoSpots)) {
#         indexOfZeroSpots <- which(spots[,iColumn_transFNNumCoSpots[i]] == 0)
#         interestedSpots <- spots[indexOfZeroSpots,]
#         closestSpotsInFN <- closestSpot(interestedSpots, transFNSpots[[i]], countOnlyUnusedSpots = TRUE)
#         closestSpotIndexInFN <- closestSpotIndex(interestedSpots, transFNSpots[[i]], countOnlyUnusedSpots = TRUE)        
#         
#         # test if the closest spot to closestSpotsInFN is interestedSpot, if not, disgard
#         closestSpotInSpots <- closestSpot(closestSpotsInFN, refSpots)
#         # test if the distance between the interestedSpot and closestSpotsInFN is < secondaryCospotDist, if not, disgar
#         iDistanceSmall <- which(distance(interestedSpots, closestSpotsInFN) < secondaryCospotDist &
#                                     closestSpotInSpots$Pos_X == interestedSpots$Pos_X &
#                                     closestSpotInSpots$Pos_Y == interestedSpots$Pos_Y)        
#         
#         interestedSpots[iDistanceSmall, iColumn_transFNNumCoSpots[i]] <- 1
#         interestedSpots[iDistanceSmall, iColumn_transFNCoSpotsIndice[i]] <- closestSpotIndexInFN[iDistanceSmall]
#         spots[indexOfZeroSpots, ] <- interestedSpots
#     }
   
}


# Fill in one more coSpot if it's not far from the refSpot and the refSpot has exisiting x coSpot
dao_fillInOneCospotToXCospotsWithUnusedSpots = function (ref_spots, x, secondaryCospotDist = 1, 
                                                         comparedRefSpots = ref_spots) {
    for (mergedName in mergedNames) {
        indexOfXSpots <- which(ref_spots[paste("transMergedNumCoSpots", mergedName, sep = "")] == x)
        interestedSpots <- ref_spots[indexOfXSpots,]
        closestSpotsInFN <- closestSpot(interestedSpots, transMergedSpots[[mergedName]], 
                                        countOnlyUnusedSpots = TRUE)
        closestSpotIndexInFN <- closestSpotIndex(interestedSpots, transMergedSpots[[mergedName]], 
                                                 countOnlyUnusedSpots = TRUE)        
        # Test if the closest spot to closestSpotsInFN is interestedSpot, if not, disgard
        closestSpotInRefSpots <- closestSpot(closestSpotsInFN, comparedRefSpots)
        
        # Test if the distance between the interestedSpot and closestSpotsInFN is < secondaryCospotDist, if not, disgard
        iDistanceSmall <- which((distance(interestedSpots, closestSpotsInFN) < secondaryCospotDist) &
                                    (closestSpotInRefSpots$Pos_X == interestedSpots$Pos_X) &
                                    (closestSpotInRefSpots$Pos_Y == interestedSpots$Pos_Y)) 
        
        # Update the columns of transMergedNumCoSpots and transMergedCoSpotsIndice of ref_spots
        if (length(iDistanceSmall) != 0) {
            interestedSpots[iDistanceSmall, paste("transMergedNumCoSpots", mergedName, sep = "")] <- 
                mapply(sum, interestedSpots[iDistanceSmall, paste("transMergedNumCoSpots", mergedName, sep = "")],
                       1, SIMPLIFY=TRUE)
            interestedSpots[iDistanceSmall,][[paste("transMergedCoSpotsIndice", mergedName, sep = "")]] <- 
                mapply(c, interestedSpots[iDistanceSmall, paste("transMergedCoSpotsIndice", mergedName, sep = "")], 
                       closestSpotIndexInFN[iDistanceSmall], SIMPLIFY=FALSE)
            ref_spots[indexOfXSpots, ] <- interestedSpots
        }   
    }
    ref_spots
}


# Debug
# Fill in the more coSpot if there are close unused coSpot within 1 um distance.
# Fill in by changing the column of transFNNumCoSpots and transFNCoSpotsIndice
# This might not only add 1 coSpot but several coSpots within the secondaryCospotDist
dao_fillInMoreCospotsWithUnusedSpots = function (refSpots, secondaryCospotDist = 1) {
    for (mergedName in mergedNames) {
        numOfClosestSpotsInFN <- num_neighbors_circle(refSpots, transMergedSpots[[mergedName]], 
                                                 secondaryCospotDist, countOnlyUnusedSpots = TRUE)
        closestSpotIndexInFN <- index_neighbors_circle(refSpots, transMergedSpots[[mergedName]], 
                                                       secondaryCospotDist, countOnlyUnusedSpots = TRUE)
        
        # Update the columns of transMergedNumCoSpots and transMergedCoSpotsIndice of refSpots
        refSpots[[paste("transMergedNumCoSpots", mergedName, sep = "")]] <- 
            refSpots[[paste("transMergedNumCoSpots", mergedName, sep = "")]] + numOfClosestSpotsInFN
        refSpots[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]] <- mapply(c,
            refSpots[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]], closestSpotIndexInFN)
    }
    refSpots
}



# For singlylLoneRefSpots (lonely_ref_spots that has 0 or 1 coSpot in each time), assign
# 1 color at each time point, either with FQ if numCoSpots = 1 or GaussianFit if numCoSpots = 0 
dao_assignAllColorsFromFQOrGaussianFit = function(lonely_ref_spots, coSpot_dist = 1, 
                                            minNumOfSpotFQFound = 1, isToPrintFit = FALSE) {
    allColors <- NULL
    for (i in seq_along(lonely_ref_spots[,1])) {
        print(i)
        refSpot <- lonely_ref_spots[i,]
        separation <- ""
        allColor <- NULL
        for(mergedName in mergedNames) {
            if (refSpot[paste("transMergedNumCoSpots", mergedName, sep = "")] >= minNumOfSpotFQFound) {
                transMergedSpotsIndice <- refSpot[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]][[1]]
                # Choose the first color, so it there are two coSpots, only the one with brighter index is chosen
                # assuming the index is ordered by cospot brightness
                color <- transMergedSpots[[mergedName]][transMergedSpotsIndice[1],]$Color                
            } else {
                color <- dao_chooseAColorFromGaussianFit(refSpot, mergedName, coSpot_dist, isToPrintFit)
                                                      
            }
            allColor <- paste(allColor, color, sep = separation)
            separation = " "
        }
        allColors[i] <- allColor
    }
    allColors
}



# Given a refSpot, fit Gaussian at all 4 images of that time point.
# Choose the color with highest Gaussian amp that is within coSpot_dist (1um) from the refSpot
# Fit a Gaussian first. 
# If a Gaussian fit return NA, substract nearby spots within 2um radius as Gaussians and refit.  
# If the Gaussian fit returns NA again, subtract nearby spots within 3um radius as Gaussians and refit.
# Note coSpotDistColor <- sqrt((loneColorGaussianX0 - refSpot$Pos_X)^2 + (loneColorGaussianY0 - refSpot$Pos_Y)^2)
# or coSpotDistColor <- sqrt((loneColorGaussianX0 - transRefSpotOnColor$Pos_X)^2 + (loneColorGaussianY0 - transRefSpotOnColor$Pos_Y)^2)
dao_chooseAColorFromGaussianFit = function(refSpot, mergedName, coSpot_dist = 1, isToPrintFit = FALSE, 
                                           SNRThreshold = c(0, 0, 0, 0), numOfColors = 1) {
    
    # Gaussian Amps & Bgs
    loneGaussianFits <- dao_4GaussianFitsWithMergedName(refSpot, mergedName, isToPrintFit)    
    loneGaussianAmps <- mapply(getAmpFromGaussianFit, loneGaussianFits)
    loneGaussianBgs <- mapply(getBackgroundFromGaussianFit, loneGaussianFits)
    
    # Normalized Gaussian Amps
    mergedMaxMinusBgTiffAmpPerColor4 <- mergedMaxMinusBgTiffAmpPerColor[[mergedName]]
    normalizedLoneGaussianAmps <- loneGaussianAmps / mergedMaxMinusBgTiffAmpPerColor4
    
    # Get transRefSpots
    currentFluorNames <- fluorNames[substr(fluorNames,2,2) == substr(mergedName,2,2)]
    transRefSpotOnFluors <- list()
    for(currentFluorName in currentFluorNames) {
        transRefSpotOnFluors[[currentFluorName]] <- moveRefSpotsByTransXYPerPolygon(
            refSpot, polygonMatrix, transFluorXYPerPolygon[[currentFluorName]])
    }
    
    # Only wants fitted Gaussian centers to be close to the transRefSpot
    coSpotDists <- mapply(dao_spot_fit_distance, loneGaussianFits, transRefSpotOnFluors)
    
    # Debug, because I already subtract nearby bright spots when necessary, might not need this.
    normalizedLoneGaussianAmps[coSpotDists > coSpot_dist] = 0
    
    snr <- loneGaussianAmps / loneGaussianBgs
    normalizedLoneGaussianAmps[snr < SNRThreshold] = 0
    colorOrderIndex <- order(normalizedLoneGaussianAmps, decreasing = TRUE)
    color <- colorNames[colorOrderIndex[c(1:numOfColors)]]
    color
}


# Get the colorAmps of 4 colors at the time point (mergedName), given a refSpot
# The fitted Gaussian center must be within coSpot_dist of the refSpot
dao_4NormalizedGaussianAmp = function(refSpot, mergedName, coSpot_dist = 1, isToPrintFit = FALSE) {  
    # Gaussian Amps
    loneGaussianFits <- dao_4GaussianFitsWithMergedName(refSpot, mergedName, isToPrintFit)    
    loneGaussianAmps <- mapply(getAmpFromGaussianFit, loneGaussianFits)
    
    # Normalized Gaussian Amps
    mergedMaxMinusBgTiffAmpPerColor4 <- mergedMaxMinusBgTiffAmpPerColor[[mergedName]]
    normalizedLoneGaussianAmps <- loneGaussianAmps / mergedMaxMinusBgTiffAmpPerColor4
    
    # Only wants fitted Gaussian centers to be close to the transRefSpot
    currentFluorNames <- fluorNames[substr(fluorNames,2,2) == substr(mergedName,2,2)]
    transRefSpotOnFluors <- list()
    for(currentFluorName in currentFluorNames) {
        transRefSpotOnFluors[[currentFluorName]] <- moveRefSpotsByTransXYPerPolygon(
            refSpot, polygonMatrix, transFluorXYPerPolygon[[currentFluorName]])
    }
    
    # Debug: I used to set coSpotDists = c(0,0,0,0)
    coSpotDists <- mapply(dao_spot_fit_distance, loneGaussianFits, transRefSpotOnFluors)
    normalizedLoneGaussianAmps[coSpotDists > coSpot_dist] = 0
    as.numeric(normalizedLoneGaussianAmps)
}

# get fluorNames from mergedName
dao_mergedNameToFluorNames = function (mergedName) {
    fluorNamesIndice <- substr(fluorNames,2,2) == substr(mergedName,2,2)
    currentFluorNames <- fluorNames[fluorNamesIndice]
    currentFluorNames
}


# Given a refSpot and mergeName (indicating time point), fit the 4 Gaussians for 4 colors
# at that time point. The Gaussisans is a named list; names are g1, y1, o1, r1, etc.
dao_4GaussianFitsWithMergedName = function(refSpot, mergedName, isToPrintFit = FALSE) {
    fluorNamesIndice <- substr(fluorNames,2,2) == substr(mergedName,2,2)
    transFluorXYPerPolygons4 <- transFluorXYPerPolygon[fluorNamesIndice]
    fluorTiffMatrices4 <- fluorTiffMatrices[fluorNamesIndice]
    fluorSpots4 <- fluorSpots[fluorNamesIndice]
    currentFluorNames <-  fluorNames[fluorNamesIndice]
    mergedMaxMinusBgTiffAmpPerColor4 <- mergedMaxMinusBgTiffAmpPerColor[[mergedName]]
    names(mergedMaxMinusBgTiffAmpPerColor4) <- names(transFluorXYPerPolygons4)
    loneGaussianFits <- dao_4GaussiansFits(refSpot, transFluorXYPerPolygons4, 
                                    fluorTiffMatrices4, fluorSpots4, isToPrintFit,
                                    currentFluorNames, mergedMaxMinusBgTiffAmpPerColor4)
    loneGaussianFits
}

# Given a refSpot, transFluorXYPerPolygons4, fluorTiffMatrices4, and scaled
# Return the colorName for the fluorTiffMatrix that contain the highest amplitude
dao_brightestColorName = function(refSpot, transFluorXYPerPolygons4, fluorTiffMatrices4, 
                                  mergedMaxMinusBgTiffAmpPerColor4) {
    scakedAmps4 <- numeric()
    for (fluorName in names(transFluorXYPerPolygons4)) {
        transRefNeighborsOnFluor <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, 
                                                                    transFluorXYPerPolygons4[[fluorName]])
        transRefNeighborsIndice <- spotsToIndice2D(transRefNeighborsOnFluor, pix_XY)
        amp <- fluorTiffMatrices4[[fluorName]][transRefNeighborsIndice]
        nearbyMinAmp <- getMinNearbyAmp(fluorTiffMatrices4[[fluorName]], transRefNeighborsIndice, radius = 1/pix_XY)
        scaledAmp <- (amp - nearbyMinAmp) / mergedMaxMinusBgTiffAmpPerColor4[fluorName]
        scakedAmps4[fluorName] <- scaledAmp
    }
    names(which.max(scakedAmps4))
}


# Given a refSpot, fit Gaussian at all 4 images of that time point.
# translate the refSpot to the fluorSpots coordinates and fit within diffraction limit of the TiffMatrix
# If a Gaussian fit return NA, substract nearby spots within 2um radius as Gaussians and refit.  
# If the Gaussian fit returns NA again, subtract nearby spots within 3um radius as Gaussians and refit.
dao_4GaussiansFits = function(refSpot, transFluorXYPerPolygons4, fluorTiffMatrices4,
                              fluorSpots4, isToPrintFit = FALSE, currentFluorNames,
                              mergedMaxMinusBgTiffAmpPerColor4) {
    
    loneGaussianFits <- list()
    for (i in seq_along(colorNames)) {
        loneGaussianFit <- NA
        fluorName <- names(transFluorXYPerPolygons4)[i]
        colorName <- names(diffractionLimits)[i]
        
        transRefSpotOnFluor <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, 
                                                               transFluorXYPerPolygons4[[fluorName]])
        fluorTiffMatrix <- fluorTiffMatrices4[[fluorName]]
        fitDistance <- diffractionLimits[[colorName]]
        closetFluorSpot <- closestSpot(transRefSpotOnFluor, fluorSpots4[[fluorName]])
        
        # debug
        bg0 <- closetFluorSpot$MSKY  # If the closestSpot is very far, than bg0 is not correct
        minNearbyIntensity <- getMinNearbyAmp(fluorTiffMatrix, spotsToIndice2D(transRefSpotOnFluor, pix_XY), radius = 3 / pix_XY)
        bg0 <- ifelse(spot_distance(transRefSpotOnFluor, closetFluorSpot) < 5, bg0, minNearbyIntensity)

        loneGaussianFit <- fitGaussian2D(transRefSpotOnFluor, fluorTiffMatrix, fitDistance)
        #loneGaussianFit <- dao_fitGaussian2DFixBg(transRefSpotOnFluor, fluorTiffMatrix, fitDistance, bg0 = bg0)
        
        # debug to remove or add later
        # get nearby neighbors, if the neighbor is used, fit the neighbor with fix-centered gaussian and substract it
        # neighborSpots <- neighbors_exclude_spot_circle(transRefSpotOnColor, colorSpotsList[[i]], neighbor_dist = 1)
        # usedNeighborSpots <- neighborSpots %>% filter(isUsed = TRUE)
        
        # debug for the isBadFit criteria.
        # Check if the previous Gaussian fit is a good fit
        # Since there is not a coSpot, the numNeighbor only include neighboring spot
        numNeighbors <- num_neighbors_circle(transRefSpotOnFluor, fluorSpots4[[fluorName]], 2)
        isBadFit <- FALSE
        if (!is.na(loneGaussianFit)[1]) {
            sigma <- getSigmaRFromGaussianFit(loneGaussianFit)
            numRefNeighbors <- num_neighbors_exclude_spot_circle(refSpot, refSpots, 2)
            offDistance <- sqrt((transRefSpotOnFluor$Pos_X - getX0FromGaussianFit(loneGaussianFit))^2 
                                + (transRefSpotOnFluor$Pos_Y - getY0FromGaussianFit(loneGaussianFit))^2)
            isBadFit <- (sigma >= 0.45) & (numNeighbors >= 1 | numRefNeighbors >= 1) & (offDistance > coSpotDist) 
            # debug: use to be offDistance > fitDistance
        }

        # Refit by subtrating the neighbor fluorSpots near the transRefSpotOnFluor
        if (is.na(loneGaussianFit)[1] | isBadFit) {
            # I changed from refSpot to transRefSpotOnFluor           
            numRefNeighbors <- num_neighbors_exclude_spot_circle(refSpot, refSpots, 2)
            neighborSpots <- neighbors_circle(transRefSpotOnFluor, fluorSpots4[[fluorName]], 2)
            # This loneGuassian Fit are fixed centers and bg, so getX0FromGaussianFit or getY0FromGaussianFit are NA
            if (nrow(neighborSpots) > 0) {
                loneGaussianFit <- dao_fitGuassianWithSubtraction(
                    transRefSpotOnFluor, neighborSpots, fluorTiffMatrix, fitDistance, 2, bg0)
            }   
        }
            
        # Refit by subtracting the neighbor refSpots near the transRefSpotOnFluor
        # It happens when the DAOSTORM didn't get the neighborFluorSpots but got neightborRefSpots
        # Need to get the amplifitude at the neightobrRefSpots for the four fluorplots, and see which fluorColor the 
        # neighborRefSpots is at.
        if (is.na(loneGaussianFit)[1] & numRefNeighbors > numNeighbors) {
            transRefAmpOnFluorMatrix <- fluorTiffMatrix[spotsToIndice2D(transRefSpotOnFluor, pix_XY)]
            refNeighbors <- neighbors_exclude_spot_circle(refSpot, refSpots, 2)
            transRefNeighborsToSubtract <- spotDataFrame(refNeighbors)
            for (j in seq_along(refNeighbors[,1])) {
                refNeighbor <- refNeighbors[j,]
                brightestColorName <- dao_brightestColorName(refNeighbor, transFluorXYPerPolygons4, 
                                                             fluorTiffMatrices4, mergedMaxMinusBgTiffAmpPerColor4)
                transRefNeighborOnFluor <- moveRefSpotsByTransXYPerPolygon(refNeighbor, polygonMatrix, 
                                                                            transFluorXYPerPolygons4[[fluorName]])
                closeFluorSpotToTransRefNeighbor <- closestSpot(transRefNeighborOnFluor, fluorSpots4[[fluorName]])
                spotDistance <- spot_distance(transRefNeighborOnFluor, closeFluorSpotToTransRefNeighbor)
                transRefNeighborAmpOnFluorMatrix <- fluorTiffMatrix[spotsToIndice2D(transRefNeighborOnFluor, pix_XY)]
                if (brightestColorName == fluorName & spotDistance > 0.5 & 
                        transRefNeighborAmpOnFluorMatrix > brightestColorName) {
                    transRefNeighborsToSubtract <- rbind(transRefNeighborsToSubtract, transRefNeighborOnFluor)
                }                 
            }
            #interferingRefNeighborsOnFluor <- transRefNeighborsToSubtract[
            #    fluorTiffMatrix[spotsToIndice2D(transRefNeighborsToSubtract, pix_XY)] > 
            #        fluorTiffMatrix[spotsToIndice2D(transRefSpotOnFluor, pix_XY)],]
            if (nrow(transRefNeighborsToSubtract) > 0) {
                # Might need to debug to set multiple Gaussian fit instead of subtraction
                # Otherwise the subtraction might be too strong and reduce refSpot brightness
                loneGaussianFit <- dao_fitGuassianWithSubtraction(transRefSpotOnFluor, interferingRefNeighborsOnFluor,
                                                                  fluorTiffMatrix, fitDistance, 2, bg0)
                #GaussianFitFor2 <- fitGaussian2DFor2(transRefSpotOnFluor, interferingRefNeighborsOnFluor, 
                #                                       fluorTiffMatrix, fitDistance)
            }
        }
        
        # debug to remove or not
        #if (is.na(loneGaussianFit)) {
        #    # I changed from refSpot to transRefSpotOnFluor
        #    loneGaussianFit <- dao_fitGuassianWithSubtraction(transRefSpotOnFluor, neighborSpots,
        #                                                      fluorTiffMatrix, fitDistance, 3, bg0)
        #}
        
        #        loneColorGaussianX0 <- getX0FromGaussianFit(loneColorFit)
        #        loneColorGaussianY0 <- getY0FromGaussianFit(loneColorFit)
        if (isToPrintFit) {
            print(paste(colorNames[i]))
            print(loneGaussianFit)
        }
        loneGaussianFits[[currentFluorNames[i]]] <- loneGaussianFit
    }
    loneGaussianFits
}


# distance between a spot and the center of a GaussianFit
dao_spot_fit_distance = function(loneGaussianFit, spot) {
    xFit <- getX0FromGaussianFit(loneGaussianFit)
    yFit <- getY0FromGaussianFit(loneGaussianFit)
    xSpot <- spot$Pos_X
    ySpot <- spot$Pos_Y
    distance <- sqrt((xSpot - xFit)^2 + (ySpot - yFit)^2)
    as.numeric(distance)
}


# Fit gaussian of the spot after subtracting nearby bright FQ/Dao identified spots as Gaussians
# And by substract nearby fake bright spots that has local maximum as Gaussians
# Fit Gaussian with fixed center
dao_fitGuassianWithSubtraction = function (spot, neighborSpots, fluorTiffMatrix, 
                                           fitDistance, neighborDistance, bg0) {
    
    # Subtract neighbor fluorSpots1 as Gaussians from the tiffMatrix
    #neighborSpots <- neighbors_exclude_spot_circle(spot, fluorSpots1, neighborDistance)
    subtractedTiffMatrix <- dao_subtractGaussianFitInOrder(fluorTiffMatrix, 
                                                       neighborSpots, neighborDistance, bg0)
    
    # debug
    # Subtract neighbor fake brighter local maximum spots as Guassian from the substractedTiffMatrix
    #localMaximumIndice2D <- getLocalMaximumIndice2D(subtractedTiffMatrix, spot, 
    #                                                neighborDistance)
    #extraSpots <- createSpotsFrom2DIndice(subtractedTiffMatrix, localMaximumIndice2D)
    #extraBrightSpots <- extraSpots[extraSpots$tiffAmp > subtractedTiffMatrix[
    #    spotsToIndice2D(spot, pix_XY)],]
    #subtractedTiffMatrix2 <- dao_subtractGaussianFitInOrder(subtractedTiffMatrix,
    #                                                    extraBrightSpots, fitDistance)

    # debug
    #loneColorFit <- fitGaussian2DFixCenter(spot, subtractedTiffMatrix2, fitDistance,
    #                                       trace = FALSE)
    #loneColorFit <- fitGaussian2DFixCenter(spot, subtractedTiffMatrix, fitDistance,
    #                                       trace = FALSE)
    loneColorFit <- dao_fitGaussian2DFixCenterAndBg(spot, subtractedTiffMatrix, fitDistance,
                                           trace = FALSE, bg0 = bg0)
    loneColorFit  
    #plotNearbyTiffRaster(spot, extraBrightSpots, pix_XY, subtractedTiffMatrix, isToTranslateRefSpots = FALSE, plot_area_dist = 5)
    #plotNearbyTiffRaster(spot, extraBrightSpots, pix_XY, subtractedTiffMatrix2, isToTranslateRefSpots = FALSE, plot_area_dist = 5)
}



# Mark the bleakSpot in trasMergedSpots[[mergedName]] as bleak
# Givin the crowd2RefSpots
# If the coSpots contain a bleakSpot, mark the bleakSpot in trasMergedSpots[[mergedName]] as bleak
dao_markTransMergedSpotsAsBleack = function (crowd2RefSpots, transMergedSpots, mergedName, 
                                             bleakColor, trueColor, ampThreshold) {
    for (r in seq_along(crowd2RefSpots[,1])) {
        spot <- crowd2RefSpots[r,]
        # Check that the number of coSpots is 2
        if (spot[paste("transMergedNumCoSpots", mergedName, sep = "")] == 2) {
            transSpotsIndice <- spot[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]][[1]]
            bleakSpotIndex <- transSpotsIndice[transMergedSpots[[mergedName]]
                                               [transSpotsIndice,]$Color == bleakColor]
            trueSpotIndex <- transSpotsIndice[transMergedSpots[[mergedName]]
                                              [transSpotsIndice,]$Color == trueColor]
            bleakColorSpot <- transMergedSpots[[mergedName]][bleakSpotIndex,]
            trueColorSpot <- transMergedSpots[[mergedName]][trueSpotIndex,]
            
            # Check that one spot is yellow and one spot is orange
            if (nrow(bleakColorSpot) == 1 & nrow(trueColorSpot) == 1) {
                bleackColorAmp <- bleakColorSpot$ScaledAmp
                trueColorAmp <- trueColorSpot$ScaledAmp
                
                # Check that the yellow cospot is at least 7 times dimmer than the orange cospot
                if (trueColorAmp/bleackColorAmp > ampThreshold) {
                    #crowd2RefSpots[r,][paste("transMergedNumCoSpots",mergedName, sep="")] <- 1
                    #crowd2RefSpots[r,][[paste("transMergedCoSpotsIndice",mergedName, sep="")]][[1]] <- trueSpotIndex
                    transMergedSpots[[mergedName]][bleakSpotIndex,]$isBleak <- TRUE
                }
            }
        }
    }
    transMergedSpots
}


# Remove the bleak coSpots from the refSpot's column of $transMergedNumCoSpots  and $transMergedCoSpotsIndice
# using the transMergedSpots which has the isBleak column
dao_removeBleakCoSpots = function(crowd2RefSpots, transMergedSpots, mergedName) {
    for (r in seq_along(crowd2RefSpots[,1])) {
        spot <- crowd2RefSpots[r,]
        transSpotsIndice <- spot[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]][[1]]
        for (transSpotIndex in transSpotsIndice) {
            transSpot <- transMergedSpots[[mergedName]][transSpotIndex,]
            if (transSpot$isBleak) {
                #print(r)
                #print(transSpotIndex)
                crowd2RefSpots[r,][paste("transMergedNumCoSpots",mergedName, sep="")] <- 
                    crowd2RefSpots[r,][paste("transMergedNumCoSpots",mergedName, sep="")] - 1
                crowd2RefSpots[r,][[paste("transMergedCoSpotsIndice",mergedName, sep="")]][[1]] <- 
                    crowd2RefSpots[r,][[paste("transMergedCoSpotsIndice",mergedName, sep="")]][[1]][
                        crowd2RefSpots[r,][[paste("transMergedCoSpotsIndice",mergedName, sep="")]][[1]] != transSpotIndex]
            }
        }
    }
    crowd2RefSpots
}

# Givine transMergedCoSpotsIndice1, which is a list, and each list element is a vector
# Reorder the vector so the index with brighter coSpot in the transMergedSpots is at front
dao_reorderIndiceByBrightness = function (transMergedCoSpotsIndice1, transMergedSpots1) {
    for (i in seq_along(transMergedCoSpotsIndice1)) {
        coSpotIndice <- transMergedCoSpotsIndice1[[i]]
        brightness <- transMergedSpots1[coSpotIndice,]$ScaledAmp
        transMergedCoSpotsIndice1[[i]] <- coSpotIndice[order(brightness, decreasing = TRUE)]
    }
    transMergedCoSpotsIndice1
}

# Very simple assignment: assuming daostorm doesn't miss one color of a two color coSpots.
# Debug: Need modification later....
# Assign allColors for doublet refSpots as list of vectors, with each vector contains 2 string elements
# When there is no cospots, use 2D Gaussian to find a color.
# When thers is a coSpot, use it for two color sequences.
# When there are two coSpots, use the two color of the two coSpots.
dao_assignAllColorsForDoublets = function(doubletRefSpotsInCellNoBleak, isToPrintFit = FALSE) {
    allColors <- list()
    for (i in seq_along(doubletRefSpotsInCellNoBleak[,1])) {
        refSpot <- doubletRefSpotsInCellNoBleak[i,]
        separation <- ""
        allColor <- NULL
        for (mergedName in mergedNames) {
            if (refSpot[paste("transMergedNumCoSpots", mergedName, sep = "")] != 2) {
                colorAmps <- dao_4NormalizedGaussianAmp(refSpot, mergedName)
                
                colorsInOrder <- colorNames[order(colorAmps, decreasing = TRUE)]
                colorAmpsInOrder <- colorAmps[order(colorAmps, decreasing = TRUE)]
                bestColor <- colorsInOrder[1]
                secondBestColor <- colorsInOrder[2]
                bestColorAmp <- colorAmpsInOrder[1]
                secondBestColorAmp <- colorAmpsInOrder[2]
                if (refSpot[paste("transMergedNumCoSpots", mergedName, sep = "")] == 0) {
                    if (bestColorAmp / secondBestColorAmp > 6.5) {
                        color <- c(bestColor, bestColor)
                    } else {
                        color <- c(bestColor, secondBestColor)
                    }
                } else {
                    singletSpotIndex <- refSpot[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]][[1]]
                    color <- transMergedSpots[[mergedName]][singletSpotIndex,]$Color
                    if (color == bestColor) {
                        if (bestColorAmp / secondBestColorAmp > 6.5) {
                            color <- c(bestColor, bestColor)
                        } else {
                            color <- c(bestColor, secondBestColor)
                        }
                    } else if (color == secondBestColor) {
                        color <- c(color, secondBestColor)
                    } else {
                        color <- c(color, color) # Debug
                    }
                }
            } else {
                doubletSpotIndice <- refSpot[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]][[1]]
                color <- transMergedSpots[[mergedName]][doubletSpotIndice,]$Color
            }
            allColor <- paste(allColor, color, sep = separation)
            separation <- " "
        }
        allColors[[i]] <- allColor
    }
    allColors
}

# Give doubletSpots, if any 2 doublets at a time point is too close, the spot is doublet
dao_getIndiceOfCloseDoubletSpots = function(doubletSpots, closeDist = 0.1) {
    refCloseDoubletSpotsIndice <- NULL
    for (mergedName in mergedNames) {
        fnCosposNum <- doubletSpots[[paste("transMergedNumCoSpots", mergedName, sep = "")]]
        fnCospotsIndice <- doubletSpots[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]]
        refDoubletSpotsIndice <- which(fnCosposNum == 2)
        fnCoSpotsIndiceForRefDoubletSpots <- fnCospotsIndice[refDoubletSpotsIndice]
        firstSpotIndice <- unlist(lapply(fnCoSpotsIndiceForRefDoubletSpots, '[[', 1))
        secondSpotIndice <- unlist(lapply(fnCoSpotsIndiceForRefDoubletSpots, '[[', 2))
        distance <- spot_distance(transMergedSpots[[mergedName]][firstSpotIndice, ], 
                                  transMergedSpots[[mergedName]][secondSpotIndice, ])
        refCloseDoubletSpotsIndice <- c(refCloseDoubletSpotsIndice, 
                                        refDoubletSpotsIndice[distance < closeDist]) 
    }
    unique(refCloseDoubletSpotsIndice)
}


# givin the centerSpot, the tiffMatrix, the initial parameters
# Fit Gaussian2D at the spot's indice2D
# If the nls fit doesn't work, don't print out error message, instead,
# return NA
# Note that this Gaussian doesn't check whether it's a good Gaussian fit or not ...
dao_fitGaussian2DFixBg = function(spot, tiffMatrix, fit_distance_um,
                                lowerSigmaBound = 0.25, upperSigmaBound = 0.5, trace = FALSE,
                                bg0) {
    
    # get the matrixCoordX and matrixCoordY in um
    centerIndice <- spotsToIndice2D(spot, pix_XY)
    fit_distance <- floor(fit_distance_um / pix_XY)
    nearbyIndice2D <- getNearbyIndice2D(tiffMatrix, centerIndice, fit_distance)
    matrixCoordX <- getXCoordsFromIndice2D(nearbyIndice2D)
    matrixCoordY <- getYCoordsFromIndice2D(nearbyIndice2D)
    
    # get the Gaussian data of neighbor pixels, unit is intensity
    Gaussian <- tiffMatrix[nearbyIndice2D]
    
    # guess the init parameters    
    nearbyRadius <- fit_distance_um # 2
    avgRadius <- fit_distance_um # 2
    maxAvgIndice <- getMaxNearbyAvgIndice(tiffMatrix, centerIndice, nearbyRadius, avgRadius)
    maxAvgAmp <- getMaxNearbyAvgAmp(tiffMatrix, centerIndice, nearbyRadius, avgRadius)
    amp <- maxAvgAmp - bg0
    x0 <- getXCoordsFromIndice2D(maxAvgIndice)
    y0 <- getYCoordsFromIndice2D(maxAvgIndice)
    parametersInt <- c(amp, x0, y0, 0.3)
    
    lowerParameterBound <- c(0, x0-fit_distance_um, y0-fit_distance_um, lowerSigmaBound)
    upperParameterBound <- c(Inf, x0+fit_distance_um, y0+fit_distance_um, upperSigmaBound)
    
    # fit with nls
    fit <- NA
    try(fit <- nls(Gaussian ~ Gaussian2D(amp, x0, y0, sigma_r, bg0, matrixCoordX, matrixCoordY), 
                   start = list(amp = parametersInt[1], x0 = parametersInt[2], y0 = parametersInt[3],
                                sigma_r = parametersInt[4]), 
                   lower = list(amp = lowerParameterBound[1], x0 = lowerParameterBound[2], y0 = lowerParameterBound[3],
                                sigma_r = lowerParameterBound[4]),
                   upper = list(amp = upperParameterBound[1], x0 = upperParameterBound[2], y0 = upperParameterBound[3],
                                sigma_r = upperParameterBound[4]), 
                   algorithm = "port", trace = trace),
        silent = TRUE) 
    fit
}


# givin the centerSpot, the tiffMatrix, the initial parameters
# Fit Gaussian2D at the spot's indice2D
# If the nls fit doesn't work, don't print out error message, instead,
# return NA
# Note that this Gaussian doesn't check whether it's a good Gaussian fit or not ...
dao_fitGaussian2DFixCenterAndBg = function(spot, tiffMatrix, fit_distance_um,
                                           lowerSigmaBound = 0.25, upperSigmaBound = 0.5, trace = FALSE,
                                           bg0) {
    
    # get the matrixCoordX and matrixCoordY in um
    centerIndice <- spotsToIndice2D(spot, pix_XY)
    fit_distance <- floor(fit_distance_um / pix_XY)
    nearbyIndice2D <- getNearbyIndice2D(tiffMatrix, centerIndice, fit_distance)
    matrixCoordX <- getXCoordsFromIndice2D(nearbyIndice2D)
    matrixCoordY <- getYCoordsFromIndice2D(nearbyIndice2D)
    
    # get the Gaussian data of neighbor pixels, unit is intensity
    Gaussian <- tiffMatrix[nearbyIndice2D]
    
    # guess the init parameters    
    x0 <- spot$Pos_X
    y0 <- spot$Pos_Y
    centerAmp <- tiffMatrix[centerIndice]
    #bg0 <- spot$MSKY
    amp <- centerAmp - bg0
    parametersInt <- c(amp, 0.3)
    
    lowerParameterBound <- c(0, lowerSigmaBound)
    upperParameterBound <- c(Inf, upperSigmaBound)
    
    # fit with nls
    fit <- NA
    try(fit <- nls(Gaussian ~ Gaussian2D(amp, x0, y0, sigma_r, bg0, matrixCoordX, matrixCoordY), 
                   start = list(amp = parametersInt[1], sigma_r = parametersInt[2]), 
                   lower = list(amp = lowerParameterBound[1], sigma_r = lowerParameterBound[2]),
                   upper = list(amp = upperParameterBound[1], sigma_r = upperParameterBound[2]), 
                   control = nls.control(maxiter = 100), algorithm = "port", na.action=na.exclude,trace = trace),
        silent = FALSE) 
    fit
}


dao_fitGaussian2DFixCenterAmpAndBg = function(spot, tiffMatrix, fit_distance_um,
                                           lowerSigmaBound = 0.25, upperSigmaBound = 0.5, trace = FALSE,
                                           bg0) {
    
    # get the matrixCoordX and matrixCoordY in um
    centerIndice <- spotsToIndice2D(spot, pix_XY)
    fit_distance <- floor(fit_distance_um / pix_XY)
    nearbyIndice2D <- getNearbyIndice2D(tiffMatrix, centerIndice, fit_distance)
    matrixCoordX <- getXCoordsFromIndice2D(nearbyIndice2D)
    matrixCoordY <- getYCoordsFromIndice2D(nearbyIndice2D)
    
    # get the Gaussian data of neighbor pixels, unit is intensity
    Gaussian <- tiffMatrix[nearbyIndice2D]
    
    # guess the init parameters    
    x0 <- spot$Pos_X
    y0 <- spot$Pos_Y
    centerAmp <- tiffMatrix[centerIndice]
    #bg0 <- spot$MSKY
    amp0 <- centerAmp - bg0
    parametersInt <- c(0.3)
    
    lowerParameterBound <- c(lowerSigmaBound)
    upperParameterBound <- c(upperSigmaBound)
    
    # fit with nls
    fit <- NA
    try(fit <- nls(Gaussian ~ Gaussian2D(amp0, x0, y0, sigma_r, bg0, matrixCoordX, matrixCoordY), 
                   start = list(sigma_r = parametersInt[1]), 
                   lower = list(sigma_r = lowerParameterBound[1]),
                   upper = list(sigma_r = upperParameterBound[1]), 
                   control = nls.control(maxiter = 100), algorithm = "port", na.action=na.exclude,trace = trace),
        silent = FALSE) 
    fit
}


# givin the centerSpot, the tiffMatrix, the initial parameters
# Fit Gaussian2D at the spot's indice2D
# If the nls fit doesn't work, don't print out error message, instead,
# return NA
# Note that this Gaussian doesn't check whether it's a good Gaussian fit or not ...
dao_fitGaussian2DFixCenterSigmaAndBg = function(spot, tiffMatrix, fit_distance_um, trace = FALSE) {
    
    # get the matrixCoordX and matrixCoordY in um
    centerIndice <- spotsToIndice2D(spot, pix_XY)
    fit_distance <- floor(fit_distance_um / pix_XY)
    nearbyIndice2D <- getNearbyIndice2D(tiffMatrix, centerIndice, fit_distance)
    matrixCoordX <- getXCoordsFromIndice2D(nearbyIndice2D)
    matrixCoordY <- getYCoordsFromIndice2D(nearbyIndice2D)
    
    # get the Gaussian data of neighbor pixels, unit is intensity
    Gaussian <- tiffMatrix[nearbyIndice2D]
    
    # guess the init parameters    
    x0 <- spot$Pos_X
    y0 <- spot$Pos_Y
    centerAmp <- tiffMatrix[centerIndice]
    sigma_r0 <- 0.3 # Debug to 0.21  * lambda / 0.4
    bg0 <- spot$MSKY
    amp <- centerAmp - bg0
    parametersInt <- c(amp)
    
    lowerParameterBound <- c(0)
    upperParameterBound <- c(Inf)
    
    # fit with nls
    fit <- NA
    try(fit <- nls(Gaussian ~ Gaussian2D(amp, x0, y0, sigma_r0, bg0, matrixCoordX, matrixCoordY), 
                   start = list(amp = parametersInt[1]), 
                   lower = list(amp = lowerParameterBound[1]),
                   upper = list(amp = upperParameterBound[1]), 
                   control = nls.control(maxiter = 100), algorithm = "port", na.action=na.exclude,trace = trace),
        silent = FALSE) 
    fit
}


# Given list of spots and a tiff matrix, subtract the Gaussian fit of each spot
# from strong intensity to weak intensity from the tiff matrix
# Return a subtracted tiff matrix
dao_subtractGaussianFitInOrder = function(tiffMatrix, spots, fit_distance_um, bg0) {
    currentTiffMatrix <- tiffMatrix
    if (nrow(spots) != 0) {
        orderedSpots <- spots[order(spots$tiffAmp, decreasing = TRUE),] 
        numSpots <- dim(orderedSpots)[1]
        amps <- vector(length = numSpots)
        sigma_rs <- vector(length = numSpots)
        bgs <- vector(length = numSpots)
        
        for (i in seq_along(orderedSpots[,1])) {
            spot <- orderedSpots[i,]
            
            # debug
            # Fit the Gaussian2DFixedCenter
            #curerntGaussianFit <- dao_fitGaussian2DFixCenterAndBg(spot, currentTiffMatrix, 
            #                                         fit_distance_um, 0.25, 0.5, bg0 = spot$MSKY)
            curerntGaussianFit <- dao_fitGaussian2DFixCenterAmpAndBg(spot, currentTiffMatrix, 
                                                                  fit_distance_um, 0.25, 0.5, bg0 = bg0)
            #amps[i] <- getAmpFromGaussianFit(curerntGaussianFit)
            sigma_rs[i] <- getSigmaRFromGaussianFit(curerntGaussianFit)
            amps[i] <- spot$tiffAmp - spot$MSKY            
            #sigma_rs[i]  <- 0.3
            # Subtract the Gaussian from the tiffplot
            centerIndice <- spotsToIndice2D(spot, pix_XY)
            nearbyIndice2D <- getNearbyIndice2D(currentTiffMatrix, centerIndice, 
                                                fit_distance_um / pix_XY)
            matrixCoordX <- getXCoordsFromIndice2D(nearbyIndice2D)
            matrixCoordY <- getYCoordsFromIndice2D(nearbyIndice2D)
            Gaussian <- Gaussian2D(amps[i], spot$Pos_X, spot$Pos_Y, sigma_rs[i], 
                                   0, matrixCoordX, matrixCoordY)
            currentTiffMatrix[nearbyIndice2D] <- currentTiffMatrix[nearbyIndice2D] -
                Gaussian
            #a<- matrix(0, size(currentTiffMatrix)[1], size(currentTiffMatrix)[2])
            #a[nearbyIndice2D] <- Gaussian
            #plotNearbyTiffRaster(spots, spots, pix_XY, a, isToTranslateRefSpots = FALSE, plot_area_dist = 2)
            #plotNearbyTiffRaster(spots, spots, pix_XY, tiffMatrix, isToTranslateRefSpots = FALSE, plot_area_dist = 2)        
            #plotNearbyTiffRaster(spots, spots, pix_XY, currentTiffMatrix, isToTranslateRefSpots = FALSE, plot_area_dist = 2)        
        }
    }
    currentTiffMatrix
}

# Return a indice2D for spots, this time don't round the pixel.
# Somehow for FQ spots, we need to add 1 pixel in both XY direction.
dao_spotsToIndice2D_not_integer = function (spots, pix_XY) {
    pixX <- spots$Pos_X / pix_XY
    pixY <- spots$Pos_Y / pix_XY
    indice2D <- getIndice2DFromXY(pixX, pixY)
    indice2D
}

# Plot all nearby Spots in transF1, transF2, transF3, and transF4,
# and plot all nearby Spots in 16 tiff matrix images with the refSpot pixel
# circled in red and neighbor fluor spots circled in white
dao_plotAllNearby = function (spot, plot_area_dist = 2, toTranslateRef = TRUE, fluorMatrices = fluorScaledMatrices) {
    
    zoomedMergedPlots <- list()
    for(mergedName in mergedNames) {
        zoomedMergedPlots[[mergedName]] <- plotNearbySpots(spot, transMergedSpots[[mergedName]], 
                                                          refSpots, plot_area_dist = plot_area_dist, 
                                                          paste("Time", substr(mergedName,2,2), "merged dots"))
    }
    grid.arrange(zoomedMergedPlots[["f1"]], zoomedMergedPlots[["f2"]], zoomedMergedPlots[["f3"]],
                 zoomedMergedPlots[["f4"]], ncol = 4)
    
    par(mfrow = c(1,1))
    plotNearbyTiffRaster(spot, refSpots, pix_XY, refScaledMatrix, isToTranslateRefSpots = FALSE, col = greenPalette, plot_area_dist = plot_area_dist, plotname = "nearby ref spots on tiff")
    par(mfrow = c(4,4))
    for(fluorName in fluorNames) {
        plotNearbyTiffRaster(spot, fluorSpots[[fluorName]], pix_XY, fluorMatrices[[fluorName]], 
                             transFluorXYPerPolygon[[fluorName]], isToTranslateRefSpots = toTranslateRef, 
                             col = palettes[[fluorName]], plot_area_dist, 
                             plotname = paste(fluorName, "spots on fluor image"))
    }
}

