source("/Users/jackie/Downloads/SGTC microscope/2016-08-11-10mRNAs vs 4 exon mRNAs-37C-45C Amp 3rd time/4 exon mRNAs/spotPlotFunctions.R")
library("Bessel")

# spotLength = function(spot) {
#     sqrt(spot$Pos_X^2 + spot$Pos_Y^2)
# }

# Debug to remove because of confusion
# Distance between two spots (um), could be a vector
distance = function(spot1, spot2) {
    sqrt((spot1$Pos_X - spot2$Pos_X)^2 + (spot1$Pos_Y - spot2$Pos_Y)^2)
}

# Distance between two spots (um), could be a vector
spot_distance = function(spot1, spot2) {
    sqrt((spot1$Pos_X - spot2$Pos_X)^2 + (spot1$Pos_Y - spot2$Pos_Y)^2)
}

# spot2$Pos_X - spot1$Pos_X, could be a vector
# If spots have no rows, the answer is numeric()
get_dx = function(spots1, spots2) {
    spots2$Pos_X - spots1$Pos_X
}

# spot2$Pos_Y - spot1$Pos_Y, could be a vector
# If spots have no rows, the answer is numeric()
get_dy = function(spots1, spots2) {
    spots2$Pos_Y - spots1$Pos_Y
}

# If spots have no rows, the answer is NaN
avg_dx = function(spots1, spots2) {
    mean(get_dx(spots1, spots2))
}

# If spots have no rows, the answer is NaN
avg_dy = function(spots1, spots2) {
    mean(get_dy(spots1, spots2))
}


# spot is the anchor, draw line between anchorSpot-spot2, and anchorSpot-spot1, find hte angle between the 2 lines
# angle2to1 = function(spot2, spot1, anchorSpot) {
#     angle2_xAxis = angleToXAxis(spot2, anchorSpot)
#     angle1_xAxis = angleToXAxis(spot1, anchorSpot)
#     angle2_xAxis - angle1_xAxis[,1]
# }

# if spot and anchorSpot are the same, will return NAN
# angleToXAxis = function(spot, anchorSpot) {
#     angle_xAxis <- atan((spot$Pos_Y - anchorSpot$Pos_Y) / (spot$Pos_X - anchorSpot$Pos_X))
#     toAdd <- angle_xAxis #matrix(0, dim(spot)[1])
#     toAdd[spot$Pos_X - anchorSpot$Pos_X < 0] <- pi
#     toAdd[!(spot$Pos_X - anchorSpot$Pos_X < 0)] <- 0
#     angle_xAxis = angle_xAxis + toAdd
#     angle_xAxis
# }

# can rotate multiple spots
# or can roate multiple angles
# rotate = function(spot, anchorSpot, angle) {
#     rotatedSpot <- spot
#     isSpotSameAsAnchor <- (spot$Pos_X == anchorSpot$Pos_X) & (spot$Pos_Y == anchorSpot$Pos_Y)
#     spotsSameAsAnchor <- spot[isSpotSameAsAnchor,]
#     spotsNotSameAsAnchor <- spot[!isSpotSameAsAnchor,]
#     rotatedSpot[isSpotSameAsAnchor,] <- data.frame(Pos_Y = spotsSameAsAnchor$Pos_Y,
#                                                    Pos_X = spotsSameAsAnchor$Pos_X, 
#                                                    AMP = spotsSameAsAnchor$AMP, 
#                                                    RES = spotsSameAsAnchor$RES, 
#                                                    Color = spotsSameAsAnchor$Color,
#                                                    tiffAMP = spotsSameAsAnchor$tiffAMP)
#     if (sum(isSpotSameAsAnchor) != nrow(spot)) {
#         rotatedSpot[!isSpotSameAsAnchor,] <- rotateNoSameSpot(spotsNotSameAsAnchor, 
#                                                               anchorSpot, angle)
#     }
#     rotatedSpot
# }

# angle should be a n x 1 matrix
# rotateNoSameSpot = function(spot, anchorSpot, angle) {
#     angle_xAxis <- as.vector(angleToXAxis(spot, anchorSpot))
#     newAngle <- angle_xAxis + as.vector(angle) #[,1]
#     length <- distance(spot, anchorSpot)
#     deltaX <- length * cos(newAngle)
#     deltaY <- length * sin(newAngle)
#     newX <- anchorSpot$Pos_X + deltaX
#     newY <- anchorSpot$Pos_Y + deltaY
#     data.frame(Pos_Y = newY, Pos_X = newX, AMP = spot$AMP, RES = spot$RES, Color = spot$Color,
#                tiffAMP = spot$tiffAMP)        
# }

# scale = function(spot, anchorSpot, scaler) {
#     deltaX <- spot$Pos_X - anchorSpot$Pos_X
#     deltaY <- spot$Pos_Y - anchorSpot$Pos_Y
#     deltaX <- scaler * deltaX
#     deltaY <- scaler * deltaY
#     newX <- anchorSpot$Pos_X + deltaX
#     newY <- anchorSpot$Pos_Y + deltaY
#     data.frame(Pos_Y = newY, Pos_X = newX, AMP = spot$AMP, RES = spot$RES, Color = spot$Color, 
#                tiffAMP = spot$tiffAMP)
# }


# Translate spots in x, y direction to fit to refSpots.  
# Translate spots in a cell at a time. 
# maxPreCospotDist is the max distance allowed between a initial spot and a refSpot, 
# if not within the distance, it's likely that the spot is dirt
# minSpotNumPerCell is the min number of spots in a cell, if not over the number,
# it's likely there are very few spots to have a good translation, so we might as well not translate.
# and also return a data.frame that has row number as polygon number and deltaX and deltaY column.
# The deltaX and teltaY (um) is to move spot to refSpot
translateAllCells = function (spots, refSpots, polygons_um, maxPreCospotDist, minSpotNumPerCell = 3) {
    spotsToTranslate <- spots
    toMoveXYPerCell <- data.frame(deltaX = numeric(), deltaY = numeric())

    # The minusBgTiffAmp threshold should be the dimmest 10% spots
    if (nrow(spots) != 0) {
        brightThresholdIndex <- order(spots$minusBgTiffAmp)[ceil(nrow(spots) * 0.1)]
        minusBgTiffAmpThreshold <- spots$minusBgTiffAmp[brightThresholdIndex]
    } else {
        minusBgTiffAmpThreshold <- 0
    }
    
    for (i in 1 : length(polygons_um)) {
        encSpots <- enclosedSpots(spots, polygons_um[[i]])
        encRefSpots <- enclosedSpots(refSpots, polygons_um[[i]])
        if (nrow(encSpots) >= minSpotNumPerCell) {
            # changed
            translateACellOutput <- translateACell(encSpots, encRefSpots, maxPreCospotDist, minusBgTiffAmpThreshold)
            #translateACellOutput <- translateACell(encSpots, encRefSpots, maxPreCospotDist)
            translatedEncSpots <- translateACellOutput[[1]]
            dx <- translateACellOutput[[2]]
            dy <- translateACellOutput[[3]]
            spotsToTranslate[in.out(polygons_um[[i]], coord(spots)),] <- translatedEncSpots
        } else {
            dx <- 0
            dy <- 0
        }
        toMoveXYPerCell <- rbind(toMoveXYPerCell, data.frame(deltaX = dx, deltaY = dy))
    }
    list(spotsToTranslate, toMoveXYPerCell)
}


# Translate spots in a cell to refSpots.
# maxPreCospotDist is the max distance allowed between a initial good spot and a refSpot, 
# if not within the distance, it's likely that the spot is dirt
# minusBgTiffAmpThreshold is the minimum minusBgTiffAmp of the encSpots, which are bright enough 
# to be used as ancher for translations, so we don't accidentally
# include overwhelming dim fake spots, which often appear near the side of nuclei.
# and also return the deltaX and deltaY to translate this cell.
translateACell = function (encSpots, encRefSpots, maxPreCospotDist, minusBgTiffAmpThreshold) {
    brightEncSpots <- encSpots %>% filter(minusBgTiffAmp > minusBgTiffAmpThreshold)
    closestRefSpots <- closestSpot (brightEncSpots, encRefSpots)
    distances <- distance(brightEncSpots, closestRefSpots)
    shortestDistanceOrder <- order(distances)
    numOfEncSpots <- nrow(brightEncSpots)
    numOfTopSpots <- 90  # used to be 10
    topEncSpots <- brightEncSpots[shortestDistanceOrder[1:min(numOfEncSpots, numOfTopSpots)],]
    
    currentShorestAvgDistBright <- avgDistance(brightEncSpots, closestRefSpots)
    currentDeltaX <- 0
    currentDeltaY <- 0
    currentMovedSpots <- encSpots
    for (i in seq_along(topEncSpots[,1])) {
        topEncSpot <- topEncSpots[i,]
        regSpotsAndAvgDistBright <- regSpotsAndAvgDistanceUsingASpot(topEncSpot, brightEncSpots, 
                                                               encRefSpots, maxPreCospotDist)
        
        shortestAvgDistBright <- regSpotsAndAvgDistBright[[2]]
        if (shortestAvgDistBright < currentShorestAvgDistBright) {
            currentShorestAvgDistBright <- shortestAvgDistBright
            regSpotsAndAvgDist <- regSpotsAndAvgDistanceUsingASpot(topEncSpot, encSpots, 
                                                                encRefSpots, maxPreCospotDist)
            currentMovedSpots <- regSpotsAndAvgDist[[1]]
            currentDeltaX <- regSpotsAndAvgDist[[3]]
            currentDeltaY <- regSpotsAndAvgDist[[4]]
        }
    }
    list(currentMovedSpots, currentDeltaX, currentDeltaY)
}

# get the best coSpot of a spot within maxPreCospotDist, so that 
# after translating all the spots, the average distance between all translatedSpots
# and their closest refSpots is shortest.
# and also return the deltaX and deltaY to register
regSpotsAndAvgDistanceUsingASpot = function (spot, spots, refSpots, maxPreCospotDist) {
    neighborRefSpots <- neighborhood(spot, refSpots, maxPreCospotDist)
    
    if (nrow(neighborRefSpots) == 0) {
        movedSpots <- spots
        shortestAvgDistance <- avgDistanceToClosestRefSpot(movedSpots, refSpots)
        deltaX <- 0
        deltaY <- 0
    } else {
        averageDistance <- 0
        for (i in seq_along(nrow(neighborRefSpots))) {
            neighborRefSpot <- neighborRefSpots[i,]
            movedSpots <- moveXYFromSpot1ToSpot2(spots, spot, neighborRefSpot)
            averageDistance[i] <- avgDistanceToClosestRefSpot(movedSpots, refSpots)
        }
        
        shortestDistanceIndex <- order(averageDistance)[i]
        bestRefSpotForReg <- neighborRefSpots[shortestDistanceIndex,]
        movedSpots <- moveXYFromSpot1ToSpot2(spots, spot, bestRefSpotForReg)
        shortestAvgDistance <- averageDistance[shortestDistanceIndex]
        deltaX <- bestRefSpotForReg$Pos_X - spot$Pos_X
        deltaY <- bestRefSpotForReg$Pos_Y - spot$Pos_Y
    } 
    list(movedSpots, shortestAvgDistance, deltaX, deltaY)
}


# Move all Spots with relative deltaY = spot2$Pos_Y - spot1$Pos_Y and 
# deltaX = spot2$Pos_X - spot1$Pos_X
moveXYFromSpot1ToSpot2 = function (spots, spot1, spot2) {
    deltaY <- spot2$Pos_Y - spot1$Pos_Y
    deltaX <- spot2$Pos_X - spot1$Pos_X
    movedSpots <- moveXY(spots, deltaX = deltaX, deltaY = deltaY)
    movedSpots
}


# spot could be multiple spots.
moveXY = function(spot, deltaX, deltaY) {
    movedSpot <- spot
    movedSpot$Pos_Y = spot$Pos_Y + deltaY
    movedSpot$Pos_X = spot$Pos_X + deltaX
    movedSpot
}


# spot is on plot, want to move it to refPlot. spot1 and spot are on plot while coSpot1
# anc coSpot2 are on ref plot
# moveWithCoPairs = function (spot, anchorSpots) {
#     spot1 <- anchorSpots[[1]] 
#     coSpot1 <- anchorSpots[[2]] 
#     spot2 <- anchorSpots[[3]]  
#     coSpot2 <- anchorSpots[[4]] 
#     
#     deltaX <- coSpot1$Pos_X - spot1$Pos_X
#     deltaY <- coSpot1$Pos_Y - spot1$Pos_Y
#     
#     refAngle21 <- angleToXAxis(coSpot2, coSpot1)
#     angle21 <- angleToXAxis(spot2, spot1)
#     angleDiff21 <- refAngle21 - angle21
#     
#     movedSpot <- moveXY(spot, deltaX, deltaY)
#     movedSpot <- rotate(movedSpot, coSpot1, angleDiff21)
#     
#     distance <- distance(spot1, spot2)
#     distanceOnRef <- distance(coSpot1, coSpot2)
#     movedSpot <- scale(movedSpot, coSpot1, distanceOnRef / distance)
#     movedSpot
# }

# spot is on plot, want to move it to refPlot. spot1 and spot are on plot while coSpot1
# anc coSpot2 are on ref plot
# reverseMoveWithCoPairs = function (spot, anchorSpots) {    
#     spot1 <- anchorSpots[[1]] 
#     refCoSpot1 <- anchorSpots[[2]]
#     spot2 <- anchorSpots[[3]]
#     refCoSpot2 <- anchorSpots[[4]]
#     
#     distance <- distance(spot1, spot2)
#     distanceOnRef <- distance(refCoSpot1, refCoSpot2)
#     reverseMovedSpot <- scale(spot, refCoSpot1, distance/ distanceOnRef)
#     
#     refAngle21 <- angleToXAxis(refCoSpot2, refCoSpot1)
#     angle21 <- angleToXAxis(spot2, spot1)
#     angleDiff21 <- refAngle21 - angle21
#     reverseMovedSpot <- rotate(reverseMovedSpot, refCoSpot1, -angleDiff21)
#     
#     deltaX <- refCoSpot1$Pos_X - spot1$Pos_X
#     deltaY <- refCoSpot1$Pos_Y - spot1$Pos_Y
#     reverseMovedSpot <- moveXY(reverseMovedSpot, -deltaX, -deltaY)
#     reverseMovedSpot
# }


# mapNeighbor2RelativeToNeighbor1 = function(neighbor2, neighbor1, anchorSpot, 
#                                            targetNeighbor1, targetAnchorSpot) {
#     angle2to1 = angle2to1(neighbor2, neighbor1, anchorSpot)
#     shortNeighbor2OnTarget <- rotate(targetNeighbor1, targetAnchorSpot, angle2to1)
#     neighbor2OnTarget <- scale(shortNeighbor2OnTarget, targetAnchorSpot, 
#                                distance(neighbor2, anchorSpot) / distance(targetNeighbor1, targetAnchorSpot))
#     neighbor2OnTarget
# }


# Get the closest refspots to spots, and calculate the average distance between
# closest refSpots and spots.
avgDistanceToClosestRefSpot = function (spots, refSpots) {
    closestRefSpots <- closestSpot (spots, refSpots)
    averageDistance <- avgDistance(spots, closestRefSpots)
    averageDistance
}


# get average distance between 2 sets of spots. The number of spots1 and spots2 are the same.
avgDistance = function(spots1, spots2) {
    distances <- distance(spots1, spots2)
    averageDistance <- sum(distances) / length(distances)
    averageDistance
}


# Debug, remove comments
# Get closestSpot of spots choosing from refSpots
# If spots is a vector, go through spots to find out the closest spot for each element of the vector
closestSpot = function(spots, refSpots, countOnlyUnusedSpots = FALSE) {
    closestSpotIndice <- closestSpotIndex(spots, refSpots, countOnlyUnusedSpots)  
    closestSpot <- refSpots[closestSpotIndice,]                 
#     if (countOnlyUnusedSpots) {
#         refSpots <- refSpots[refSpots$IsUsed == FALSE,]
#     }
#     closestSpot <- spotDataFrame(refSpots)
#     for (r in 1 : nrow(spots)) {
#         distances = distance(spots[r,], refSpots)
#         closestSpot <- rbind(closestSpot, refSpots[which.min(distances),]);
#     }
     closestSpot
}


# Get closestSpot of spots choosing from refSpots, exclude spot if spot is on refSpots
closestSpotExcludeSelf = function(spots, refSpots) {
    closestSpot <- spotDataFrame(refSpots)
    for (r in 1 : nrow(spots)) {
        distances = distance(spots[r,], refSpots)
        refSpotsWithoutASpot <- refSpots[distances != 0,]
        newDistances <- distances[distances != 0]
        closestSpot <- rbind(closestSpot, refSpotsWithoutASpot[which.min(newDistances),]);
    }
    closestSpot
}


# spot with lots of neighbors, but can't be composed of two spots, so RES should be lower
# than average RES, and AMP should be higher than average AMP
# it could return a null, return spots contain an extra column $ numNeighbors
# localDenseSpots = function(spots, minY = min(spots$Pos_Y),  maxY = max(spots$Pos_Y),
#                            minX = min(spots$Pos_X),  maxX = max(spots$Pos_X), neighbor_dist = 10000,
#                            dense_spot_num = 1, min_neighbor_num = 3) {
#     
#     realMinY <- minY
#     realMaxY <- maxY
#     realMinX <- minX
#     realMaxX <- maxX
#     restrictedSpots <- spots %>% filter(Pos_Y >= realMinY & Pos_Y <= realMaxY &
#                                             Pos_X >= realMinX & Pos_X <= realMaxX)
#     
#     # might not be necessary, debug, usually the FQ-detected spot near doublet has bigger RES, 
#     # try to avoid the doublet
#     #    averageRES <- mean(restrictedSpots$RES)
#     #    averageAMP <- mean(restrictedSpots$AMP)
#     #    currentAMP <- 0
#     #     densestSpot = restrictedSpots[1, ];
#     #     desnestSpotNeighbors <- neighbors_exclude_spot(densestSpot, restrictedSpots, neighbor_dist)
#     #     mostNeighborNum = dim(desnestSpotNeighbors)[1]
#     
#     denseSpots <- spotDataFrame()
#     for(r in 1 : nrow(restrictedSpots)) {
#         spot <- restrictedSpots[r, ];
#         spotNeighbors <- neighbors_exclude_spot(spot, restrictedSpots, neighbor_dist)
#         numNeighbors = nrow(spotNeighbors)
#         if (numNeighbors >= min_neighbor_num) {
#             #    closestNeighbor <- closestSpot(spot, spotNeighbors)
#             #    closestNeighborDistance <- distance(spot, closestNeighbor)
#             #    if ((numNeighbors > mostNeighborNum | (numNeighbors == mostNeighborNum & spot$AMP > 
#             #                            currentAMP)) & spot$AMP > averageAMP) {#& spot$RES < averageRES & spot$AMP > averageAMP & closestNeighborDistance <= max_neighbor_distance) {
#             #        densestSpot <- spot
#             #        currentAMP <- spot$AMP
#             #        mostNeighborNum <- numNeighbors
#             #    }
#             spot$numNeighbors <- numNeighbors
#             denseSpots <- rbind(denseSpots, spot)
#         }       
#     }
#     if (nrow(denseSpots) == 0) {
#         denseSpots
#     } else {
#         denseSpots <- denseSpots[rev(order(denseSpots$numNeighbors)),]
#         if (nrow(denseSpots) <= dense_spot_num) {
#             denseSpots
#         } else {
#             denseSpots[1:dense_spot_num, ]
#         }
#     } 
# }

# get multiple dense spots across the plot
# denseSpots = function(spots, numAcrossY, numAcrossX, neighbor_dist = 10, dense_spot_num = 1, 
#                       min_neighbor_num = 3, isToMicron = TRUE) {
#     dfDenseSpots <- spotDataFrame()
#     maxY <- max(spots$Pos_Y)
#     maxX <- max(spots$Pos_X)
#     minY <- min(spots$Pos_Y)
#     minX <- min(spots$Pos_X)
#     ySpan <- (maxY - minY) / numAcrossY
#     xSpan <- (maxX - minX) / numAcrossX
#     for (i in 1 : numAcrossY) {
#         for (j in 1 : numAcrossX) {
#             someDenseSpots <- localDenseSpots(spots, minY + (i-1) * ySpan, minY + i * ySpan, 
#                                               minX + (j-1) * xSpan, minX + j * xSpan, neighbor_dist, 
#                                               dense_spot_num, min_neighbor_num)
#             if (nrow(someDenseSpots) != 0) {
#                 #neighbors <- neighbors_exclude_spot(aDenseSpot, spots, neighbor_dist)
#                 #if (nrow(neighbors) >= 3) {
#                 someDenseSpots$yPortion <- i
#                 someDenseSpots$xPortion <- j
#                 #aDenseSpot$numNeighbors <- nrow(neighbors)
#                 dfDenseSpots <- rbind(dfDenseSpots, someDenseSpots)
#                 #}
#             }  
#         }
#     }
#     dfDenseSpots
# }

spotDataFrame = function(spot) {
    spot[0,]
    # data.frame(Pos_X = double(), Pos_Y = double(), AMP = double(), Color = character())
}

# isSublist = function (sublist, list, maxTolerableDiff, sublist2 = sublist, list2 = list, 
#                       maxTolerableDiff2 = maxTolerableDiff, minMatchedRatio) {
#     correctMatchedRatio = correctMatchedSublistRatio (sublist, list, maxTolerableDiff, 
#                                                       sublist2, list2, maxTolerableDiff2)
#     correctMatchedRatio >= minMatchedRatio
# }

# correctMatchedSublistRatio = function (sublist, list, maxTolerableDiff, sublist2 = sublist, list2 = list, 
#                                        maxTolerableDiff2 = maxTolerableDiff) {
#     count <- 0
#     for (i in 1: length(sublist)) {
#         if (isSubElement(sublist[i], list, maxTolerableDiff, sublist2[i], list2, 
#                          maxTolerableDiff2)) {         
#             count <- count + 1
#         }   
#     }
#     count / length(sublist)
# }

# isSubElement = function (element, list, maxTolerableDiff, element2 = element, list2 = list,  
#                          maxTolerableDiff2 = maxTolerableDiff) {
#     isElementInList = abs(list - element) <= maxTolerableDiff
#     isElement2InList2 = abs(list2 - element2) <= maxTolerableDiff2
#     areElemtntAnd2InLists = isElementInList & isElement2InList2
#     sum(areElemtntAnd2InLists) > 0
# }

# indiceOfSublist = function(sublist, list, maxTolerableDiff, sublist2 = sublist, 
#                            list2 = list, maxTolerableDiff2 = maxTolerableDiff) {
#     indice <- numeric()
#     for (i in 1: length(sublist)) {
#         extraIndice <- indiceOfSubElement(sublist[i], list, maxTolerableDiff,
#                                           sublist2[i], list2, maxTolerableDiff2)
#         indice <- c(indice, extraIndice)
#     }
#     unique(indice)
# }

# indiceOfSubElement = function(element, list, maxTolerableDiff, element2 = sublist, 
#                               list2 = list, maxTolerableDiff2 = maxTolerableDiff) {
#     isElementInList = abs(list - element) <= maxTolerableDiff
#     isElement2InList2 = abs(list2 - element2) <= maxTolerableDiff2
#     areElemtntAnd2InLists = isElementInList & isElement2InList2
#     which(areElemtntAnd2InLists)
# }

# 
# isSubNeighbors = function(neighborDists, neighborAngles, refSpot, refNeighbors, 
#                           minNeighborMatchRatio, maxTolerableDistDiff, maxTolerableAngleDiff) {
#     refNeighborDists <- distance(refSpot, refNeighbors)
#     refNeighborAngles<- as.vector(angleToXAxis(refNeighbors, refSpot))
#     count <- 0
#     for (iNeighbor in 1: length(neighborDists)) {
#         if (isSubElement(neighborDists[iNeighbor], refNeighborDists, maxTolerableDistDiff,
#                          neighborAngles[iNeighbor], refNeighborAngles, maxTolerableAngleDiff)) {
#             count <- count + 1
#         }   
#     }
#     count / length(neighborDists) >= minNeighborMatchRatio    
# }

# Get initial coSpotsOnRef, rank then in order of correctNeighborRatio 
# initialCoSpots = function(spot, spots, refSpots, neighbor_dist = 10000, 
#                           ini_cospot_dist = 20000, coSpot_dist = 500,
#                           angle_tolerance = pi/10, minCorrectNeighborRatio = 0.5, 
#                           min_correct_neighbor = 3) {
#     neighbors <- neighbors_exclude_spot(spot, spots, neighbor_dist)
#     neighborDists <- distance(spot, neighbors)
#     neighborAngles <- as.vector(angleToXAxis(neighbors, spot))
#     coSpotsOnRef <- neighborhood(spot, refSpots, ini_cospot_dist)
#     
#     # at least 2 out of all or 60% has to be correct
#     minCorrectNeighborRatio <- max(min_correct_neighbor/nrow(neighbors), minCorrectNeighborRatio)
#     coSpotsOnRefWithGoodNeighbors <- spotDataFrame()    
#     for (r in 1: nrow(coSpotsOnRef)) {
#         refSpot <- coSpotsOnRef[r,]
#         refNeighbors <- neighbors_exclude_spot(refSpot, refSpots, neighbor_dist)
#         refNeighborDists <- distance(refSpot, refNeighbors)
#         refNeighborAngles <- as.vector(angleToXAxis(refNeighbors, refSpot))
#         correctedMatchedNeighborRatio <- correctMatchedSublistRatio(neighborDists, 
#                                                                     refNeighborDists, coSpot_dist, neighborAngles, refNeighborAngles, angle_tolerance)
#         if (correctedMatchedNeighborRatio >= minCorrectNeighborRatio) {
#             refSpot$goodNeighborRatio <- correctedMatchedNeighborRatio
#             coSpotsOnRefWithGoodNeighbors <- rbind(coSpotsOnRefWithGoodNeighbors, refSpot)
#         }
#     }
#     if (nrow(coSpotsOnRefWithGoodNeighbors) == 0) {
#         coSpotsOnRefWithGoodNeighbors
#     } else {
#         coSpotsOnRefWithGoodNeighbors[rev(order(coSpotsOnRefWithGoodNeighbors$goodNeighborRatio)),]
#     } 
# }





# Add the coSpotOnRef to plot, with coSpotOnRef as a big red circle, refNeighbors as big 
# blue circles
# addRefCoSpotAndNeighborsToPlot = function (p, spot, spots, refSpot, refSpots, 
#                                            neighbor_dist = 10000, coSpot_dist = 500, angle_tolerance = pi/10) {
#     neighbors <- neighbors_exclude_spot(spot, spots, neighbor_dist) 
#     neighborDists <- distance(spot, neighbors)
#     neighborAngles <- as.vector(angleToXAxis(neighbors, spot))    
#     refNeighbors <- neighbors_exclude_spot(refSpot, refSpots, neighbor_dist) 
#     refNeighborDists <-  distance(refSpot, refNeighbors)
#     refNeighborAngles <- as.vector(angleToXAxis(refNeighbors, refSpot))
#     indiceOfGoodRefNeighbors <- indiceOfSublist(neighborDists, refNeighborDists, coSpot_dist,
#                                                 neighborAngles, refNeighborAngles, angle_tolerance)
#     goodNeighbors <- refNeighbors[indiceOfGoodRefNeighbors, ] 
#     x <- p + geom_point(data = refSpot, aes(x = Pos_X, y = Pos_Y), color = "red", size = I(5), shape = 21) +
#         geom_point(data = goodNeighbors, aes(x = Pos_X, y = Pos_Y), color = "blue", size = I(5), shape = 21) 
#     x
# }






# This might no be useful because it doesn't fit FSH
# for every refSpots, transform it to the spots on g1Plot, and get the normTiffAMP for the spot
# as the maxNormTiffAMP with radius = 3
# nearbyMaxTiffAmp = function(refSpots, anchorSpots, Pix_XY, tiffMatrix, radius = 3) {   
#     maxTiffAmpIndex <- maxTiffAmpXYIndexNearSpot(refSpots, anchorSpots, Pix_XY, tiffMatrix, radius)
#     minValueInMatrix <- min(min(tiffMatrix))
#     maxTiffAmp <- rep(minValueInMatrix, nrow(maxTiffAmpIndex))
#     maxTiffAmp[maxTiffAmpIndex[, 1] != 0] <- tiffMatrix[maxTiffAmpIndex[maxTiffAmpIndex[, 1] != 0,]]
#     maxTiffAmp
# }


# This might no be useful because it doesn't fit FSH
# nearbyMaxTiffAmpOverSurround = function(refSpots, anchorSpots, Pix_XY, tiffMatrix, radius = 3) {
#     maxTiffAmpIndex <- maxTiffAmpYXIndexNearSpot(refSpots, anchorSpots, Pix_XY, tiffMatrix, radius)
#     #localMinAmp <- 
#     
#     isIndexInbound <- maxTiffAmpIndex[, 1] != 0
#     deltaSurroundIndice <- list(c(-1,-1), c(-1,0), c(-1,1), c(0,-1), c(0,1), c(1,-1),
#                                 c(1,0), c(1,1))
#     sumSurroundAmp <- rep(0, nrow(maxTiffAmpIndex))
#     for (c in 1 : length(deltaSurroundIndice)) {
#         surroundIndex <- maxTiffAmpIndex + deltaSurroundIndice[[c]]
#         inboundSuroundIndex <- surroundIndex[isIndexInbound,]
#         
#         if (length(inboundSuroundIndex) == 2) {   # when there is only 1 index
#             inboundSuroundIndexY <- inboundSuroundIndex[1]
#             inboundSuroundIndexX <- inboundSuroundIndex[2]
#             inboundSuroundIndex <- cbind(inboundSuroundIndexY, inboundSuroundIndexX)
#         }
#         
#         sumSurroundAmp[isIndexInbound] <- sumSurroundAmp[isIndexInbound] + 
#             tiffMatrix[inboundSuroundIndex] 
#         #print(surroundIndex[isIndexInbound,])
#         #print(tiffMatrix[inboundSuroundIndex])
#     }
#     averageSurroundAmp <- sumSurroundAmp / length(deltaSurroundIndice)
#     
#     inboundMaxAmpIndex <- surroundIndex[isIndexInbound,]
#     if (length(inboundMaxAmpIndex) == 2) {   # when there is only 1 index
#         inboundMaxAmpIndexY <- inboundMaxAmpIndex[1]
#         inboundMaxAmpIndexX <- inboundMaxAmpIndex[2]
#         inboundMaxAmpIndex <- cbind(inboundMaxAmpIndexY, inboundMaxAmpIndexX)
#     }
#     maxTiffAmp <- rep(0, nrow(maxTiffAmpIndex))
#     maxTiffAmp[isIndexInbound] <- tiffMatrix[inboundMaxAmpIndex] 
#     
#     ampOverSurround <- rep(1, nrow(maxTiffAmpIndex))
#     ampOverSurround[isIndexInbound] <- maxTiffAmp[isIndexInbound] / averageSurroundAmp[isIndexInbound]
#     ampOverSurround
# }

# This might no be useful because it doesn't fit FSH
# transform the refspots to coSpots, and map them to pixel, and get the maxium indice
# of spots within radius.  Return (Xs and Ys).  If the spot radius is outside 
# the pixel region, x, y are both 0. 
# maxTiffAmpYXIndexNearSpot = function(refSpots, anchorSpots, Pix_XY, normTiffMatrix, radius) {   
#     
#     pixXY <- spotToXY(refSpots, pix_XY, anchorSpots)
#     pixY <- pixXY[,2]
#     pixX <- pixXY[,1]
#     trueRadius <- radius - 1
#     minX <- pixX - trueRadius
#     maxX <- pixX + trueRadius
#     minY <- pixY - trueRadius
#     maxY <- pixY + trueRadius
#     
#     maxAMPXYIndice <- maxTiffAmpYXIndiceWithinPixel(minY, maxY, minX, maxX, normTiffMatrix)
#     maxAMPXYIndice
# }



# Plot the raster plot, and circle the refSpot revesreMoved to the plot 
# (g1Plot, for example) as a white circle, while the (g1) neighbors in the plot 
# (g1Plot, for example) as blue circle
# plotNearbyTiffRaster = function (refSpot, spots, pix_XY, normTiffMatrix, 
#                                  anchorSpots = list(data.frame(Pos_X = 1, Pos_Y = 1), 
#                                                     data.frame(Pos_X = 1, Pos_Y = 1), 
#                                                     data.frame(Pos_X = 2, Pos_Y = 2), 
#                                                     data.frame(Pos_X = 2, Pos_Y = 2)),
#                                  plot_area_dist = 2000, 
#                                  color = colorRampPalette(c("black", "green")) (1000),
#                                  name = "") {   
#     
#     resX <- dim(normTiffMatrix)[2]
#     resY <- dim(normTiffMatrix)[1]
#     xy <- spotToXY(refSpot, pix_XY, anchorSpots)
#     pixX <- xy[1]
#     pixY <- xy[2]  
#     pixPlotArea <- round(plot_area_dist / pix_XY)
#     minPixX <- max(1, pixX - pixPlotArea)
#     maxPixX <- min(pixX + pixPlotArea, resX)
#     minPixY <- max(1, pixY - pixPlotArea)
#     maxPixY <- min(pixY + pixPlotArea, resY)
#     
#     nearbyMatrix <- normTiffMatrix[minPixY : maxPixY, minPixX : maxPixX]
#     newResX <- dim(nearbyMatrix)[2]
#     newResY <- dim(nearbyMatrix)[1]
#     rasterImage <- raster(nearbyMatrix, ymn = 0, ymx = newResY, xmn = 0, xmx = newResX) 
#     #plot(rasterImage, col = color, xlim=c(minPixX, maxPixX),ylim=c(minPixY, maxPixY), main = name)
#     plot(rasterImage, col = color, xlim=c(0, newResX),ylim=c(0, newResY), main = name)
#     pixXOnNewMatrix <- pixX - minPixX + 1 + 0.5
#     pixYOnNewMatrix <- newResY - (pixY - minPixY) - 1 - 0.5
#     points(pixXOnNewMatrix, pixYOnNewMatrix, col = "red", cex = 2)
#     
#     # plot neighbors spots on the tiff plot
#     refSpotOnPlot <- reverseMoveWithCoPairs(refSpot, anchorSpots)
#     neighborSpots <- neighborhood(refSpotOnPlot, spots[,1: ncol(refSpotOnPlot)], plot_area_dist)
#     neighborXY <- spotToXY(neighborSpots, pix_XY)
#     neighborPixX <- neighborXY[,1]
#     neighborPixY <- neighborXY[,2]
#     neighborPixXOnNewMatrix <- neighborPixX - minPixX + 1 + 0.5
#     neighborPixYOnNewMatrix <- newResY - (neighborPixY - minPixY) - 1 - 0.5
#     points(neighborPixXOnNewMatrix, neighborPixYOnNewMatrix, col = "white", cex = 2)
# }



# from refSpot, to get coSpot in a slice, and then get its pixX and pixY as cbind()
# Note that the spots from FQ has y-axis from top to bottom while the pixel from readTIFF
# has y-axis from bottom to top, I have to change the y-axis as maxY - yaxis
# spotToXY = function (refSpot, pix_XY, anchorSpots = list(data.frame(Pos_X = 1, Pos_Y = 1), 
#                                                          data.frame(Pos_X = 1, Pos_Y = 1), 
#                                                          data.frame(Pos_X = 2, Pos_Y = 2), 
#                                                          data.frame(Pos_X = 2, Pos_Y = 2))) {
#     coSpot <- reverseMoveWithCoPairs(refSpot, anchorSpots)
#     pixPlotArea <- round(plot_area_dist / pix_XY)
#     pixX <- round(coSpot$Pos_X / pix_XY)
#     pixY <- round(coSpot$Pos_Y / pix_XY)
#     cbind(pixX, pixY)
# }

# Remove refSpots that when reverse moved to any single Tiff image (such as g1) will be outside of the 
# Tiff image.
# removeOutSpots <- function(refSpots, pix_XY, allAnchorSpots, bufferPixel = 1) {
#     isInsideRange <- rep(TRUE, nrow(refSpots))
#     for (i in 1 : length(allAnchorSpots)) {
#         anchorSpotsForOneImage <- allAnchorSpots[[i]]
#         xyPixelOnOneImage <- spotToXY(refSpots, pix_XY, anchorSpotsForOneImage)
#         # note I also remove spots that's exactly at the edge in case FQ coudn't find them.
#         isOutsideRangeOfOneImage <- xyPixelOnOneImage[,1] < 1 + bufferPixel | 
#             xyPixelOnOneImage[,1] > resX - bufferPixel | 
#             xyPixelOnOneImage[,2] < 1 + bufferPixel | 
#             xyPixelOnOneImage[,2] > resY - bufferPixel
#         isInsideRange <- isInsideRange & !isOutsideRangeOfOneImage
#     }
#     refSpots[isInsideRange,]
# }








