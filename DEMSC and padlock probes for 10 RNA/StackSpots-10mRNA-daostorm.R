rm(list = ls())
setwd("/Users/jackie/Downloads/SGTC microscope/2016-08-11-10mRNAs vs 4 exon mRNAs-37C-45C Amp 3rd time/10 mRNAs/registered/");
source("/Users/jackie/Downloads/SGTC microscope/2016-08-11-10mRNAs vs 4 exon mRNAs-37C-45C Amp 3rd time/4 exon mRNAs/spotPlotFunctions.R")
source("/Users/jackie/Downloads/SGTC microscope/2016-08-11-10mRNAs vs 4 exon mRNAs-37C-45C Amp 3rd time/4 exon mRNAs/spotRegisterFunctions.R")
source("/Users/jackie/Downloads/SGTC microscope/2016-08-11-10mRNAs vs 4 exon mRNAs-37C-45C Amp 3rd time/4 exon mRNAs/spotStackFunctions.R")
source("/Users/jackie/Downloads/SGTC microscope/2016-08-11-10mRNAs vs 4 exon mRNAs-37C-45C Amp 3rd time/4 exon mRNAs/dao_spotFunctions.R")

# 1. Plot the tiffPlots, note that the plot has lower left as (0, 0)
fluorNames <- c("g1", "y1", "o1", "r1", "g2", "y2", "o2", "r2", "g3", "y3", "o3", "r3", "g4", "y4", "o4", "r4")

preRefFilename <- "H1 10mRNA vigorous DI wash old G25-G4 Gfilter 2s registered"
preFluorFilenames <- list(g1 = "I1 10mRNA F1-NF234 Gfilter 2s registered",
                          y1 = "I1 10mRNA F1-NF234 Yfilter 1s registered",
                          o1 = "I1 10mRNA F1-NF234 Ofilter 2s registered",
                          r1 = "I1 10mRNA F1-NF234 Rfilter 2s registered",
                          g2 = "J1 10mRNA F2-NF134 Gfilter 2s registered",
                          y2 = "J1 10mRNA F2-NF134 Yfilter 1s registered",
                          o2 = "J1 10mRNA F2-NF134 Ofilter 2s registered",
                          r2 = "J1 10mRNA F2-NF134 Rfilter 2s registered",
                          g3 = "L1 10mRNA F3-NF124 45C Gfilter 2s registered",
                          y3 = "L1 10mRNA F3-NF124 45C Yfilter 1s registered",
                          o3 = "L1 10mRNA F3-NF124 45C Ofilter 2s registered",
                          r3 = "L1 10mRNA F3-NF124 45C Rfilter 2s registered",
                          g4 = "M1 10mRNA F4-NF123 Gfilter 2s registered",
                          y4 = "M1 10mRNA F4-NF123 Yfilter 1s registered",
                          o4 = "M1 10mRNA F4-NF123 Ofilter 2s registered",
                          r4 = "M1 10mRNA F4-NF123 Rfilter 2s registered")

refTiffFilename <- paste(preRefFilename, ".tif", sep = "")
fluorTiffFilenames <- list()
for (fluorName in fluorNames) {
    fluorTiffFilenames[[fluorName]] <- paste(preFluorFilenames[[fluorName]], ".tif", sep = "")
}

# 1.1 Read the tif image's intensity values into matrix
refTiffMatrix <- readTIFF(refTiffFilename, as.is = TRUE)
fluorTiffMatrices <- list()
for (fluorName in fluorNames) {
    fluorTiffMatrices[[fluorName]] <- readTIFF(fluorTiffFilenames[[fluorName]], as.is = TRUE)
}

# 1.2 Scale the tiffMatrix to be between 0 and 1 with contrast
refScaledMatrix <- scaled_matrix(refTiffMatrix, 5)
fluorScaledMatrices <- list()
scaledFactors <- c(5, 3, 10, 10, 5, 3, 10, 10, 5, 3, 10, 10, 5, 3, 10, 10)
for (i in seq_along(fluorNames)) {
    fluorName <- fluorNames[i]
    scaledFactor <- scaledFactors[i]
    fluorScaledMatrices[[fluorName]] <- scaled_matrix(fluorTiffMatrices[[fluorName]], scaledFactor)
}

# 1.3 Set up the color palette to draw the tiff image
greenPalette = colorRampPalette(c("black", "green")) (1000)
yellowPalette = colorRampPalette(c("black", "gold")) (1000)
orangePalette = colorRampPalette(c("black", "orange")) (1000)
redPalette = colorRampPalette(c("black", "red")) (1000)
palettes <- list()
for (fluorName in fluorNames) {
    palettes[[fluorName]] <- switch(substr(fluorName,1,1), "g" = greenPalette, "y" = yellowPalette, "o" = orangePalette, "r" = redPalette)
}

# 1.4 Get pix_XY size = 0.3225 um, resY = height of image, resX = width of image
refOldSpotsFilename <- "../spots/H1 10mRNA vigorous DI wash old G25-G4 Gfilter 2s stacked__spots.txt"
pix_XY <- read_pixel_xy(refOldSpotsFilename)
resY = dim(refScaledMatrix)[1]
resX = dim(refScaledMatrix)[2]

# 1.5 Get the cell outline polygons
outlineFilename <- "../outline/H1 10mRNA vigorous DI wash old G25-G4 Gfilter 2s stacked__outline.txt"
polygons <- readOutlines(outlineFilename, 13)
polygons_revY <- polygons_reverseY (polygons, resY)
polygons_um <- polygons_toMicron (polygons, pix_XY)
polygonMatrix <- getPolygonMatrix(polygons, resX, resY)

# 1.6 Get the cell nuclei polygons
nucleiOutlineFilename <- "../outline/G1 10mRNA new F2-NF134 Gfilter 2s outline of nuclei.txt"
nucleiPolygons <- readOutlines(nucleiOutlineFilename, 13)
nucleiPolygons_revY <- polygons_reverseY (nucleiPolygons, resY)
nucleiPolygons_um <- polygons_toMicron (nucleiPolygons, pix_XY)

# 1.7 Plot the tiff images
isToDrawOutlines <- TRUE
withOutlineNumber <- TRUE
par(mfrow=c(1,1))
plotMatrix(refScaledMatrix, greenPalette, "ref plot", polygons_revY, isToDrawOutlines, withOutlineNumber, nucleiPolygons_revY)

par(mfrow=c(2,2))
for (fluorName in fluorNames) {
    plotMatrix(fluorScaledMatrices[[fluorName]], palettes[[fluorName]], paste(fluorName, "plot"), polygons_revY, 
               isToDrawOutlines, withOutlineNumber)
}

# 2. Read all the spots, Pos_X and Pos_Y are in um
refSpotsFilename <- paste("fits images/", gsub(" ", "_", preRefFilename), ".pos.out", sep = "")
refSpots <- read_daostorm_spots(refSpotsFilename, "gray", refTiffFilename)

fluorSpotsFilenames <- list()
for (fluorName in fluorNames) {
    fluorSpotsFilenames[[fluorName]] <- paste("fits images/", gsub(" ", "_", preFluorFilenames[[fluorName]]), ".pos.out", sep = "")
}

fluorColors <- list()
fluorSpots <- list()
for (fluorName in fluorNames) {
    fluorColors[[fluorName]] <- switch(substr(fluorName,1,1), "g" = "green", "y" = "yellow", "o" = "orange", "r" = "red")
    fluorSpots[[fluorName]] <- read_daostorm_spots(fluorSpotsFilenames[[fluorName]], 
                                                   fluorColors[[fluorName]], fluorTiffFilenames[[fluorName]])
}

# 2.1 Plot the spots against refSpot, the unit is in um
spotSize <- 1.5
fluorPlots <- list()
for (fluorName in fluorNames) {
    fluorPlots[[fluorName]] <- plotSpots(fluorSpots[[fluorName]], refSpots, spotSize = spotSize, 
                                         polygons_um = polygons_um, plotname = paste(fluorName, "Plot"))
    print(fluorPlots[[fluorName]])
}

# 2.2 Choose the refSpots that's within cell. Therefore, the outside dirt spots could be removed
refSpots <- enclosedSpotsWithinBoundaries(refSpots, polygons_um)


# 3. For every cell outline (polygon), regid translate the spots to the refSpots.
maxPreCospotDist <- 1.5

transFluorSpotsInfo <- list()
transFluorSpots <- list()
transFluorXYPerPolygon <- list()
for (fluorName in fluorNames) {
    print(fluorName)
    transFluorSpotsInfo[[fluorName]] <- translateAllCells(fluorSpots[[fluorName]], refSpots, 
                                                          polygons_um, maxPreCospotDist)
    transFluorSpots[[fluorName]] <- transFluorSpotsInfo[[fluorName]][[1]]
    transFluorXYPerPolygon[[fluorName]] <- transFluorSpotsInfo[[fluorName]][[2]]     
}


# 3.1 Plot the translated spots on refSpots
transFluorPlots <- list()
for (fluorName in fluorNames) {
    transFluorPlots[[fluorName]] <- plotSpots(transFluorSpots[[fluorName]], refSpots, 
                                              spotSize = spotSize, polygons_um = polygons_um, plotname = "transG1Plot")   
    print(transFluorPlots[[fluorName]])
}


# 4. Combine f1Spots, f2Spots, f3Spots, f4Spots and plot
mergedNames <- c("f1", "f2", "f3", "f4")
mergedSpots <- list()
mergedSpots[["f1"]] <- rbind(fluorSpots[["g1"]], fluorSpots[["y1"]], fluorSpots[["o1"]], fluorSpots[["r1"]])
mergedSpots[["f2"]] <- rbind(fluorSpots[["g2"]], fluorSpots[["y2"]], fluorSpots[["o2"]], fluorSpots[["r2"]])
mergedSpots[["f3"]] <- rbind(fluorSpots[["g3"]], fluorSpots[["y3"]], fluorSpots[["o3"]], fluorSpots[["r3"]])
mergedSpots[["f4"]] <- rbind(fluorSpots[["g4"]], fluorSpots[["y4"]], fluorSpots[["o4"]], fluorSpots[["r4"]])

mergedPlots <- list()
for (mergedName in mergedNames) {
    mergedPlots[[mergedName]] <- plotSpots(mergedSpots[[mergedName]], refSpots, spotSize = 2, 
                                           polygons_um = polygons_um, plotname = paste(mergedName, "Plot"))
    print(mergedPlots[[mergedName]])
}

# 4.1 Combine transF1Spots, transF2Spots, transF3Spots, transF4Spots and plot
transMergedSpots <- list()
transMergedSpots[["f1"]] <- rbind(transFluorSpots[["g1"]], transFluorSpots[["y1"]], transFluorSpots[["o1"]], transFluorSpots[["r1"]])
transMergedSpots[["f2"]] <- rbind(transFluorSpots[["g2"]], transFluorSpots[["y2"]], transFluorSpots[["o2"]], transFluorSpots[["r2"]])
transMergedSpots[["f3"]] <- rbind(transFluorSpots[["g3"]], transFluorSpots[["y3"]], transFluorSpots[["o3"]], transFluorSpots[["r3"]])
transMergedSpots[["f4"]] <- rbind(transFluorSpots[["g4"]], transFluorSpots[["y4"]], transFluorSpots[["o4"]], transFluorSpots[["r4"]])

transMergedPlots <- list()
for (mergedName in mergedNames) {
    transMergedPlots[[mergedName]] <- plotSpots(transMergedSpots[[mergedName]], refSpots, spotSize = 2, 
                                           polygons_um = polygons_um, plotname = paste("trans", mergedName, "Plot"))
    print(transMergedPlots[[mergedName]])
}

# 4.2 Check number of spots in each slice. Map each spot from each slice to refSpots to find 
# refCoSpots, and record the coSpots of each refSpot. They are 1520, 1792, 1925, 1763, 1715
nrow(refSpots)
for (mergedName in mergedNames) {
    print(nrow(transMergedSpots[[mergedName]]))
}


# 5. Get max intensity of each color at each time point, in order to do normalization in later step
# Choose spots that maps to a refence spot, so it's less likely to be a dirt.
# f1MaxMinusBgTiffAmpPerColor = 883  127  575 -Inf
# f2MaxMinusBgTiffAmpPerColor = 767 1977 -Inf -Inf
# f3MaxMinusBgTiffAmpPerColor = 642 -Inf -Inf  323
# f4MaxMinusBgTiffAmpPerColor = -Inf 1443 -Inf  339
# Use weighted avarage instead
# f1MaxMinusBgTiffAmpPerColor = 883, (883*1977/767 = 2275.999), 575, (883*323/642 = 444.2508)
# f2MaxMinusBgTiffAmpPerColor = 767, 1977, (767*575/883 = 499.4621), (767*323/642 + 1977*339/1443)/2 = 425.1703
# f3MaxMinusBgTiffAmpPerColor = 642, (642*1977/767 + 323*1443/339)/2 = 1514.848, (642*575/883 = 418.0634), 323
# f4MaxMinusBgTiffAmpPerColor = (1443*767/1977 + 339*642/323)/2 = 616.8152, 1443, 575, 339
coSpotDist <- 0.5

for (mergedName in mergedNames) {
    refSpots[[paste("transMergedNumCoSpots", mergedName, sep = "")]] <- num_neighbors_circle(
        refSpots, transMergedSpots[[mergedName]], coSpotDist)
    refSpots[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]] <- index_neighbors_circle(
        refSpots, transMergedSpots[[mergedName]], coSpotDist)
    
    # Remove the transMergedNumCoSpots and transMergedCoSpotsIndice of the refSpots if 
    # the coSpot is not a coSpot of the refSpot
    allNeighborMergedIndice <- unlist(refSpots[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]])
    duplicatedMergedIndice <- allNeighborMergedIndice[duplicated(allNeighborMergedIndice)]
    for(duplicatedMergedIndex in duplicatedMergedIndice) {
        duplicatedRefIndice <- which(sapply(refSpots[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]], 
                                            function(x) {is.element(duplicatedMergedIndex,x)}))
        distances <- spot_distance(refSpots[duplicatedRefIndice,],transMergedSpots[[mergedName]][duplicatedMergedIndex,])
        numberOfCoSpots <-  refSpots[duplicatedRefIndice,][[paste("transMergedNumCoSpots", mergedName, sep = "")]]
        distances[numberOfCoSpots != 1] <- distances[numberOfCoSpots != 1] + 10 # Random large number to set the distance larget
        refSpots[duplicatedRefIndice,][-which.min(distances),][[paste("transMergedCoSpotsIndice", mergedName, sep = "")]] <-
            lapply(refSpots[duplicatedRefIndice,][-which.min(distances),][[paste("transMergedCoSpotsIndice", mergedName, sep = "")]],
                   function(x) {x[x != duplicatedMergedIndex]})
        refSpots[duplicatedRefIndice,][-which.min(distances),][[paste("transMergedNumCoSpots", mergedName, sep = "")]] <-
            sapply(refSpots[duplicatedRefIndice,][-which.min(distances),][[paste("transMergedCoSpotsIndice", mergedName, sep = "")]],
                   length) 
    }
}

closestTransMergedCospots <- list()
for (mergedName in mergedNames) {
    closestTransMergedCospots[[mergedName]] <- dao_getClosestTransCospots(
        refSpots, transMergedSpots, mergedName)
}

refMaxMinusBgTiffAmp <- max(cleanRefSpots$minusBgTiffAmp)

fluorMaxMinusBgTiffAmp <- list()
for (fluorName in fluorNames) {
    mergedName <- paste("f", substr(fluorName, 2,2), sep = "")
    fluorMaxMinusBgTiffAmp[[fluorName]] <- max(closestTransMergedCospots[[mergedName]][
        closestTransMergedCospots[[mergedName]]$Color == fluorColors[[fluorName]],]$minusBgTiffAmp)
}

colorNames <- c("green", "yellow", "orange", "red")
mergedMaxMinusBgTiffAmpPerColor <- list()
for (mergedName in mergedNames) {
    for (colorName in colorNames) {
        fluorName <- paste(substr(colorName, 1, 1), substr(mergedName, 2, 2), sep = "")
        mergedMaxMinusBgTiffAmpPerColor[[mergedName]] <- c(mergedMaxMinusBgTiffAmpPerColor[[mergedName]],
                                                           fluorMaxMinusBgTiffAmp[[fluorName]]) 
    }   
}
mergedMaxMinusBgTiffAmpPerColor

# 5.1 set the minimum threshold of maximum intensity per color (green, yellow, orange, red)
# yellow color should be brightest while red color should be dimmest.
# If the maxItensityOfThatColor in that time is below the mins threshold, set it to be the 
# average of other time.
minMaxIntensityPerColor <- c(400, 1000, 400, 250)
countMaxMinusBgTiffAmpPerColor <- numeric(4)
sumMaxMinusBgTiffAmpPerColor <- numeric(4)
for(mergedName in mergedNames) {
    mergedMaxMinusBgTiffAmpPerColor[[mergedName]][
        mergedMaxMinusBgTiffAmpPerColor[[mergedName]] < minMaxIntensityPerColor] = 0
    countMaxMinusBgTiffAmpPerColor <- countMaxMinusBgTiffAmpPerColor + 
        (mergedMaxMinusBgTiffAmpPerColor[[mergedName]] > minMaxIntensityPerColor)
    sumMaxMinusBgTiffAmpPerColor <- sumMaxMinusBgTiffAmpPerColor + mergedMaxMinusBgTiffAmpPerColor[[mergedName]]
}
averageMaxMinusBgTiffAmpPerColor <- sumMaxMinusBgTiffAmpPerColor / countMaxMinusBgTiffAmpPerColor

# This averages over the same color over different time point.
# I think it is good enough to esitimate the color brightness for no fluor
for (mergedName in mergedNames) {
    tooDimColorIndex <- mergedMaxMinusBgTiffAmpPerColor[[mergedName]] < minMaxIntensityPerColor
    mergedMaxMinusBgTiffAmpPerColor[[mergedName]][tooDimColorIndex] <-
        averageMaxMinusBgTiffAmpPerColor[tooDimColorIndex] 
}
mergedMaxMinusBgTiffAmpPerColor

# Change to mergedMaxMinusBgTiffAmp: from list(FN) of vector to list of list(color) for future convenience
mergedMaxMinusBgTiffAmp <- list()
for (mergedName in mergedNames) {
    mergedMaxMinusBgTiffAmp[[mergedName]] <- list()
    for (colorName in colorNames) {
        mergedMaxMinusBgTiffAmp[[mergedName]][[colorName]] <-
            mergedMaxMinusBgTiffAmpPerColor[[mergedName]][colorNames == colorName]
    }    
}

# 5.2 Add the ScaledAmp to transMergedSpots
for (mergedName in mergedNames) {
    transMergedSpots[[mergedName]]$ScaledAmp <- 
        transMergedSpots[[mergedName]]$AMP / 
        unlist(mergedMaxMinusBgTiffAmp[[mergedName]][transMergedSpots[[mergedName]]$Color])
}


# 6. Separate refSpots into spots with close spots within 0.5 um and spots without close spots
# All 1520 cleanSpots have no closeNeighbors and is within cell boundaries
neighborDist <- 0.5
numNeighborRefSpotsExcludingSelf <- num_neighbors_exclude_spot_circle(refSpots, 
                                                                      refSpots, neighborDist)
neighborCountToFrequency <- table(numNeighborRefSpotsExcludingSelf)
neighborCountToFrequency

# 6.1 Check the frequency of transFNNumCoSpots for refSpots
# f1 (0, 1, 2, 3) = (193, 1221, 112, 6)
# f2 (0, 1, 2, 3) = (98, 1325, 105, 4)
# f3 (0, 1, 2, 3) = (236, 1204, 88, 4)
# f4 (0, 1, 2, 3) = (186, 1295, 48, 3)

# After remove dupicated transMergedSpots indice from refSpots
# f1 (0, 1, 2, 3) = (190, 1215, 110, 5)
# f2 (0, 1, 2, 3) = (95, 1321, 100, 4)
# f3 (0, 1, 2, 3) = (236, 1195, 85, 4)
# f4 (0, 1, 2, 3) = (183, 1289, 45, 3)
for (mergedName in mergedNames) {
    print(table(refSpots[paste("transMergedNumCoSpots",mergedName, sep = "")]))
}


# 6.3 for all refSpots$transMergedNumCoSpotsFN that is 0, use its neighbor usedSpots to translate the transMergedSpots
# Hmmm seem to not help at all. 
# for (mergedName in mergedNames) {
#     fnZeroCospotIndice <- setdiff(c(1:nrow(transMergedSpots[[mergedName]])), 
#                                   unlist(refSpots[[paste(transMergedCoSpotsIndice,mergedName, sep = "")]]))
#     
#     fnSpotsIsUsed <- transMergedSpots[[mergedName]] %>% filter(isUsed)    
#     for (fnZeroCospotIndex in fnZeroCospotIndice) {
#         fnWithZeroCospot <- transMergedSpots[[mergedName]][fnZeroCospotIndex,]
#         fnSpotsIsUsedSameColor <- fnSpotsIsUsed %>% filter(Color == fnWithZeroCospot$Color)
#         neighborUsedSameColorSpots <- neighbors_circle(fnWithZeroCospot, fnSpotsIsUsedSameColor, 3)
#         neighborRefCospots <- refSpots[neighborUsedSameColorSpots$CoRefSpotIndex,]
#         averageDx <- avg_dx(neighborRefCospots, neighborUsedSameColorSpots)
#         averageDy <- avg_dy(neighborRefCospots, neighborUsedSameColorSpots)
#         if (!is.nan(averageDx)) {
#             fnWithZeroCospot$Pos_X <- fnWithZeroCospot$Pos_X + averageDx
#         }
#         if (!is.nan(averageDy)) {
#             fnWithZeroCospot$Pos_Y <- fnWithZeroCospot$Pos_Y + averageDy
#         }
#         transMergedSpots[[mergedName]][fnZeroCospotIndex,] <- fnWithZeroCospot
#     }
# }


# 7. Prepare to stack for the refSpots that has no closeNeighbors, in this case, all refSpots
# Set up the lambda and diffractionLimit for all 4 colors
lonelyRefSpots <- refSpots#[numNeighborRefSpotsExcludingSelf == 0,]
tempTransMergedSpots <- transMergedSpots

numericAperture <- 0.4
diffractionFactor <- 0.61

lambda_um <- c(0.519, 0.573, 0.617, 0.664)
lambdas <- list()
diffractionLimits <- list() # (0.791, 0.874, 0.941, 1.013)
for (colorName in colorNames) {
    lambdas[[colorName]] <- lambda_um[which(colorNames == colorName)]
    diffractionLimits[[colorName]] <- diffractionFactor * lambdas[[colorName]] / numericAperture
}

mRNAToColor <- data.frame(name = c("ACTB", "ACTG", "RN7SK", "FTH1", "RPLP0", "EEF2", 
                                   "GAPDH", "FN1", "COL1A2", "TGFBI"),
                          allColors = c("green green red yellow", "green yellow green red",
                                        "green green green yellow", "orange yellow red red",
                                        "orange green red red", "orange green green yellow",
                                        "orange yellow red yellow", "orange yellow green yellow", 
                                        "orange green red yellow", "orange yellow green red"),
                          plotColor = c("red", "orange", "yellow", "darkgreen", "blue", "purple", 
                                        "pink", "green", "skyblue", "mediumpurple1"))

par(mfrow=c(1,1))
interesgedRefSpots <- rbind.fill(fourOffspringRefSpots, 
                                 threeOffspringRefSpots, 
                                 twoOffspringRefSpots,
                                 singletCrowd2RefSpotsSpot,
                                 closeDoubletCrowd2RefSpots,
                                 farDoubletCrowdRefSpots,
                                 verycrowdRefSpots)
interesgedRefSpots <- merge(interesgedRefSpots, mRNAToColor, by = "allColors", all = TRUE)
table(interesgedRefSpots$name) 

plotNearbyTiffRasterWithMRNAs(spot, interesgedRefSpots, 
                              mRNAToColor, pix_XY, refScaledMatrix, 2000, 0.3, 
                              "10 mRNAs on reference image", maxNumY = 1)


# 8. Stack for coSpots are either 0 or 1
# Conclusion: choose the number of coSpots as 2, 3, 4, but not 0, 1

# 8.1 Mark the FNSpots that has a 1-on-1 or more coSpot on cleanRefSpots as IsUsed
# Note that it's possible that one transMergedSpot is used for multiple refSpots....
transMergedSpots <- dao_updateIsUsedForTransMergedSpots(transMergedSpots, refSpots)

# Gradually increase secondaryCospotDist to add unused coSpot for refSpots that had zero coSpots
for (secondaryCospotDist in 0.1*c(6:11)) {
    refSpots <- dao_fillInOneCospotToXCospotsWithUnusedSpots(refSpots, 0, secondaryCospotDist, refSpots)
    transMergedSpots <- dao_updateIsUsedForTransMergedSpots(transMergedSpots, refSpots)
}



# 8.1 Stack a refSpots with 1 coSpot in each wash using the daostorm spots
# Make the cospot distance gradually expand to 1.1 um, so total fourOffspringRefSpots: 1004
#[1] "      green green red yellow,   249,       ACTB"
#[1] "   orange green green yellow,     7,       EEF2"
#[1] "     orange green red yellow,    40,     COL1A2"
#[1] "       orange yellow red red,    34,       FTH1"
#[1] "  orange yellow green yellow,    67,        FN1"
#[1] "      green yellow green red,    39,       ACTG"
#[1] "    orange yellow red yellow,    71,      GAPDH"
#[1] "    green green green yellow,   481,      RN7SK"
#[1] "     orange yellow green red,    15,      TGFBI"
#[1] "        orange green red red,     1,      RPLP0"
fourOffspringRefSpots <- refSpots %>% filter(transMergedNumCoSpotsf1 == 1 
                                             & transMergedNumCoSpotsf2 == 1
                                             & transMergedNumCoSpotsf3 == 1 
                                             & transMergedNumCoSpotsf4 == 1)
nrow(fourOffspringRefSpots)

fourOffspringRefSpots$allColors <- dao_assignAllColorsFromFQOrGaussianFit(
    fourOffspringRefSpots, coSpot_dist)
printCountOfMRNA(fourOffspringRefSpots$allColors, mRNAToColor)


which(fourOffspringRefSpots$allColors == "green yellow green yellow")
spot <- fourOffspringRefSpots[815,]
plotAllNearby(spot, 5)
plotNearbyTiffRaster(refSpot = spot, fourOffspringLonelyRefSpots[c(6, 9, 31, 49),], pix_XY, refScaledMatrix, isToTranslateRefSpots = FALSE, plot_area_dist = 2000)
dao_assignAllColorsFromFQOrGaussianFit(spot, coSpot_dist, 1, isToPrintFit = TRUE)

# Plot the nearby refMatrixPlot along with 10 mRNA spos
plotNearbyTiffRaster(spot, a, pix_XY, refTiffMatrix, isToTranslateRefSpots = FALSE, 
                     plot_area_dist = 50, color = colorRampPalette(c("black", "white")) (1000))


# 8.2 Check a lonely RefSpot with 1 coSpot in 3 washes and 0 cospot in rest 1 washes using the ds spots
# Make the cospot distance gradually expand to 1.1 um, so total threeOffspringRefSpots: 115
#[1] "    orange yellow red yellow,    28,      GAPDH"
#[1] "      green green red yellow,    27,       ACTB"
#[1] "    green green green yellow,    31,      RN7SK"
#[1] "       orange yellow red red,     7,       FTH1"
#[1] "     green green green green,     1, "                => dirt
#[1] "  orange yellow green yellow,     7,        FN1"
#[1] "     orange green red yellow,     7,     COL1A2"
#[1] "       green green green red,     2, "                => GGGY, Y4 not translated well in nuclei
#[1] "      green yellow green red,     2,       ACTG"
#[1] "     orange yellow green red,     1,      TGFBI"
#[1] "   orange green green yellow,     2,       EEF2"
# spot 7 is "green green green green", but g4 shouldn't have any color, so it should be a dirt.
# spot 38 & 68 are "green green green red", but they should be GGGY by eye. Y4 spots at nuclei are not translated well.
threeOffspringRefSpots <- refSpots %>% filter ((transMergedNumCoSpotsf1 <= 1 
                                                 & transMergedNumCoSpotsf2 <= 1 
                                                 & transMergedNumCoSpotsf3 <= 1 
                                                 & transMergedNumCoSpotsf4 <= 1)
                                                & (transMergedNumCoSpotsf1
                                                   + transMergedNumCoSpotsf2
                                                   + transMergedNumCoSpotsf3
                                                   + transMergedNumCoSpotsf4 == 3))
nrow(threeOffspringRefSpots)

threeOffspringRefSpots$allColors <- dao_assignAllColorsFromFQOrGaussianFit(
    threeOffspringRefSpots, coSpot_dist = 0.5)
printCountOfMRNA(threeOffspringRefSpots$allColors, mRNAToColor)
threeOffspringRefSpots[38,]$allColors <- "green green green yellow"
threeOffspringRefSpots[68,]$allColors <- "green green green yellow"

which(threeOffspringRefSpots$allColors == "green green green red")
spot <- threeOffspringRefSpots[38,]; spot; dao_plotAllNearby(spot)
color <- dao_chooseAColorFromGaussianFit(spot, "f4", coSpot_dist, isToPrintFit = TRUE)



# 8.3 Check a lonely RefSpot with 1 coSpot in 2 washes and 0 cospot in rest 2 washes using the ds spots
# Make the cospot distance gradually expand to 1.1 um, so total twoOffspringRefSpots: 71
#[1] "    orange yellow red yellow,    35,      GAPDH"
#[1] "    green green green yellow,    10,      RN7SK"
#[1] "      green green red yellow,    10,       ACTB"
#[1] "       orange yellow red red,     5,       FTH1"
#[1] "     orange green red yellow,     3,     COL1A2"
#[1] "        orange green red red,     2,      RPLP0"
#[1] "     green yellow red yellow,     1, "             => mix of two dim spots
#[1] "       green green green red,     1, "             => GGGY, Y4 is not traslated well near nuclei
#[1] "      green yellow green red,     2,       ACTG"
#[1] "   orange green green yellow,     1,       EEF2"
#[1] "      red green green yellow,     1, "             => OGGY
# spot 25 is "green yellow red yellow", but it should be a mix of two dim spots, can't tell the exact colors of the two spots 
# spot 35 is "green green green red", but should be GGGY, but y4 is not registered well and a bright R4 spot is nearby
# spot 63 is "red green green yellow", but should be OGGY of EEF2, but the o1 is dim doublets,
        # so the Gaussian fit of O1 is too far from the transRefSpot
# Should I consider disregarding any empty refSpots that is too dim or close to bright spots?
twoOffspringRefSpots <- refSpots%>% filter((transMergedNumCoSpotsf1 <= 1 
                                            & transMergedNumCoSpotsf2 <= 1 
                                            & transMergedNumCoSpotsf3 <= 1 
                                            & transMergedNumCoSpotsf4 <= 1) 
                                           & (transMergedNumCoSpotsf1 
                                              + transMergedNumCoSpotsf2 
                                              + transMergedNumCoSpotsf3 
                                              + transMergedNumCoSpotsf4 == 2))
nrow(twoOffspringRefSpots)

twoOffspringRefSpots$allColors <- dao_assignAllColorsFromFQOrGaussianFit(
    twoOffspringRefSpots, coSpot_dist = 0.6)
printCountOfMRNA(twoOffspringRefSpots$allColors, mRNAToColor)
twoOffspringRefSpots[35,]$allColors <- "green green green yellow"
twoOffspringRefSpots[63,]$allColors <- "orange green green yellow"

which(twoOffspringRefSpots$allColors == "green red green yellow")
ispot <- twoOffspringRefSpots[25,];ispot;dao_plotAllNearby(ispot, 5)
color <- dao_chooseAColorFromGaussianFit(ispot, "f2", coSpotDist, isToPrintFit = TRUE)



# 8.4 Check a lonely RefSpot with 1 coSpot in a wash and 0 cospot in rest 3 washes using the FQ spots
# Make the cospot distance gradually expand to 1.1 um, so total oneOffspringRefSpots: 46
#[1] "  orange yellow green yellow,     2,        FN1"
#[1] "    orange yellow red yellow,    13,      GAPDH"
#[1] "    green green green yellow,     6,      RN7SK"
#[1] "       orange yellow red red,     7,       FTH1"
#[1] "     orange green red yellow,     1,     COL1A2"
#[1] "     green yellow red yellow,     1, "
#[1] "      green green red yellow,     5,       ACTB"
#[1] "      green red green yellow,     3, "
#[1] "       green green green red,     1, "
#[1] "  green yellow yellow yellow,     1, "
#[1] "   green yellow green yellow,     1, "
#[1] "        orange green red red,     1,      RPLP0"
#[1] "        orange red red green,     1, "
# The spots were dim, or were composed of two dim spots, so the merged dot color were weak.
# too lazy for other spots, for now just disregard the refspots with just 1 coSpots
# Can't use single Guassian to fit nuclei spots...
# The red background from Atto-647N could cause a problem, maybe need to use dim Alexa 647 instead.
oneOffspringRefSpots <- refSpots%>% filter((transMergedNumCoSpotsf1 <= 1 
                                                  & transMergedNumCoSpotsf2 <= 1 
                                                  & transMergedNumCoSpotsf3 <= 1 
                                                  & transMergedNumCoSpotsf4 <= 1) 
                                                 & (transMergedNumCoSpotsf1 
                                                    + transMergedNumCoSpotsf2 
                                                    + transMergedNumCoSpotsf3 
                                                    + transMergedNumCoSpotsf4 == 1))
nrow(oneOffspringRefSpots)

oneOffspringRefSpots$allColors <- dao_assignAllColorsFromFQOrGaussianFit(
    oneOffspringRefSpots, coSpot_dist)
printCountOfMRNA(oneOffspringRefSpots$allColors, mRNAToColor)

# It's easier to get coSpots when the spot is bright
a <- merge(oneOffspringRefSpots, mRNAToColor, by = "allColors", all = TRUE)
test <- a[order(a$AMP),]

# It's easier to get coSpots when there is no nearby refSpots: 15, 14 are correct
#[1] "  orange yellow green yellow,     1,        FN1"
#[1] "    orange yellow red yellow,     3,      GAPDH"
#[1] "       orange yellow red red,     4,       FTH1"
#[1] "     orange green red yellow,     1,     COL1A2"
#[1] "      green green red yellow,     3,       ACTB"
#[1] "    green green green yellow,     2,      RN7SK"
#[1] "        green red red yellow,     1, "
numNeighborRefSpotsExcludingSelf <- num_neighbors_exclude_spot_circle(
    oneOffspringRefSpots, refSpots, 2)
oneOffspringLonelyRefSpotsLonely <- oneOffspringRefSpots[numNeighborRefSpotsExcludingSelf == 0,]
nrow(oneOffspringLonelyRefSpotsLonely)


oneOffspringLonelyRefSpotsLonely$allColors <- dao_assignAllColorsFromFQOrGaussianFit(
    oneOffspringLonelyRefSpotsLonely, coSpot_dist)
printCountOfMRNA(oneOffspringLonelyRefSpotsLonely$allColors, mRNAToColor)

which(oneOffspringRefSpots$allColors == "green red green yellow")
spot <- oneOffspringRefSpots[35,]; spot; dao_plotAllNearby(spot)

plotNearbyTiffRaster(refSpot = spot, oneOffspringLonelyRefSpots[c(43,44),], pix_XY, refScaledMatrix, isToTranslateRefSpots = FALSE, plot_area_dist = 2000)
dao_assignAllColorsFromFQOrGaussianFit(spot, coSpot_dist, isToPrintFit = TRUE)

# 8.5 Check a lonely RefSpot with 0 coSpot in each wash using the FQ spots: 48
# should disregard this spots
#[1] "      green green red yellow,     6,       ACTB"
#[1] "    green green green yellow,     9,      RN7SK"
#[1] "       green green red green,     2, "
#[1] "        red red green yellow,     1, "
#[1] "          orange red red red,     2, "
#[1] "     orange yellow green red,     1,      TGFBI"
#[1] "           green red red red,     1, "
#[1] "        green red red yellow,     3, "
#[1] "        red red green orange,     1, "
#[1] "     orange green red yellow,     4,     COL1A2"
#[1] "     green yellow red yellow,     1, "
#[1] "       green red green green,     1, "
#[1] "     orange red green yellow,     1, "
#[1] "   green yellow green yellow,     2, "
#[1] "       green green green red,     2, "
#[1] "    orange yellow red yellow,     3,      GAPDH"
#[1] "           red red green red,     1, "
#[1] "  orange yellow green yellow,     1,        FN1"
#[1] "       orange red red yellow,     2, "
#[1] "      red green green yellow,     1, "
#[1] "         green green red red,     1, "
#[1] "        green red red orange,     1, "
#[1] "      green red green yellow,     1, "
noOffspringRefSpots <- refSpots%>% filter(transMergedNumCoSpotsf1 == 0 
                                           & transMergedNumCoSpotsf2 == 0 
                                           & transMergedNumCoSpotsf3 == 0 
                                           & transMergedNumCoSpotsf4 == 0) 
nrow(noOffspringRefSpots)

noOffspringRefSpots$allColors <- dao_assignAllColorsFromFQOrGaussianFit(
    noOffspringRefSpots, coSpot_dist)
printCountOfMRNA(noOffspringRefSpots$allColors, mRNAToColor)

# It's easier to get coSpots when there is no nearby refSpots: 16
#[1] "      green green red yellow,     2,       ACTB"
#[1] "        red red green orange,     1, "
#[1] "       green red green green,     1, "
#[1] "        green red red yellow,     1, "
#[1] "     orange green red yellow,     3,     COL1A2"
#[1] "    orange yellow red yellow,     1,      GAPDH"
#[1] "           red red green red,     1, "
#[1] "    green green green yellow,     1,      RN7SK"
#[1] "       orange red red yellow,     1, "
#[1] "      red green green yellow,     1, "
#[1] "         green green red red,     1, "
#[1] "      green red green yellow,     1, "
#[1] "       green green green red,     1, "
numNeighborRefSpotsExcludingSelf <- num_neighbors_exclude_spot_circle(
    noOffspringRefSpots, refSpots, 2)
noOffspringLonelyRefSpotsLonely <- noOffspringRefSpots[numNeighborRefSpotsExcludingSelf == 0,]
nrow(noOffspringLonelyRefSpotsLonely)

noOffspringLonelyRefSpotsLonely$allColors <- dao_assignAllColorsFromFQOrGaussianFit(
    noOffspringLonelyRefSpotsLonely, coSpot_dist)
printCountOfMRNA(noOffspringLonelyRefSpotsLonely$allColors, mRNAToColor)

which(noOffspringRefSpots$allColors == "orange yellow red yellow")
spot <- noOffspringLonelyRefSpots[2,]
plotAllNearby(spot, 2)


# 9. Check refSpots that has up to 2 coSpots: 222
# Algorithm:
# 1. If the two spots are very close to each other with 0.1uM, could be a dominent color 
#       and a false-positive color => mark the false-positive color from transCoSpots
#                                  => change refSpots$transMergedNumCoSpots to 1
#                                  => change refSpot$transMergedCoSpotsIndice to one index
# 2. If the two spots are within 0.3uM could be the two spots are too close to each other to be distingquished
# 3. If the two spots have a Pos_X and Pos_Y difference, need to split the two spots.
crowd2RefSpotsIndice <- ((refSpots$transMergedNumCoSpotsf1 >= 2 
                          | refSpots$transMergedNumCoSpotsf2 >= 2 
                          | refSpots$transMergedNumCoSpotsf3 >= 2 
                          | refSpots$transMergedNumCoSpotsf4 >= 2) 
                         & !(refSpots$transMergedNumCoSpotsf1 >= 3 
                             | refSpots$transMergedNumCoSpotsf2 >= 3
                             | refSpots$transMergedNumCoSpotsf3 >= 3 
                             | refSpots$transMergedNumCoSpotsf4 >= 3))
crowd2RefSpots <- refSpots[crowd2RefSpotsIndice,]
nrow(crowd2RefSpots)

# Gradually add coSpots for crowdRefSpots that has only one coSpots
tempRefSpots2 <- refSpots
tempTransMergedSpots2<- transMergedSpots

for (secondaryCospotDist in 0.1*c(6:8)) {
    refSpots[crowd2RefSpotsIndice, ] <- dao_fillInOneCospotToXCospotsWithUnusedSpots(
        refSpots[crowd2RefSpotsIndice, ], 1, secondaryCospotDist, refSpots)
    transMergedSpots <- dao_updateIsUsedForTransMergedSpots(transMergedSpots, refSpots)
}

# Change the order of transMergedCoSpotsIndice so that brighter CoSpotIndice is at front and 
# dimmer at back
for (mergedName in mergedNames) {
    refSpots[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]] <- 
        dao_reorderIndiceByBrightness(
            refSpots[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]],
            transMergedSpots[[mergedName]])
}

crowd2RefSpots <- refSpots[crowd2RefSpotsIndice, ]





# Stragety:
1. Deal with doublet:
    a. find the extra dot when there is only one detected dot using Gaussian:
        i. Still only has the one detected dot
        ii. An extra deteced dot
    b. Get color sequence
        i. For one dot, use the dot color
        ii.For two dots, use the direction to match doublets


# 9.1 Deal with y1 spots (or any future yellow spots) that's bleakthrough from o1: 25 of them
# Mark these transMergedSpots[["f1"]] of yellow color as isBleak
for (mergedName in mergedNames) {
    transMergedSpots[[mergedName]]$isBleak <- FALSE    
}
transMergedSpots <- dao_markTransMergedSpotsAsBleack(crowd2RefSpots, transMergedSpots, mergedName = "f1", 
                                                 bleakColor = "yellow", trueColor = "orange", 
                                                 ampThreshold = 6.5) 
crowd2RefSpotsNoBleak <- dao_removeBleakCoSpots(crowd2RefSpots, transMergedSpots, mergedName = "f1")

# 9.3 Deal with single spots again: 19 spots, all are correct
#[1] "       orange yellow red red,     3,       FTH1"
#[1] "  orange yellow green yellow,     7,        FN1"
#[1] "    orange yellow red yellow,     4,      GAPDH"
#[1] "     orange yellow green red,     4,      TGFBI"
#[1] "     orange green red yellow,     1,     COL1A2"
singletCrowd2RefSpotsSpot <- crowd2RefSpotsNoBleak %>% filter((transMergedNumCoSpotsf1 <= 1  
                                                              & transMergedNumCoSpotsf2 <= 1 
                                                              & transMergedNumCoSpotsf3 <= 1 
                                                              & transMergedNumCoSpotsf4 <= 1)) 

nrow(singletCrowd2RefSpotsSpot)

singletCrowd2RefSpotsSpot$allColors <- dao_assignAllColorsFromFQOrGaussianFit(
    singletCrowd2RefSpotsSpot, coSpot_dist)
printCountOfMRNA(singletCrowd2RefSpotsSpot$allColors, mRNAToColor)

# 9.4 Deal with doubelt spots: 203 spots
# a. Try to find the 2 spots with 2D-Gaussian: 16 was wrong
#[1] "     orange green red yellow,    19,     COL1A2"
#[1] "      green green red yellow,    44,       ACTB"
#[1] "    orange yellow red yellow,    60,      GAPDH"
#[1] "  orange yellow green yellow,    30,        FN1"
#[1] "    green green green yellow,   203,      RN7SK"
#[1] "      green yellow green red,     9,       ACTG"
#[1] "     orange yellow green red,     7,      TGFBI"
#[1] "   orange green green yellow,    12,       EEF2"
#[1] "       orange yellow red red,    16,       FTH1"
#[1] "   green yellow green yellow,    11, "
#[1] "     green yellow red yellow,     2, "
#[1] "      green red green yellow,     1, "
#[1] "      orange green green red,     2, "
#[1] "        orange green red red,     1,      RPLP0"
#[1] "        green yellow red red,     1, "
# spot 48 has "yellow green green yellow"
doubletCrowd2RefSpots <- crowd2RefSpotsNoBleak %>% filter(!(transMergedNumCoSpotsf1 <= 1  
                                                                   & transMergedNumCoSpotsf2 <= 1 
                                                                   & transMergedNumCoSpotsf3 <= 1 
                                                                   & transMergedNumCoSpotsf4 <= 1)) 

nrow(doubletCrowd2RefSpots)

# 9.4 Seperate the doublet spot to closeDoublet and farDoublet
# When at one time point the doublets are too close to distinguish, give up the dimmer spot
# Number of close doublet: 7 spots
#[1] "     orange green red yellow,     5,     COL1A2"
#[1] "  orange yellow green yellow,     1,        FN1"
#[1] "      green green red yellow,     1,       ACTB"
indiceOfCloseDoubletSpots <- dao_getIndiceOfCloseDoubletSpots(doubletCrowd2RefSpots, 0.1) 
closeDoubletCrowd2RefSpots <- doubletCrowd2RefSpots[indiceOfCloseDoubletSpots, ]
nrow(closeDoubletCrowd2RefSpots)

closeDoubletCrowd2RefSpots$allColors <- dao_assignAllColorsFromFQOrGaussianFit(
    closeDoubletCrowd2RefSpots)
printCountOfMRNA(closeDoubletCrowd2RefSpots$allColors, mRNAToColor)

# 9.5 Number of farDoublet: 196 spots
# Try to also just get the brightest spot of the 2
#[1] "     orange green red yellow,     7,     COL1A2"
#[1] "    orange yellow red yellow,    28,      GAPDH"
#[1] "  orange yellow green yellow,    12,        FN1"
#[1] "    green green green yellow,   102,      RN7SK"
#[1] "      green green red yellow,    22,       ACTB"
#[1] "      green yellow green red,     6,       ACTG"
#[1] "     orange yellow green red,     6,      TGFBI"
#[1] "       orange yellow red red,     9,       FTH1"
#[1] "      green red green yellow,     1, "
#[1] "   green yellow green yellow,     3, "
# spot 94 is "green red green yellow": not sure what's wrong, there shouldn't be any R2 spots.
# spot 139 is "green yellow green yellow", should be GGGY but G2 should include 3 spots: 2G and 1Y, 
        #a further G brighter than Y is not within coSpotDist
# spot 162 is "green yellow green yellow", should be 2 correct spots of similar brightness, GGRY & OYGR
        #Therefore position is more important than brightness
# spot 189 is "green yellow green yellow", sould be 2 correct spots of similar brightness, OYGR & GGGY
        # However, daostorm miss a Y4 spot in the crowd space
farDoubletCrowdRefSpots <- doubletCrowd2RefSpots[-indiceOfCloseDoubletSpots, ]
nrow(farDoubletCrowdRefSpots)

farDoubletCrowdRefSpots$allColors <- dao_assignAllColorsFromFQOrGaussianFit(
    farDoubletCrowdRefSpots)
printCountOfMRNA(farDoubletCrowdRefSpots$allColors, mRNAToColor)
farDoubletCrowdRefSpots[139,]$allColors <- "green green green yellow"
farDoubletCrowdRefSpots[162,]$allColors <- "green green red yellow"
farDoubletCrowdRefSpots[189,]$allColors <- "orange yellow green red"

indice <- grep("green yellow green yellow", farDoubletCrowdRefSpots$allColors);indice
spot <- farDoubletCrowdRefSpots[189,]; spot; dao_plotAllNearby(spot)



doubletCrowd2RefSpots$allColors <- dao_assignAllColorsForDoublets(doubletCrowd2RefSpots)
printCountOfMRNA(doubletCrowd2RefSpots$allColors, mRNAToColor)



indice <- grep("green red green yellow", doubletCrowd2RefSpotsSpot$allColors);indice
spot <- doubletCrowd2RefSpotsSpot[1,]; spot; plotAllNearby(spot)
spot$allColors <- dao_assignAllColorsForDoublets(spot, TRUE)
printCountOfMRNA(spot$allColors, mRNAToColor)



# 10. refSpots that have 3 coSpots at a time point: 14 spots
#[1] "       green green green red,     1, "
#[1] "      green green red yellow,     2,       ACTB"
#[1] "     orange green red yellow,     3,     COL1A2"
#[1] "    green green green yellow,     4,      RN7SK"
#[1] "    orange yellow red yellow,     1,      GAPDH"
#[1] "      green yellow green red,     1,       ACTG"
#[1] "     green yellow red yellow,     1, "
#[1] "   orange green green yellow,     1,       EEF2"
# spot 1 was green green green red, but it was 2 GGGY sandwiching a OYGR FN1
    # It seems that scaledAMP is not as good as normalTiffAmp when the red color had a strong bg
# spot 13 was green yellow red yellow. but it was 2 spots of similar brightness: OYRY GAPDH & GYGR ACTG
verycrowdRefSpots <- refSpots%>% filter((transMergedNumCoSpotsf1 >= 3 
                                         | transMergedNumCoSpotsf2 >= 3
                                         | transMergedNumCoSpotsf3 >= 3 
                                         | transMergedNumCoSpotsf4 >= 3))
nrow(verycrowdRefSpots)

verycrowdRefSpots$allColors <- dao_assignAllColorsFromFQOrGaussianFit(
    verycrowdRefSpots)
printCountOfMRNA(verycrowdRefSpots$allColors, mRNAToColor)
verycrowdRefSpots[1,]$allColors <- "green green green yellow"
verycrowdRefSpots[13,]$allColors <- "orange yellow red yellow"




indice <- grep("green yellow red yellow", verycrowdRefSpots$allColors);indice
spot <- verycrowdRefSpots[13,]; spot; dao_plotAllNearby(spot)




1. Go through all the doublet spots, and get the relative position of the spot and amp 
2. The easiest strategy: assume that daostorm get all the spots, so 1 spot a composite of 2 spots,
so allColors is a list of vectors.



for (mergedName in mergedNames) {
    if (spot[paste("transMergedNumCoSpots", mergedName, sep = "")] == 1) {
        print(mergedName)
        singletSpotIndex <- spot[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]][[1]]
        singletSpot <- transMergedSpots[[mergedName]][singletSpotIndex,]
        singletScaledAmp <- singletSpot$ScaledAmp
        print(paste(singletScaledAmp))
    }
    if (spot[paste("transMergedNumCoSpots", mergedName, sep = "")] == 2) {
        print(mergedName)
        transCospotsIndice <- spot[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]][[1]]
        brightIndex <- transCospotsIndice[which.max(transMergedSpots[[mergedName]][transCospotsIndice, ]$ScaledAmp)]
        dimIndex <- transCospotsIndice[transCospotsIndice != strongIndex]
        brightSpot <- transMergedSpots[[mergedName]][brightIndex,]
        dimSpot <- transMergedSpots[[mergedName]][dimIndex,]
        brightScaledAmp <- brightSpot$ScaledAmp
        dimScaledAmp <- dimSpot$ScaledAmp
        relativeX <- brightSpot$Pos_X - dimSpot$Pos_X
        relativeY <- brightSpot$Pos_Y - dimSpot$Pos_Y
        print(paste(brightScaledAmp, dimScaledAmp, relativeX, relativeY))
    }
}





# 9.1 Check refSpots that has 2 coSpots in a wash, and rest 0 or 1 coSpots: 154
crowd2RefSpots <- refSpots%>% filter((transMergedNumCoSpotsf1 == 2 
                                         & transMergedNumCoSpotsf2 < 2 
                                         & transMergedNumCoSpotsf3 < 2 
                                         & transMergedNumCoSpotsf4 < 2) 
                                        | (transMergedNumCoSpotsf1 < 2 
                                           & transMergedNumCoSpotsf2 == 2 
                                           & transMergedNumCoSpotsf3 < 2 
                                           & transMergedNumCoSpotsf4 < 2) 
                                        | (transMergedNumCoSpotsf1 < 2 
                                           & transMergedNumCoSpotsf2 < 2 
                                           & transMergedNumCoSpotsf3 == 2 
                                           & transMergedNumCoSpotsf4 < 2) 
                                        | (transMergedNumCoSpotsf1 < 2 
                                           & transMergedNumCoSpotsf2 < 2 
                                           & transMergedNumCoSpotsf3 < 2 
                                           & transMergedNumCoSpotsf4 == 2))
nrow(crowd2RefSpots)

# 9.2 Separate into crowd2SameColorRefSpots (2 coSpots in a time have the same color): 95
# and crowd2DiffColorRefSpots (2 coSpots in a time have different color): 59
crowd2SameColorRefSpots <- spotDataFrame(crowd2RefSpots)
crowd2DiffColorRefSpots <- spotDataFrame(crowd2RefSpots)
for(i in seq_along(crowd2RefSpots[,1])) {
    spot <- crowd2RefSpots[i,]
    for (mergedName in mergedNames) {
        if (spot[[paste("transMergedNumCoSpots", mergedName, sep = "")]] == 2) {
            transCospot1 <- transMergedSpots[[mergedName]][
                spot[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]][[1]][1],]
            transCospot2 <- transMergedSpots[[mergedName]][
                spot[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]][[1]][2],]
        }
    }
    if (transCospot1$Color == transCospot2$Color) {
        crowd2SameColorRefSpots <- rbind(crowd2SameColorRefSpots, spot)
    } else {
        crowd2DiffColorRefSpots <- rbind(crowd2DiffColorRefSpots, spot)
    }
}
nrow(crowd2SameColorRefSpots)
nrow(crowd2DiffColorRefSpots)


# 9.3 Separate crowd2SameColorRefSpots into 
# crowd2CloseSameColorRefSpots (2 same color spots distance within 1 um): 95 and 
# crowd2FarSameColorRefSpots (2 same color spot distance above 1 um): 0
crowd2CloseSameColorRefSpots <- spotDataFrame(crowd2SameColorRefSpots)
crowd2FarSameColorRefSpots <- spotDataFrame(crowd2SameColorRefSpots)
for(i in seq_along(crowd2SameColorRefSpots[,1])) {
    spot <- crowd2SameColorRefSpots[i,]
    for (mergedName in mergedNames) {
        if (spot[[paste("transMergedNumCoSpots", mergedName, sep = "")]] == 2) {
            transCospot1 <- transMergedSpots[[mergedName]][
                spot[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]][[1]][1],]
            transCospot2 <- transMergedSpots[[mergedName]][
                spot[[paste("transMergedCoSpotsIndice", mergedName, sep = "")]][[1]][2],]
        }
    }
    if (distance(transCospot1, transCospot2) < 1) {
        crowd2CloseSameColorRefSpots <- rbind(crowd2CloseSameColorRefSpots, spot)   
    } else {
        crowd2FarSameColorRefSpots <- rbind(crowd2FarSameColorRefSpots, spot)
    }
}
nrow(crowd2CloseSameColorRefSpots)
nrow(crowd2FarSameColorRefSpots)


# 9.4 crowd2CloseSameColorRefSpots observations:
# crowd2CloseSameColorLonelyRefSpots_1Offsprings (with 1 offspring): 1 RN7SK, correct => 0
# crowd2CloseSameColorLonelyRefSpots_234Offsprings (with 2, 3, or 4 offsprings): 95
#[1] "       orange yellow red red,     7,       FTH1" x2
#[1] "  orange yellow green yellow,     5,        FN1" x2
#[1] "      green green red yellow,     7,       ACTB" x2
#[1] "    green green green yellow,    48,      RN7SK" x2
#[1] "    orange yellow red yellow,    24,      GAPDH" x2
#[1] "      green yellow green red,     3,       ACTG" x2
#[1] "     orange green red yellow,     1,     COL1A2" x2
# should times the number of spots by 2 => 95 * 2 = 190
# => keep crowd2CloseSameColorLonelyRefSpots_1Offsprings as long as other images have visible bright spots
# => keep crowd2CloseSameColorLonelyRefSpots_234Offsprings 
crowd2CloseSameColorRefSpots_1Offsprings <- crowd2CloseSameColorRefSpots %>% filter(
    (transMergedNumCoSpotsf1 + transMergedNumCoSpotsf2 + transMergedNumCoSpotsf3 + transMergedNumCoSpotsf4) == 2)

crowd2CloseSameColorRefSpots_1Offsprings$allColors <- dao_assignAllColorsFromFQOrGaussianFit(
    crowd2CloseSameColorRefSpots_1Offsprings, coSpot_dist)
printCountOfMRNA(crowd2CloseSameColorRefSpots_1Offsprings$allColors, mRNAToColor)
nrow(crowd2CloseSameColorRefSpots_1Offsprings)

crowd2CloseSameColorRefSpots_234Offsprings <- crowd2CloseSameColorRefSpots %>% filter(
    (transMergedNumCoSpotsf1 + transMergedNumCoSpotsf2 + transMergedNumCoSpotsf3 + transMergedNumCoSpotsf4) > 2)

crowd2CloseSameColorRefSpots_234Offsprings$allColors <- dao_assignAllColorsFromFQOrGaussianFit(
    crowd2CloseSameColorRefSpots_234Offsprings, diffractionLimit, coSpot_dist)
printCountOfMRNA(crowd2CloseSameColorRefSpots_234Offsprings$allColors, mRNAToColor)
nrow(crowd2CloseSameColorRefSpots_234Offsprings)


# 9.5 crowd2FarSameColorLonelyRefSpots observations (2 same color spot distance above 1 um)
# 0 spots, should times 2 as well
crowd2FarSameColorRefSpots$allColors <- dao_assignAllColorsFromFQOrGaussianFit(
    crowd2FarSameColorRefSpots, coSpot_dist)
nrow(crowd2FarSameColorRefSpots)
printCountOfMRNA(crowd2FarSameColorRefSpots$allColors, mRNAToColor)


# 9.6 crowd2DiffColorRefSpots observations: 59, which takes too long to analyze, give up
# crowd2DiffColorRefSpots_nuclei: 25 were witin nuclei, which could be too hard to analyze
# crowd2DiffColorRefSpots_not_nuclei: 34 were outside of nuclei
#[1] "      green green red yellow,    37,       ACTB"
#[1] "     orange green red yellow,    31,     COL1A2"
#[1] "    orange yellow red yellow,     4,      GAPDH"
#[1] "     green yellow red yellow,     4, "
#[1] "    green green green yellow,     7,      RN7SK"
#[1] "   green yellow green yellow,     2, "
#[1] "  orange yellow green yellow,     3,        FN1"
#[1] "   orange green green yellow,     2,       EEF2"
#[1] "      orange green red green,     1, "
#[1] "     orange yellow red green,     1, "
#[1] "        green yellow red red,     1, "
#[1] "       orange yellow red red,     1,       FTH1"
# Also, the red fluorescent background is really annoying. Need a better image.
# spot 9 69 75 90 are GYRY, which are bright spot along with a dim spot, so in other time points that are either
# 0 or 1 spot, I should look for another dim spot or a doublet spot of the same color.
crowd2DiffColorRefSpots_nuclei <- enclosedSpotsWithinBoundaries(crowd2DiffColorRefSpots, 
                                                                     nucleiPolygons_um)
nrow(crowd2DiffColorRefSpots_nuclei)
crowd2DiffColorRefSpots_not_nuclei <- crowd2DiffColorRefSpots[!rownames(crowd2DiffColorRefSpots) %in% 
                                                                              (rownames(crowd2DiffColorRefSpots_nuclei)),]
nrow(crowd2DiffColorRefSpots_not_nuclei)

crowd2DiffColorRefSpots_not_nuclei_splitSpots <- assign4ColorsFromFQOrGaussianFit_splitSpot(
    crowd2DiffColorRefSpots_not_nuclei, diffractionLimit, coSpot_dist)
printCountOfMRNA(crowd2DiffColorLonelyRefSpots_not_nuclei_splitSpots$allColors, mRNAToColor)
crowd2DiffColorLonelyRefSpots_not_nuclei_splitSpots <- 
    merge(crowd2DiffColorLonelyRefSpots_not_nuclei_splitSpots, mRNAToColor, by = "allColors", all.x = TRUE)

which(crowd2DiffColorLonelyRefSpots_not_nuclei_splitSpots$allColors == "green yellow red yellow")
spot <- crowd2DiffColorLonelyRefSpots_not_nuclei_splitSpots[90,]
plotAllNearby(spot)


printCountOfMRNA(crowd2FarSameColorLonelyRefSpots$allColors, mRNAToColor)



spot <- crowd2DiffColorLonelyRefSpots_not_nuclei[2,]
plotAllNearby(spot)


# change
# Assume it is 2 spots, later will assume it is 3 spots
# Get the first time that has 2 spots.
# Use the perfectly fit Gaussian as the location
transSpot1 <- transF2Spots[33, ]
transSpot2 <- transF2Spots[777,]
refSpot1 <- spot
refSpot1[,1:2] <- transSpot1[,1:2]
refSpot1[,11:14] <- c(0,1,0,0)
refSpot1$transF2CoSpotsIndice <- 33
refSpot1$allColors <- assign4ColorsFromFQOrGaussianFit(
    refSpot1, diffractionLimit, coSpot_dist)
plotAllNearby(refSpot1)
# y2 is a perfect GaussianFit
# g2 has large enough amplitube, but a fat radius, so it could be composed of two spots
g2Fit <- fitGaussian2D(spot, g2TiffMatrix,diffractionLimit)
y2Fit <- fitGaussian2D(spot, y2TiffMatrix,diffractionLimit)

fit <- fitGaussian2DFor2From1Spot(spot, refTiffMatrix, nearby_distance_um = 1.5)
spot1 <- spot
spot1$Pos_X <- 127.19
spot1$Pox_Y <- 39.03
spot2 <- spot 
spot2$Pos_X <- 125.79
spot2$Pox_Y <- 38.04
plotAllNearby(spot1)
getNearbyMatrix(refTiffMatrix, spotsToIndice2D(spot,pix_XY), 3)



4: YRY
3: GGG
2: GYG
1: GGG


# 9.3 rest spots: 83, which takes too long to analyze, give up
restRefSpots <- refSpots %>% filter(!(transMergedNumCoSpotsf1 <= 1 
                                      & transMergedNumCoSpotsf2 <= 1 
                                      & transMergedNumCoSpotsf3 <= 1 
                                      & transMergedNumCoSpotsf4 <= 1)
                                    &!((transMergedNumCoSpotsf1 == 2 
                                        & transMergedNumCoSpotsf2 < 2 
                                        & transMergedNumCoSpotsf3 < 2 
                                        & transMergedNumCoSpotsf4 < 2) 
                                       |(transMergedNumCoSpotsf1 < 2 
                                         & transMergedNumCoSpotsf2 == 2 
                                         & transMergedNumCoSpotsf3 < 2 
                                         & transMergedNumCoSpotsf4 < 2) 
                                       |(transMergedNumCoSpotsf1 < 2 
                                         & transMergedNumCoSpotsf2 < 2 
                                         & transMergedNumCoSpotsf3 == 2 
                                         & transMergedNumCoSpotsf4 < 2) 
                                       |(transMergedNumCoSpotsf1 < 2 
                                         & transMergedNumCoSpotsf2 < 2 
                                         & transMergedNumCoSpotsf3 < 2 
                                         & transMergedNumCoSpotsf4 == 2)))



