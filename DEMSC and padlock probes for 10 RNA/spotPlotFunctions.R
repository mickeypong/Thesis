library(plyr)
library(dplyr)
library(ggplot2)
library(tiff)
library(gridExtra)
library(raster)  # for raster images
library(pracma)  # for lsqcurvefit
library(sp)      # for create shape file   
library(rgdal)   # for create shape file
library(mgcv)    # for checking if sets of spots are within coords of a polygon
library(AtmRay)  # for meshgrid to generate index X and Y of a matrix


# Scale the matrix to minus the min (should be bg), and the set value from 0 to 1 for plot
scaled_matrix = function (tiffMatrix, curve = 1) {
    minIntensity <- min(min(tiffMatrix))
    maxIntensity <- max(max(tiffMatrix))
    scaledTiffMatrix <- (tiffMatrix - minIntensity) / (maxIntensity - minIntensity)
    scaledTiffMatrix <- 1 - (1 - scaledTiffMatrix)^curve 
}


# Load outline file from FQ to create shape file.
# Input: outlineFilename
# Output: list of polygons, each as set of indice2D 
## dataframe of xIndice and yIndice
readOutlines = function (outlineFilename, first_important_line) {
    all_content = readLines(outlineFilename);
    outline_content = all_content[-1:-(first_important_line - 1)];
    
    polygons <- list()
    for (i in 0 : (length(outline_content) / 5 - 1)) {
        cellName <- outline_content[i * 5 + 1]
        cellName <- strsplit(cellName, "\t")[[1]][2]
        xIndice <- strsplit(outline_content[i * 5 + 2], "\t")[[1]][-1]
        xIndice <- as.numeric(c(xIndice, xIndice[1]))
        yIndice <- strsplit(outline_content[i * 5 + 3], "\t")[[1]][-1]
        yIndice <- as.numeric(c(yIndice, yIndice[1]))
        indice2D <- getIndice2DFromXY(xIndice, yIndice)
        # indice2D <- cbind(xIndice, yIndice)
        polygons <- c(polygons, list(indice2D))  
    }
    polygons
}


# Input: 2 vectors, or 2 named vectors xIndice, yIndice
# Output: indice in 2D, actually a matrix with 2 columns with yIndice in the 1st column 
# and xIndice in the 2nd column
getIndice2DFromXY = function (xIndice, yIndice) {
    indice2D <- cbind(yIndice[], xIndice[])
    indice2D
}


# Get the xIndice from indice2D as the 2nd column
getXIndiceFromIndice2D = function(indice2D) {
    indice2D[,2]
}


# Get the yIndice from indice2D as the 1st column
getYIndiceFromIndice2D = function(indice2D) {
    indice2D[,1]
} 


# Get the x coordinate (um) from xIndice
getXCoordsFromIndice2D = function(indice2D, pix_XY. = pix_XY) {
    xIndice <- getXIndiceFromIndice2D(indice2D)
    x <- (xIndice - 1) * pix_XY.
    x
}


# Get the y coordinate (um) from yIndice
getYCoordsFromIndice2D = function(indice2D, pix_XY. = pix_XY) {
    yIndice <- getYIndiceFromIndice2D(indice2D)
    y <- (yIndice - 1) * pix_XY.
    y
}


# Return whether matrixIndice as cbind (y, x) are in a polygon.
# isInPolygon = function (matrixIndice, polygon) {
#     coords <- cbind(matrixIndice[,2], matrixIndice[,1])
#     isInThePolygon <- in.out(polygon, coords)
#     isInThePolygon
# }
isInPolygon = function (indice2D, polygon) {
    isInThePolygon <- in.out(polygon, indice2D)
    isInThePolygon
}



# Return a matrix that has value 0 at a Indice of it's not in any polygon
# and return 1 if it's at polygon number 1, etc.
getPolygonMatrix = function (polygons, resX, resY) {
    matrix <- matrix(0, resY, resX)
    allIndice1D <- 1 : (resX * resY)
    allIndice2D <- getIndice2DFrom1D(allIndice1D, resY)
    for (i in seq_along(polygons)) {
        polygon <- polygons[[i]]
        matrix[in.out(polygon, allIndice2D)] <- i
    }
    matrix
}


# return 2d indice (y, x) of a matrix using 1d vector index
# Note that 1d vector index for a 2D matrix is from up to down, from left to right
getIndice2DFrom1D = function (indice1D, ncolMatrix) {
    y2D <- indice1D %% ncolMatrix
    x2D <- floor((indice1D - 1) / ncolMatrix)
    indice2D <- getIndice2DFromXY(x2D, y2D)
    indice2D
}


# get index of polygon for spots
getPolygonIndex = function(spots, pix_XY, polygonMatrix) {
    indice2D <- spotsToIndice2D(spots, pix_XY)
    indice <- polygonMatrix[indice2D]
    indice
}


#Get new polygons with Y reverse coordinates as resY - Y, for plotting.
polygons_reverseY = function (polygons, resY) {
    polygons_revY <- list()
    for (i in 1 : length(polygons)) {
        polygon <- polygons[[i]]
        xIndice <- getXIndiceFromIndice2D(polygon)
        yIndice <- getYIndiceFromIndice2D(polygon)
        polygon_revY<- getIndice2DFromXY(xIndice, resY - yIndice)
            # cbind(polygon[,1], resY - polygon[,2])
        polygons_revY <- c(polygons_revY, list(polygon_revY))
    }
    polygons_revY
}


# Change the polygon units from pixel to um, this time 1st column is X_Pos, 
# and 2nd column is Y_Pos
polygons_toMicron = function (polygons, pix_XY) {
    polygons_um <- list()
    for (i in 1 : length(polygons)) {
        polygon <- polygons[[i]]
        xIndice <- getXIndiceFromIndice2D(polygon)
        yIndice <- getYIndiceFromIndice2D(polygon)
        polygon_um <- cbind(xIndice, yIndice) * pix_XY
        polygons_um <- c(polygons_um, list(polygon_um))
    }    
    polygons_um
}


# Note that the matrix index is (Y, X), and the plot has y bottom as 0 and top as resY
# while x left as 0 and right as resX
# plot the matrix as raster image. Can deside whether to draw outline of cells or not.
# Note the matrix plot row is intrisically from top to down, on the other hand,
# spot plot intrinsically y axis is from down to up, so need reverseY for spot plot
plotMatrix = function (matrix, palette = greenPalette, title = "", polygons_revY, isDrawOutlines = TRUE, withOutlineNumber = FALSE,
                       nuclei_polygons_revY = NULL) {
    resY = dim(matrix)[1]
    resX = dim(matrix)[2]
    plot(raster(matrix, ymn = 0, ymx = resY, xmn = 0, xmx = resX), 
         col = palette, xlim=c(0,resX),ylim=c(0,resY), main = title)
    if (isDrawOutlines) {
        for (i in 1 : length(polygons_revY)) {
            polygon_revY <- polygons_revY[[i]]
            lines(x = getXIndiceFromIndice2D(polygon_revY), y = getYIndiceFromIndice2D(polygon_revY), col = "white")
        }
    }
    if (withOutlineNumber) {
        for (i in 1 : length(polygons_revY)) {
            polygon_revY <- polygons_revY[[i]]
            xIndice <- getXIndiceFromIndice2D(polygon_revY)
            yIndice <- getYIndiceFromIndice2D(polygon_revY)
            text(xIndice[1] + 10, yIndice[1], labels = i, col = "white")
        }
    }
    if (!is.null(nuclei_polygons_revY)) {
        for (i in 1 : length(nuclei_polygons_revY)) {
            nuclei_polygon_revY <- nuclei_polygons_revY[[i]]
            xIndice <- getXIndiceFromIndice2D(nuclei_polygon_revY)
            yIndice <- getYIndiceFromIndice2D(nuclei_polygon_revY)
            lines(x = xIndice, y = yIndice, col = "yellow")
            text(xIndice[1] + 10, yIndice[1], labels = i, col = "yellow")
        }
    }
}


# Plot bounded matrix by the polygon_revY, using the pixel bounded by the min and 
# max polygon_revY coords.
plotBoundedMatrix = function (matrix, palette = greenPalette, title = "", polygon_revY, isDrawOutlines = TRUE) {
    resY = dim(refScaledMatrix)[1]
    resX = dim(refScaledMatrix)[2]
    minX <- min(getXIndiceFromIndice2D(polygon_revY))
    maxX <- max(getXIndiceFromIndice2D(polygon_revY))
    minY <- min(getYIndiceFromIndice2D(polygon_revY))
    maxY <- max(getYIndiceFromIndice2D(polygon_revY))
    plot(raster(matrix, ymn = 0, ymx = resY, xmn = 0, xmx = resX), 
         col = palette, xlim=c(minX,maxX),ylim=c(minY,maxY), main = title)
    if (isDrawOutlines) {
        lines(x = getXIndiceFromIndice2D(polygon_revY), y = getYIndiceFromIndice2D(polygon_revY), col = "white")
    }
}



# Plot all nearby Spots in transF1, transF2, transF3, and transF4,
# and plot all nearby Spots in 16 tiff matrix images with the refSpot pixel
# circled in red and neighbor fluor spots circled in white
plotAllNearby = function (spot, plot_area_dist = 2, toTranslateRef = TRUE) {
    zoomedF1Plot <- plotNearbySpots(spot, transF1Spots, refSpots, plot_area_dist = plot_area_dist, "nearby transF1")
    zoomedF2Plot <- plotNearbySpots(spot, transF2Spots, refSpots, plot_area_dist = plot_area_dist, "nearby transF2")
    zoomedF3Plot <- plotNearbySpots(spot, transF3Spots, refSpots, plot_area_dist = plot_area_dist, "nearby transF3")
    zoomedF4Plot <- plotNearbySpots(spot, transF4Spots, refSpots, plot_area_dist = plot_area_dist, "nearby transF4")
    grid.arrange(zoomedF1Plot, zoomedF2Plot, zoomedF3Plot, zoomedF4Plot, ncol = 2)
    
    par(mfrow = c(1,1))
    plotNearbyTiffRaster(spot, refSpots, pix_XY, refScaledMatrix, isToTranslateRefSpots = FALSE, col = greenPalette, plot_area_dist = plot_area_dist, plotname = "nearby ref spots on tiff")
    par(mfrow = c(2,2))
    plotNearbyTiffRaster(spot, g1Spots, pix_XY, g1ScaledMatrix, transG1XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = greenPalette, plot_area_dist, plotname = "nearby g1 spots on tiff")
    plotNearbyTiffRaster(spot, y1Spots, pix_XY, y1ScaledMatrix, transY1XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = yellowPalette, plot_area_dist, plotname = "nearby y1 spots on tiff")
    plotNearbyTiffRaster(spot, o1Spots, pix_XY, o1ScaledMatrix, transO1XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = orangePalette, plot_area_dist, plotname = "nearby o1 spots on tiff")
    plotNearbyTiffRaster(spot, r1Spots, pix_XY, r1ScaledMatrix, transR1XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = redPalette, plot_area_dist, plotname = "nearby r1 spots on tiff")
    plotNearbyTiffRaster(spot, g2Spots, pix_XY, g2ScaledMatrix, transG2XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = greenPalette, plot_area_dist, plotname = "nearby g2 spots on tiff")
    plotNearbyTiffRaster(spot, y2Spots, pix_XY, y2ScaledMatrix, transY2XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = yellowPalette, plot_area_dist, plotname = "nearby y2 spots on tiff")
    plotNearbyTiffRaster(spot, o2Spots, pix_XY, o2ScaledMatrix, transO2XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = orangePalette, plot_area_dist, plotname = "nearby o2 spots on tiff")
    plotNearbyTiffRaster(spot, r2Spots, pix_XY, r2ScaledMatrix, transR2XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = redPalette, plot_area_dist, plotname = "nearby r2 spots on tiff")
    plotNearbyTiffRaster(spot, g3Spots, pix_XY, g3ScaledMatrix, transG3XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = greenPalette, plot_area_dist, plotname = "nearby g3 spots on tiff")
    plotNearbyTiffRaster(spot, y3Spots, pix_XY, y3ScaledMatrix, transY3XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = yellowPalette, plot_area_dist, plotname = "nearby y3 spots on tiff")
    plotNearbyTiffRaster(spot, o3Spots, pix_XY, o3ScaledMatrix, transO3XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = orangePalette, plot_area_dist, plotname = "nearby o3 spots on tiff")
    plotNearbyTiffRaster(spot, r3Spots, pix_XY, r3ScaledMatrix, transR3XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = redPalette, plot_area_dist, plotname = "nearby r3 spots on tiff")
    plotNearbyTiffRaster(spot, g4Spots, pix_XY, g4ScaledMatrix, transG4XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = greenPalette, plot_area_dist, plotname = "nearby g4 spots on tiff")
    plotNearbyTiffRaster(spot, y4Spots, pix_XY, y4ScaledMatrix, transY4XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = yellowPalette, plot_area_dist, plotname = "nearby y4 spots on tiff")
    plotNearbyTiffRaster(spot, o4Spots, pix_XY, o4ScaledMatrix, transO4XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = orangePalette, plot_area_dist, plotname = "nearby o4 spots on tiff")
    plotNearbyTiffRaster(spot, r4Spots, pix_XY, r4ScaledMatrix, transR4XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = redPalette, plot_area_dist, plotname = "nearby r4 spots on tiff")
}


# Plot all nearby Spots in transF1, transF2, transF3, and transF4,
# and plot all nearby Spots in 16 tiff matrix images with the refSpot pixel
# circled in red and neighbor fluor spots circled in white
plotAllNearbyFor6mRNAs = function (spot, plot_area_dist = 2, toTranslateRef = TRUE) {
    zoomedF1Plot <- plotNearbySpots(spot, transF1Spots, refSpots, plot_area_dist = plot_area_dist, "nearby transF1")
    zoomedF2Plot <- plotNearbySpots(spot, transF2Spots, refSpots, plot_area_dist = plot_area_dist, "nearby transF2")
    zoomedF3Plot <- plotNearbySpots(spot, transF3Spots, refSpots, plot_area_dist = plot_area_dist, "nearby transF3")
    zoomedF4Plot <- plotNearbySpots(spot, transF4Spots, refSpots, plot_area_dist = plot_area_dist, "nearby transF4")
    grid.arrange(zoomedF1Plot, zoomedF2Plot, zoomedF3Plot, zoomedF4Plot, ncol = 2)
    
    par(mfrow = c(1,1))
    plotNearbyTiffRaster(spot, refSpots, pix_XY, refScaledMatrix, isToTranslateRefSpots = FALSE, col = greenPalette, plot_area_dist = plot_area_dist, plotname = "nearby ref spots on tiff")
    par(mfrow = c(1,2))
    plotNearbyTiffRaster(spot, g1Spots, pix_XY, g1ScaledMatrix, transG1XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = greenPalette, plot_area_dist, plotname = "nearby g1 spots on tiff")
    #plotNearbyTiffRaster(spot, y1Spots, pix_XY, y1ScaledMatrix, transY1XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = yellowPalette, plot_area_dist, plotname = "nearby y1 spots on tiff")
    plotNearbyTiffRaster(spot, o1Spots, pix_XY, o1ScaledMatrix, transO1XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = orangePalette, plot_area_dist, plotname = "nearby o1 spots on tiff")
    #plotNearbyTiffRaster(spot, r1Spots, pix_XY, r1ScaledMatrix, transR1XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = redPalette, plot_area_dist, plotname = "nearby r1 spots on tiff")
    plotNearbyTiffRaster(spot, g2Spots, pix_XY, g2ScaledMatrix, transG2XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = greenPalette, plot_area_dist, plotname = "nearby g2 spots on tiff")
    plotNearbyTiffRaster(spot, y2Spots, pix_XY, y2ScaledMatrix, transY2XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = yellowPalette, plot_area_dist, plotname = "nearby y2 spots on tiff")
    #plotNearbyTiffRaster(spot, o2Spots, pix_XY, o2ScaledMatrix, transO2XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = orangePalette, plot_area_dist, plotname = "nearby o2 spots on tiff")
    #plotNearbyTiffRaster(spot, r2Spots, pix_XY, r2ScaledMatrix, transR2XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = redPalette, plot_area_dist, plotname = "nearby r2 spots on tiff")
    plotNearbyTiffRaster(spot, g3Spots, pix_XY, g3ScaledMatrix, transG3XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = greenPalette, plot_area_dist, plotname = "nearby g3 spots on tiff")
    #plotNearbyTiffRaster(spot, y3Spots, pix_XY, y3ScaledMatrix, transY3XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = yellowPalette, plot_area_dist, plotname = "nearby y3 spots on tiff")
    #plotNearbyTiffRaster(spot, o3Spots, pix_XY, o3ScaledMatrix, transO3XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = orangePalette, plot_area_dist, plotname = "nearby o3 spots on tiff")
    plotNearbyTiffRaster(spot, r3Spots, pix_XY, r3ScaledMatrix, transR3XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = redPalette, plot_area_dist, plotname = "nearby r3 spots on tiff")
    #plotNearbyTiffRaster(spot, g4Spots, pix_XY, g4ScaledMatrix, transG4XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = greenPalette, plot_area_dist, plotname = "nearby g4 spots on tiff")
    plotNearbyTiffRaster(spot, y4Spots, pix_XY, y4ScaledMatrix, transY4XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = yellowPalette, plot_area_dist, plotname = "nearby y4 spots on tiff")
    #plotNearbyTiffRaster(spot, o4Spots, pix_XY, o4ScaledMatrix, transO4XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = orangePalette, plot_area_dist, plotname = "nearby o4 spots on tiff")
    plotNearbyTiffRaster(spot, r4Spots, pix_XY, r4ScaledMatrix, transR4XYPerPolygon, isToTranslateRefSpots = toTranslateRef, col = redPalette, plot_area_dist, plotname = "nearby r4 spots on tiff")
}


# Plot spots on tiff matrix, with each kind of mRNAs of different color
# depending the mRNAToColor's plotColor value
plotNearbyTiffRasterWithMRNAs = function (refSpot, spots, mRNAToColor, pix_XY, normTiffMatrix,
                                          plot_area_dist = 20, spotSize = 2, plotname = "", maxNumY = 7, 
                                          dx = 120, dy = -40, pointTextXDistance = 60) {
    
    spots <- merge(spots, mRNAToColor, by = "allColors", all.x = TRUE)
    resX <- dim(normTiffMatrix)[2]
    resY <- dim(normTiffMatrix)[1]
    refIndice2D_int <- spotsToIndice2D(refSpot, pix_XY)
    pixX_int <- getXIndiceFromIndice2D(refIndice2D_int)
    pixY_int <- getYIndiceFromIndice2D(refIndice2D_int)
    pixPlotArea <- round(plot_area_dist / pix_XY)
    minPixX <- max(1, pixX_int - pixPlotArea)
    maxPixX <- min(pixX_int + pixPlotArea, resX)
    minPixY <- max(1, pixY_int - pixPlotArea)
    maxPixY <- min(pixY_int + pixPlotArea, resY)
    
    nearbyMatrix <- normTiffMatrix[minPixY : maxPixY, minPixX : maxPixX]
    rasterImage <- raster(nearbyMatrix, ymn = minPixY-0.5, ymx = maxPixY+0.5, xmn = minPixX-0.5, xmx = maxPixX+0.5) 
    plot(rasterImage, col = colorRampPalette(c("black", "white")) (1000), 
         xlim=c(minPixX-0.5, maxPixX+0.5),ylim=c(minPixY-0.5, maxPixY+0.5), 
         main = plotname)
    
    neighborSpots <- neighbors_exclude_spot(refSpot, spots, plot_area_dist)
    neightborIndice2D <- spotsToIndice2D_not_integer(neighborSpots, pix_XY)
    neighborPixXOnNewMatrix <- getXIndiceFromIndice2D(neightborIndice2D)
    neighborPixYOnNewMatrix <- maxPixY - (getYIndiceFromIndice2D(neightborIndice2D) - minPixY)
    points(neighborPixXOnNewMatrix, neighborPixYOnNewMatrix, 
           col = as.character(neighborSpots$plotColor), 
           cex = spotSize, pch = 16)
    
    for(i in 1:nrow(mRNAToColor)) {
        xIndex <- ceil(i / maxNumY)
        yIndex <- i %% maxNumY
        if (yIndex == 0) {yIndex <- maxNumY}
        text(dx * xIndex, dy * yIndex, mRNAToColor$name[i], cex = 0.7)
        points(dx * xIndex - pointTextXDistance, dy * yIndex, 
               col = as.character(mRNAToColor$plotColor[i]), 
               cex = 0.7, pch = 16)
    }
}


# Plot the raster plot, and circle the refSpot revesreMoved to the plot 
# (g1Plot, for example) as a white circle, while the (g1) neighbors in the plot 
# (g1Plot, for example) as blue circle
# Use the normalTiffMatrix to plot, and add the spots, and add the transRefSpots on top
# as -deltaX, and -deltaY
plotNearbyTiffRaster = function (refSpot, spots, pix_XY, normTiffMatrix, transXYPerPolygon,
                                 isToTranslateRefSpots = TRUE, plot_area_dist = 20, 
                                 color = colorRampPalette(c("black", "green")) (1000),
                                 plotname = "", spotSize = 2, neighborSpotSize = 2) {   
    resX <- dim(normTiffMatrix)[2]
    resY <- dim(normTiffMatrix)[1]
    refIndice2D_int <- spotsToIndice2D(refSpot, pix_XY)
    pixX_int <- getXIndiceFromIndice2D(refIndice2D_int)
    pixY_int <- getYIndiceFromIndice2D(refIndice2D_int)
    pixPlotArea <- round(plot_area_dist / pix_XY)
    minPixX <- max(1, pixX_int - pixPlotArea)
    maxPixX <- min(pixX_int + pixPlotArea, resX)
    minPixY <- max(1, pixY_int - pixPlotArea)
    maxPixY <- min(pixY_int + pixPlotArea, resY)
    
    nearbyMatrix <- normTiffMatrix[minPixY : maxPixY, minPixX : maxPixX]
    rasterImage <- raster(nearbyMatrix, ymn = minPixY-0.5, ymx = maxPixY+0.5, xmn = minPixX-0.5, xmx = maxPixX+0.5) 
    plot(rasterImage, col = color, xlim=c(minPixX-0.5, maxPixX+0.5),ylim=c(minPixY-0.5, maxPixY+0.5), main = plotname)
    
    # plot the ref spots on the tiff plot
    newRefSpot <- refSpot
    if (isToTranslateRefSpots == TRUE) {
        newRefSpot <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transXYPerPolygon)
    }

    refIndice2D <- spotsToIndice2D_not_integer(newRefSpot, pix_XY)
    pixXOnNewMatrix <- getXIndiceFromIndice2D(refIndice2D)
    pixYOnNewMatrix <- maxPixY - (getYIndiceFromIndice2D(refIndice2D) - minPixY)
    points(pixXOnNewMatrix, pixYOnNewMatrix, col = "purple", cex = spotSize)
    points(pixXOnNewMatrix, pixYOnNewMatrix, col = "yellow", cex = spotSize + 0.5)
    
    # plot neighbors spots on the tiff plot
    neighborSpots <- neighbors_exclude_spot(refSpot, spots, plot_area_dist)
    neightborIndice2D <- spotsToIndice2D_not_integer(neighborSpots, pix_XY)
    neighborPixXOnNewMatrix <- getXIndiceFromIndice2D(neightborIndice2D)
    neighborPixYOnNewMatrix <- maxPixY - (getYIndiceFromIndice2D(neightborIndice2D) - minPixY)
    points(neighborPixXOnNewMatrix, neighborPixYOnNewMatrix, col = "white", cex = neighborSpotSize)
    #points(neighborPixXOnNewMatrix, neighborPixYOnNewMatrix, col = as.character(neighborSpots$plotColor), cex = neighborSpotSize, pch = 16) 
}


# Might not be useful. Get the nearby Matrix around indice_2D, and avoid out of the matrix area.
getNearbyMatrix = function(normTiffMatrix, indice2D_int, pix_dist) {
    pixX_int <- getXIndiceFromIndice2D(indice2D_int)
    pixY_int <- getYIndiceFromIndice2D(indice2D_int)
    minPixX <- max(1, pixX_int - pix_dist)
    maxPixX <- min(pixX_int + pix_dist, resX)
    minPixY <- max(1, pixY_int - pix_dist)
    maxPixY <- min(pixY_int + pix_dist, resY)
    nearbyMatrix <- normTiffMatrix[minPixY : maxPixY, minPixX : maxPixX]
    nearbyMatrix
}


# Might not be useful. Get the nearby neighbor indice_2D around indice_2D_int, and 
# avoid out of the matrix area.
getNearbyIndice2D = function(normTiffMatrix, indice2D_int, pix_dist) {
    pixX_int <- getXIndiceFromIndice2D(indice2D_int)
    pixY_int <- getYIndiceFromIndice2D(indice2D_int)
    minPixX <- max(1, pixX_int - pix_dist)
    maxPixX <- min(pixX_int + pix_dist, resX)
    minPixY <- max(1, pixY_int - pix_dist)
    maxPixY <- min(pixY_int + pix_dist, resY)
    grid <- meshgrid(minPixX : maxPixX, minPixY: maxPixY)
    indiceX <- as.vector(grid$x)
    indiceY <- as.vector(grid$y)
    indice2D <- getIndice2DFromXY(indiceX,indiceY)
    indice2D
}


# Move the reference spots by -deltaX and -deltaY according to the transXYPerPolygon index
moveRefSpotsByTransXYPerPolygon = function(refSpots, polygonMatrix, transXYPerPolygon) {
    newRefSpots <- refSpots
    polygonIndexOfRefSpot <- polygonMatrix[spotsToIndice2D(newRefSpots, pix_XY)]
    deltaX <- ifelse(polygonIndexOfRefSpot == 0, 0, transXYPerPolygon$deltaX[polygonIndexOfRefSpot])
    deltaY <- ifelse(polygonIndexOfRefSpot == 0, 0, transXYPerPolygon$deltaY[polygonIndexOfRefSpot])
    newRefSpots <- moveXY(refSpots, - deltaX, - deltaY)
    newRefSpots
}




# Return a indice2D for spots
# Somehow for FQ spots, we need to add 1 pixel in both XY direction.
# Note that both the indice and spots are from top to down and left to right.
# Only when it comes to tiff images labeling is pixel up to down but y_label down to up.
spotsToIndice2D = function (spots, pix_XY) {
    pixX <- round(spots$Pos_X / pix_XY) + 1
    pixY <- round(spots$Pos_Y / pix_XY) + 1
    indice2D <- getIndice2DFromXY(pixX, pixY)
    indice2D
}


# Return a indice 2D for Pox_Xs and Pos_Ys
# Somehow for FQ spots, we need to add 1 pixel in both XY direction.
# Note that both the indice and spots are from top to down and left to right.
# Only when it comes to tiff images labeling is pixel up to down but y_label down to up.
getIndice2DFromPosXY = function(Pos_Xs, Pos_Ys, pix_XY) {
    pixX <- round(Pos_Xs / pix_XY) + 1
    pixY <- round(Pos_Ys / pix_XY) + 1
    indice2D <- getIndice2DFromXY(pixX, pixY)
    indice2D
}


# Return a indice2D for spots, this time don't round the pixel.
# Somehow for FQ spots, we need to add 1 pixel in both XY direction.
spotsToIndice2D_not_integer = function (spots, pix_XY) {
    pixX <- spots$Pos_X / pix_XY + 1
    pixY <- spots$Pos_Y / pix_XY + 1
    indice2D <- getIndice2DFromXY(pixX, pixY)
    indice2D
}


# From the spotFile, get the pixel_XY in uM
read_pixel_xy = function(filename, isToMicron = TRUE) {
    all_content = readLines(filename);
    line_11_to_12 = all_content[11:12];
    skip_1_to_18 = all_content[-1:-18];
    pix_XY = read.csv(textConnection(line_11_to_12), sep = "\t")$Pix.XY;
    if (isToMicron) {
        pix_XY <- pix_XY / 1000
    }
    pix_XY
}


# Read spots from a spotFile, and use tiffFile to add intensity info
# (debug for accurate information)
# Note that the spots Y axis is from top to down, and x axis is from left to right
# The corresponding pixel of Y is y_um / pix_XY + 1
# The corresponding pixel of X is x_um / pix_XY + 1
read_spots = function (filename, color, tifFilename, isToMicron = TRUE, isMore = FALSE) {
    # Load spots of each graph with X-cord, Y-cord, Amp, and Color
    all_content = readLines(filename);
    line_11_to_12 = all_content[11:12];
    skip_1_to_18 = all_content[-1:-18];
    skip_1_to_18 = skip_1_to_18[-length(skip_1_to_18)]  # remove the last line
    pix_XY = read.csv(textConnection(line_11_to_12), sep = "\t")$Pix.XY;
    
    if (length(skip_1_to_18) == 0) {
        spots <- data.frame(Pos_Y = numeric(), 
                           Pos_X = numeric(),
                           AMP = numeric(),
                           RES = numeric(),
                           INT_filt = numeric())
        #maxMinusBgTiffAmp <- NA
        #averageTiffAmp <- NA
        standardDeviationTiffAmp <- NA
    } else {
        spots.raw <- read.csv(textConnection(skip_1_to_18), sep = "\t");
        # remove rows with camera defects of bright single pixels, and NA
        spots.raw <- spots.raw[as.numeric(paste(spots.raw$SigmaX)) > pix_XY, ];
        spots.raw <- na.omit(spots.raw)
        if (nrow(spots.raw) == 0) {
            spots <- data.frame(Pos_Y = numeric(), 
                                Pos_X = numeric(),
                                AMP = numeric(),
                                RES = numeric(),
                                INT_filt = numeric())
        } else {
            # remove columns other than Pos_Y, Pos_X and AMP
            # change Pos_Y and Pox_X unit from nm to um
            spots <- spots.raw[,c("Pos_Y","Pos_X","AMP","RES","INT_filt")]
            spots$Pos_Y <- as.numeric(paste(spots$Pos_Y))
            spots$Pos_X <- as.numeric(paste(spots$Pos_X))
            
            if (isToMicron) {
                pix_XY <- pix_XY / 1000
                spots$Pos_Y <- spots$Pos_Y / 1000
                spots$Pos_X <- spots$Pos_X / 1000
            }
            spots$AMP <- as.numeric(paste(spots$AMP))
            
            
            # averageAMP <- mean(spots$AMP)
            # standardDeviationAMP <- sd(spots$AMP)
            # spots$AMP <- (spots$AMP - averageAMP) / standardDeviationAMP
            
            spots$RES <- as.numeric(paste(spots$RES))
            # averageRES <- mean(spots$RES)
            # standardDeviationRES <- sd(spots$RES)
            # spots$RES <- (spots$RES - averageRES) / standardDeviationRES
            
            spots$INT_filt <- as.numeric(paste(spots$INT_filt))
            
            # add the Color column
            spots$Color = color
            spots
            
            # use the tiff plot to get exact pixel brightness
            tiff = readTIFF(tifFilename, as.is = TRUE)
            #spotPixelX <- round(spots$Pos_X / pix_XY)
            #spotPixelY <- round(spots$Pos_Y / pix_XY)
            #matrixIndex <- cbind(spotPixelY, spotPixelX)
            matrixIndex <- spotsToIndice2D(spots, pix_XY)
            spots$tiffAmp <- tiff[matrixIndex] 
            spots$bgTiffAmp <- getMinNearbyAmp(tiff, matrixIndex, 5)
            spots$minusBgTiffAmp <- spots$tiffAmp - spots$bgTiffAmp
            maxMinusBgTiffAmp <- max(spots$minusBgTiffAmp)
            spots$normalTiffAmp <-  spots$minusBgTiffAmp / maxMinusBgTiffAmp
            
            #averageTiffAmp <- mean(spots$tiffAmp)
            #standardDeviationTiffAmp <- sd(spots$tiffAmp)
            #spots$tiffAMP <- (spots$tiffAMP - averageTiffAMP) / standardDeviationTiffAMP
            
            if (isMore) {
                spots.extra <- spots.raw[,c("Pos_Z", "SigmaX","SigmaY","SigmaZ","BGD")]
                spots$Pos_Z <- as.numeric(paste(spots.extra$Pos_Z))
                spots$SigmaX <- as.numeric(paste(spots.extra$SigmaX))   
                spots$SigmaY <- as.numeric(paste(spots.extra$SigmaY))  
                spots$SigmaZ <- as.numeric(paste(spots.extra$SigmaZ))
                spots$BGD <- as.numeric(paste(spots.extra$BGD))  
                
                if (isToMicron) {
                    spots$Pos_Z <- spots$Pos_Z / 1000
                    spots$SigmaX <- spots$SigmaX / 1000
                    spots$SigmaY <- spots$SigmaY / 1000
                    spots$SigmaZ <- spots$SigmaZ / 1000
                }
            }
        }
        spots$isUsed <- FALSE
    }
    spots
}


# Get the min amplitube within certain radius around many indice of a matrix
# radius is in pixel and the min amplitube can't be 0 in case of sensor of a pixel broken
getMinNearbyAmp = function(tiffMatrix, centerIndice, radius) {
    resY <- dim(tiffMatrix)[1]
    resX <- dim(tiffMatrix)[2]
    centerYs <- getYIndiceFromIndice2D(centerIndice)
    centerXs <- getXIndiceFromIndice2D(centerIndice)
    minAmp <- tiffMatrix[centerIndice] 
    for (dx in -radius : radius) {
        for (dy in -radius : radius) {
            distance <- sqrt(dx^2 + dy^2)
            x <- centerXs + dx
            y <- centerYs + dy
            indice <- centerIndice
            goodIndiceLogic <- distance <= radius & x >= 1 & x <= resX & y >= 1 & y <= resY
            indice[goodIndiceLogic,] <- getIndice2DFromXY(x[goodIndiceLogic], y[goodIndiceLogic])
            amp <- tiffMatrix[indice]
            minAmp <- ifelse(amp < minAmp & amp != 0, amp, minAmp)
        }
    }
    minAmp
}


# Get the indice that have max amplitube within certain radius around certain indice of a matrix
getMaxNearbyIndice = function(tiffMatrix, centerIndice, radius) {
    resY <- dim(tiffMatrix)[1]
    resX <- dim(tiffMatrix)[2]
    centerYs <- getYIndiceFromIndice2D(centerIndice)
    centerXs <- getXIndiceFromIndice2D(centerIndice)
    maxAmp <- tiffMatrix[centerIndice] 
    maxIndice <- centerIndice
    for (dx in -radius : radius) {
        for (dy in -radius : radius) {
            distance <- sqrt(dx^2 + dy^2)
            x <- centerXs + dx
            y <- centerYs + dy
            neighborIndice <- centerIndice
            goodIndiceLogic <- distance <= radius & x >= 1 & x <= resX & y >= 1 & y <= resY
            neighborIndice[goodIndiceLogic,] <- getIndice2DFromXY(x[goodIndiceLogic], y[goodIndiceLogic])
            amp <- tiffMatrix[neighborIndice]
            maxXIndice <- ifelse(amp > maxAmp, getXIndiceFromIndice2D(neighborIndice), getXIndiceFromIndice2D(maxIndice))
            maxYIndice <- ifelse(amp > maxAmp, getYIndiceFromIndice2D(neighborIndice), getYIndiceFromIndice2D(maxIndice))
            maxIndice <- getIndice2DFromXY(maxXIndice, maxYIndice)
            maxAmp <- ifelse(amp > maxAmp, amp, maxAmp)
        }
    }
    maxIndice
}


# Get the max amplitube within certain radius around certain indice of a matrix
getMaxNearbyAmp = function(tiffMatrix, centerIndice, radius) {
    maxIndice <- getMaxNearbyIndice(tiffMatrix, centerIndice, radius)
    maxAmp <- tiffMatrix[maxIndice]
    maxAmp
}


# Get the indice that have max avg amplitube where centers are within certain nearbyRadius 
# around many indice of a matrix and avg is over an avgRadius (in pixel)
getMaxNearbyAvgIndice = function(tiffMatrix, centerIndice, nearbyRadius, avgRadius) {
    resY <- dim(tiffMatrix)[1]
    resX <- dim(tiffMatrix)[2]
    centerYs <- getYIndiceFromIndice2D(centerIndice)
    centerXs <- getXIndiceFromIndice2D(centerIndice)
    maxAvgAmp <- getAvgAmpSurroundingIndice2D(tiffMatrix, centerIndice, avgRadius)
    maxAvgIndice <- centerIndice
    for (dx in -nearbyRadius : nearbyRadius) {
        for (dy in -nearbyRadius : nearbyRadius) {
            distance <- sqrt(dx^2 + dy^2)
            x <- centerXs + dx
            y <- centerYs + dy
            neighborIndice <- centerIndice
            goodIndiceLogic <- distance <= nearbyRadius & x >= 1 & x <= resX & y >= 1 & y <= resY
            neighborIndice[goodIndiceLogic,] <- getIndice2DFromXY(x[goodIndiceLogic], y[goodIndiceLogic])
            avgAmp <- getAvgAmpSurroundingIndice2D(tiffMatrix, neighborIndice, avgRadius)
            maxAvgXIndice <- ifelse(avgAmp > maxAvgAmp, getXIndiceFromIndice2D(neighborIndice), getXIndiceFromIndice2D(maxAvgIndice))
            maxAvgYIndice <- ifelse(avgAmp > maxAvgAmp, getYIndiceFromIndice2D(neighborIndice), getYIndiceFromIndice2D(maxAvgIndice))
            maxAvgIndice <- getIndice2DFromXY(maxAvgXIndice, maxAvgYIndice)
            maxAvgAmp <- ifelse(avgAmp > maxAvgAmp, avgAmp, maxAvgAmp)
        }
    }
    maxAvgIndice
}


# Get the max avg amplitube where centers are within certain nearbyRadius 
# around many indice of a matrix and avg is over an avgRadius
getMaxNearbyAvgAmp = function(tiffMatrix, centerIndice, nearbyRadius, avgRadius) {
    maxAvgIndice <- getMaxNearbyAvgIndice(tiffMatrix, centerIndice, nearbyRadius, avgRadius)
    maxAvgAmp <- getAvgAmpSurroundingIndice2D(tiffMatrix, maxAvgIndice, avgRadius)
    maxAvgAmp
}


# Get the average amp of a vector of indice2D, average over surrounding pixels within avgRadius
# When radius is zero, it just cound the pixel itself, when 1 it counts all pixels that is less
# than 1 pixel distance from the the pixel, etc. 
getAvgAmpSurroundingIndice2D = function(tiffMatrix, indice2D, avgRadius) {
    resY <- dim(tiffMatrix)[1]
    resX <- dim(tiffMatrix)[2]
    centerYs <- getYIndiceFromIndice2D(indice2D)
    centerXs <- getXIndiceFromIndice2D(indice2D)
    count <- numeric(nrow(indice2D))
    amp <- numeric(nrow(indice2D))
    for (dx in -avgRadius : avgRadius) {
        for (dy in -avgRadius : avgRadius) {
            distance <- sqrt(dx^2 + dy^2)
            x <- centerXs + dx
            y <- centerYs + dy
            indice <- getIndice2DFromXY(x, y)
            goodIndiceLogic <- distance <= avgRadius & x >= 1 & x <= resX & y >= 1 & y <= resY
            count <- count + goodIndiceLogic
            amp[goodIndiceLogic] <- amp[goodIndiceLogic] + tiffMatrix[indice[goodIndiceLogic,, drop = FALSE]]
        }
    }
    avgAmp <- amp / count
    avgAmp
}




# Note that the spot has x, y in mum, and y has top as 0 and x has left as 0
# plot the spots, which we get from the spot file.
# plot the polygons_um as well if polygons_um is provided.
# Dont' plot the polygons by set isDrawOutlines to be FALSE
plotSpots = function (spots, spots2 = spotDataFrame(spots), spotSize = 1, spot2Size = 1, 
                      polygons_um, isDrawOutlines = TRUE, isSpotsBlack = FALSE, isSpots2Black = TRUE, 
                      isBackgroundBlack = FALSE, plotname = ""){
    # name <- deparse(substitute(spots))
    p <- ggplot() + ggtitle(plotname) + scale_y_reverse() + labs(x = "Pos_X (um)", y = "Pos_Y (um)")
    # p <- p + xlim(0, max(spots2$Pos_X) + 10)
    if (nrow(spots) != 0) {
        if (isSpotsBlack) {
            p <- p + geom_point(data = spots, aes(x = Pos_X, y = Pos_Y), color = "black", size = I(spotSize))
        } else {
            p <- p + geom_point(data = spots[spots$Color == "red",], aes(x = Pos_X, y = Pos_Y), color = "red", size = I(spotSize))
            p <- p + geom_point(data = spots[spots$Color == "orange",], aes(x = Pos_X, y = Pos_Y), color = "orange", size = I(spotSize))
            p <- p + geom_point(data = spots[spots$Color == "yellow",], aes(x = Pos_X, y = Pos_Y), color = "gold", size = I(spotSize))
            p <- p + geom_point(data = spots[spots$Color == "green",], aes(x = Pos_X, y = Pos_Y), color = "green", size = I(spotSize))
            p <- p + geom_point(data = spots[spots$Color == "gray",], aes(x = Pos_X, y = Pos_Y), color = "green", size = I(spotSize))
        }
    }
    
    if (isSpots2Black) {
        p <- p + geom_point(data = spots2, aes(x = Pos_X, y = Pos_Y), color = "black", size = I(spot2Size))
    } else {
        p <- p + geom_point(data = spots2[spots2$Color == "red",], aes(x = Pos_X, y = Pos_Y), color = "red", size = I(spot2Size))
        p <- p + geom_point(data = spots2[spots2$Color == "orange",], aes(x = Pos_X, y = Pos_Y), color = "orange", size = I(spot2Size))
        p <- p + geom_point(data = spots2[spots2$Color == "yellow",], aes(x = Pos_X, y = Pos_Y), color = "gold", size = I(spot2Size))
        p <- p + geom_point(data = spots2[spots2$Color == "green",], aes(x = Pos_X, y = Pos_Y), color = "green", size = I(spot2Size))
        p <- p + geom_point(data = spots2[spots2$Color == "gray",], aes(x = Pos_X, y = Pos_Y), color = "green", size = I(spot2Size))
    }
    
    if (isBackgroundBlack) {
        p <- p + theme(
            panel.background = element_rect(fill = "black",
                                            colour = "black",
                                            size = 0.5, linetype = "solid"),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                            colour = "black"), 
            panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                            colour = "black")
        )
    }
    
    if (isDrawOutlines == TRUE) {
        for (i in seq_along(polygons_um)) {
            polygon_um <- polygons_um[[i]]
            polygon.frame <- data.frame(x = polygon_um[,1], y = polygon_um[,2])
            p <- p + geom_polygon(data = polygon.frame, aes(x = x, y = y), alpha = 0.5)
        }
    }
    p
}

# given spots, return their coords as Pos_X, and Pos_Y in um
coord = function (spots) {
    coords <- cbind(spots$Pos_X, spots$Pos_Y)
    coords
}

# choose spots within a boundary 
enclosedSpots = function (spots, polygon_um) {
    coords <- coord(spots)
    spotsWithinPolygon <- spots[in.out(polygon_um, coords),]
    spotsWithinPolygon
}

# get the indice of spots within a boundary
enclosedSpotsIndice = function(spots, polygon_um) {
    coords <- coord(spots)
    indice <- which(in.out(polygon_um, coords))
    indice
}

# choose spots within list of boundaries
enclosedSpotsWithinBoundaries = function (spots, polygons_um) {
    spot_enclosed <- spotDataFrame(spots)
    for (i in seq_along(polygons_um)) {
        spot_enclosed <- rbind(spot_enclosed, enclosedSpots(spots, polygons_um[[i]]))
    }
    spot_enclosed
}

# choose the index of the spot within list of boundaries
enclosedSpotsIndiceWithinBoundaries = function (spots, polygons_um) {
    enclosedIndice <- NULL
    for (i in seq_along(polygons_um)) {
        enclosedIndice <- c(enclosedIndice, enclosedSpotsIndice(spots, polygons_um[[i]]))
    }
    enclosedIndice
}


# Plot spots with refSpots within a boundary polygon_um
plotBoundedSpots = function(spots, spots2, polygon_um, spotSize = 1, plotname = "") {
    boundedSpots <- spots[in.out(polygon_um, coord(spots)),]
    boundedSpots2 <- spots2[in.out(polygon_um, coord(spots2)),]
    p <- plotSpots(boundedSpots, boundedSpots2, spotSize = spotSize, isDrawOutlines = FALSE, plotname = plotname)
    p
}


# Plot a spot and its nearby spots within 2 plots, within +- plot_area_dist (um)
# Can give the plot a title.
plotNearbySpots = function(spot, spots, spots2, plot_area_dist = 200, plotname = "") {
    subSpots <- neighborhood(spot, spots, plot_area_dist)
    subSpots2 <- neighborhood(spot, spots2, plot_area_dist)
    p <- plotSpots(subSpots, subSpots2, spotSize = 2, spot2Size = 1, isDrawOutlines = FALSE)
    p <- p + geom_point(data = spot, aes(x = Pos_X, y = Pos_Y), color = "red", size = I(5), shape = 1)
    p <- p + ggtitle(plotname)
    p
}


# Plot the area near the spot with spot as a big red dot, neighbors as big blue dots,
# potential RefSpots as big green dots.
# can label indice of each spots with rowname
plotNearbySpotsWithNeighbors = function(spot, spots, refSpots, plot_area_dist = 20, 
                           neighbor_dist = 10, ini_cospot_dist = 20, withLabel = TRUE) {
    subSpots <- neighborhood(spot, spots, plot_area_dist)
    subRefSpots <- neighborhood(spot, refSpots, plot_area_dist)
    neighbors <- neighbors_exclude_spot(spot, subSpots, neighbor_dist)
    coSpotsOnRef <- neighborhood(spot, subRefSpots, ini_cospot_dist)
    neighbors$neighborIndex <- rownames(neighbors)
    coSpotsOnRef$coRefSpotIndex <- rownames(coSpotsOnRef)
    p <- plotSpots(subSpots, subRefSpots)
    p <- p + geom_point(data = coSpotsOnRef, aes(x = Pos_X, y = Pos_Y), color = "green", size = I(3))
    p <- p + geom_point(data = neighbors, aes(x = Pos_X, y = Pos_Y), color = "blue", size = I(3)) 
    p <- p + geom_point(data = spot, aes(x = Pos_X, y = Pos_Y), color = "red", size = I(3))
    if (withLabel) {
        p <- p + geom_text(data = coSpotsOnRef, aes(x = Pos_X, y = Pos_Y, label = coRefSpotIndex), 
                           color = "black", check_overlap = TRUE, hjust = -0.5, size = 3)
        p <- p + geom_text(data = neighbors, aes(x = Pos_X, y = Pos_Y, label = neighborIndex), 
                           color = "blue", check_overlap = TRUE, hjust = -0.5, size = 3)
    }
    p
}


# Get the spots within +- neighbor_dist (um) in x, y direction of the spot 
# (not excluding the spot itself.)
neighborhood = function (spot, spots, neighbor_dist = 10) {
    neighborhood <- spots %>% filter(Pos_X > spot$Pos_X - neighbor_dist &
                                          Pos_X < spot$Pos_X + neighbor_dist &
                                          Pos_Y > spot$Pos_Y - neighbor_dist &
                                          Pos_Y < spot$Pos_Y + neighbor_dist)
    neighborhood
}


# Get the spots within +- neighbor_dist (um) in x, y direction of the spot 
# (excluding the spot itself.)
neighbors_exclude_spot = function (spot, spots, neighbor_dist) {
    neighbors <- neighborhood(spot, spots, neighbor_dist)
    neighbors <- neighbors %>% filter(!(Pos_X == spot$Pos_X & Pos_Y == spot$Pos_Y))
    neighbors
}


