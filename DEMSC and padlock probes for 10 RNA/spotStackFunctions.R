source("/Users/jackie/Downloads/SGTC microscope/2016-08-11-10mRNAs vs 4 exon mRNAs-37C-45C Amp 3rd time/4 exon mRNAs/spotPlotFunctions.R")
source("/Users/jackie/Downloads/SGTC microscope/2016-08-11-10mRNAs vs 4 exon mRNAs-37C-45C Amp 3rd time/4 exon mRNAs/spotRegisterFunctions.R")

# Gaussian fit spots
# The point spread function (PSF) describes the response of an image to a point source.
# The point source forms an Airy pattern at the focus plane of a lens.
# Airy radius = Rayleigh diffraction limit: r_airy = 0.61 * lambda / NA =
# 0.61 * lambda / (n * sin(theta)) =  0.61 * lambda / (n * lens_radius / focal_length)
# numerical aperture: NA = n * sin(theta) ~ n * theta = n * lens_radius / focal_length
# Airy function can be approximated by Gaussian function g = I * exp( (-(r-ro)^2 / (2*sigma_r^2) )  + bg
# where Gaussian standard deviation: sigma_r ~ 0.42 * lambda * N
# N is f-number of the lens = focal_length / lens_diameter = 1 / (2 * NA) 
# sigma_r ~ 0.42 * lambda * focal_length / lens_diameter = 0.42 * lambda / (2 * NA) ~ 1/3 r_airy
# My 20X lens has NA = 0.4, lambda ~ 519nm, 573nm, 617nm and 669nm for Alexa-488, Alexa 546, Alex 594 and Atto-647,
# so sigma_r ~ 0.42 * (519 | 573 | 617 | 669) / 0.8 ~ 272nm, 301nm, 324nm, 351nm
# each pixel = 330nm, so sigma_r ~ 0.85 - 1.09 pixels and r_airy ~ 2.45 - 3.16 pixels

AiryDisk2D = function(amp, x0, y0, nu, bg, matrixCoordX, matrixCoordY) {
    x <- sqrt((matrixCoordX - x0)^2 + (matrixCoordY - y0)^2)
    AiryDisk <- amp * (2 * besselJ(x, nu) / x)^2 + bg
    AiryDisk[is.nan(AiryDisk)] <- amp + bg 
    AiryDisk
}


# givin the initial centerSpot, the tiffMatrix, the initial parameters
# Fit AiryDisk2D around the spot's indice2D, within nearby_distance (um)
# If the nls fit doesn't work, don't print out error message, instead,
# return NA
# Note that this Gaussian doesn't check whether it's a good Gaussian fit or not ...
fitAiryDisk2D = function(spot, tiffMatrix, fit_distance_um, nearby_distance_um = 0.3, 
                         lowerRadiusBound = 0.5, upperRadiusBound = 1.5, trace = FALSE) {
    # get the matrixCoordX and matrixCoordY in um
    centerIndice <- spotsToIndice2D(spot, pix_XY)
    fit_distance <- floor(fit_distance_um / pix_XY)
    nearbyIndice2D <- getNearbyIndice2D(tiffMatrix, centerIndice, fit_distance)
    matrixCoordX <- getXCoordsFromIndice2D(nearbyIndice2D)
    matrixCoordY <- getYCoordsFromIndice2D(nearbyIndice2D)
    
    # get the AiryDisk data of neighbor pixels, unit is intensity
    AiryDisk <- tiffMatrix[nearbyIndice2D]
    
    # guess the init parameters
    maxAvgIndice <- getMaxNearbyAvgIndice(tiffMatrix, centerIndice, nearbyRadius = 2, avgRadius = 2)
    maxAvgAmp <- getMaxNearbyAvgAmp(tiffMatrix, centerIndice, nearbyRadius = 2, avgRadius = 2)
    bg <- getMinNearbyAmp(tiffMatrix, maxAvgIndice, 5)
    amp <- maxAvgAmp - bg
    x0 <- getXCoordsFromIndice2D(maxAvgIndice)
    y0 <- getYCoordsFromIndice2D(maxAvgIndice)
    parametersInt <- c(amp, x0, y0, 1, bg)
    
    lowerParameterBound <- c(0, x0-nearby_distance_um, y0-nearby_distance_um, lowerRadiusBound, 0)
    upperParameterBound <- c(Inf, x0+nearby_distance_um, y0+nearby_distance_um, upperRadiusBound, Inf)
    
    # fit with nls
    fit <- NA
    try(fit <- nls(AiryDisk ~ AiryDisk2D(amp, x0, y0, nu, bg, matrixCoordX, matrixCoordY), 
                   start = list(amp = parametersInt[1], x0 = parametersInt[2], y0 = parametersInt[3],
                                nu = parametersInt[4], bg = parametersInt[5]), 
                   lower = list(amp = lowerParameterBound[1], x0 = lowerParameterBound[2], y0 = lowerParameterBound[3],
                                nu = lowerParameterBound[4], bg = lowerParameterBound[5]),
                   upper = list(amp = upperParameterBound[1], x0 = upperParameterBound[2], y0 = upperParameterBound[3],
                                nu = upperParameterBound[4], bg = upperParameterBound[5]), 
                   algorithm = "port", trace = trace),
        silent = TRUE) 
    fit
}

# g = amp * exp( -((x-x0)^2 + (y-y0)^2) / (2*sigma_r^2) )  + bg
# where parametar = [amp, x0, y0, sigma_r, bg]
# and matrixCoordX(um) are the xCoords matrix indicating the xCoords of each pixel
# and matrixCoordY(um) are the yCoords matrix indicating the yCoords of each pixel

# Gaussian2D = function(parameters, matrixCoordX, matrixCoordY) {
#     Gaussian <- parameters[1] * exp(-((matrixCoordX - parameters[2])^2 + 
#                 (matrixCoordY - parameters[3])^2) / (2 * parameters[4]^2)) +
#                 parameters[5]
#     Gaussian
# }
Gaussian2D = function(amp, x0, y0, sigma_r, bg, matrixCoordX, matrixCoordY) {
    Gaussian <- amp * exp(-((matrixCoordX - x0)^2 + 
                    (matrixCoordY - y0)^2) / (2 * sigma_r^2)) + bg
    Gaussian
}

# Gaussian3D, only used for FishQuant fit data
Gaussian3D = function(amp, x0, y0, z0, sigma_r, sigma_z, bg, matrixCoordX, matrixCoordY) {
    Gaussian <- amp * exp(-((matrixCoordX - x0)^2 + 
                                (matrixCoordY - y0)^2) / (2 * sigma_r^2) - 
                              z0^2 / (2 * sigma_z^2)) + bg
    Gaussian
}
    
# g = amp1 * exp( -((x-x01)^2 + (y-y01)^2) / (2*sigma_r1^2) ) + 
# amp2 * exp( -((x-x02)^2 + (y-y02)^2) / (2*sigma_r2^2) ) + bg
# where parametar = [amp1, x01, y01, sigma_r1, amp2, x02, y02, sigma_r2, bg]
# and matrixCoordX(um) are the xCoords matrix indicating the xCoords of each pixel
# and matrixCoordY(um) are the yCoords matrix indicating the yCoords of each pixel
Gaussian2DFor2 = function(parameters, matrixCoordX, matrixCoordY) {
    Gaussian <- parameters[1] * exp(-((matrixCoordX - parameters[2])^2 + 
                (matrixCoordY - parameters[3])^2) / (2 * parameters[4]^2)) +
                parameters[5] * exp(-((matrixCoordX - parameters[6])^2 + 
                (matrixCoordY - parameters[7])^2) / (2 * parameters[8]^2)) +
                parameters[9]
    Gaussian
}


# Gaussian2DFor2 = function(amp1, x1, y1, sigma1, amp2, x2, y2, sigma2, bg, matrixCoordX, matrixCoordY) {
#     Gaussian <- amp1 * exp(-((matrixCoordX - x1)^2 + 
#                                           (matrixCoordY - y1)^2) / (2 * sigma1^2)) +
#         amp2 * exp(-((matrixCoordX - x2)^2 + 
#                                   (matrixCoordY - y2)^2) / (2 * sigma2^2)) +
#         bg
#     Gaussian
# }




# g = amp1 * exp( -((x-x01)^2 + (y-y01)^2) / (2*sigma_r1^2) ) + 
# amp2 * exp( -((x-x02)^2 + (y-y02)^2) / (2*sigma_r2^2) ) + 
# amp3 * exp( -((x-x03)^2 + (y-y03)^2) / (2*sigma_r3^2) ) + bg
# where parametar = [amp1, x01, y01, sigma_r1, amp2, x02, y02, sigma_r2, amp3, x03, y03, sigma_r3, bg]
# and matrixCoordX(um) are the xCoords matrix indicating the xCoords of each pixel
# and matrixCoordY(um) are the yCoords matrix indicating the yCoords of each pixel
Gaussian2DFor3 = function(parameters, matrixCoordX, matrixCoordY) {
    Gaussian <- parameters[1] * exp(-((matrixCoordX - parameters[2])^2 + 
                (matrixCoordY - parameters[3])^2) / (2 * parameters[4]^2)) +
        parameters[5] * exp(-((matrixCoordX - parameters[6])^2 + 
                (matrixCoordY - parameters[7])^2) / (2 * parameters[8]^2)) +
        parameters[9] * exp(-((matrixCoordX - parameters[10])^2 + 
                (matrixCoordY - parameters[11])^2) / (2 * parameters[12]^2)) +
        parameters[13]
    Gaussian
}


# The sum of resdiual squares (RSS) where residue is the Gaussian2D data minus
# the mathematical value calculated from parameters and x, y coordinates
# Where Gaussian, matrixCoordX and matrixCoordY are predefined number values
# I use Gaussian2DResidualSqrSum for optim fit 
# i.e. optim(par = parametersInt, fn = Gaussian2DResidualSqrSum, method = "BFGS")$par
# instead of Gaussian2D for nls fit
# i.e. fit = nls(Gaussian ~ Gaussian2D(parameters, matrixCoordX, matrixCoordY), start = list(parameters = parametersInt), trace = TRUE)
# because nls can't fit artificial clean signal and it needs to specify the
# initial parameters closer to the real value.
Gaussian2DResidualSqrSum <- function(data, parameters) {
    with(data, sum((Gaussian - Gaussian2D(parameters, matrixCoordX, matrixCoordY))^2))
}


# givin the initial centerSpot, the tiffMatrix, the initial parameters
# Fit Gaussian2D around the spot's indice2D, within nearby_distance (um)
# If the nls fit doesn't work, don't print out error message, instead,
# return NA
# Note that this Gaussian doesn't check whether it's a good Gaussian fit or not ...
fitGaussian2D = function(spot, tiffMatrix, fit_distance_um, nearby_distance_um = 0.3, 
                         lowerSigmaBound = 0.25, upperSigmaBound = 0.5, trace = FALSE) {
    # get the matrixCoordX and matrixCoordY in um
    centerIndice <- spotsToIndice2D(spot, pix_XY)
    fit_distance <- floor(fit_distance_um / pix_XY)
    nearbyIndice2D <- getNearbyIndice2D(tiffMatrix, centerIndice, fit_distance)
    matrixCoordX <- getXCoordsFromIndice2D(nearbyIndice2D)
    matrixCoordY <- getYCoordsFromIndice2D(nearbyIndice2D)

    # get the Gaussian data of neighbor pixels, unit is intensity
    Gaussian <- tiffMatrix[nearbyIndice2D]
    
    # guess the init parameters
    nearbyRadius <- nearby_distance_um # 2
    avgRadius <- fit_distance_um # 2
    maxAvgIndice <- getMaxNearbyAvgIndice(tiffMatrix, centerIndice, nearbyRadius, avgRadius)
    maxAvgAmp <- getMaxNearbyAvgAmp(tiffMatrix, centerIndice, nearbyRadius, avgRadius)
    bg <- getMinNearbyAmp(tiffMatrix, maxAvgIndice, 5)
    amp <- maxAvgAmp - bg
    x0 <- getXCoordsFromIndice2D(maxAvgIndice)
    y0 <- getYCoordsFromIndice2D(maxAvgIndice)
    parametersInt <- c(amp, x0, y0, 0.3, bg)
    
    lowerParameterBound <- c(0, x0-nearby_distance_um, y0-nearby_distance_um, lowerSigmaBound, 0)
    upperParameterBound <- c(Inf, x0+nearby_distance_um, y0+nearby_distance_um, upperSigmaBound, Inf)
    
    # fit with nls
    fit <- NA
    try(fit <- nls(Gaussian ~ Gaussian2D(amp, x0, y0, sigma_r, bg, matrixCoordX, matrixCoordY), 
                    start = list(amp = parametersInt[1], x0 = parametersInt[2], y0 = parametersInt[3],
                                 sigma_r = parametersInt[4], bg = parametersInt[5]), 
                    lower = list(amp = lowerParameterBound[1], x0 = lowerParameterBound[2], y0 = lowerParameterBound[3],
                                 sigma_r = lowerParameterBound[4], bg = lowerParameterBound[5]),
                    upper = list(amp = upperParameterBound[1], x0 = upperParameterBound[2], y0 = upperParameterBound[3],
                                 sigma_r = upperParameterBound[4], bg = upperParameterBound[5]), 
                    algorithm = "port", trace = trace),
        silent = TRUE) 
    fit
}


# givin the centerSpot, the tiffMatrix, the initial parameters
# Fit Gaussian2D at the spot's indice2D
# If the nls fit doesn't work, don't print out error message, instead,
# return NA
# Note that this Gaussian doesn't check whether it's a good Gaussian fit or not ...
fitGaussian2DFixCenter = function(spot, tiffMatrix, fit_distance_um,
                         lowerSigmaBound = 0.25, upperSigmaBound = 0.5, trace = FALSE) {
    
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
    bg <- getMinNearbyAmp(tiffMatrix, centerIndice, 5)
    amp <- centerAmp - bg
    parametersInt <- c(amp, 0.3, bg)
    
    lowerParameterBound <- c(0, lowerSigmaBound, 0)
    upperParameterBound <- c(Inf, upperSigmaBound, Inf)
    
    # fit with nls
    fit <- NA
    try(fit <- nls(Gaussian ~ Gaussian2D(amp, x0, y0, sigma_r, bg, matrixCoordX, matrixCoordY), 
                   start = list(amp = parametersInt[1], sigma_r = parametersInt[2], 
                                bg = parametersInt[3]), 
                   lower = list(amp = lowerParameterBound[1], sigma_r = lowerParameterBound[2], 
                                bg = lowerParameterBound[3]),
                   upper = list(amp = upperParameterBound[1], sigma_r = upperParameterBound[2], 
                                bg = upperParameterBound[3]), 
                   control = nls.control(maxiter = 100), algorithm = "port", na.action=na.exclude,trace = trace),
        silent = FALSE) 
    fit
}


# givin the initial centerSpot, the tiffMatrix, the initial parameters
# Fit Gaussian2D around the spot's indice2D, within nearby_distance (um)
# If the nls fit doesn't work, don't print out error message, instead,
# return NA
# Note that this Gaussian doesn't check whether it's a good Gaussian fit or not ...
fitGaussian2DFixSigma = function(spot, tiffMatrix, fit_distance_um, sigma_r,
                                 nearby_distance_um = 0.3, trace = FALSE) {
    # get the matrixCoordX and matrixCoordY in um
    centerIndice <- spotsToIndice2D(spot, pix_XY)
    fit_distance <- floor(fit_distance_um / pix_XY)
    nearbyIndice2D <- getNearbyIndice2D(tiffMatrix, centerIndice, fit_distance)
    matrixCoordX <- getXCoordsFromIndice2D(nearbyIndice2D)
    matrixCoordY <- getYCoordsFromIndice2D(nearbyIndice2D)
    
    # get the Gaussian data of neighbor pixels, unit is intensity
    Gaussian <- tiffMatrix[nearbyIndice2D]
    
    # guess the init parameters
    maxAvgIndice <- getMaxNearbyAvgIndice(tiffMatrix, centerIndice, nearbyRadius = 2, avgRadius = 2)
    maxAvgAmp <- getMaxNearbyAvgAmp(tiffMatrix, centerIndice, nearbyRadius = 2, avgRadius = 2)
    bg <- getMinNearbyAmp(tiffMatrix, maxAvgIndice, 5)
    amp <- maxAvgAmp - bg
    x0 <- getXCoordsFromIndice2D(maxAvgIndice)
    y0 <- getYCoordsFromIndice2D(maxAvgIndice)
    parametersInt <- c(amp, x0, y0, bg)
    
    lowerParameterBound <- c(0, x0-nearby_distance_um, y0-nearby_distance_um, 0)
    upperParameterBound <- c(Inf, x0+nearby_distance_um, y0+nearby_distance_um, Inf)
    
    # fit with nls
    fit <- NA
    try(fit <- nls(Gaussian ~ Gaussian2D(amp, x0, y0, sigma_r, bg, matrixCoordX, matrixCoordY), 
                   start = list(amp = parametersInt[1], x0 = parametersInt[2], y0 = parametersInt[3],
                                bg = parametersInt[4]), 
                   lower = list(amp = lowerParameterBound[1], x0 = lowerParameterBound[2], 
                                y0 = lowerParameterBound[3], bg = lowerParameterBound[4]),
                   upper = list(amp = upperParameterBound[1], x0 = upperParameterBound[2], 
                                y0 = upperParameterBound[3], bg = upperParameterBound[4]), 
                   algorithm = "port", trace = trace),
        silent = TRUE) 
    fit
}


# givin the centerSpot, the tiffMatrix, the initial parameters
# Fit Gaussian2D at the spot's indice2D
# If the nls fit doesn't work, don't print out error message, instead,
# return NA
# Note that this Gaussian doesn't check whether it's a good Gaussian fit or not ...
fitGaussian2DFixCenterAndSigma = function(spot, tiffMatrix, fit_distance_um,
                                  sigma_r, trace = FALSE) {
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
    bg <- getMinNearbyAmp(tiffMatrix, centerIndice, 5)
    amp <- centerAmp - bg
    parametersInt <- c(amp, bg)
    
    lowerParameterBound <- c(0, 0)
    upperParameterBound <- c(Inf, Inf)
    
    # fit with nls
    fit <- NA
    try(fit <- nls(Gaussian ~ Gaussian2D(amp, x0, y0, sigma_r, bg, matrixCoordX, matrixCoordY), 
                   start = list(amp = parametersInt[1], bg = parametersInt[2]), 
                   lower = list(amp = lowerParameterBound[1], bg = lowerParameterBound[2]),
                   upper = list(amp = upperParameterBound[1], bg = upperParameterBound[2]), 
                   control = nls.control(maxiter = 100), algorithm = "port", na.action=na.exclude,trace = trace),
        silent = FALSE) 
    fit
}


# Givin the initial centerSpot, the tiffMatrix, the initial parameters
# Fit Gaussian2DFor2 around the spot's indice2D, within nearby_distance (um)
# If the nls fit doesn't work, don't print out error message, instead,
# return NA
fitGaussian2DFor2From1Spot = function(spot, tiffMatrix, nearby_distance_um = 2, trace = FALSE) {
    # get the matrixCoordX and matrixCoordY in um
    centerIndice <- spotsToIndice2D(spot, pix_XY)
    nearby_distance <- floor(nearby_distance_um / pix_XY)
    nearbyIndice2D <- getNearbyIndice2D(tiffMatrix, centerIndice, nearby_distance)
    matrixCoordX <- getXCoordsFromIndice2D(nearbyIndice2D)
    matrixCoordY <- getYCoordsFromIndice2D(nearbyIndice2D)
    
    # get the Gaussian data of neighbor pixels, unit is intensity
    Gaussian <- tiffMatrix[nearbyIndice2D]
    
    # guess the init parameters
    maxAvgIndice <- getMaxNearbyAvgIndice(tiffMatrix, centerIndice, nearbyRadius = 2, avgRadius = 2)
    maxAvgAmp <- getMaxNearbyAvgAmp(tiffMatrix, centerIndice, nearbyRadius = 2, avgRadius = 2)
    bg <- getMinNearbyAmp(tiffMatrix, maxAvgIndice, 5)
    amp <- maxAvgAmp - bg
    
    # change 
    #x0 <- getXCoordsFromIndice2D(maxAvgIndice)
    #y0 <- getYCoordsFromIndice2D(maxAvgIndice)
    #parametersInt <- c(amp/2, x0, y0, 0.3, amp/2, x0-1, y0-0.5, 0.3, bg)
    twoSpotsDistance <- 1
    twoSpotPositions <- getPositionsTwoGaussianSpotsPositionsFromOneSpot(spot, tiffMatrix, 
                                                                         twoSpotsDistance, pix_XY)
    x1 <- twoSpotPositions$pos_X1
    y1 <- twoSpotPositions$pos_Y1
    x2 <- twoSpotPositions$pos_X2
    y2 <- twoSpotPositions$pos_Y2
    parametersInt <- c(amp/2, x1, y1, 0.3, amp/2, x2, y2, 0.3, bg)
    lowerParameterBound <- c(0, x1-nearby_distance_um, y1-nearby_distance_um, 0.25, 
                             0, x2-nearby_distance_um, y2-nearby_distance_um, 0.25, 0)
    upperParameterBound <- c(Inf, x1+nearby_distance_um, y1+nearby_distance_um, 0.5, 
                             Inf, x2+nearby_distance_um, y2+nearby_distance_um, 0.5, Inf)

    # fit with nls
    fit <- NA
    try(fit <- nls(Gaussian ~ Gaussian2DFor2(amp1, x1, y1, sigma1, amp2, x2, y2, sigma2, bg, matrixCoordX, matrixCoordY), 
                   start = list(amp1 = parametersInt[1], x1 = parametersInt[2], y1 = parametersInt[3],
                                sigma1 = parametersInt[4], amp2 = parametersInt[5], x2 = parametersInt[6],
                                y2 = parametersInt[7], sigma2 = parametersInt[8], bg = parametersInt[9]),
                   lower = list(amp1 = lowerParameterBound[1], x1 = lowerParameterBound[2], y1 = lowerParameterBound[3],
                                sigma1 = lowerParameterBound[4], amp2 = lowerParameterBound[5], x2 = lowerParameterBound[6],
                                y2 = lowerParameterBound[7], sigma2 = lowerParameterBound[8], bg = lowerParameterBound[9]), 
                   upper = list(amp1 = upperParameterBound[1], x1 = upperParameterBound[2], y1 = upperParameterBound[3],
                                sigma1 = upperParameterBound[4], amp2 = upperParameterBound[5], x2 = upperParameterBound[6],
                                y2 = upperParameterBound[7], sigma2 = upperParameterBound[8], bg = upperParameterBound[9]), 
                   control = nls.control(maxiter = 100), algorithm = "port", na.action=na.exclude,
                   trace = trace))#, silent = TRUE) 
    fit
}


# FQ detected one spot, but it is actually composed of two spots, whose centers are
# twoSpotDistance away from each other, we want to get the x, y positions of the two spots
# with unit of pixel
getPositionsTwoGaussianSpotsPositionsFromOneSpot = function(spot, tiffMatrix, twoSpotDistance_um, pix_XY) {
    pos_X <- spot$Pos_X
    pos_Y <- spot$Pos_Y
    
    # direction of of two spots: vertical, right-tile 30 degree, right tilt 45 degree, right-filt 60 degree, 
    #                            horizontal, right-tilt 120 degree, right-tilt 135 degree, right-tilt 150 degree)
    directionsX1Y1X2Y2 <- list(c(0,1,0,-1), c(0.5,sqrt(3)/2,-0.5,-sqrt(3)),
                               c(sqrt(2),sqrt(2),-sqrt(2),-sqrt(2)), c(sqrt(3)/2, 0.5,-sqrt(3)/2,-0.5),
                               c(1,0,-1,0), c(sqrt(3)/2,-0.5,-sqrt(3)/2,0.5),
                               c(sqrt(2),-sqrt(2),-sqrt(2),sqrt(2)),  c(0.5,-sqrt(3)/2,-0.5,sqrt(3)))
    
    current.Pos_X1 <- 0
    current.Pos_Y1 <- 0
    current.Pos_X2 <- 0
    current.Pos_Y2 <- 0
    current.avgAmp <- 0
    for (i in seq_along(directionsX1Y1X2Y2)) {
        direction <- directionsX1Y1X2Y2[[i]]
        pos_X1 <- pos_X + direction[1] * twoSpotDistance_um/2
        pos_Y1 <- pos_Y + direction[2] * twoSpotDistance_um/2
        pos_X2 <- pos_X + direction[3] * twoSpotDistance_um/2
        pos_Y2 <- pos_Y + direction[4] * twoSpotDistance_um/2
        amp1 <- tiffMatrix[getIndice2DFromPosXY(pos_X1, pos_Y1, pix_XY)]
        amp2 <- tiffMatrix[getIndice2DFromPosXY(pos_X2, pos_Y2, pix_XY)]
        avgAmp <- (amp1 + amp2)/2
        if (avgAmp > current.avgAmp) {
            current.Pos_X1 <- pos_X1
            current.Pos_Y1 <- pos_Y1
            current.Pos_X2 <- pos_X2 
            current.Pos_Y2 <- pos_Y2
            current.avgAmp <- avgAmp
        }
    }
    twoSpotPositions <- data.frame(pos_X1 = current.Pos_X1, pos_Y1 = current.Pos_Y1, 
                                   pos_X2 = current.Pos_X2, pos_Y2 = current.Pos_Y2)
    twoSpotPositions
}

# Givin the initial 2 Spots, the tiffMatrix, the initial parameters
# Fit Gaussian2DFor2 around the spot's indice2D, within nearby_distance (um)
# If the nls fit doesn't work, don't print out error message, instead,
# return NA
fitGaussian2DFor2 = function(spot1, spot2, tiffMatrix, nearby_distance_um = 0.3, trace = FALSE) {
    # get the matrixCoordX and matrixCoordY in um
    spot1Indice <- spotsToIndice2D(spot1, pix_XY)
    spot2Indice <- spotsToIndice2D(spot2, pix_XY)
    centerIndice <- (spot1Indice + spot2Indice) / 2
    
    nearby_distance <- floor(nearby_distance_um / pix_XY)
    nearbyIndice2D <- getNearbyIndice2D(tiffMatrix, centerIndice, nearby_distance)
    matrixCoordX <- getXCoordsFromIndice2D(nearbyIndice2D)
    matrixCoordY <- getYCoordsFromIndice2D(nearbyIndice2D)
    
    # get the Gaussian data of neighbor pixels, unit is intensity
    Gaussian <- tiffMatrix[nearbyIndice2D]
    
    # guess the init parameters
    maxAvgIndice <- getMaxNearbyAvgIndice(tiffMatrix, centerIndice, nearbyRadius = 2, avgRadius = 2)
    maxAvgAmp <- getMaxNearbyAvgAmp(tiffMatrix, centerIndice, nearbyRadius = 2, avgRadius = 2)
    bg <- getMinNearbyAmp(tiffMatrix, maxAvgIndice, 5)
    amp <- maxAvgAmp - bg
    x1 <- getXCoordsFromIndice2D(spot1Indice)
    y1 <- getYCoordsFromIndice2D(spot1Indice)
    x2 <- getXCoordsFromIndice2D(spot2Indice)
    y2 <- getYCoordsFromIndice2D(spot2Indice)
    parametersInt <- c(amp/2, x1, y1, 0.3, amp/2, x2, y2, 0.3, bg)
    å
    # fit with nls
    fit <- NA
    try(fit <- nls(Gaussian ~ Gaussian2DFor2(parameters, matrixCoordX, matrixCoordY), 
                   start = list(parameters = parametersInt), trace = trace))#, silent = TRUE) 
    fit
}


# Givin the initial centerSpot, the tiffMatrix, the initial parameters
# Fit Gaussian2DFor3 around the spot's indice2D, within nearby_distance (um)
# If the nls fit doesn't work, don't print out error message, instead,
# return NA
fitGaussian2DFor3 = function(spot, tiffMatrix, nearby_distance_um, trace = FALSE) {
    # get the matrixCoordX and matrixCoordY in um
    centerIndice <- spotsToIndice2D(spot, pix_XY)
    nearby_distance <- floor(nearby_distance_um / pix_XY)
    nearbyIndice2D <- getNearbyIndice2D(tiffMatrix, centerIndice, nearby_distance)
    matrixCoordX <- getXCoordsFromIndice2D(nearbyIndice2D)
    matrixCoordY <- getYCoordsFromIndice2D(nearbyIndice2D)
    
    # get the Gaussian data of neighbor pixels, unit is intensity
    Gaussian <- tiffMatrix[nearbyIndice2D]
    
    # guess the init parameters
    maxAvgIndice <- getMaxNearbyAvgIndice(tiffMatrix, centerIndice, nearbyRadius = 2, avgRadius = 2)
    maxAvgAmp <- getMaxNearbyAvgAmp(tiffMatrix, centerIndice, nearbyRadius = 2, avgRadius = 2)
    bg <- getMinNearbyAmp(tiffMatrix, maxAvgIndice, 5)
    amp <- maxAvgAmp - bg
    x0 <- getXCoordsFromIndice2D(maxAvgIndice)
    y0 <- getYCoordsFromIndice2D(maxAvgIndice)
    parametersInt <- c(amp, x0, y0, 0.3, amp, x0-1, y0-0.5, 0.3, amp, x0, y0+1, 0.3, bg)
    
    # fit with nls
    fit <- NA
    try(fit <- nls(Gaussian ~ Gaussian2DFor2(parameters, matrixCoordX, matrixCoordY), 
                   start = list(parameters = parametersInt), trace = trace), silent = TRUE) 
    fit
}


# from the Gaussian2D fit, get whether it's a good fit
# i.e. whether if has good Pr(>|t|) for every parameters
# Use legnth(fit) == 1 to represent fit as NA. Otherwise, length(fit) == 8
isLowPValueGaussianFit = function(fit, pValueLimit = 0.1) {    
    isLowPValue <- TRUE
    if (length(fit) == 1) {
        isLowPValue <- FALSE
    } else {
        parameters <- summary(fit)$parameters
        pValues <- parameters[,4]
        isLowPValue <- sum(pValues > pValueLimit) == 0
    }
    isLowPValue
} 


# is not null, and has low pValue, and is sigma is not at lowerBound or upperBound
isGoodGaussianFit = function(fit, pValueLimit = 0.1, lowerSigmaBound = 0.25, upperSigmaBound = 0.5) {
    sigmaR <- getSigmaRFromGaussianFit(fit)
    isGood <- isLowPValueGaussianFit(fit, pValueLimit) & 
        (sigmaR > lowerSigmaBound) # & (sigmaR < upperSigmaBound)
    isGood
}


# from the Gaussin2D fit, get amp
getAmpFromGaussianFit = function(fit) {
    amp <- 0
    if (length(fit) != 1) {
        parameters <- coef(fit)
        amp <- parameters["amp"]
    }
    amp
}


# from the Gaussin2D fit, get x0
getX0FromGaussianFit = function(fit) {
    x0 <- 0
    if (length(fit) != 1) {
        parameters <- coef(fit)
        x0 <- parameters["x0"]
    }
    x0
}


# from the Gaussin2D fit, get y0
getY0FromGaussianFit = function(fit) {
    y0 <- 0
    if (length(fit) != 1) {
        parameters <- coef(fit)
        y0 <- parameters["y0"]  
    }
    y0
}


# from the Gaussin2D fit, get sigma_r
getSigmaRFromGaussianFit = function(fit) {
    sigma_r <- 0
    if (length(fit) != 1) {
        parameters <- coef(fit)
        sigma_r <- parameters["sigma_r"] 
    }
    sigma_r
}


# from the Gaussin2D fit, get bg
getBackgroundFromGaussianFit = function(fit) {
    bg <- 0
    if (length(fit) != 1) {
        parameters <- coef(fit)
        bg <- parameters["bg"]
    }
    bg
}


# from the Gaussin2DFor2 fit, get amp1
getAmp1FromGaussianFor2Fit = function(fit) {
    amp <- 0
    if (length(fit) != 1) {
        parameters <- summary(fit)$parameters
        amp <- parameters[1,1]
    }
    amp
}


# from the Gaussin2DFor2 fit, get amp2
getAmp2FromGaussianFor2Fit = function(fit) {
    amp <- 0
    if (length(fit) != 1) {
        parameters <- summary(fit)$parameters
        amp <- parameters[5,1]
    }
    amp
}


# from the Gaussin2DFor2 fit, get x1 (x position of the 1st Gaussian)
getX1FromGaussianFor2Fit = function(fit) {
    x <- 0
    if (length(fit) != 1) {
        parameters <- summary(fit)$parameters
        x <- parameters[2,1]
    }
    x
}


# from the Gaussin2DFor2 fit, get x2 (x position of the 2nd Gaussian)
getX2FromGaussianFor2Fit = function(fit) {
    x <- 0
    if (length(fit) != 1) {
        parameters <- summary(fit)$parameters
        x <- parameters[6,1]
    }
    x
}


# from the Gaussin2DFor2 fit, get y1 (y position of the 1st Gaussian)
getY1FromGaussianFor2Fit = function(fit) {
    y <- 0
    if (length(fit) != 1) {
        parameters <- summary(fit)$parameters
        y <- parameters[3,1]
    }
    y
}


# from the Gaussin2DFor2 fit, get y2 (y position of the 2nd Gaussian)
getY2FromGaussianFor2Fit = function(fit) {
    y <- 0
    if (length(fit) != 1) {
        parameters <- summary(fit)$parameters
        y <- parameters[7,1]
    }
    y
}


# from the Gaussin2DFor2 fit, get sigma1
getSigma1FromGaussianFor2Fit = function(fit) {
    sigma <- 0
    if (length(fit) != 1) {
        parameters <- summary(fit)$parameters
        sigma <- parameters[4,1]
    }
    sigma
}


# from the Gaussin2DFor2 fit, get sigma2
getSigma2FromGaussianFor2Fit = function(fit) {
    sigma <- 0
    if (length(fit) != 1) {
        parameters <- summary(fit)$parameters
        sigma <- parameters[8,1]
    }
    sigma
}


# choose the best colors from a Gaussian normalized amp data.frame
# Input: data.frame 1st col as green norm amp, 2nd col as yellow, 3rd as orange, 4th as red, 5th a ref
# Output: a character vector indication colors, i.e. 3 as orange, 31 and orange and then green.
chooseBestColorFromAGaussianAmp = function(GaussianAmpNormal, firstColorAmpLimit = 0.3, 
                                        secondColorRatioLimit = 0.3, 
                                        thirdColorRatioLimit = 0.3,
                                        fourthColorRatioLimit = 0.3) {
    ampInOrder <- t(apply(GaussianAmpNormal[,1:4], 1, order, decreasing = TRUE))
    ampHighToLow <- t(apply(GaussianAmpNormal[,1:4], 1, sort, decreasing = TRUE))
    ampRatio <- cbind(ampHighToLow[,2] / ampHighToLow[,1], 
                        ampHighToLow[,3] / ampHighToLow[,2], 
                        ampHighToLow[,4] / ampHighToLow[,3])
    isFirtColorGood <- ampHighToLow[,1] > firstColorAmpLimit
    isSecondColorGood <- isFirtColorGood & ampRatio[,1] > secondColorRatioLimit
    isThirdColorGood <- isSecondColorGood & ampRatio[,2] > thirdColorRatioLimit
    isFourthColorGood <- isThirdColorGood & ampRatio[,3] > fourthColorRatioLimit
    
    colors <- ifelse(isFirtColorGood, as.character(ampInOrder[,1]), "")
    colors <- paste(colors, ifelse(isSecondColorGood, ampInOrder[,2], ""), sep = "")
    colors <- paste(colors, ifelse(isThirdColorGood, ampInOrder[,3], ""), sep = "")
    colors <- paste(colors, ifelse(isFourthColorGood, ampInOrder[,4], ""), sep = "")
    colors
}


# for spots1, find for each spot1 the closest spot2 to it, and if the distance between them is
# less than coSpot_dist, add the color of spot2
correspondColorOf2On1 = function(spots1, spots2, coSpot_dist = 1) {
    spots2ClosestSpotTo1 <- closestSpot(spots1, spots2)
    closestDistances <- distance(spots1, spots2ClosestSpotTo1)
    color <- spots1$Color
    iDistanceSmall <- closestDistances < coSpot_dist
    color[iDistanceSmall] <- as.character(spots2ClosestSpotTo1$Color[iDistanceSmall])
    color
}


# Get the index of the spots that's closest to the spot. If the spot is a vector, the index
# is a vector too.
closestSpotIndex = function(spot, spots, countOnlyUnusedSpots = FALSE) {
    index <- numeric()
    for (r in 1 : nrow(spot)) {
        distances = spot_distance(spot[r,], spots)
        if (countOnlyUnusedSpots) {
            distances[spots$IsUsed] = .Machine$integer.max
        }
        index <- c(index, which.min(distances));
    }
    index
}

# Assume every spot can have only 1 refCoSpot but each refSpot could have multiple coSpots
# because ref plot has only 1 color and the FQ cannot distinguish two close points well.
# Problem with the function, might need to discard.
# because there could be two very close-by refSpots sharing the same coSpot, but only 1 cospot
# is 4gned to each refSpot
indiceOfCoSpots = function(refSpots, spots, coSpot_dist = 1) {
    refCoSpot <- closestSpot(spots, refSpots)
    closestDistances <- distance(spots, refCoSpot)
    refCoSpotIndex <- closestSpotIndex(spots, refSpots)
    iDistanceSmall <- closestDistances < coSpot_dist
    coSpotsIndice <- list()
    for(i in 1 : nrow(refSpots)) {
        coSpotsIndice[[i]] <- which(((refCoSpotIndex == i) & iDistanceSmall) == TRUE)
        #numCoSpotsForTheRefSpot <- sum((refCoSpotIndex == i) & iDistanceSmall)
        #numCoSpots <- c(numCoSpots, numCoSpotsForTheRefSpot)
    }
    coSpotsIndice
}


# Problem with the function, might need to discard.
# because there could be two very close-by refSpots sharing the same coSpot, 
# but only 1 cospot is 4gned to each refSpot
numCoSpots = function(refSpots, spots, coSpot_dist = 0.3) {
    refCoSpot <- closestSpot(spots, refSpots)   # dim of spots
    closestDistances <- distance(spots, refCoSpot)  # dim of spots
    refCoSpotIndex <- closestSpotIndex(spots, refSpots)   # dim of spots
    iDistanceSmall <- closestDistances < coSpot_dist
    numCoSpots <- numeric()
    for(i in 1 : nrow(refSpots)) {
        numCoSpotsForTheRefSpot <- sum((refCoSpotIndex == i) & iDistanceSmall)
        numCoSpots <- c(numCoSpots, numCoSpotsForTheRefSpot)
    }
    numCoSpots
}


# If a refSpot has only 1 slice that has 2 cospots and the rest 1 coSpot, make the refSpot 
# to be 2 refSpots map the 2 coSpots back to 2 refSpots and remove the original refSpot.
splitADoubletCospots = function(refSpots, f1Spots, f2Spots, f3Spots, f4Spots) {
    splitRefSpots <- spotDataFrame()
    for (i in 1 : nrow(refSpots)) {
        refSpot <- refSpots[i,]
        if (refSpot$f1NumCoSpots == 2 & refSpot$f2NumCoSpots == 1 & refSpot$f3NumCoSpots == 1 & 
                refSpot$f4NumCoSpots == 1) {
            f1CoSpots <- f1Spots[refSpot$f1CoSpotsIndice[[1]],]
            for (j in 1 : nrow(f1CoSpots)) {
                newRefSpot <- refSpot
                newRefSpot$Pos_X <- f1CoSpots[j, ]$Pos_X
                newRefSpot$Pos_Y <- f1CoSpots[j, ]$Pos_Y
                splitRefSpots <- rbind(splitRefSpots, newRefSpot)
            }
        } else if (refSpot$f1NumCoSpots == 1 & refSpot$f2NumCoSpots == 2 & refSpot$f3NumCoSpots == 1 & 
                       refSpot$f4NumCoSpots == 1) {
            f2CoSpots <- f2Spots[refSpot$f2CoSpotsIndice[[1]],]
            for (j in 1 : nrow(f2CoSpots)) {
                newRefSpot <- refSpot
                newRefSpot$Pos_X <- f2CoSpots[j, ]$Pos_X
                newRefSpot$Pos_Y <- f2CoSpots[j, ]$Pos_Y
                splitRefSpots <- rbind(splitRefSpots, newRefSpot)
            }
        } else if (refSpot$f1NumCoSpots == 1 & refSpot$f2NumCoSpots == 1 & refSpot$f3NumCoSpots == 2 & 
                       refSpot$f4NumCoSpots == 1) {
            f3CoSpots <- f3Spots[refSpot$f3CoSpotsIndice[[1]],]
            for (j in 1 : nrow(f3CoSpots)) {
                newRefSpot <- refSpot
                newRefSpot$Pos_X <- f3CoSpots[j, ]$Pos_X
                newRefSpot$Pos_Y <- f3CoSpots[j, ]$Pos_Y
                splitRefSpots <- rbind(splitRefSpots, newRefSpot)
            }
        } else if (refSpot$f1NumCoSpots == 1 & refSpot$f2NumCoSpots == 1 & refSpot$f3NumCoSpots == 1 & 
                       refSpot$f4NumCoSpots == 2) {
            f4CoSpots <- f4Spots[refSpot$f4CoSpotsIndice[[1]],]
            for (j in 1 : nrow(f4CoSpots)) {
                newRefSpot <- refSpot
                newRefSpot$Pos_X <- f4CoSpots[j, ]$Pos_X
                newRefSpot$Pos_Y <- f4CoSpots[j, ]$Pos_Y
                splitRefSpots <- rbind(splitRefSpots, newRefSpot)
            }
        } else {
            splitRefSpots <- rbind(splitRefSpots, refSpot)
        }    
    }
    splitRefSpots
}



# This might no be useful because it doesn't fit FSH
# minY could be a vector, so could cbind(y's, x's), x and y are both 0 if 
# the minY, maxY, minX, maxX are outside the matrix index region
maxTiffAmpYXIndiceWithinPixel = function(minY, maxY, minX, maxX, normTiffMatrix) {
    resY <- dim(normTiffMatrix)[1]
    resX <- dim(normTiffMatrix)[2]
    #minNormTiffAMP <- min(min(normTiffMatrix))
    
    #     isOutsideBoundary <- maxX <= 1 | minX >= resX | maxY <= 1 | minY >= resY
    #     trueMinX <- ifelse(minX < 1, 1, minX)
    #     trueMaxX <- ifelse(maxX > resX, resX, maxX)
    #     trueMinY <- ifelse(minY < 1, 1, minY)
    #     trueMaxY <- ifelse(maxY > resY, resY, maxY) 
    isInsideBoundary <- maxX < resX & minX > 1 & maxY < resY & minY > 1
    
    size <- length(minY)
    Xs <- rep(0, size)
    Ys <- rep(0, size)
    maxIndexInMatrix_X <- rep(0, size)
    maxIndexInMatrix_Y <- rep(0, size)
    for (i in 1 : size) {
        if (isInsideBoundary[i]) {
            subMatrix <- normTiffMatrix[minY[i] : maxY[i] , minX[i] : maxX[i]] 
            maxIndexInSubMatrix <- which(subMatrix == max(subMatrix), arr.ind = TRUE)
            Xs[i] <- maxIndexInSubMatrix[1,][2] # there could be multiple pixel that have max value
            Ys[i] <- maxIndexInSubMatrix[1,][1]
            maxIndexInMatrix_Y[i] <- minY[i] + maxIndexInSubMatrix[1] - 1
            maxIndexInMatrix_X[i] <- minX[i] + maxIndexInSubMatrix[2] - 1
        }
    }
    cbind(maxIndexInMatrix_Y, maxIndexInMatrix_X)
}

# This might no be useful because it doesn't fit FSH
ampOverSurroundingPixels = function(pixel, matrix) {
    
}

# This might no be useful because it doesn't fit FSH
# given a set of indice, get only the part that's withint cbind (1:resY, 1:resX)
getMaxAMPIndice = function(minY, maxY, minX, maxX, normTiffMatrix) {
    resY <- dim(normTiffMatrix)[1]
    resX <- dim(normTiffMatrix)[2]
    #minNormTiffAMP <- min(min(normTiffMatrix))
    
    isOutsideBoundary <- maxX <= 1 | minX >= resX | maxY <= 1 | minY >= resY
    trueMinX <- ifelse(minX < 1, 1, minX)
    trueMaxX <- ifelse(maxX > resX, resX, maxX)
    trueMinY <- ifelse(minY < 1, 1, minY)
    trueMaxY <- ifelse(maxY > resY, resY, maxY) 
    maxAMP <- rep(0, length(minY))
     
    #maxToSurroundAMP <- rep(0, length(minY))
    for (i in 1 : length(minY)) {
        if (isOutsideBoundary[i]) {
            maxAMP[i] <- -Inf
        } else {
            subMatrix <- normTiffMatrix[trueMinY[i] : trueMaxY[i] , trueMinX[i] : trueMaxX[i]] 
            maxIndexInSubMatrix <- which(subMatrix == max(subMatrix), arr.ind = TRUE)
            maxIndexXInSubMatrix <- maxIndexInSubMatrix[1,][2] # there could be multiple pixel that have max value
            maxIndexYInSubMatrix <- maxIndexInSubMatrix[1,][1]
            maxAMP[i] <- subMatrix[maxIndexYInSubMatrix,maxIndexXInSubMatrix]  
            
            maxIndexInMatrix_Y <- trueMinY[i] + maxIndexInSubMatrix[1] - 1
            maxIndexInMatrix_X <- trueMinX[i] + maxIndexInSubMatrix[2] - 1
            surroundMaxMatrix <- normTiffMatrix[(maxIndexInMatrix_Y - 1) : (maxIndexInMatrix_Y + 1),
                                                (maxIndexInMatrix_X - 1) : (maxIndexInMatrix_X + 1)]
            
            averageAMPOfSurroundedPixel <- (((sum(sum(surroundMaxMatrix))) - maxAMP[i]) / 8)
            maxToSurroundAMP[i] <- maxAMP[i] / averageAMPOfSurroundedPixel
        }    
    }
    cbind(maxAMP, maxToSurroundAMP)
}


# Get the spots within neighbor_dist (um) of the spot (not exlcuding the spot itself)
neighbors_circle = function (spot, spots, neighbor_dist = 1, countOnlyUnusedSpots = FALSE) {
    neighbors <- spots %>% filter(((Pos_X - spot$Pos_X)^2 + (Pos_Y - spot$Pos_Y)^2) < 
                                      neighbor_dist^2)
    if (countOnlyUnusedSpots) {
        neighbors <- neighbors[neighbors$IsUsed == FALSE,]
    }
    neighbors    
}


# Get the spots within neighbor_dist (um) of the spot (exlcuding the spot itself)
neighbors_exclude_spot_circle = function (spot, spots, neighbor_dist = 1) {
    neighbors <- neighbors_circle(spot, spots, neighbor_dist)
    neighbors <- neighbors %>% filter(!(Pos_X == spot$Pos_X & Pos_Y == spot$Pos_Y))
    neighbors
}


# Get the number of closeby neighbors within neighbor_dist (um) 
# (not excluding the spots within neighborSpots)
num_neighbors_circle = function(spots, neighborSpots, neighbor_dist = 1, countOnlyUnusedSpots = FALSE) {
    count <- integer(nrow(spots))
    for (r in seq_along(spots[,1])) {
        spot <- spots[r,]
        neighborsCircle <- neighbors_circle(spot, neighborSpots, neighbor_dist, countOnlyUnusedSpots)
        count[r] <- nrow(neighborsCircle)
    }
    count
} 


# Get the number of closeby neighbors within neighbor_dist (um) 
# (excluding the spots within neighborSpots)
num_neighbors_exclude_spot_circle = function(spots, neighborSpots, neighbor_dist = 1) {
    count <- integer(nrow(spots))
    for (r in seq_along(spots[,1])) {
        spot <- spots[r,]
        neighborsExcludingSelfCircle <- neighbors_exclude_spot_circle(spot, neighborSpots, neighbor_dist)
        count[r] <- nrow(neighborsExcludingSelfCircle)
    }
    count
} 


# Get the indice (respoect to the neighbor spots dataframe) of closeby neighbors 
# within neighbor_dist (um) as a list of vectors (not excluding the spots within neighborSpots) 
index_neighbors_circle = function(spots, neighborSpots, neighbor_dist = 1, countOnlyUnusedSpots = FALSE) {
    coSpotsIndice <- list()
    for (r in seq_along(spots[,1])) {
        spot <- spots[r,]
        neighborsCircle <- neighbors_circle(spot, neighborSpots, neighbor_dist, countOnlyUnusedSpots)
        neighborsIndice <- integer()
        for (rCircle in seq_along(neighborsCircle[,1])) {
            neighborCircle <- neighborsCircle[rCircle,]
            neighborIndex <- which(neighborSpots$Pos_X == neighborCircle$Pos_X &
                                         neighborSpots$Pos_Y == neighborCircle$Pos_Y &
                                         neighborSpots$AMP == neighborCircle$AMP)
            neighborsIndice <- c(neighborsIndice, neighborIndex)
        }
        coSpotsIndice[[r]] <- neighborsIndice
    }
    coSpotsIndice
}

# transColorXYPerPolygonList has registration parameters for the four colors at a time stamp
# colorTiffMatrixList has the colorTiffMatrix for the four colors at a time stamp
# colorSpotsList has the colorSpots  for the four colors at a time stamp
# return a color out of the four colors
# Fit a Gaussian first. 
# If a Gaussian fit return NA, substract nearby spots within 2um radius as Gaussians and refit.  
# If the Gaussian fit returns NA again, subtract nearby spots within 3um radius as Gaussians and refit.
# Note coSpotDistColor <- sqrt((loneColorGaussianX0 - refSpot$Pos_X)^2 + (loneColorGaussianY0 - refSpot$Pos_Y)^2)
# or coSpotDistColor <- sqrt((loneColorGaussianX0 - transRefSpotOnColor$Pos_X)^2 + (loneColorGaussianY0 - transRefSpotOnColor$Pos_Y)^2)
chooseAColorFromGaussianFitFromList = function (transColorXYPerPolygonList, colorTiffMatrixList, 
                                                colorSpotsList,
                                                maxMinusBgTiffAmpPerColor, colors, refSpot,  
                                                diffractionLimits, coSpot_dist = 1, 
                                                isToPrintFit = FALSE, SNRThreshold = c(0, 0, 0, 0)) {
    loneGaussianAmps <- NULL
    loneGaussianBgs <- NULL
    coSpotDists <- NULL
    for (i in 1: length(transColorXYPerPolygonList)) {
        transRefSpotOnColor <- moveRefSpotsByTransXYPerPolygon(refSpot, 
                                polygonMatrix, transColorXYPerPolygonList[[i]])
        
        # debug to remove or add later
        # get nearby neighbors, if the neighbor is used, fit the neighbor with fix-centered gaussian and substract it
        # neighborSpots <- neighbors_exclude_spot_circle(transRefSpotOnColor, colorSpotsList[[i]], neighbor_dist = 1)
        # usedNeighborSpots <- neighborSpots %>% filter(IsUsed = TRUE)
        
        
        loneColorFit <- fitGaussian2D(transRefSpotOnColor, colorTiffMatrixList[[i]], 
                                      diffractionLimits[i])
        fitDistance <- diffractionLimits[i]
#        loneColorFit <- fitGaussian2DFixCenter(transRefSpotOnColor, 
#                                colorTiffMatrixList[[i]], fitDistance, trace = FALSE)
        
        # If the fit is not good, subtract nearby strong spots from the tiffMatrix and  refit
        # However, this doesn't seem to help at all. loneColorFit is still NA
        # However, subtraction can 
        subtractDistance <- diffractionLimits[i]

        if (is.na(loneColorFit)) { #| getSigmaRFromGaussianFit(loneColorFit) == 0.5) {    
            neighborDistance <- 2
            neighborSpots <- neighbors_exclude_spot_circle(transRefSpotOnColor, 
                                                colorSpotsList[[i]], neighborDistance)
            subtractedTiffMatrix <- subtractGaussianFitInOrder(colorTiffMatrixList[[i]], 
                                            neighborSpots, subtractDistance)
            localMaximumIndice2D <- getLocalMaximumIndice2D(subtractedTiffMatrix, refSpot, neighborDistance)
            extraSpots <- createSpotsFrom2DIndice(subtractedTiffMatrix, localMaximumIndice2D)
            extraBrightSpots <- extraSpots[extraSpots$tiffAmp > subtractedTiffMatrix[spotsToIndice2D(refSpot, pix_XY)],]
            subtractedTiffMatrix2 <- subtractGaussianFitInOrder(subtractedTiffMatrix,
                                            extraBrightSpots, subtractDistance)

            loneColorFit <- fitGaussian2DFixCenter(transRefSpotOnColor, 
                    subtractedTiffMatrix2, fitDistance, trace = FALSE)
            
            #plotNearbyTiffRaster(refSpot, extraBrightSpots, pix_XY, subtractedTiffMatrix, isToTranslateRefSpots = FALSE, plot_area_dist = 5)
            #plotNearbyTiffRaster(refSpot, extraBrightSpots, pix_XY, subtractedTiffMatrix2, isToTranslateRefSpots = FALSE, plot_area_dist = 5)
        } 
        if (is.na(loneColorFit)) {
            neighborDistance <- 3
            neighborSpots <- neighbors_exclude_spot_circle(transRefSpotOnColor, 
                                                           colorSpotsList[[i]], neighborDistance)
            subtractedTiffMatrix <- subtractGaussianFitInOrder(colorTiffMatrixList[[i]], 
                                                               neighborSpots, subtractDistance)
            localMaximumIndice2D <- getLocalMaximumIndice2D(subtractedTiffMatrix, refSpot, neighborDistance)
            extraSpots <- createSpotsFrom2DIndice(subtractedTiffMatrix, localMaximumIndice2D)
            extraBrightSpots <- extraSpots[extraSpots$tiffAmp > subtractedTiffMatrix[spotsToIndice2D(refSpot, pix_XY)],]
            subtractedTiffMatrix2 <- subtractGaussianFitInOrder(subtractedTiffMatrix,
                                                                extraBrightSpots, subtractDistance)
            
            loneColorFit <- fitGaussian2DFixCenter(transRefSpotOnColor, subtractedTiffMatrix2, 
                                                   fitDistance, trace = FALSE)
            
        }

#        loneColorGaussianX0 <- getX0FromGaussianFit(loneColorFit)
#        loneColorGaussianY0 <- getY0FromGaussianFit(loneColorFit)
        coSpotDistColor <- 0#sqrt((loneColorGaussianX0 - refSpot$Pos_X)^2 + 
                             #   (loneColorGaussianY0 - refSpot$Pos_Y)^2)
        loneGaussianAmps <- c(loneGaussianAmps, getAmpFromGaussianFit(loneColorFit))
        loneGaussianBgs <- c(loneGaussianBgs, getBackgroundFromGaussianFit(loneColorFit))
        coSpotDists <- c(coSpotDists, coSpotDistColor)
        if (isToPrintFit) {
            print(paste(colors[i]))
            print(loneColorFit)
        }
    }
    normalizedLoneGaussianAmps <- loneGaussianAmps / maxMinusBgTiffAmpPerColor
    normalizedLoneGaussianAmps[coSpotDists > coSpot_dist] = 0
    snr <- loneGaussianAmps / loneGaussianBgs
    normalizedLoneGaussianAmps[snr < SNRThreshold] = 0
    colorIndex <- which.max(normalizedLoneGaussianAmps)
    color <- colors[colorIndex]
    color
}

# Given list of spots and a tiff matrix, subtract the Gaussian fit of each spot
# from strong intensity to weak intensity from the tiff matrix
# Return a subtracted tiff matrix
subtractGaussianFitInOrder = function(tiffMatrix, spots, fit_distance_um) {
    currentTiffMatrix <- tiffMatrix
    if (nrow(spots) != 0) {
        orderedSpots <- spots[order(spots$tiffAmp, decreasing = TRUE),] 
        numSpots <- dim(orderedSpots)[1]
        amps <- vector(length = numSpots)
        sigma_rs <- vector(length = numSpots)
        bgs <- vector(length = numSpots)
        
        for (i in seq_along(orderedSpots[,1])) {
            spot <- orderedSpots[i,]
            
            # Fit the Gaussian2DFixedCenter
            curerntGaussianFit <- fitGaussian2DFixCenter(spot, currentTiffMatrix, 
                                                         fit_distance_um, 0.25, 0.5)
            amps[i] <- getAmpFromGaussianFit(curerntGaussianFit)
            sigma_rs[i] <- getSigmaRFromGaussianFit(curerntGaussianFit)
            bgs[i] <- getBackgroundFromGaussianFit(curerntGaussianFit)
            
            # Subtract the Gaussian from the tiffplot
            centerIndice <- spotsToIndice2D(spot, pix_XY)
            nearbyIndice2D <- getNearbyIndice2D(manySpotsTiffMatrix, centerIndice, 
                                                fit_distance_um / pix_XY)
            matrixCoordX <- getXCoordsFromIndice2D(nearbyIndice2D)
            matrixCoordY <- getYCoordsFromIndice2D(nearbyIndice2D)
            Gaussian <- Gaussian2D(amps[i], spot$Pos_X, spot$Pos_Y, sigma_rs[i], 
                                   0, matrixCoordX, matrixCoordY)
            currentTiffMatrix[nearbyIndice2D] <- currentTiffMatrix[nearbyIndice2D] -
                Gaussian
        }
    }
    currentTiffMatrix
}

# Note coSpotDistColor <- sqrt((loneColorGaussianX0 - refSpot$Pos_X)^2 + (loneColorGaussianY0 - refSpot$Pos_Y)^2)
# or coSpotDistColor <- sqrt((loneColorGaussianX0 - transRefSpotOnColor$Pos_X)^2 + (loneColorGaussianY0 - transRefSpotOnColor$Pos_Y)^2)
# Might need to implement diffractionLimits
chooseAColorFromTiffAmpFromList = function (transColorXYPerPolygonList, colorTiffMatrixList, 
                                                maxMinusBgTiffAmpPerColor, colors, refSpot, 
                                                diffractionLimits, coSpot_dist) {
    amps <- NULL
    coSpotDists <- NULL
    for (i in 1: length(transColorXYPerPolygonList)) {
        transRefSpotOnColor <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transColorXYPerPolygonList[[i]])
        centerIndice <- spotsToIndice2D(transRefSpotOnColor, pix_XY)
        tiffMatrix <- colorTiffMatrixList[[i]]
        centerAmp <- tiffMatrix[centerIndice]
        bg <- getMinNearbyAmp(tiffMatrix, centerIndice, 5)
        amp <- centerAmp - bg
        coSpotDistColor <- 0#sqrt((loneColorGaussianX0 - refSpot$Pos_X)^2 + 
        #        (loneColorGaussianY0 - refSpot$Pos_Y)^2)
        amps <- c(amps, amp)
        coSpotDists <- c(coSpotDists, coSpotDistColor)
    }
    normalizedLoneGaussianAmps <- amps / maxMinusBgTiffAmpPerColor
    normalizedLoneGaussianAmps[coSpotDists > coSpot_dist] = 0 
    colorIndex <- which.max(normalizedLoneGaussianAmps)
    color <- colors[colorIndex]
    color
}


# Givine a refSpot, and a time stamp, fit Gaussian at all 4 images of that time stamp
# and choose the color with highest Gaussian amp that is within coSpot_dist (1um) from 
# the refSpot
chooseAColorFromGaussianFit = function(refSpot, time, coSpot_dist = 1, 
                                       isToPrintFit = FALSE) {
    if (time == 1) {
        transColorXYPerPolygonList <- list(transG1XYPerPolygon, transY1XYPerPolygon, transO1XYPerPolygon, transR1XYPerPolygon)
        colorTiffMatrixList <- list(g1TiffMatrix, y1TiffMatrix, o1TiffMatrix, r1TiffMatrix)
        colorSpotsList <- list(g1Spots, y1Spots, o1Spots, r1Spots)
        maxMinusBgTiffAmpPerColor <- f1MaxMinusBgTiffAmpPerColor
    } else if (time == 2) {
        transColorXYPerPolygonList <- list(transG2XYPerPolygon, transY2XYPerPolygon, transO2XYPerPolygon, transR2XYPerPolygon)
        colorTiffMatrixList <- list(g2TiffMatrix, y2TiffMatrix, o2TiffMatrix, r2TiffMatrix)
        colorSpotsList <- list(g2Spots, y2Spots, o2Spots, r2Spots)
        maxMinusBgTiffAmpPerColor <- f2MaxMinusBgTiffAmpPerColor
    } else if (time == 3) {
        transColorXYPerPolygonList <- list(transG3XYPerPolygon, transY3XYPerPolygon, transO3XYPerPolygon, transR3XYPerPolygon)
        colorTiffMatrixList <- list(g3TiffMatrix, y3TiffMatrix, o3TiffMatrix, r3TiffMatrix)
        colorSpotsList <- list(g3Spots, y3Spots, o3Spots, r3Spots)
        maxMinusBgTiffAmpPerColor <- f3MaxMinusBgTiffAmpPerColor
    } else if (time == 4) {
        transColorXYPerPolygonList <- list(transG4XYPerPolygon, transY4XYPerPolygon, transO4XYPerPolygon, transR4XYPerPolygon)
        colorTiffMatrixList <- list(g4TiffMatrix, y4TiffMatrix, o4TiffMatrix, r4TiffMatrix)
        colorSpotsList <- list(g4Spots, y4Spots, o4Spots, r4Spots)
        maxMinusBgTiffAmpPerColor <- f4MaxMinusBgTiffAmpPerColor
    } 
    if (time != 1 && time != 2 && time != 3 && time != 4) {
        color <- "grey"
    } else {
        colors <- c("green", "yellow", "orange", "red")
        diffractionLimits <- c(diffractionLimitG, diffractionLimitY, 
                               diffractionLimitO, diffractionLimitR)
        color <- chooseAColorFromGaussianFitFromList(transColorXYPerPolygonList, colorTiffMatrixList,
                                                     colorSpotsList, 
                                                     maxMinusBgTiffAmpPerColor, colors, refSpot, 
                                                     diffractionLimits, coSpot_dist,
                                                     isToPrintFit)
    }
    color
}



# Givine a refSpot, and a time stamp, fit Gaussian at all 4 images of that time stamp
# and choose the color with highest Gu4na amp that is within coSpot_dist (1um) from 
# the refSpot
# Note: this is for 6mRNAs on 3/3/16 only, not a general function
chooseAColorFromGaussianFit6mRNAs = function(refSpot, time, coSpot_dist = 1, 
                                             isToPrintFit = FALSE) {
    if (time == 1) {
        transColorXYPerPolygonList <- list(transG1XYPerPolygon, transO1XYPerPolygon)
        colorTiffMatrixList <- list(g1TiffMatrix, o1TiffMatrix)
        colorSpotsList <- list(g1Spots, o1Spots)    
        maxMinusBgTiffAmpPerColor <- f1MaxMinusBgTiffAmpPerColor
        colors <- c("green", "orange")
        diffractionLimits <- c(diffractionLimitG, diffractionLimitO)
        #sigmaRs <- c(sigmaR_G, sigmaR_O)
    } else if (time == 2) {
        transColorXYPerPolygonList <- list(transG2XYPerPolygon, transY2XYPerPolygon)
        colorTiffMatrixList <- list(g2TiffMatrix, y2TiffMatrix)
        colorSpotsList <- list(g2Spots, y2Spots)    
        maxMinusBgTiffAmpPerColor <- f2MaxMinusBgTiffAmpPerColor
        colors <- c("green", "yellow")
        diffractionLimits <- c(diffractionLimitG, diffractionLimitY)
        #sigmaRs <- c(sigmaR_G, sigmaR_Y)
    } else if (time == 3) {
        transColorXYPerPolygonList <- list(transG3XYPerPolygon, transR3XYPerPolygon)
        colorTiffMatrixList <- list(g3TiffMatrix, r3TiffMatrix)
        colorSpotsList <- list(g3Spots, r3Spots)    
        maxMinusBgTiffAmpPerColor <- f3MaxMinusBgTiffAmpPerColor
        colors <- c("green", "red")
        diffractionLimits <- c(diffractionLimitG, diffractionLimitR)
        #sigmaRs <- c(sigmaR_G, sigmaR_R)
    } else if (time == 4) {
        transColorXYPerPolygonList <- list(transY4XYPerPolygon, transR4XYPerPolygon)
        colorTiffMatrixList <- list(y4TiffMatrix, r4TiffMatrix)
        colorSpotsList <- list(y4Spots, r4Spots)    
        maxMinusBgTiffAmpPerColor <- f4MaxMinusBgTiffAmpPerColor
        colors <- c("yellow", "red")
        diffractionLimits <- c(diffractionLimitY, diffractionLimitR)
        #sigmaRs <- c(sigmaR_Y, sigmaR_R)
    } 
    if (time != 1 && time != 2 && time != 3 && time != 4) {
        color <- "grey"
    } else {
        color <- chooseAColorFromGaussianFitFromList(transColorXYPerPolygonList, colorTiffMatrixList,
                                                     colorSpotsList,    
                                                     maxMinusBgTiffAmpPerColor, colors, refSpot, 
                                                     diffractionLimits, coSpot_dist, 
                                                     isToPrintFit)
    }
    color
}


# Note: might be able to get rid of this function.
# Givine a refSpot, and a time stamp, fit Gaussian at all 4 images of that time stamp
# and choose the color with highest Gu4na amp that is within coSpot_dist (1um) from 
# the refSpot
# chooseAColorFromGaussianFit = function(refSpot, time, diffractionLimit = 1.5 , coSpot_dist = 1) {
#     colors <- c("green", "yellow", "orange", "red")
#     if (time == 1) {
#         transRefSpotOnG1 <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transG1XYPerPolygon)
#         transRefSpotOnY1 <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transY1XYPerPolygon)
#         transRefSpotOnO1 <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transO1XYPerPolygon)
#         transRefSpotOnR1 <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transR1XYPerPolygon)    
#         
#         loneG1Fit <- fitGaussian2D(transRefSpotOnG1, g1TiffMatrix, diffractionLimit)
#         loneY1Fit <- fitGaussian2D(transRefSpotOnY1, y1TiffMatrix, diffractionLimit)
#         loneO1Fit <- fitGaussian2D(transRefSpotOnO1, o1TiffMatrix, diffractionLimit)
#         loneR1Fit <- fitGaussian2D(transRefSpotOnR1, r1TiffMatrix, diffractionLimit)
#         
#         loneG1GaussianAmp <- getAmpFromGaussianFit(loneG1Fit)
#         loneY1GaussianAmp <- getAmpFromGaussianFit(loneY1Fit)
#         loneO1GaussianAmp <- getAmpFromGaussianFit(loneO1Fit)
#         loneR1GaussianAmp <- getAmpFromGaussianFit(loneR1Fit)
#         
#         loneG1GaussianX0 <- getX0FromGaussianFit(loneG1Fit)
#         loneY1GaussianX0 <- getX0FromGaussianFit(loneY1Fit)
#         loneO1GaussianX0 <- getX0FromGaussianFit(loneO1Fit)
#         loneR1GaussianX0 <- getX0FromGaussianFit(loneR1Fit)
#         
#         loneG1GaussianY0 <- getY0FromGaussianFit(loneG1Fit)
#         loneY1GaussianY0 <- getY0FromGaussianFit(loneY1Fit)
#         loneO1GaussianY0 <- getY0FromGaussianFit(loneO1Fit)
#         loneR1GaussianY0 <- getY0FromGaussianFit(loneR1Fit)
#         
#         coSpotDistG1 <- sqrt((loneG1GaussianX0 - refSpot$Pos_X)^2 + (loneG1GaussianY0 - refSpot$Pos_Y)^2)
#         coSpotDistY1 <- sqrt((loneY1GaussianX0 - refSpot$Pos_X)^2 + (loneY1GaussianY0 - refSpot$Pos_Y)^2)
#         coSpotDistO1 <- sqrt((loneO1GaussianX0 - refSpot$Pos_X)^2 + (loneO1GaussianY0 - refSpot$Pos_Y)^2)
#         coSpotDistR1 <- sqrt((loneR1GaussianX0 - refSpot$Pos_X)^2 + (loneR1GaussianY0 - refSpot$Pos_Y)^2)
#         
#         loneGaussianAmps <- c(loneG1GaussianAmp, loneY1GaussianAmp, loneO1GaussianAmp, loneR1GaussianAmp)
#         normalizedLoneGaussianAmps <- loneGaussianAmps / f1MaxMinusBgTiffAmpPerColor
#         coSpotDists <- c(coSpotDistG1, coSpotDistY1, coSpotDistO1, coSpotDistR1)
#         normalizedLoneGaussianAmps[coSpotDists > 1] = 0 
# 
#         colorIndex <- which.max(normalizedLoneGaussianAmps)
#         color <- colors[colorIndex]
#     } else if (time == 2) {
#         transRefSpotOnG2 <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transG2XYPerPolygon)
#         transRefSpotOnY2 <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transY2XYPerPolygon)
#         transRefSpotOnO2 <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transO2XYPerPolygon)
#         transRefSpotOnR2 <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transR2XYPerPolygon)    
#         
#         loneG2Fit <- fitGaussian2D(transRefSpotOnG2, g2TiffMatrix, diffractionLimit)
#         loneY2Fit <- fitGaussian2D(transRefSpotOnY2, y2TiffMatrix, diffractionLimit)
#         loneO2Fit <- fitGaussian2D(transRefSpotOnO2, o2TiffMatrix, diffractionLimit)
#         loneR2Fit <- fitGaussian2D(transRefSpotOnR2, r2TiffMatrix, diffractionLimit)
#         
#         loneG2GaussianAmp <- getAmpFromGaussianFit(loneG2Fit)
#         loneY2GaussianAmp <- getAmpFromGaussianFit(loneY2Fit)
#         loneO2GaussianAmp <- getAmpFromGaussianFit(loneO2Fit)
#         loneR2GaussianAmp <- getAmpFromGaussianFit(loneR2Fit)
#         
#         loneG2GaussianX0 <- getX0FromGaussianFit(loneG2Fit)
#         loneY2GaussianX0 <- getX0FromGaussianFit(loneY2Fit)
#         loneO2GaussianX0 <- getX0FromGaussianFit(loneO2Fit)
#         loneR2GaussianX0 <- getX0FromGaussianFit(loneR2Fit)
#         
#         loneG2GaussianY0 <- getY0FromGaussianFit(loneG2Fit)
#         loneY2GaussianY0 <- getY0FromGaussianFit(loneY2Fit)
#         loneO2GaussianY0 <- getY0FromGaussianFit(loneO2Fit)
#         loneR2GaussianY0 <- getY0FromGaussianFit(loneR2Fit)
#         
#         coSpotDistG2 <- sqrt((loneG2GaussianX0 - refSpot$Pos_X)^2 + (loneG2GaussianY0 - refSpot$Pos_Y)^2)
#         coSpotDistY2 <- sqrt((loneY2GaussianX0 - refSpot$Pos_X)^2 + (loneY2GaussianY0 - refSpot$Pos_Y)^2)
#         coSpotDistO2 <- sqrt((loneO2GaussianX0 - refSpot$Pos_X)^2 + (loneO2GaussianY0 - refSpot$Pos_Y)^2)
#         coSpotDistR2 <- sqrt((loneR2GaussianX0 - refSpot$Pos_X)^2 + (loneR2GaussianY0 - refSpot$Pos_Y)^2)
#         
#         loneGaussianAmps <- c(loneG2GaussianAmp, loneY2GaussianAmp, loneO2GaussianAmp, loneR2GaussianAmp)
#         normalizedLoneGaussianAmps <- loneGaussianAmps / f2MaxMinusBgTiffAmpPerColor
#         coSpotDists <- c(coSpotDistG2, coSpotDistY2, coSpotDistO2, coSpotDistR2)
#         normalizedLoneGaussianAmps[coSpotDists > 1] = 0 
#         
#         colorIndex <- which.max(normalizedLoneGaussianAmps)
#         color <- colors[colorIndex]
#     } else if (time == 3) {
#         transRefSpotOnG3 <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transG3XYPerPolygon)
#         transRefSpotOnY3 <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transY3XYPerPolygon)
#         transRefSpotOnO3 <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transO3XYPerPolygon)
#         transRefSpotOnR3 <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transR3XYPerPolygon)    
#         
#         loneG3Fit <- fitGaussian2D(transRefSpotOnG3, g3TiffMatrix, diffractionLimit)
#         loneY3Fit <- fitGaussian2D(transRefSpotOnY3, y3TiffMatrix, diffractionLimit)
#         loneO3Fit <- fitGaussian2D(transRefSpotOnO3, o3TiffMatrix, diffractionLimit)
#         loneR3Fit <- fitGaussian2D(transRefSpotOnR3, r3TiffMatrix, diffractionLimit)
#         
#         loneG3GaussianAmp <- getAmpFromGaussianFit(loneG3Fit)
#         loneY3GaussianAmp <- getAmpFromGaussianFit(loneY3Fit)
#         loneO3GaussianAmp <- getAmpFromGaussianFit(loneO3Fit)
#         loneR3GaussianAmp <- getAmpFromGaussianFit(loneR3Fit)
#         
#         loneG3GaussianX0 <- getX0FromGaussianFit(loneG3Fit)
#         loneY3GaussianX0 <- getX0FromGaussianFit(loneY3Fit)
#         loneO3GaussianX0 <- getX0FromGaussianFit(loneO3Fit)
#         loneR3GaussianX0 <- getX0FromGaussianFit(loneR3Fit)
#         
#         loneG3GaussianY0 <- getY0FromGaussianFit(loneG3Fit)
#         loneY3GaussianY0 <- getY0FromGaussianFit(loneY3Fit)
#         loneO3GaussianY0 <- getY0FromGaussianFit(loneO3Fit)
#         loneR3GaussianY0 <- getY0FromGaussianFit(loneR3Fit)
#         
#         coSpotDistG3 <- sqrt((loneG3GaussianX0 - refSpot$Pos_X)^2 + (loneG3GaussianY0 - refSpot$Pos_Y)^2)
#         coSpotDistY3 <- sqrt((loneY3GaussianX0 - refSpot$Pos_X)^2 + (loneY3GaussianY0 - refSpot$Pos_Y)^2)
#         coSpotDistO3 <- sqrt((loneO3GaussianX0 - refSpot$Pos_X)^2 + (loneO3GaussianY0 - refSpot$Pos_Y)^2)
#         coSpotDistR3 <- sqrt((loneR3GaussianX0 - refSpot$Pos_X)^2 + (loneR3GaussianY0 - refSpot$Pos_Y)^2)
#         
#         loneGaussianAmps <- c(loneG3GaussianAmp, loneY3GaussianAmp, loneO3GaussianAmp, loneR3GaussianAmp)
#         normalizedLoneGaussianAmps <- loneGaussianAmps / f3MaxMinusBgTiffAmpPerColor
#         coSpotDists <- c(coSpotDistG3, coSpotDistY3, coSpotDistO3, coSpotDistR3)
#         normalizedLoneGaussianAmps[coSpotDists > 1] = 0 
#         
#         colorIndex <- which.max(normalizedLoneGaussianAmps)
#         color <- colors[colorIndex]
#     } else if (time == 4) {
#         transRefSpotOnG4 <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transG4XYPerPolygon)
#         transRefSpotOnY4 <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transY4XYPerPolygon)
#         transRefSpotOnO4 <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transO4XYPerPolygon)
#         transRefSpotOnR4 <- moveRefSpotsByTransXYPerPolygon(refSpot, polygonMatrix, transR4XYPerPolygon)    
#         
#         loneG4Fit <- fitGaussian2D(transRefSpotOnG4, g4TiffMatrix, diffractionLimit)
#         loneY4Fit <- fitGaussian2D(transRefSpotOnY4, y4TiffMatrix, diffractionLimit)
#         loneO4Fit <- fitGaussian2D(transRefSpotOnO4, o4TiffMatrix, diffractionLimit)
#         loneR4Fit <- fitGaussian2D(transRefSpotOnR4, r4TiffMatrix, diffractionLimit)
#         
#         loneG4GaussianAmp <- getAmpFromGaussianFit(loneG4Fit)
#         loneY4GaussianAmp <- getAmpFromGaussianFit(loneY4Fit)
#         loneO4GaussianAmp <- getAmpFromGaussianFit(loneO4Fit)
#         loneR4GaussianAmp <- getAmpFromGaussianFit(loneR4Fit)
#         
#         loneG4GaussianX0 <- getX0FromGaussianFit(loneG4Fit)
#         loneY4GaussianX0 <- getX0FromGaussianFit(loneY4Fit)
#         loneO4GaussianX0 <- getX0FromGaussianFit(loneO4Fit)
#         loneR4GaussianX0 <- getX0FromGaussianFit(loneR4Fit)
#         
#         loneG4GaussianY0 <- getY0FromGaussianFit(loneG4Fit)
#         loneY4GaussianY0 <- getY0FromGaussianFit(loneY4Fit)
#         loneO4GaussianY0 <- getY0FromGaussianFit(loneO4Fit)
#         loneR4GaussianY0 <- getY0FromGaussianFit(loneR4Fit)
#         
#         coSpotDistG4 <- sqrt((loneG4GaussianX0 - refSpot$Pos_X)^2 + (loneG4GaussianY0 - refSpot$Pos_Y)^2)
#         coSpotDistY4 <- sqrt((loneY4GaussianX0 - refSpot$Pos_X)^2 + (loneY4GaussianY0 - refSpot$Pos_Y)^2)
#         coSpotDistO4 <- sqrt((loneO4GaussianX0 - refSpot$Pos_X)^2 + (loneO4GaussianY0 - refSpot$Pos_Y)^2)
#         coSpotDistR4 <- sqrt((loneR4GaussianX0 - refSpot$Pos_X)^2 + (loneR4GaussianY0 - refSpot$Pos_Y)^2)
#         
#         loneGaussianAmps <- c(loneG4GaussianAmp, loneY4GaussianAmp, loneO4GaussianAmp, loneR4GaussianAmp)
#         normalizedLoneGaussianAmps <- loneGaussianAmps / f4MaxMinusBgTiffAmpPerColor
#         coSpotDists <- c(coSpotDistG4, coSpotDistY4, coSpotDistO4, coSpotDistR4)
#         normalizedLoneGaussianAmps[coSpotDists > 1] = 0 
#         
#         colorIndex <- which.max(normalizedLoneGaussianAmps)
#         color <- colors[colorIndex]
#     } else {
#         color <- "grey"
#     }
#     color
# }


# change: a repeated equation, could be deleted
# For singlylLoneRefSpots (lonelyRefSpots that has 0 or 1 coSpot in each time), assign
# 1 color at each time point, either with FQ if numCoSpots = 1 or GaussianFit if numCoSpots = 0 
assign4ColorsFromFQOrGaussianFit = function(singlyLonelyRefSpots, coSpot_dist = 1, 
                                            minNumOfSpotFQFound = 1,
                                            isFor6mRNAs = FALSE, isToPrintFit = FALSE) {
    
    allColors <- NULL
    for (i in seq_along(singlyLonelyRefSpots[,1])) {
        refSpot <- singlyLonelyRefSpots[i,]
        if (refSpot$transF1NumCoSpots >= minNumOfSpotFQFound) {
            color1 <- transF1Spots[refSpot$transF1CoSpotsIndice[[1]][1],]$Color
        } else {
            if (isFor6mRNAs) {
                color1 <- chooseAColorFromGaussianFit6mRNAs(refSpot, time = 1, #diffractionLimit, 
                                                      coSpot_dist, isToPrintFit)
            } else {
                color1 <- chooseAColorFromGaussianFit(refSpot, time = 1, #diffractionLimit, 
                                                      coSpot_dist, isToPrintFit)
            }   
        }

        if (refSpot$transF2NumCoSpots >= minNumOfSpotFQFound) {
            color2 <- transF2Spots[refSpot$transF2CoSpotsIndice[[1]][1],]$Color
        } else {
            if (isFor6mRNAs) {
                color2 <- chooseAColorFromGaussianFit6mRNAs(refSpot, 2, 
                                                      coSpot_dist, isToPrintFit)
            } else {
                color2 <- chooseAColorFromGaussianFit(refSpot, 2, 
                                                      coSpot_dist, isToPrintFit)
            } 
        } 
        
        if (refSpot$transF3NumCoSpots >= minNumOfSpotFQFound) {
            color3 <- transF3Spots[refSpot$transF3CoSpotsIndice[[1]][1],]$Color
        } else {
            if (isFor6mRNAs) {
                color3 <- chooseAColorFromGaussianFit6mRNAs(refSpot, 3, 
                                                      coSpot_dist, isToPrintFit)
            } else {               
                color3 <- chooseAColorFromGaussianFit(refSpot, 3, 
                                                      coSpot_dist, isToPrintFit)
            }
        } 
        
        if (refSpot$transF4NumCoSpots >= minNumOfSpotFQFound) {
            color4 <- transF4Spots[refSpot$transF4CoSpotsIndice[[1]][1],]$Color
        } else {
            if (isFor6mRNAs) {
                color4 <- chooseAColorFromGaussianFit6mRNAs(refSpot, 4, 
                                                      coSpot_dist, isToPrintFit)
            } else {
                color4 <- chooseAColorFromGaussianFit(refSpot, 4, 
                                                      coSpot_dist, isToPrintFit)
            }
        } 
        allColors[i] <- paste(color1, color2, color3, color4)
    }
    allColors
}


# For singlylLoneRefSpots (lonelyRefSpots that has 0 or 1 coSpot in each time), assign
# 1 color at each time point, either with FQ if numCoSpots = 1 or GaussianFit if numCoSpots = 0 
assign4ColorsFromFQOrGaussianFit = function(singlyLonelyRefSpots, coSpot_dist = 1, 
                                            minNumOfSpotFQFound = 1,
                                            isFor6mRNAs = FALSE, isToPrintFit = FALSE) {
    
    transF1234NumCoSpots_columnIndice <- which(names(singlyLonelyRefSpots) %in% 
                                    c("transF1NumCoSpots", "transF2NumCoSpots", "transF3NumCoSpots", "transF4NumCoSpots"))
    transF1234CoSpotsIndice_columnIndice <-  which(names(singlyLonelyRefSpots) %in% 
                                    c("transF1CoSpotsIndice", "transF2CoSpotsIndice", "transF3CoSpotsIndice", "transF4CoSpotsIndice"))
    transF1234Spots_list <- list(transF1Spots, transF2Spots, transF3Spots, transF4Spots)

    # Go throuhg every spot
    allColors <- NULL
    for (iSpot in seq_along(singlyLonelyRefSpots[,1])) {
        refSpot <- singlyLonelyRefSpots[iSpot,]
        
        # Go throuhg every time point
        aSpotColors <- NULL
        for (iTime in seq_along(transF1234Spots_list)) {            
            transF1234NumCoSpots <- refSpot[,transF1234NumCoSpots_columnIndice[iTime]]
            if (transF1234NumCoSpots >= minNumOfSpotFQFound) {
                transF1234Spots <- transF1234Spots_list[[iTime]]
                transF1234CoSpotsIndice <- refSpot[,transF1234CoSpotsIndice_columnIndice[iTime]]
                aSpotColors[iTime] <- transF1234Spots[transF1234CoSpotsIndice[[1]][1],]$Color
            } else {
                if (isFor6mRNAs) {
                    aSpotColors[iTime] <- chooseAColorFromGaussianFit6mRNAs(refSpot, time = iTime, #diffractionLimit, 
                                                                coSpot_dist, isToPrintFit)
                } else {
                    aSpotColors[iTime] <- chooseAColorFromGaussianFit(refSpot, time = iTime, #diffractionLimit, 
                                                          coSpot_dist, isToPrintFit)
                }
            }
        }
        allColors[iSpot] <- paste(aSpotColors, collapse = ' ')
    }
    allColors
}


# For up to one time point with 2 coSpots, other time points has 0 or 1 coSpots
# split spots when 2 coSpots happend, assign colors to each of the split spots
# GaussianFit if numCoSpots = 0 
# return all the spots including the split spots
assign4ColorsFromFQOrGaussianFit_splitSpot = function(singlyLonelyRefSpots, coSpot_dist = 1, 
                                            minNumOfSpotFQFound = 1,
                                            isFor6mRNAs = FALSE, isToPrintFit = FALSE) {
    
    transF1234NumCoSpots_columnIndice <- which(names(singlyLonelyRefSpots) %in% 
                                                   c("transF1NumCoSpots", "transF2NumCoSpots", "transF3NumCoSpots", "transF4NumCoSpots"))
    transF1234CoSpotsIndice_columnIndice <-  which(names(singlyLonelyRefSpots) %in% 
                                                       c("transF1CoSpotsIndice", "transF2CoSpotsIndice", "transF3CoSpotsIndice", "transF4CoSpotsIndice"))
    transF1234Spots_list <- list(transF1Spots, transF2Spots, transF3Spots, transF4Spots)

    # Go throuhg every spot
    newSpots <- spotDataFrame(singlyLonelyRefSpots)
    for (iSpot in seq_along(singlyLonelyRefSpots[,1])) {
        refSpot <- singlyLonelyRefSpots[iSpot,]
        
        # Go throuhg every time point
        split_spots <- refSpot
        num_split_spots <- 1
        split_spot_colors <- NULL
        
        for (iTime in seq_along(transF1234Spots_list)) {            
            transF1234NumCoSpots <- refSpot[,transF1234NumCoSpots_columnIndice[iTime]]
            if (transF1234NumCoSpots >= minNumOfSpotFQFound) {
                transF1234Spots <- transF1234Spots_list[[iTime]]
                transF1234CoSpotsIndice <- refSpot[,transF1234CoSpotsIndice_columnIndice[iTime]]
                coTransF1234CoSpots <- transF1234Spots[transF1234CoSpotsIndice[[1]],]
                currentColor <- coTransF1234CoSpots$Color
                if (transF1234NumCoSpots == 2) {
                    # Decide whether to include both spots (brihgt + dim), usually < 0.06 is considered dim,
                    # yet it still possibel to get two colors even when it's dim by using 
                    # Gaussin fit to find not only one spot but two spots.
                    distanceTwoSpots <- distance(coTransF1234CoSpots[1,], coTransF1234CoSpots[2,])
                    
                    split_spots <- rbind(refSpot, refSpot)
                    split_spots[,1:2] <- coTransF1234CoSpots[,1:2]
                    # change
                    split_spots$coBrightnes <- coTransF1234CoSpots$normalTiffAmp
                    num_split_spots <- 2
                    split_spot_colors <- c(split_spot_colors, split_spot_colors)
                }
            } else {
                # so far assume Gaussian only return one spot
                if (isFor6mRNAs) {
                    currentColor <- chooseAColorFromGaussianFit6mRNAs(refSpot, time = iTime, #diffractionLimit, 
                                                                            coSpot_dist, isToPrintFit)
                } else {
                    currentColor <- chooseAColorFromGaussianFit(refSpot, time = iTime, #diffractionLimit, 
                                                                      coSpot_dist, isToPrintFit)
                }
            }
            sep <- ifelse(iTime == 1, "", " ")
            split_spot_colors <- paste(split_spot_colors, currentColor, sep = sep)
        }
        split_spots$allColors <- split_spot_colors
        newSpots <- rbind(newSpots, split_spots)
        
    }
    newSpots
}


# Given the stacked color of spots as a vector ("orange yellow red gree" ... for example)
# and a dataframe to mRNA as colume 1 and corresponding color as colume 2
# print out the color, count, and mRNA
printCountOfMRNA = function(allColors, mRNAsToColor) {
    allColors <- unlist(allColors)
    uniqueGoodColors <- unique(allColors)
    for (color in uniqueGoodColors) {
        num <- sum(allColors == color)
        mRNAName <- as.character(mRNAToColor$name[mRNAToColor$allColors == color]) 
        print(paste(formatC(color, width = 28), formatC(num, width = 5), 
                    formatC(mRNAName, width = 10), sep = ", "))
    }
}


# Given refSpots and the time point, use transF1CospotIndice (if time is 1) to get 
# the closest transF1Spots for each refspot and choose the truely closest one.
# return dataFrame of transF1Spots that might truely overlap with refSpots.
getClosestTransCospots = function(refSpots, time) {
    if (time == 1) {
        coSpotIndiceList <- refSpots$transF1CoSpotsIndice
        transSpots <- transF1Spots
    } else if (time == 2) {
        coSpotIndiceList <- refSpots$transF2CoSpotsIndice 
        transSpots <- transF2Spots
    } else if (time == 3) {
        coSpotIndiceList <- refSpots$transF3CoSpotsIndice
        transSpots <- transF3Spots
    } else if (time == 4) {
        coSpotIndiceList <- refSpots$transF4CoSpotsIndice 
        transSpots <- transF4Spots
    }
    closestTransCospots <- spotDataFrame(transSpots)
    for (r in seq_along(refSpots[,1])) {
        refSpot <- refSpots[r, ]
        coSpotIndice <- coSpotIndiceList[[r]]
        if (length(coSpotIndice) >= 1) {
            closestTransCospot <- closestSpot(refSpot, transSpots[coSpotIndice,])
            closestTransCospots <- rbind(closestTransCospots, closestTransCospot)
        }
    }
    closestTransCospots
}


# Get the maximum minusBgTiffAmp for each color at time x
# the value of each color can't exist maxAllowedMinusBgTiffAmp
maxMinusBgTiffAmpPerColor = function(maxAllowedMinusBgTiffAmp, time) {
    if (time == 1) {
        g1MinusBgTiffAmp <- closestTransF1Cospots[closestTransF1Cospots$Color == "green",]$minusBgTiffAmp
        y1MinusBgTiffAmp <- closestTransF1Cospots[closestTransF1Cospots$Color == "yellow",]$minusBgTiffAmp
        o1MinusBgTiffAmp <- closestTransF1Cospots[closestTransF1Cospots$Color == "orange",]$minusBgTiffAmp
        r1MinusBgTiffAmp <- closestTransF1Cospots[closestTransF1Cospots$Color == "red",]$minusBgTiffAmp
        
        g1MaxMinusBgTiffAmp <- max(g1MinusBgTiffAmp[g1MinusBgTiffAmp < maxAllowedMinusBgTiffAmp[1]])
        y1MaxMinusBgTiffAmp <- max(y1MinusBgTiffAmp[y1MinusBgTiffAmp < maxAllowedMinusBgTiffAmp[2]])
        o1MaxMinusBgTiffAmp <- max(o1MinusBgTiffAmp[o1MinusBgTiffAmp < maxAllowedMinusBgTiffAmp[3]])
        r1MaxMinusBgTiffAmp <- max(r1MinusBgTiffAmp[r1MinusBgTiffAmp < maxAllowedMinusBgTiffAmp[4]])
        maxMinusBgTiffAmpPerColor <- c(g1MaxMinusBgTiffAmp, y1MaxMinusBgTiffAmp, o1MaxMinusBgTiffAmp, r1MaxMinusBgTiffAmp)
    
    } else if (time == 2) {
        g2MinusBgTiffAmp <- closestTransF2Cospots[closestTransF2Cospots$Color == "green",]$minusBgTiffAmp
        y2MinusBgTiffAmp <- closestTransF2Cospots[closestTransF2Cospots$Color == "yellow",]$minusBgTiffAmp
        o2MinusBgTiffAmp <- closestTransF2Cospots[closestTransF2Cospots$Color == "orange",]$minusBgTiffAmp
        r2MinusBgTiffAmp <- closestTransF2Cospots[closestTransF2Cospots$Color == "red",]$minusBgTiffAmp
        
        g2MaxMinusBgTiffAmp <- max(g2MinusBgTiffAmp[g2MinusBgTiffAmp < maxAllowedMinusBgTiffAmp[1]])
        y2MaxMinusBgTiffAmp <- max(y2MinusBgTiffAmp[y2MinusBgTiffAmp < maxAllowedMinusBgTiffAmp[2]])
        o2MaxMinusBgTiffAmp <- max(o2MinusBgTiffAmp[o2MinusBgTiffAmp < maxAllowedMinusBgTiffAmp[3]])
        r2MaxMinusBgTiffAmp <- max(r2MinusBgTiffAmp[r2MinusBgTiffAmp < maxAllowedMinusBgTiffAmp[4]])   
        maxMinusBgTiffAmpPerColor <- c(g2MaxMinusBgTiffAmp, y2MaxMinusBgTiffAmp, o2MaxMinusBgTiffAmp, r2MaxMinusBgTiffAmp)
        
    } else if (time == 3) {
        g3MinusBgTiffAmp <- closestTransF3Cospots[closestTransF3Cospots$Color == "green",]$minusBgTiffAmp
        y3MinusBgTiffAmp <- closestTransF3Cospots[closestTransF3Cospots$Color == "yellow",]$minusBgTiffAmp
        o3MinusBgTiffAmp <- closestTransF3Cospots[closestTransF3Cospots$Color == "orange",]$minusBgTiffAmp
        r3MinusBgTiffAmp <- closestTransF3Cospots[closestTransF3Cospots$Color == "red",]$minusBgTiffAmp
        
        g3MaxMinusBgTiffAmp <- max(g3MinusBgTiffAmp[g3MinusBgTiffAmp < maxAllowedMinusBgTiffAmp[1]])
        y3MaxMinusBgTiffAmp <- max(y3MinusBgTiffAmp[y3MinusBgTiffAmp < maxAllowedMinusBgTiffAmp[2]])
        o3MaxMinusBgTiffAmp <- max(o3MinusBgTiffAmp[o3MinusBgTiffAmp < maxAllowedMinusBgTiffAmp[3]])
        r3MaxMinusBgTiffAmp <- max(r3MinusBgTiffAmp[r3MinusBgTiffAmp < maxAllowedMinusBgTiffAmp[4]])
        maxMinusBgTiffAmpPerColor <- c(g3MaxMinusBgTiffAmp, y3MaxMinusBgTiffAmp, o3MaxMinusBgTiffAmp, r3MaxMinusBgTiffAmp)

    } else if (time == 4) {
        g4MinusBgTiffAmp <- closestTransF4Cospots[closestTransF4Cospots$Color == "green",]$minusBgTiffAmp
        y4MinusBgTiffAmp <- closestTransF4Cospots[closestTransF4Cospots$Color == "yellow",]$minusBgTiffAmp
        o4MinusBgTiffAmp <- closestTransF4Cospots[closestTransF4Cospots$Color == "orange",]$minusBgTiffAmp
        r4MinusBgTiffAmp <- closestTransF4Cospots[closestTransF4Cospots$Color == "red",]$minusBgTiffAmp
        
        g4MaxMinusBgTiffAmp <- max(g4MinusBgTiffAmp[g4MinusBgTiffAmp < maxAllowedMinusBgTiffAmp[1]])
        y4MaxMinusBgTiffAmp <- max(y4MinusBgTiffAmp[y4MinusBgTiffAmp < maxAllowedMinusBgTiffAmp[2]])
        o4MaxMinusBgTiffAmp <- max(o4MinusBgTiffAmp[o4MinusBgTiffAmp < maxAllowedMinusBgTiffAmp[3]])
        r4MaxMinusBgTiffAmp <- max(r4MinusBgTiffAmp[r4MinusBgTiffAmp < maxAllowedMinusBgTiffAmp[4]]) 
        maxMinusBgTiffAmpPerColor <- c(g4MaxMinusBgTiffAmp, y4MaxMinusBgTiffAmp, o4MaxMinusBgTiffAmp, r4MaxMinusBgTiffAmp)
    }
    maxMinusBgTiffAmpPerColor
}



# Get the local maxium of a tiff Image (which is a pixel with higher intensity
# than all it's 3 neighbors)
# Return 2D indice which are local maximum.
# Note, we are no counting all the side pixel since they don't have 8 neighbors.
getLocalMaximumIndice2D = function(tiffMatrix, spot = NA, nearDistance_um = 2000) {
    indice2D <- NULL
    surroundingIndexDiffJ <- c(-1, 0, 1, 1, 1, 0, -1, -1)
    surroundingIndexDiffI <- c(1, 1, 1, 0, -1, -1, -1, 0)
    minI <- 2
    maxI <- dim(tiffMatrix)[1] - 1
    minJ <- 2
    maxJ <- dim(tiffMatrix)[2] - 1
    if (!is.na(spot)[1]) {
        indice2D <- spotsToIndice2D(spot, pix_XY)
        midI <- getYIndiceFromIndice2D(indice2D)
        midJ <- getXIndiceFromIndice2D(indice2D)
        nearDistance_pix <- floor(nearDistance_um / pix_XY)
        minI <- max(minI, midI - nearDistance_pix)
        maxI <- min(maxI, midI + nearDistance_pix)
        minJ <- max(minJ, midJ - nearDistance_pix)
        maxJ <- min(maxJ, midJ + nearDistance_pix)
    }
    
    for(i in minI : maxI) {
        for (j in minJ : maxJ) {
            centerIntensity <- tiffMatrix[i, j]
            isCenterMaxium <- TRUE
            for(k in 1: length(surroundingIndexDiffJ)) {
                nearbyIntensity <- tiffMatrix[i + surroundingIndexDiffI[k],
                                              j + surroundingIndexDiffJ[k]]
                if (centerIntensity <  nearbyIntensity) {
                    isCenterMaxium <- FALSE
                    break
                }
            }
            if (isCenterMaxium) {
                indice2D <- rbind(indice2D, c(i, j))
            }
        }
    }
    indice2D
}


# from indice2D, creat fake spot with Pos_X, Pos_Y as center of indice,
# and tiffAmp from tiffMatrix
createSpotsFrom2DIndice = function(tiffMatrix, indice2D) {
    xIndice <- getXIndiceFromIndice2D(indice2D)
    yIndice <- getYIndiceFromIndice2D(indice2D)
    Pos_X <- (xIndice - 0.5) * pix_XY
    Pos_Y <- (yIndice - 0.5) * pix_XY
    tiffAmp <- tiffMatrix[indice2D]
    spots <- data.frame(Pos_Y = Pos_Y, Pos_X = Pos_X, tiffAmp = tiffAmp)
    spots
}


# Giving a vector of values, give the single top X percentile tile value
# i.e. if there are 100 values, get the 10th largest value.
# Do this instead of the max value to avoid outliers
getTopPercentValue = function (values, percentage) {
    if (length(values) == 0) {
        highValue <- max(values)
    } else {
        orderedValues <- values[order(values, decreasing = TRUE)]
        ind <- ceil(percentage / 100 * length(values))
        highValue <- orderedValues[ind]
    }
    highValue
}


# Fill in the empty coSpot if there is a close unused coSpot within 1 um distance.
# Fill in by changing the column of transFNNumCoSpots and transFNCoSpotsIndice
fillInEmptyCospotsWithUnusedSpots = function (spots, secondaryCospotDist, refSpots) {
    iColumn_transFNNumCoSpots <- which(grepl("NumCoSpots", names(spots)))
    iColumn_transF1CoSpotsIndice <- which(grepl("CoSpotsIndice", names(spots)))
    transFNSpots <- list(transF1Spots, transF2Spots, transF3Spots, transF4Spots)
    for (i in 1:length(iColumn_transFNNumCoSpots)) {
        indexOfZeroSpots <- which(spots[,iColumn_transFNNumCoSpots[i]] == 0)
        interestedSpots <- spots[indexOfZeroSpots,]
        closestSpotsInFN <- closestSpot(interestedSpots, transFNSpots[[i]], countOnlyUnusedSpots = TRUE)
        closestSpotIndexInFN <- closestSpotIndex(interestedSpots, transFNSpots[[i]], countOnlyUnusedSpots = TRUE)        
        
        # test if the closest spot to closestSpotsInFN is interestedSpot, if not, disgard
        closestSpotInSpots <- closestSpot(closestSpotsInFN, refSpots)
        # test if the distance between the interestedSpot and closestSpotsInFN is < secondaryCospotDist, if not, disgar
        iDistanceSmall <- which(distance(interestedSpots, closestSpotsInFN) < secondaryCospotDist &
                                    closestSpotInSpots$Pos_X == interestedSpots$Pos_X &
                                    closestSpotInSpots$Pos_Y == interestedSpots$Pos_Y)        
        
        interestedSpots[iDistanceSmall, iColumn_transFNNumCoSpots[i]] <- 1
        interestedSpots[iDistanceSmall, iColumn_transF1CoSpotsIndice[i]] <- closestSpotIndexInFN[iDistanceSmall]
        spots[indexOfZeroSpots, ] <- interestedSpots
    }
    spots
}