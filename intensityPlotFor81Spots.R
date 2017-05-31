rm(list = ls())
setwd("/Users/jackie/Downloads/SGTC microscope/2016-07-22-81 strands Atto sequence color/Analyze/");
source("intensityPlotFor81SpotsFunctions.R")

# 1. read a) colorSpotColor for 81 spots which is a matrix 
# b) reference colorSpots, and
# c) all 16 color spots to get pureIntensiy as avgIntensity - bgIntensity
# For whole spot, sigma of pureIntensity is 0.7066744
# For 1 pixel,  sigma of pureIntensity is 0.9329983
correctSpotColors <- readCorrectSpotColors("correctSpotColors.csv")

G25 <- readColorSpots("intensityA2 Amine G25-81strands Gfilter 2s.csv");

G1 <- readColorSpots("intensityD2 Amine NG25-F1-NF234 again Gfilter 2s.csv");
O1 <- readColorSpots("intensityD2 Amine NG25-F1-NF234 again Ofilter 2s.csv");
R1 <- readColorSpots("intensityD2 Amine NG25-F1-NF234 again Rfilter 2s.csv");
Y1 <- readColorSpots("intensityD2 Amine NG25-F1-NF234 again Yfilter 1s.csv");

G2 <- readColorSpots("intensityE2 Amine F2-NF134 Gfilter 2s.csv");
O2 <- readColorSpots("intensityE2 Amine F2-NF134 Ofilter 2s.csv");
R2 <- readColorSpots("intensityE2 Amine F2-NF134 Rfilter 2s.csv");
Y2 <- readColorSpots("intensityE2 Amine F2-NF134 Yfilter 1s.csv");

G3 <- readColorSpots("intensityF2 Amine F3-NF124 Gfilter 2s.csv");
O3 <- readColorSpots("intensityF2 Amine F3-NF124 Ofilter 2s.csv");
R3 <- readColorSpots("intensityF2 Amine F3-NF124 Rfilter 2s.csv");
Y3 <- readColorSpots("intensityF2 Amine F3-NF124 Yfilter 1s.csv");

G4 <- readColorSpots("intensityG2 Amine F4-NF123 Gfilter 2s.csv");
O4 <- readColorSpots("intensityG2 Amine F4-NF123 Ofilter 2s.csv");
R4 <- readColorSpots("intensityG2 Amine F4-NF123 Rfilter 2s.csv");
Y4 <- readColorSpots("intensityG2 Amine F4-NF123 Yfilter 1s.csv");


G25 <- readColorSpots("intensityA2 Amine G25-81strands Gfilter 2s small.csv");

G1 <- readColorSpots("intensityD2 Amine NG25-F1-NF234 again Gfilter 2s small.csv");
O1 <- readColorSpots("intensityD2 Amine NG25-F1-NF234 again Ofilter 2s small.csv");
R1 <- readColorSpots("intensityD2 Amine NG25-F1-NF234 again Rfilter 2s small.csv");
Y1 <- readColorSpots("intensityD2 Amine NG25-F1-NF234 again Yfilter 1s small.csv");

G2 <- readColorSpots("intensityE2 Amine F2-NF134 Gfilter 2s small.csv");
O2 <- readColorSpots("intensityE2 Amine F2-NF134 Ofilter 2s small.csv");
R2 <- readColorSpots("intensityE2 Amine F2-NF134 Rfilter 2s small.csv");
Y2 <- readColorSpots("intensityE2 Amine F2-NF134 Yfilter 1s small.csv");

G3 <- readColorSpots("intensityF2 Amine F3-NF124 Gfilter 2s small.csv");
O3 <- readColorSpots("intensityF2 Amine F3-NF124 Ofilter 2s small.csv");
R3 <- readColorSpots("intensityF2 Amine F3-NF124 Rfilter 2s small.csv");
Y3 <- readColorSpots("intensityF2 Amine F3-NF124 Yfilter 1s small.csv");

G4 <- readColorSpots("intensityG2 Amine F4-NF123 Gfilter 2s small.csv");
O4 <- readColorSpots("intensityG2 Amine F4-NF123 Ofilter 2s small.csv");
R4 <- readColorSpots("intensityG2 Amine F4-NF123 Rfilter 2s small.csv");
Y4 <- readColorSpots("intensityG2 Amine F4-NF123 Yfilter 1s small.csv");

sigmasOfPureIntensity <- getSigmasForPositiveColors(listOfColorSpots, "pureIntensity", correctSpotColors)
mean(sigmasOfPureIntensity, na.rm = TRUE)

# 2. Add relativeToRefIntensity to 16 spots as pureIntensity / refSpot's pureIntensity
# Note: relativeToRefIntensity is suppose to correct spot difference, but it's
# not perfect because ref DNA and encoder DNAs have different binding efficiency.
# For whole spots, sigma of relativeToRefIntensity is 0.3556208.
# For 1 pixel, sigma of relativeToRefIntensity is 0.4918928
listOfColorSpots <- list(G1, O1, R1, Y1, G2, O2, R2, Y2, G3, O3, R3, Y3, G4, O4, R4, Y4)
listOfColorSpots <- addRelativeToRefIntensity(listOfColorSpots, G25)

plotGORY2(listOfColorSpots, "relativeToRefIntensity")
plotGORY2(listOfColorSpots, "relativeToRefIntensity", washPoint = 2)
plotGORY2(listOfColorSpots, "relativeToRefIntensity", washPoint = 3)
plotGORY2(listOfColorSpots, "relativeToRefIntensity", washPoint = 4)
plotGORYForAll(listOfColorSpots, "relativeToRefIntensity")

sigmasOfRelativeToRefIntensity <- getSigmasForPositiveColors(listOfColorSpots, "relativeToRefIntensity", correctSpotColors)
mean(sigmasOfRelativeToRefIntensity, na.rm = TRUE)


# 3. add colorCalibratedRelativeToRefIntensity to 16 colorSpots as
# relativeToRefIntensity / normalizationFactorForColor
# where the normalizationFactorForColor is the avg of top 5-10% relativeToRefIntensity 
# of a color if the color is fluor. Non-fluor color's nf is estimated from the fluor color.
# Note: normalizationFactorForColor is supposed to correct color difference
# For whole spots, sigma of colorCalibratedRelativeToRefIntensity from relativeToRefIntensity is 0.3556208.
# For 1 pixel, sigma of colorCalibratedRelativeToRefIntensity from relativeToRefIntensity is 0.4918928.
# For whole spots,  sigma of colorCalibratedRelativeToRefIntensity from pureIntensity is 0.7066744.
# For 1 pixel,  sigma of colorCalibratedRelativeToRefIntensity from pureIntensity is 0.9329983
normalizationFactorsForColor <- getNormalizationFactorsForColor(listOfColorSpots, 
                                                                50, "pureIntensity")
listOfColorSpots <- addColorCalibratedRelativeToRefIntensity(listOfColorSpots, 
                                                             normalizationFactorsForColor,
                                                             "pureIntensity")

plotGORY2(listOfColorSpots, "colorCalibratedRelativeToRefIntensity")
plotGORY2(listOfColorSpots, "colorCalibratedRelativeToRefIntensity", washPoint = 2)
plotGORY2(listOfColorSpots, "colorCalibratedRelativeToRefIntensity", washPoint = 3)
plotGORY2(listOfColorSpots, "colorCalibratedRelativeToRefIntensity", washPoint = 4)
plotGORYForAll(listOfColorSpots, "colorCalibratedRelativeToRefIntensity")

sigmasOfColorCalibratedRelativeToRefIntensity <- getSigmasForPositiveColors(listOfColorSpots, "colorCalibratedRelativeToRefIntensity", correctSpotColors)
mean(sigmasOfColorCalibratedRelativeToRefIntensity, na.rm = TRUE)


# 4. add spotNormalizedIntensity to 16 colorSpots so the 16 colors of a spot
# is normalized to be 1, i.e. color1^2 + .... color16^12 = 1
# Note: this almost correct spot difference perfectly
# For whole spots, sigmas of spotNormalizedIntensity from relativeToRefIntensity is 0.1270563.
# for 1 pixel, sigmas of spotNormalizedIntensity from relativeToRefIntensity is 0.2186318.
# For whole spots, sigmas of spotNormalizedIntensity from pureIntensity is 0.1389537.
# For 1 pixel, sigmas of spotNormalizedIntensity from pureIntensity is 0.2323793
listOfColorSpots <- addSpotNormalizedIntensities(listOfColorSpots)
plotGORY2(listOfColorSpots, "spotNormalizedIntensity")
plotGORY2(listOfColorSpots, "spotNormalizedIntensity", washPoint = 2)
plotGORY2(listOfColorSpots, "spotNormalizedIntensity", washPoint = 3)
plotGORY2(listOfColorSpots, "spotNormalizedIntensity", washPoint = 4)
plotGORYForAll(listOfColorSpots, "spotNormalizedIntensity")

sigmasOfSpotNormalizedIntensity <- getSigmasForPositiveColors(listOfColorSpots, "spotNormalizedIntensity", correctSpotColors)
mean(sigmasOfSpotNormalizedIntensity, na.rm = TRUE)


# 5. Is there a good threshold for spotNormalizedIntensity? Maybe 0.15
# Assume the best scenario is 0.5 for the 4 perfect colors and 0 for rest 12 colors
# The minimin 4th top colorIntenties of 81 spots is 0.18 for whole spots
# The maximum 5th top colorIntenties of 81 spots is 0.11 for whole spots
spotToColorIntensities <- getMatrixOfSpotToColorIntensities(
    listOfColorSpots, "spotNormalizedIntensity")

sortedIntensities <- t(apply(spotToColorIntensities, 1, sort, decreasing = TRUE))
min(sortedIntensities[,4])
max(sortedIntensities[,5])


# 6. Get the choosenSpotColors where the color is above 0.15 threshold
# and find out spots that don't match the correctSpotColors
# For 1 pixel spots using relativeToRefIntensity for color normalization, 
# spot 23, 71, 76 has 1 false postive color, spot 12 has 1 false negative color
#   bit 1 -> 0 error rate = 1 / (81*4) = 0.31%
#   bit 0 -> 1 error rate = 3 / (81*8) = 0.46%, not 81*12 because only have 3*4 = 12 fluors.
# For 1 pixel spots using pureIntensity for color normalization, 
# spot 64, 71 has 1 false postive color, spot 12 has 1 false negative color
#   bit 1 -> 0 error rate = 1 / (81*4) = 0.31%
#   bit 0 -> 1 error rate = 2 / (81*8) = 0.31%
#   they are all correctable if using the HD4 sequence
# for whole spots, the bit error rate is 0%.
# Note, we can get good color normalization factors from pureIntensity, 
# so we don't need to use reference spots to get relativeToRefIntensity 
# for color normalization factors. 
intensityThreshold <- 0.15
spotWithExtraColor <- getSpotWithExtraColor(spotToColorIntensities, 
                                            correctSpotColor, intensityThreshold)
spotWithDeficitColor <- getSpotWithDeficitColor(spotToColorIntensities, 
                                               correctSpotColor, intensityThreshold)
spotWithExtraColor
spotWithDeficitColor
nrow(spotWithExtraColor) / (81*8)
nrow(spotWithDeficitColor) / (81*4)


# -------------------------------------------------------------------
# Below is for no error corrrection i.e. directly pick strongest color per wash
# and for 80-20 mix of spots

# Check out the first image with green color (27 spots), then check out the 
# second image for other colors, get relative intensity of other color over
# the green color to see if the dominant color is correct.
G1Index <- which(correctSpotColors[,"firstIndex"] == 1);
G2$relativeToG1 <- G2$pureIntensity / G1$pureIntensity;
O2$relativeToG1 <- O2$pureIntensity / G1$pureIntensity;
R2$relativeToG1 <- R2$pureIntensity / G1$pureIntensity;
Y2$relativeToG1 <- Y2$pureIntensity / G1$pureIntensity;

barPlotData <- rbind(G2$relativeToG1, O2$relativeToG1, 
                     R2$relativeToG1, Y2$relativeToG1);
portion <- G1Index;
maxYLimit <- max(barPlotData[,portion]);
par(mfrow = c(1,1));
barplot(barPlotData[,portion], col = c("green", "darkorange", "red", "yellow"), 
        beside = TRUE, names.arg = portion, ylim = c(-0.2, maxYLimit));


# 2. Calculate intensity ratio of the same color between two time points
# It's interesting that intensity ratio goes up and up as time goes by!
# Maybe it's because the laser is brighter?  or more because the adapter / fluor hybridization better?
# True for amine slides too.
colors <- c("g", "o", "r", "y")
colors1 <- c("g", "g", "g", "y", "y", "y", "o", "o", "o",
             "g", "g", "g", "y", "y", "y", "o", "o", "o",
             "g", "g", "g", "y", "y", "y", "o", "o", "o",
             "g", "g", "g", "y", "y", "y", "o", "o", "o",
             "g", "g", "g", "y", "y", "y", "o", "o", "o",
             "g", "g", "g", "y", "y", "y", "o", "o", "o",
             "g", "g", "g", "y", "y", "y", "o", "o", "o",
             "g", "g", "g", "y", "y", "y", "o", "o", "o",
             "g", "g", "g", "y", "y", "y", "o", "o", "o")

colors2 <- c("g", "y", "r", "g", "y", "r", "g", "y", "r",
             "g", "y", "r", "g", "y", "r", "g", "y", "r",
             "g", "y", "r", "g", "y", "r", "g", "y", "r",
             "g", "y", "r", "g", "y", "r", "g", "y", "r",
             "g", "y", "r", "g", "y", "r", "g", "y", "r",
             "g", "y", "r", "g", "y", "r", "g", "y", "r",
             "g", "y", "r", "g", "y", "r", "g", "y", "r",
             "g", "y", "r", "g", "y", "r", "g", "y", "r",
             "g", "y", "r", "g", "y", "r", "g", "y", "r")

colors3 <- c("r", "r", "r", "r", "r", "r", "r", "r", "r",
             "r", "r", "r", "r", "r", "r", "r", "r", "r",
             "r", "r", "r", "r", "r", "r", "r", "r", "r",
             "o", "o", "o", "o", "o", "o", "o", "o", "o",
             "o", "o", "o", "o", "o", "o", "o", "o", "o",
             "o", "o", "o", "o", "o", "o", "o", "o", "o",
             "g", "g", "g", "g", "g", "g", "g", "g", "g",
             "g", "g", "g", "g", "g", "g", "g", "g", "g",
             "g", "g", "g", "g", "g", "g", "g", "g", "g")

colors4 <- c("r", "r", "r", "r", "r", "r", "r", "r", "r",
             "o", "o", "o", "o", "o", "o", "o", "o", "o",
             "y", "y", "y", "y", "y", "y", "y", "y", "y", 
             "r", "r", "r", "r", "r", "r", "r", "r", "r",
             "o", "o", "o", "o", "o", "o", "o", "o", "o",
             "y", "y", "y", "y", "y", "y", "y", "y", "y", 
             "r", "r", "r", "r", "r", "r", "r", "r", "r",
             "o", "o", "o", "o", "o", "o", "o", "o", "o",
             "y", "y", "y", "y", "y", "y", "y", "y", "y")

wash1 <- data.frame(color = colors1, gAvg = G1$avgIntensity, gBg = G1$bgIntensity, 
                 yAvg = Y1$avgIntensity, yBg = Y1$bgIntensity,
                 oAvg = O1$avgIntensity, oBg = O1$bgIntensity,
                 rAvg = R1$avgIntensity, rBg = R1$bgIntensity)
wash2 <- data.frame(color = colors2, gAvg = G2$avgIntensity, gBg = G2$bgIntensity, 
                 yAvg = Y2$avgIntensity, yBg = Y2$bgIntensity,
                 oAvg = O2$avgIntensity, oBg = O2$bgIntensity,
                 rAvg = R2$avgIntensity, rBg = R2$bgIntensity)
wash3 <- data.frame(color = colors3, gAvg = G3$avgIntensity, gBg = G3$bgIntensity, 
                 yAvg = Y3$avgIntensity, yBg = Y3$bgIntensity,
                 oAvg = O3$avgIntensity, oBg = O3$bgIntensity,
                 rAvg = R3$avgIntensity, rBg = R3$bgIntensity)
wash4 <- data.frame(color = colors4, gAvg = G4$avgIntensity, gBg = G4$bgIntensity, 
                 yAvg = Y4$avgIntensity, yBg = Y4$bgIntensity,
                 oAvg = O4$avgIntensity, oBg = O4$bgIntensity,
                 rAvg = R4$avgIntensity, rBg = R4$bgIntensity)

printSequentialIntensityRatio = function(wash1, wash2) {
    intensityRatio <- data.frame()
    for (i in seq_along(colors)) {
        color <- colors[i]
        commonIndex <- which(wash1$color == color & wash2$color == color)
        avgIntensity1 <- wash1[commonIndex,i*2]
        bgIntensity1 <- wash1[commonIndex, i*2 + 1]
        avgIntensity2 <- wash2[commonIndex,i*2]
        bgIntensity2 <- wash2[commonIndex, i*2 + 1]
        ratio <- (avgIntensity2 - bgIntensity2) / (avgIntensity1 - bgIntensity1)
        print(color)
        print(ratio)
    }
}

printSequentialIntensityRatio(wash1, wash2)
printSequentialIntensityRatio(wash2, wash3)
printSequentialIntensityRatio(wash3, wash4)
printSequentialIntensityRatio(wash1, wash3)
printSequentialIntensityRatio(wash2, wash4)
printSequentialIntensityRatio(wash1, wash4)


# 3. Calculate whether 80%, 20% mixture will work
# Yes, they work except mixing wash14 and wash 24 because fluorescent intensities are brighter
# at later wash (mistery...), so the mix is not 80% 20% but might be more like 60% 40%.
percentage = 0.8
isToNormal = TRUE
time1ColorSpots <- listOfColorSpots[1:4]
time2ColorSpots <- listOfColorSpots[5:8]
time3ColorSpots <- listOfColorSpots[9:12]
time4ColorSpots <- listOfColorSpots[13:16]

plotGORY_mixer(time1ColorSpots, time2ColorSpots, plotnameSuffix = "12")
plotGORY_mixer(wash1, wash3, "13")
# plotGORY_mixer(wash1, wash4, "14") 
# won't work since the color become brighter and brighter, so it's like 80% of weak color + 20% of strong color
# for example, the y brightness become several times bigger at the same spot

plotGORY_mixer(wash2, wash1, "21")
plotGORY_mixer(wash2, wash3, "23")
# plotGORY_mixer(wash2, wash4, "24") also won't work because wash 4 is much brighter than wash 2

plotGORY_mixer(wash3, wash1, "31")
plotGORY_mixer(wash3, wash2, "32")
plotGORY_mixer(wash3, wash4, "34")

plotGORY_mixer(wash4, wash1, "41")
plotGORY_mixer(wash4, wash2, "42")
plotGORY_mixer(wash4, wash3, "43")









    
    
