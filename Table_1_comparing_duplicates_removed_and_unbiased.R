###############################################################################
# This includes the code for generating figures for the CA04 illumina manuscript
###############################################################################
require(ggplot2)
require(grid)
require(stringr)
require(reshape2)
require(format)
require(plyr)

# blue/red pallette
color1 <- rgb(0,148,194, maxColorValue=255)
color2 <- rgb(74,204,245, maxColorValue=255)
color3 <- rgb(0,90,117, maxColorValue=255)
color4 <- rgb(117,11,0, maxColorValue=255)
color5 <- rgb(207,19,0, maxColorValue=255)
color6 <- rgb(245,89,74, maxColorValue=255)
color7 <- rgb(255,255,255, maxColorValue=255)

###############################################################################
# compare unbiased and duplicate reads removed
###############################################################################

# 1. run Figure_X_number_SNPs_detected_amplicon_vs_unbiased

unbiased_only_SNPs <- sort(unique(when_detected_df$reference_position))

# 2. run Figure_X_number_SNPs_detected_amplicon_vs_unbiased, except keep only duplicates removed = yes


