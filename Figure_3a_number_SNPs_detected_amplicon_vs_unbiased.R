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
# Determine how many SNPs are deteced in each dilution by amplicon vs. unbiased
###############################################################################

# read in the combined dataframe. 
data=read.table("~/Documents/CA04_SNP_detection/data/unbiased_data/amplicon_vs_unbiased.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("reference_position","type","length","reference","allele","linkage","zygosity","count","coverage","frequency","hyper-allelic","fr_balance","average_quality","input","replicate","minimum_coverage_required","method","dups_removed"))

# remove indels:
data <- data[data$type != "Deletion",]
data <- data[data$type != "Insertion",]

# remove SNPs called for minimum coverage required = 1 and dups removed
data <- data[data$minimum_coverage_required == "100",]
data <- data[data$dups_removed == "no",]

# or, to compare unbiased vs. duplicates removed, do this: 
# data2 <- data[data$dups_removed == "yes",]
# data3 <- data[data$method == "unbiased",]
# data <- rbind(data2, data3)

# for this, need to convert the input to a character
data$input = as.character(as.numeric(data$input))
data$input <- gsub("10000","1e+04",data$input)
data$input <- gsub("1000","1e+03",data$input)
data$input <- gsub("100","1e+02",data$input)

# reformat dataframe so that SNP frequencies can be compared between amplicon and unbiased methods
# select only the columns: reference_position, frequency, input, replicate and method
df <- data[c(1,10, 14:18)]

# separate out unbiased and amplicon into 2 separate dataframes, then merge back together and convert all NAs to 0s
# method.x is unbiased, method.y is ampilcon
df1 = df[df$method =="unbiased",]
df2 = df[df$method =="amplicon",]

df = merge(df1, df2, by=c("input","replicate", "reference_position"), all=TRUE)
df$frequency.x[is.na(df$frequency.x)] <- 0
df$frequency.y[is.na(df$frequency.y)] <- 0
df$method.x[is.na(df$method.x)] <- "unbiased"
df$method.y[is.na(df$method.y)] <- "amplicon"


# remove all dilutions 10^4 and less
df = df[df$input == "1e+07" | df$input == "1e+06" | df$input == "1e+05",]


# calculate the number of total SNPs in each dilution, the fraction of SNPs that are detected by both amplicon and unbiased, and the fraction only found in either 1 or the other

# run only on variants > 10% frequency in the population
#df = df[df$frequency.x >= 10 | df$frequency.y >= 10,]


dilution_list = unique(df$input)
replicate_list = unique(df$replicate)
sites_list = unique(df$reference_position)
total_SNPs = lengthunique(df$reference_position)
count = 0
count2 = 0

for (d in dilution_list)
{
  df1 <- df[df$input == d,]
  total_SNPs = nrow(df1)
  
  # define sites_list as the sites that were detected in any replicates in that dilution
  sites_list = unique(df1$reference_position)
  
  # define counts for detected_in_both, unbiasesd_only and amplicon_only
  detected_in_both = 0
  unbiased_only = 0
  amplicon_only = 0
  
  for (s in sites_list)
  {
    df2 <- df1[df1$reference_position == s,]
    
    # determine whether that site was detected in ampicon, unbiased or both
    if (max(df2$frequency.x) > 0 & max(df2$frequency.y) > 0){
      detected_in_both = detected_in_both + 1
      when_detected = "both"
    
      } else if (max(df2$frequency.x) > 0 & max(df2$frequency.y) == 0) {
      unbiased_only = unbiased_only + 1
      when_detected = "unbiased_only"
      frequency = mean(df2$frequency.x)
      
    } else if (max(df2$frequency.x) == 0 & max(df2$frequency.y) > 0) {
      amplicon_only = amplicon_only + 1
      when_detected = "amplicon_only"
    }
    
    df4 = df[df$input == d & df$reference_position == s, ]
    df4["when_detected"] <- when_detected
    
    if (count2 == 0){
      when_detected_df = df4
      count2 = count2 + 1
    } else {
      when_detected_df= rbind(df4, when_detected_df)
    }
    
  }
  
  total_sites = amplicon_only + unbiased_only + detected_in_both
  df3 = data.frame(d,amplicon_only,unbiased_only,detected_in_both, total_sites)
  
  if (count == 0){
    total_counts_df = df3
    count = count + 1
  } else {
    total_counts_df = rbind(total_counts_df,df3)
  }
}     


# now add in columns with the percentages
total_counts_df$percent_amplicon_only = (total_counts_df$amplicon_only/total_counts_df$total_sites)*100
total_counts_df$percent_unbiased_only = (total_counts_df$unbiased_only/total_counts_df$total_sites)*100
total_counts_df$percent_both = (total_counts_df$detected_in_both/total_counts_df$total_sites)*100

write.table(total_counts_df, file = "~/Documents/CA04_SNP_detection/Figures/unbiased_data_figures/amplicon_vs_unbiased_total_counts_170418.txt", sep = "\t")

# write.table(total_counts_df, file = "~/Documents/CA04_SNP_detection/Figures/unbiased_data_figures/amplicon_vs_unbiased_dups_removed_total_counts_170511.txt", sep = "\t")


# I now want to take this output and only plot the frequencies of SNPs that are detected in only unbiased or only amplicon

not_detected_in_both = when_detected_df[when_detected_df$when_detected != "both",]
not_detected_in_both = not_detected_in_both[c(1:4,8,12)]

# melt the dataframe
not_detected_in_both_melted <- melt(not_detected_in_both, id=c("input","replicate","reference_position","when_detected"))

df_melted = not_detected_in_both_melted[not_detected_in_both_melted$value != 0,]

ggplot(df_melted, aes(y=value, x=reference_position, color=input, alpha=0.9, shape=when_detected))+geom_point(size=2)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=13))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title=element_text(size=13, vjust=0.2))+theme(axis.text=element_text(size=11, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=11))+theme(legend.title=element_text(size=13, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="SNP frequency (%) (unbiased)",y="SNP frequency (%) (amplicon)")+scale_x_continuous(breaks=seq(0,1700,100), limits=c(0,1700))+scale_y_continuous(breaks=seq(0,20,5), limits=c(0,20))+scale_colour_manual(values=c(color4,color5,color6))+scale_shape_manual(values=c(1,4))
