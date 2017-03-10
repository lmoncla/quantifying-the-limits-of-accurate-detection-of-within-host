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
# Figure 4: SNPs detected more than 1 time means and standard deviations
###############################################################################

#### FIGURE 4A: MEAN AND STANDARD DEVIATION OF SNPS DETECTED >1 TIME ##########

# read in the combined dataframe. 
data=read.table("~/Documents/CA04_SNP_detection/data/combined_SNPs.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("reference_position","type","length","reference","allele","linkage","zygosity","count","coverage","frequency","hyper-allelic","fr_balance","average_quality","overlapping_annotations","coding_region_change","amino_acid_change","dups_removed","input","replicate"))

# remove indels:
data <- data[data$type != "Deletion",]
data <- data[data$type != "Insertion",]

# for this, need to convert the input to a character
data$input = as.character(as.numeric(data$input))
data$input <- gsub("10000","1e+04",data$input)
data$input <- gsub("1000","1e+03",data$input)
data$input <- gsub("100","1e+02",data$input)

# Subset the dataframe to remove polymorphic sites detected only once.  
# 1st step creates a list of sites that are duplicated; 2nd step subsets data to include only rows with those nonunique sites
nonunique <- subset(data$reference_position, duplicated(data$reference_position), )
nonunique = sort(nonunique)
nonunique = unique(nonunique)
data_nonunique <- data[data$reference_position %in% nonunique, ]

# Now calculate the mean and standard deviation for the nonunique data
dilution_list = unique(data_nonunique$input)
count = 0

data_nonunique$frequency = as.numeric(as.character(data_nonunique$frequency))

for (d in dilution_list)
{
  sites <- data_nonunique[data_nonunique$input == d,]
  site_list = unique(sites$reference_position)
  for (i in site_list)
  {
    df <- data_nonunique[data_nonunique$reference_position == i & data_nonunique$input == d,]
    if (nrow(df) == 1){
      all_values = list(df$frequency, 0, 0)
      times_detected = 1
      non_mean_frequency = df$frequency
    } else if (nrow(df) == 2) {
      all_values = list(df$frequency[1], df$frequency[2], 0)
      times_detected=2
      non_mean_frequency = NA
    }
    else if (nrow(df) == 3) {
      all_values = list(df$frequency[1], df$frequency[2], df$frequency[3])
      times_detected = 3
      non_mean_frequency = NA
    }
    # now calculate the mean, variance, and predicted variance
    mean = mean(as.numeric(all_values))
    st_dev = sd(as.numeric(all_values))
    df2 = data.frame(i, mean, st_dev, times_detected, non_mean_frequency, d)
    
    if (count == 0){
      mean_variance_df = df2
      count = count + 1
    } else {
      mean_variance_df = rbind(mean_variance_df, df2)
    }
  }
}


# loop through and add lines to dataframe so that every dilution gets plotted for every site
dilution_list = unique(mean_variance_df$d)
site_list = unique(mean_variance_df$i)
count = 0


for (d in dilution_list)
{
  for (i in site_list)
  {
    df <- mean_variance_df[mean_variance_df$i == i & mean_variance_df$d == d,]
    if (nrow(df) == 0){
      df2 = data.frame(i, 0, 0, 0, NA, d)
    }
    if (count == 0){
      df_to_append = df2
      count = count + 1
    } else {
      df_to_append = rbind(df_to_append, df2)
    }
  }
}

# rename columsn in df_to_append and append it to data_nonunique
colnames(df_to_append) = c("i","mean","st_dev","times_detected","non_mean_frequency","d")
data_nonunique2 = rbind(mean_variance_df, df_to_append)

# convert sites to character and factor them
data_nonunique2$i <- as.character(as.numeric(data_nonunique2$i)) 
data_nonunique2$i <- factor(data_nonunique2$i, levels=c("123","246","406","408","424","432","436","470","511","515","531","540","596","598","622","635","695","1502","1691"))

data_nonunique2$d <- factor(data_nonunique2$d, levels=c("1e+02","1e+03","1e+04",1e+05, 1e+06, 1e+07))

# plot
ggplot(data_nonunique2, aes(i, mean, fill=d))+geom_bar(aes(fill=d), position="dodge", stat="identity", width=0.7)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=13))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title=element_text(size=13, vjust=0.2))+theme(axis.text=element_text(size=11, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=11))+theme(legend.title=element_text(size=13, face="plain"))+theme(panel.margin=unit(0.9, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="nucleotide site",y="SNP frequency (%)")+scale_colour_manual(values=c(color1,color2,color3,color4,color5,color6))+scale_fill_manual(values=c(color1,color2,color3,color4,color5,color6))+scale_x_discrete(drop=F)+scale_y_continuous(breaks=seq(-30,110,10), limits=c(-30,110), minor_breaks=c(10,30,50,70,90))+geom_errorbar(data=data_nonunique2,aes(x=i, ymin=mean-st_dev, ymax=mean+st_dev),position="dodge", stat="identity",width=0.7, alpha=0.6)


#### FIGURE 4B: STANDARD DEVIATION VS INPUT DNA COPIES #################

ggplot(mean_variance_df, aes(x=d,y=log10(st_dev),color=d, alpha=0.8))+geom_point(size=2)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=13))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title=element_text(size=13, vjust=0.2))+theme(axis.text=element_text(size=11, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=11))+theme(legend.title=element_text(size=13, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="log10(dilution group)",y="log10(standard deviation)")+scale_y_continuous(breaks=seq(-2,2,1), limits=c(-2,2), minor_breaks=c(-2,-1,0,1,2))+scale_colour_manual(values=c(color1,color2,color3,color4,color5,color6))

