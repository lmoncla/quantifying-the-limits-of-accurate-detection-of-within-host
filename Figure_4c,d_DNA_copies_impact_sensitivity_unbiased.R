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
# Figure 3: Input copy number determines reproducibility
###############################################################################

#### FIGURE 3A: BAR PLOT OF FRACTION DETECTED 1, 2, OR 3 TIMES ################

# read in the combined dataframe. 
data=read.table("~/Documents/CA04_SNP_detection/data/unbiased_data/combined_unbiased_SNPs.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("reference_position","type","length","reference","allele","linkage","zygosity","count","coverage","frequency","hyper-allelic","fr_balance","average_quality","input","replicate","minimum_coverage_required"))

# remove indels:
data <- data[data$type != "Deletion",]
data <- data[data$type != "Insertion",]

# remove SNPs called for minimum coverage required = 1
data <- data[data$minimum_coverage_required == "100",]
#data <- data[data$minimum_coverage_required == "1",]


# for this, need to convert the input to a character
data$input = as.character(as.numeric(data$input))
data$input <- gsub("10000","1e+04",data$input)
data$input <- gsub("1000","1e+03",data$input)
data$input <- gsub("100","1e+02",data$input)

# remove all dilutions 10^4 and less
data = data[data$input == "1e+07" | data$input == "1e+06" | data$input == "1e+05",]

# merge the reference allele, reference position, and SNP allele into 1 column called "change" that will look like "A 72 T". This will make sure that if there are multiple changes at a particular site, they will all be treated separately. 
data$change = paste(data$reference, data$reference_position, data$allele)

# calculate the number of times each SNP was detected. First, set up a list and a counter. Then loop through each dilution. Subset data to include only the dilution d, then use the table method to output the count of each reference position within the dataframe. Make the output a dataframe of 4 columns, filling in d for input . If it is the first loop, then make the result a new table. Otherwise, append it to the existing table. 

dilution_list = unique(data$input)
count = 0

for (d in dilution_list)
{
  df <- data[data$input == d,]
  table <- as.data.frame(table(df$change))
  table$input = d
  
  number_sites = nrow(table)
  fraction_in_all_3 = (sum(table$Freq == 3)/number_sites)*100
  fraction_in_2 = (sum(table$Freq == 2)/number_sites)*100
  fraction_in_1 = (sum(table$Freq == 1)/number_sites)*100
  fraction = c(fraction_in_all_3, fraction_in_2, fraction_in_1)
  times_detected = c(3,2,1)
  input = c(d,d,d)
  df2 = data.frame(input,fraction,times_detected)
  
  if (count == 0){
    counts_df = table
    percentage_df = df2
    count = count + 1
  } else {
    counts_df = rbind(counts_df,table)
    percentage_df = rbind(df2,percentage_df)
  }
}


# Plot the fraction of SNPs detected 1, 2 or 3 times for each method
percentage_df$times_detected = as.character(as.numeric(percentage_df$times_detected))
percentage_df$input = factor(percentage_df$input, levels=c("1e+02", "1e+03", "1e+04", 1e+05, 1e+06, 1e+07))

ggplot(percentage_df, aes(x=input, y=fraction, color=times_detected, fill=times_detected))+geom_bar(position="stack", stat = "identity")+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+labs(x="input DNA copies",y="fraction of SNPs")+theme(plot.title=element_text(size=16))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title.y=element_text(size=16, vjust=0.5))+theme(axis.title.x=element_text(size=16, vjust=0.5))+theme(axis.text=element_text(size=16, colour="black"))+theme(axis.text.x=element_text(hjust=1, angle=45))+theme(legend.text=element_text(size=16))+theme(legend.title=element_text(size=16, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(legend.key.size=unit(0.7, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+scale_y_continuous(breaks=c(0,20,40,60,80,100), limits=c(0,100))+scale_colour_manual(values=c("grey80","black","grey40"))+scale_fill_manual(values=c("grey80","black","grey40"))


#### FIGURE 3B: HISTOGRAM OF SNPS DETECTED 1,2, OR 3 TIMES ##################

# first, calculate the mean and standard deviation frequency for each SNP. loop through each dilution and each site. Calculate the mean and variance for each site as well as the number of times it was detected. 

dilution_list = unique(data$input)
count = 0

data$frequency = as.numeric(as.character(data$frequency))

for (d in dilution_list)
{
  sites <- data[data$input == d,]
  site_list = unique(sites$change)
  for (i in site_list)
  {
    df <- data[data$change == i & data$input == d,]
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
    # now calculate the mean and variance
    mean = mean(as.numeric(all_values))
    st_dev = sd(as.numeric(all_values))
    replicate1 = all_values[[1]]
    replicate2 = all_values[[2]]
    replicate3 = all_values[[3]]
    df2 = data.frame(i, replicate1, replicate2, replicate3, mean, st_dev, times_detected, non_mean_frequency, d)
    
    if (count == 0){
      mean_variance_df = df2
      count = count + 1
    } else {
      mean_variance_df = rbind(mean_variance_df, df2)
    }
  }
}

#write.table(mean_variance_df, file = "~/Documents/CA04_SNP_detection/data/All_SNPs_170321.csv", sep = "\t")

# plot a histogram of the number of SNPs in each frequency bin detected only once
mean_variance_df$d <- factor(mean_variance_df$d, levels = c("1e+02", "1e+03", "1e+04", 1e+05, 1e+06, 1e+07))
mean_variance_df$times_detected = as.character(as.numeric(mean_variance_df$times_detected))

mean_variance_df2 = mean_variance_df[mean_variance_df$times_detected == 1,]

# make a series of fake data rows that will be white on the histogram
mean_variance_df2 = rbind(mean_variance_df2,"fake"=c(0,0,0,0,1,1e+08))
mean_variance_df2 = rbind(mean_variance_df2,"fake"=c(0,0,0,0,11,1e+08))
mean_variance_df2 = rbind(mean_variance_df2,"fake"=c(0,0,0,0,21,1e+08))
mean_variance_df2 = rbind(mean_variance_df2,"fake"=c(0,0,0,0,31,1e+08))
mean_variance_df2 = rbind(mean_variance_df2,"fake"=c(0,0,0,0,41,1e+08))
mean_variance_df2 = rbind(mean_variance_df2,"fake"=c(0,0,0,0,51,1e+08))
mean_variance_df2 = rbind(mean_variance_df2,"fake"=c(0,0,0,0,61,1e+08))
mean_variance_df2 = rbind(mean_variance_df2,"fake"=c(0,0,0,0,71,1e+08))
mean_variance_df2 = rbind(mean_variance_df2,"fake"=c(0,0,0,0,81,1e+08))
mean_variance_df2 = rbind(mean_variance_df2,"fake"=c(0,0,0,0,91,1e+08))

ggplot(mean_variance_df2, aes(x=non_mean_frequency, color=d, fill=d))+geom_histogram(stat="bin", position="dodge", bins=10, breaks=c(0,10,20,30,40,50,60,70,80,90,100))+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+labs(x="frequency in population",y="number of SNPs")+theme(plot.title=element_text(size=16))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title.y=element_text(size=16, vjust=0.5))+theme(axis.title.x=element_text(size=16, vjust=0.5))+theme(axis.text=element_text(size=16, colour="black"))+theme(axis.text.x=element_text(hjust=1, angle=45))+theme(legend.text=element_text(size=16))+theme(legend.title=element_text(size=16, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(legend.key.size=unit(0.7, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+scale_colour_manual(values=c(color4, color5, color6))+scale_fill_manual(values=c(color4, color5, color6))+scale_y_continuous(breaks=seq(0,20,2), limits=c(0,18))+scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100))
