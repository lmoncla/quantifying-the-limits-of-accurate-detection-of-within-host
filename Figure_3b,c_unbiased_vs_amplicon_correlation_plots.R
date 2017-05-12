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
# Figure X: Correlation between the mean frequency of SNPs as detected in unbiased vs. amplicon sequence data
###############################################################################

# read in the combined dataframe 
data=read.table("~/Documents/CA04_SNP_detection/data/unbiased_data/amplicon_vs_unbiased.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("reference_position","type","length","reference","allele","linkage","zygosity","count","coverage","frequency","hyper-allelic","fr_balance","average_quality","input","replicate","minimum_coverage_required","method","dups_removed"))

# remove indels:
data <- data[data$type != "Deletion",]
data <- data[data$type != "Insertion",]

# remove SNPs called for minimum coverage required = 1 and dups removed
data <- data[data$minimum_coverage_required == "100",]
#data <- data[data$minimum_coverage_required == "1" & data$method == "unbiased" | data$minimum_coverage_required == "100" & data$method == "amplicon",]
data <- data[data$dups_removed == "no",]

# for this, need to convert the input to a character
data$input = as.character(as.numeric(data$input))
data$input <- gsub("10000","1e+04",data$input)
data$input <- gsub("1000","1e+03",data$input)
data$input <- gsub("100","1e+02",data$input)

# I would like to answer the following question: if we consider every single SNP detected in this study, including those detected by both methods, what is the mean frequency of that SNP for each dilution in each method? I will then plot the correlation. 

method_list = unique(data$method)
dilution_list = unique(data$input)
site_list = unique(data$reference_position)
count = 0

data$frequency = as.numeric(as.character(data$frequency))

for (d in dilution_list)
{
  for (m in method_list)
  {
    for (i in site_list)
    {
      df <- data[data$reference_position == i & data$input == d & data$method == m,]
      if (nrow(df) == 1){
        all_values = list(df$frequency, 0, 0)
        times_detected = 1
        non_mean_frequency = df$frequency
      } else if (nrow(df) == 2) {
        all_values = list(df$frequency[1], df$frequency[2], 0)
        times_detected=2
        non_mean_frequency = NA
      } else if (nrow(df) == 3) {
        all_values = list(df$frequency[1], df$frequency[2], df$frequency[3])
        times_detected = 3
        non_mean_frequency = NA
      } else if (nrow(df) == 0) {
        all_values = list(0, 0, 0)
        times_detected = 0
        non_mean_frequency = 0
      }
      # now calculate the mean and variance
      mean = mean(as.numeric(all_values))
      st_dev = sd(as.numeric(all_values))
      replicate1 = all_values[[1]]
      replicate2 = all_values[[2]]
      replicate3 = all_values[[3]]
      df2 = data.frame(i, replicate1, replicate2, replicate3, mean, st_dev,times_detected, non_mean_frequency, d, m)
      
      if (count == 0){
        mean_variance_df = df2
        count = count + 1
      } else {
        mean_variance_df = rbind(mean_variance_df, df2)
      }
    }
  }
}


# reformat dataframe so that SNP frequencies can be compared between amplicon and unbiased methods
# select only the columns: reference_position, frequency, input, replicate and method
df <- mean_variance_df[c(1,5:10)]

# separate out unbiased and amplicon into 2 separate dataframes, then merge back together and convert all NAs to 0s
# method.x is unbiased, method.y is ampilcon
df1 = df[df$m =="unbiased",]
df2 = df[df$m =="amplicon",]

df = merge(df1, df2, by=c("i","d"), all=TRUE)
df$mean.x[is.na(df$mean.x)] <- 0
df$mean.y[is.na(df$mean.y)] <- 0

# remove all dilutions 10^4 and less
df = df[df$d == "1e+07" | df$d == "1e+06" | df$d == "1e+05",]

# compare amplicon vs. unbiased as a scatter plot
ggplot(df, aes(x=mean.x,y=mean.y, alpha=0.8, colour=d))+geom_point(size=2)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=13))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title=element_text(size=13, vjust=0.2))+theme(axis.text=element_text(size=11, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=11))+theme(legend.title=element_text(size=13, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="SNP frequency (%) (unbiased)",y="SNP frequency (%) (amplicon)")+scale_x_continuous(breaks=seq(0,100,10), limits=c(0,100))+scale_y_continuous(breaks=seq(0,100,10), limits=c(0,100))+scale_colour_manual(values=c(color4,color5,color6))

# calculate regression between amplicon and unbiased
r2.lm = lm(mean.x~mean.y, data=df)
r2.lm$residuals #get residuals
summary(r2.lm)

# perform a paired t-test to compare the mean frequency of a SNP in unbiased vs. amplicon for SNPs that were detected in both methods
df_nozeroes <- df[df$mean.x != 0 & df$mean.y != 0,]
t.test(df_nozeroes$mean.x, df_nozeroes$mean.y, paired=T)

# plot only those SNPs between 1-10% frequency
df2 = df[df$mean.x <= 10 | df$mean.y <= 10,]

# compare amplicon vs. unbiased as a scatter plot
ggplot(df2, aes(x=mean.x,y=mean.y, alpha=0.8, colour=d))+geom_point(size=2)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=13))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title=element_text(size=13, vjust=0.2))+theme(axis.text=element_text(size=11, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=11))+theme(legend.title=element_text(size=13, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="SNP frequency (%) (unbiased)",y="SNP frequency (%) (amplicon)")+scale_x_continuous(breaks=seq(0,10,1), limits=c(0,10))+scale_y_continuous(breaks=seq(0,10,1), limits=c(0,10))+scale_colour_manual(values=c(color4,color5,color6))

# calculate regression between amplicon and unbiased
r2.lm = lm(mean.y~mean.x, data=df2)
r2.lm$residuals #get residuals
summary(r2.lm)
