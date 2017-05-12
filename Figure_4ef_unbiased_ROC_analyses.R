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

# merge the reference allele, reference position, and SNP allele into 1 column called "change" that will look like "A 72 T". This will make sure that if there are multiple changes at a particular site, they will all be treated separately. 
data$change = paste(data$reference, data$reference_position, data$allele)

# remove SNPs called for minimum coverage required = 1 and dups removed
data <- data[data$minimum_coverage_required == "100",]
#data <- data[data$minimum_coverage_required == "1" & data$method == "unbiased" | data$minimum_coverage_required == "100" & data$method == "amplicon",]
data <- data[data$dups_removed == "no",]
# remove all dilutions 10^4 and less
data = data[data$input == "1e+07" | data$input == "1e+06" | data$input == "1e+05" | data$method == "amplicon",]


# I also want to try out only using SNPs that are >10% in frequency here: 
#data = data[data$frequency >= 20,]

# for this, need to convert the input to a character
data$input = as.character(as.numeric(data$input))
data$input <- gsub("10000","1e+04",data$input)
data$input <- gsub("1000","1e+03",data$input)
data$input <- gsub("100","1e+02",data$input)


###############################################################################
# loop through each method, dilution, replicate and site in the nonunique list; make a df, calculate sensitivity and specificity. True positives are defined as those that are detected at least twice (either in 2 replicates or 2 dilutions) and are therefore calculated separately for each method (not subsampled, picard, bbtools)
###############################################################################

methods_list = unique(data$method)
dilution_list = unique(data$input)
replicate_list = unique(data$replicate)
count = 0

# find all nonunique SNPs
nonunique <- subset(data$reference_position, duplicated(data$reference_position), )
nonunique = sort(nonunique)
nonunique = unique(nonunique)

# make true negatives all sites in the genome not contained in nonunique; set true positives as nonunique list
true_negatives = c(1:1701)
true_negatives <- true_negatives[! true_negatives %in% nonunique]
true_positives = nonunique

for (m in methods_list)
{
  for (d in dilution_list)
  {
    for (r in replicate_list)
    {
      df <- data[data$input == d & data$replicate == r & data$method == m,]
      site_list = unique(df$reference_position)
      true_pos = 0
      true_neg = 0
      false_pos = 0
      false_neg = 0
      
      # calculate the number of true positives and false positives
      for (s in site_list)
      {
        if (s %in% true_positives){
          true_pos = true_pos + 1
        }
        if (s %in% true_negatives){
          false_pos = false_pos + 1
        }  
      }
      
      # calculate false negatives
      for (i in true_positives)
      {
        if (!(i %in% site_list)){
          false_neg = false_neg + 1
        }
      }
      
      # calculate true negatives
      for (i in true_negatives)
      {
        if (!(i %in% site_list)){
          true_neg = true_neg + 1
        }
      }
      
      sensitivity = (true_pos/(true_pos + false_neg))
      specificity = (true_neg/(true_neg + false_pos))
      one_minus_specificity = 1 - specificity
      
      # write these out to a dataframe
      df1 = data.frame(m,d,r,true_pos, false_pos, false_neg, true_neg, sensitivity, specificity, one_minus_specificity)
      
      if (count == 0){
        ROC_df = df1
        count = count + 1
      } else {
        ROC_df = rbind(ROC_df, df1)
      }
    }
  }
}

# factor input levels
ROC_df$f = ROC_df$d
ROC_df$f = factor(ROC_df$d, levels = c(1e+02, 1e+03, 1e+04, 1e+05, 1e+06, 1e+07))
ROC_df$m = factor(ROC_df$m, levels=c("amplicon", "unbiased"))

# calculate the mean sensitivity and specificity for each dilution:

count = 0
for (m in methods_list)
{
  for (d in dilution_list)
  {
    df <- ROC_df[ROC_df$d == d & ROC_df$m ==m,]
    mean_sensitivity = mean(df$sensitivity)
    sd_sensitivity = sd(df$sensitivity)
    mean_specificity = mean(df$specificity)
    sd_specificity = sd(df$specificity)
    mean_one_minus_specificity = mean(df$one_minus_specificity)
    sd_one_minus_specificity = sd(df$one_minus_specificity)
    
    df1 = data.frame(m, d, mean_sensitivity,sd_sensitivity,mean_specificity,sd_specificity, mean_one_minus_specificity,sd_one_minus_specificity)
    
    if (count == 0){
      means_df = df1
      count = count + 1
    } else {
      means_df = rbind(means_df,df1)
    }
  }
}

means_df$f = means_df$d
means_df$f = factor(means_df$d, levels = c(1e+02, 1e+03, 1e+04, 1e+05, 1e+06, 1e+07))
means_df$m = factor(means_df$m, levels=c("amplicon", "unbiased"))


# calculate regression between group and standard deviation
df <- ROC_df[ROC_df$sensitivity != 0,]
#df = df[df$m == "amplicon",]
df = df[df$m == "unbiased",]

df$d = as.numeric(as.character(df$d))

r2.lm = lm(log10(d)~sensitivity, data=df)
r2.lm$residuals #get residuals
summary(r2.lm)

ggplot(df, aes(x=log10(d), y=sensitivity, shape=m, linetype=m))+geom_point(size=3)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+labs(x="log10(input DNA copies)",y="sensitivity")+theme(plot.title=element_text(size=16))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=16))+theme(axis.title.y=element_text(size=16, vjust=0.5))+theme(axis.title.x=element_text(size=16, vjust=0.5))+theme(axis.text=element_text(size=16, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=16))+theme(legend.title=element_text(size=16, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(legend.key.size=unit(0.7, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+scale_colour_manual(values=c(color1, color2, color3, color4, color5, color6))+scale_y_continuous(limits=c(0,1))+geom_smooth(method="lm", colour="grey65", se=FALSE)+scale_shape_manual(values=c(19,1))+scale_linetype_manual(values=c(1,5))

