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
# Figure 2: coverage plots and quality score plots
###############################################################################

####### FIGURE 2A: COVERAGE PLOT ##############################################

# read in combined pileup files
data=read.table("~/Documents/CA04_SNP_detection/coverage_depth_and_Qscores/unbiased_sequencing_sam_files/combined2.pileup", header=F, sep="\t", col.names=c("dilution","replicate","reference","site","ref_base","coverage"))

# the sites where there is 0 coverage are just omitted in the pileup file and therefore will give skewed values for the average coverage. So I need to loop through each site, replicate, and dilution and add lines with 0 coverage for the missing ones
## NOTE THIS TAKES FOREVER BUT IT DOES WORK

sites_list = (1:1701)
dilution_list = unique(data$dilution)
replicate_list = unique(data$replicate)
count = 0

for (dilution in dilution_list)
{
  for (replicate in replicate_list)
  {
    for (site in sites_list)
    {
      df <- data[data$site == site & data$dilution == dilution & data$replicate == replicate,]
      
      if (nrow(df) == 0){
        coverage = 0
        row_to_append = data.frame(dilution,replicate,"CA04_HA_GQ117044",site,"N",coverage)
      } else {
        coverage = df$coverage
        row_to_append = data.frame(dilution,replicate,"CA04_HA_GQ117044",site,"N",coverage)
      }
      
      if (count == 0){
        df2 = row_to_append
        count = count + 1
      } else {
        df2 = rbind(df2,row_to_append)
      }
    }
  }
}


mean(df2$coverage)
sd(df2$coverage)

df_1e7 = df2[df2$dilution == 1e7,]
df_1e6 = df2[df2$dilution == 1e6,]
df_1e5 = df2[df2$dilution == 1e5,]
df_1e4 = df2[df2$dilution == 1e4,]
df_1e3 = df2[df2$dilution == 1000,]
df_1e2 = df2[df2$dilution == 100,]

mean(df_1e7$coverage)
sd(df_1e7$coverage)
mean(df_1e6$coverage)
sd(df_1e6$coverage)
mean(df_1e5$coverage)
sd(df_1e5$coverage)
mean(df_1e4$coverage)
sd(df_1e4$coverage)
mean(df_1e3$coverage)
sd(df_1e3$coverage)
mean(df_1e2$coverage)
sd(df_1e2$coverage)


## plot the mean and standard deviation for each dilution on 1 plot (averaging across replicates)

data = df2
site_list = unique(data$site)
dilution_list = unique(data$dilution)
count = 0

for (d in dilution_list)
{
  for (s in site_list)
  {
    df <- data[data$site == s & data$dilution == d,]
    mean_coverage = mean(df$coverage)
    sd_coverage = sd(df$coverage)
    df2 = data.frame(d,s,mean_coverage,sd_coverage)
  
    if (count == 0){
      coverage_df = df2
      count = count + 1
    } else {
      coverage_df = rbind(coverage_df,df2)
    }
  }
}

coverage_df["lower"] = coverage_df$mean_coverage - coverage_df$sd_coverage
coverage_df["upper"] = coverage_df$mean_coverage + coverage_df$sd_coverage
coverage_df[is.na(coverage_df)] <- 0

ggplot(coverage_df, aes(x=s, y=mean_coverage))+
  facet_wrap(~d, scales="free")+
  theme(panel.grid.major=element_line(colour=NA,size=NA))+
  theme(panel.grid.minor=element_line(colour=NA,size=NA))+
  labs(x="nucleotide site",y="coverage depth")+
  theme(plot.title=element_text(size=13))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(axis.line.x=element_line(colour="black"))+
  theme(axis.line.y=element_line(colour="black"))+
  theme(strip.text.x=element_text(size=11))+
  theme(axis.title=element_text(size=13, vjust=0.2))+
  theme(axis.text=element_text(size=11, colour="black"))+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(legend.text=element_text(size=11))+
  theme(legend.title=element_text(size=13, face="plain"))+
  theme(panel.margin=unit(0.9, "lines"))+
  theme(legend.key.size=unit(0.5, "cm"))+
  theme(panel.background=element_rect(fill=NA))+
  theme(legend.key=element_rect(fill=NA))+
  scale_y_log10(limits=c(1,10000), breaks=c(1,10,100,1000,10000))+
  scale_x_continuous(breaks=seq(0,1800,300), limits=c(0,1800))+
  geom_ribbon(aes(x=s, ymin=lower, ymax=upper), fill="grey80", linetype=0, alpha=0.6)+
  geom_line(size=0.7)


####### FIGURE 2B: QSCORES PLOT ##############################################

data = read.table("~/Documents/CA04_SNP_detection/coverage_depth_and_Qscores/quality_calculations/combined.unbiased.aqhist.txt", header=T, sep="\t", col.names=c("id","quality","count","fraction"))

total_reads = sum(data$count)
data$total_fraction = (data$count/total_reads)*100

# cast the data into long format
data = dcast(data, quality~id, value.var="count", fun.aggregate=sum)

data$total = rowSums(data)
data$corrected = data$total - data$quality
data$corrected_fraction = (data$corrected/total_reads)*100

# make a bar plot
ggplot(data, aes(x=quality, y=corrected_fraction))+geom_bar(position="dodge", stat = "identity")+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+labs(x="average q-score",y="fraction of reads")+theme(plot.title=element_text(size=16))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title.y=element_text(size=16, vjust=0.5))+theme(axis.title.x=element_text(size=16, vjust=0.5))+theme(axis.text=element_text(size=16, colour="black"))+theme(axis.text.x=element_text(hjust=1, angle=45))+theme(legend.text=element_text(size=16))+theme(legend.title=element_text(size=16, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(legend.key.size=unit(0.7, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+scale_x_continuous(limits=c(25,40), breaks=seq(25,40,1))+scale_y_continuous(breaks=seq(0,35,5), limits=c(0,35))



###############################################################################
# Figure X: Comparing ampicon to unbiased SNPs detected
###############################################################################

# read in the combined dataframe. 
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

# plot only those SNPs between 1-10% frequency
df2 = df[df$mean.x <= 10 | df$mean.y <= 10,]

# compare amplicon vs. unbiased as a scatter plot
ggplot(df2, aes(x=mean.x,y=mean.y, alpha=0.8, colour=d))+geom_point(size=2)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=13))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title=element_text(size=13, vjust=0.2))+theme(axis.text=element_text(size=11, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=11))+theme(legend.title=element_text(size=13, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="SNP frequency (%) (unbiased)",y="SNP frequency (%) (amplicon)")+scale_x_continuous(breaks=seq(0,10,1), limits=c(0,10))+scale_y_continuous(breaks=seq(0,10,1), limits=c(0,10))+scale_colour_manual(values=c(color4,color5,color6))

# calculate regression between amplicon and unbiased
r2.lm = lm(mean.y~mean.x, data=df2)
r2.lm$residuals #get residuals
summary(r2.lm)


# calculate the number of total SNPs in each dilution, the fraction of SNPs that are detected by both amplicon and unbiased, and the fraction only found in either 1 or the other

df = data

# run only on variants > 10% frequency in the population
#df = df[df$frequency.x >= 10 | df$frequency.y >= 10,]
dilution_list = unique(df$input)
replicate_list = unique(df$replicate)
sites_list = unique(df$reference_position)
total_SNPs = lengthunique(df$reference_position)
count = 0

for (d in dilution_list)
{
  df1 <- df[df$input == d,]
  total_SNPs = nrow(df1)
  
  # define sites_list as the sites that were detected in any replicates in that dilution
  sites_list = unique(df1$i)
  
  # define counts for detected_in_both, unbiasesd_only and amplicon_only
  detected_in_both = 0
  unbiased_only = 0
  amplicon_only = 0
   
  for (s in sites_list)
  {
    df2 <- df1[df1$i == s,]
    
    # determine whether that site was detected in ampicon, unbiased or both
    if (max(df2$frequency.x) > 0 & max(df2$frequency.y) > 0){
      detected_in_both = detected_in_both + 1
    } else if (max(df2$frequency.x) > 0 & max(df2$frequency.y) == 0) {
      unbiased_only = unbiased_only + 1     
    } else if (max(df2$frequency.x) == 0 & max(df2$frequency.y) > 0) {
      amplicon_only = amplicon_only + 1
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


#################### ANOTHER METHOD #######

# run only on variants > 10% frequency in the population
#df = df[df$frequency.x >= 10 | df$frequency.y >= 10,]
dilution_list = unique(df$input)
replicate_list = unique(df$replicate)
method_list = unique(df$method)
sites_list = unique(df$reference_position)
total_SNPs = length(unique(df$reference_position))
count = 0

for (d in dilution_list)
{
  for (r in replicate_list)
  {
    for (m in method_list)
      df1 <- df[df$input == d & df$replicate == r & df$method == m,]
    
      # define counts for detected_in_both, unbiasesd_only and amplicon_only
      SNPs_detected = 0
  
      for (s in sites_list)
      {
        df2 <- df1[df1$reference_position == s,]
    
        # determine whether that site was detected in ampicon, unbiased or both
        if (nrow(df2) > 0){
          SNPs_detected = SNPs_detected + 1 
        } 
      }
  
    percent_detected = (SNPs_detected/total_SNPs) * 100
    df3 = data.frame(d, r, m, SNPs_detected, percent_detected)
  
    if (count == 0){
      output_df = df3
      count = count + 1
    } else {
       output_df = rbind(output_df,df3)
    }
         
  }
}



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

ggplot(mean_variance_df2, aes(x=non_mean_frequency, color=d, fill=d))+geom_histogram(stat="bin", position="dodge", bins=10)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+labs(x="frequency in population",y="number of SNPs")+theme(plot.title=element_text(size=16))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title.y=element_text(size=16, vjust=0.5))+theme(axis.title.x=element_text(size=16, vjust=0.5))+theme(axis.text=element_text(size=16, colour="black"))+theme(axis.text.x=element_text(hjust=1, angle=45))+theme(legend.text=element_text(size=16))+theme(legend.title=element_text(size=16, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(legend.key.size=unit(0.7, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+scale_colour_manual(values=c(color4, color5, color6))+scale_fill_manual(values=c(color4, color5, color6))+scale_y_continuous(breaks=seq(0,12,2), limits=c(0,12))+scale_x_discrete(breaks=c(0,10,20,30,40,50,60,70,80,90,100), limits=c(0,100))



###############################################################################
# Figure 4: SNPs detected more than 1 time means, sd
###############################################################################

#### FIGURE 4A: MEAN AND STANDARD DEVIATION OF SNPS DETECTED >1 TIME ##########

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

# Subset the dataframe to remove nonduplicated sites. 
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
    replicate1 = all_values[[1]]
    replicate2 = all_values[[2]]
    replicate3 = all_values[[3]]
    df2 = data.frame(i, replicate1,replicate2,replicate3,mean, st_dev, times_detected, non_mean_frequency, d)
    
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
      df2 = data.frame(i, 0, 0, 0, 0, 0, 0, NA, d)
    }
    if (count == 0){
      df_to_append = df2
      count = count + 1
    } else {
      df_to_append = rbind(df_to_append, df2)
    }
  }
}

# rename columns in df_to_append and append it to data_nonunique
colnames(df_to_append) = c("i","replicate1", "replicate2", "replicate3","mean","st_dev","times_detected","non_mean_frequency","d")
data_nonunique2 = rbind(mean_variance_df, df_to_append)

# convert sites to character and factor them
data_nonunique2$i <- as.character(as.numeric(data_nonunique2$i)) 
#data_nonunique2$i <- factor(data_nonunique2$i, levels=c("123","246","406","408","424","432","436","470","511","515","531","540","596","598","622","635","695","1502","1691"))

data_nonunique2$d <- factor(data_nonunique2$d, levels=c("1e+02","1e+03","1e+04",1e+05, 1e+06, 1e+07))

# plot
ggplot(data_nonunique2, aes(i, mean, fill=d))+geom_bar(aes(fill=d), position="dodge", stat="identity", width=0.7)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=13))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title=element_text(size=13, vjust=0.2))+theme(axis.text=element_text(size=11, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=11))+theme(legend.title=element_text(size=13, face="plain"))+theme(panel.margin=unit(0.9, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="nucleotide site",y="SNP frequency (%)")+scale_colour_manual(values=c(color1,color2,color3,color4,color5,color6))+scale_fill_manual(values=c(color1,color2,color3,color4,color5,color6))+scale_x_discrete(drop=F)+scale_y_continuous(breaks=seq(-30,110,10), limits=c(-30,110), minor_breaks=c(10,30,50,70,90))+geom_errorbar(data=data_nonunique2,aes(x=i, ymin=mean-st_dev, ymax=mean+st_dev),position="dodge", stat="identity",width=0.7, alpha=0.6)


#### FIGURE 4B: STANDARD DEVIATION VS INPUT DNA COPIES #################

ggplot(mean_variance_df, aes(x=d,y=log10(st_dev),color=d, alpha=0.8))+geom_point(size=2)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=13))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title=element_text(size=13, vjust=0.2))+theme(axis.text=element_text(size=11, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=11))+theme(legend.title=element_text(size=13, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="log10(dilution group)",y="log10(standard deviation)")+scale_y_continuous(breaks=seq(-2,2,1), limits=c(-2,2), minor_breaks=c(-2,-1,0,1,2))+scale_colour_manual(values=c(color1,color2,color3,color4,color5,color6))

