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
# Figure 6: Duplicate read removal has variable effects on the data
###############################################################################

### FIGURE 6A: REMOVING DUPLICATE READS RESULTS IN A SLIGHT SKEW IN SPECTRUM ##

data=read.table("~/Documents/CA04_SNP_detection/duplicate_read_removal_comparisons/combined_bbtools_picard_CLC_original4.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("reference_position","type","length","reference","allele","linkage","zygosity","count","coverage","frequency","hyper-allelic","fr_balance","average_quality","overlapping_annotations","coding_region_change","amino_acid_change","input","replicate","method"))

# remove indels:
data <- data[data$type != "Deletion",]
data <- data[data$type != "Insertion",]

# calculate the number of times a specific SNP is detected in each dilution among replicates (1,2, or 3 times) and see if that differs depending on duplicate removal method. 

method_list = unique(data$method)
dilution_list = unique(data$input)
count = 0

for (i in method_list)
{
  for (d in dilution_list)
  {
    df <- data[data$method == i & data$input == d,]
    table <- as.data.frame(table(df$reference_position))
    table$input = d
    table$method = i
    
    number_sites = nrow(table)
    fraction_in_all_3 = (sum(table$Freq == 3)/number_sites)*100
    fraction_in_2 = (sum(table$Freq == 2)/number_sites)*100
    fraction_in_1 = (sum(table$Freq == 1)/number_sites)*100
    fraction = c(fraction_in_all_3, fraction_in_2, fraction_in_1)
    times_detected = c(3,2,1)
    input = c(d,d,d)
    method = c(i,i,i)
    df2 = data.frame(input,method,fraction,times_detected)
    
    if (count == 0){
      counts_df = table
      percentage_df = df2
      count = count + 1
    } else {
      counts_df = rbind(counts_df,table)
      percentage_df = rbind(df2,percentage_df)
    }
  }
}

# Plot the fraction of SNPs detected 1, 2 or 3 times for each method
percentage_df$times_detected = as.character(as.numeric(percentage_df$times_detected))
percentage_df$input = as.character(as.numeric(percentage_df$input))
percentage_df$method <- gsub("no","not subsampled",percentage_df$method)
percentage_df$method <- gsub("0_bbtools_absorbcont","bbtools",percentage_df$method)
percentage_df$method = factor(percentage_df$method, levels=c("not subsampled","bbtools","original_picard"))
percentage_df$d <- factor(percentage_df$input, levels = c(100, 1000, 10000, 1e+05, 1e+06, 1e+07))

ggplot(percentage_df, aes(x=input, y=fraction, color=times_detected, fill=times_detected))+geom_bar(position="stack", stat = "identity")+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+labs(x="input RNA copies",y="fraction of SNPs")+theme(plot.title=element_text(size=16))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title.y=element_text(size=16, vjust=0.5))+theme(axis.title.x=element_text(size=16, vjust=0.5))+theme(axis.text=element_text(size=16, colour="black"))+theme(axis.text.x=element_text(hjust=1, angle=45))+theme(legend.text=element_text(size=16))+theme(legend.title=element_text(size=16, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(legend.key.size=unit(0.7, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+facet_wrap(~method, scales="free")+scale_y_continuous(breaks=c(0,20,40,60,80,100), limits=c(0,100))+scale_colour_manual(values=c("grey80","grey40","black"))+scale_fill_manual(values=c("grey80","grey40","black"))


#### FIGURE 6B: SNP FREQUENCY SPECTRUM  ##############################

# plot a frequency spectrum 
df <- data[c(1,10,17:19)]
df$input <- factor(df$input, levels = c(1e+02, 1e+03, 1e+04, 1e+05, 1e+06, 1e+07))
df$method <- gsub("no","not subsampled",df$method)
df$method <- gsub("0_bbtools_absorbcont","bbtools",df$method)
df$method <- gsub("original_picard","picard",df$method)

# make a series of fake data rows that will be white on the histogram
count = 0
method_list = unique(df$method)
number_list = list(1,11,21,31,41,51,61,71,81,91)

for (i in method_list)
{
  for (n in number_list)
  {
    fake_df = data.frame(0,n,1e+08,0,i)
    
    if (count == 0){
      fake_df2 = fake_df
      count = count + 1
    } else {
      fake_df2 = rbind(fake_df2,fake_df)
    }
  }
}

colnames(fake_df2) = c("reference_position","frequency","input","replicate","method")
df = rbind(df,fake_df2)
df$method = factor(df$method, levels=c("not subsampled","bbtools","picard"))

ggplot(df, aes(x=frequency, color=input, fill=input))+geom_histogram(stat="bin", position="dodge", bins=10)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+labs(x="frequency in population",y="number of SNPs")+theme(plot.title=element_text(size=16))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title.y=element_text(size=16, vjust=0.5))+theme(axis.title.x=element_text(size=16, vjust=0.5))+theme(axis.text=element_text(size=16, colour="black"))+theme(axis.text.x=element_text(hjust=1, angle=45))+theme(legend.text=element_text(size=16))+theme(legend.title=element_text(size=16, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(legend.key.size=unit(0.7, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+facet_wrap(~method, scales="free")+scale_colour_manual(values=c(color1, color2, color3, color4, color5, color6, color7))+scale_fill_manual(values=c(color1, color2, color3, color4, color5, color6, color7))+scale_x_discrete(breaks=c(0,10,20,30,40,50,60,70,80,90,100), limits=c(0,100))+scale_y_continuous(breaks=seq(0,60,10), limits=c(0,60))



#### FIGURE 6C: CORRELATION PLOTS  ########################################

data=read.table("~/Documents/CA04_SNP_detection/duplicate_read_removal_comparisons/combined_bbtools_picard_CLC_original4.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("reference_position","type","length","reference","allele","linkage","zygosity","count","coverage","frequency","hyper-allelic","fr_balance","average_quality","overlapping_annotations","coding_region_change","amino_acid_change","input","replicate","method"))

# remove indels:
data <- data[data$type != "Deletion",]
data <- data[data$type != "Insertion",]

data$method <- gsub("no","not subsampled",data$method)
data$method <- gsub("0_bbtools_absorbcont","bbtools",data$method)
data$method <- gsub("original_picard","picard",data$method)

data <- data[c(1,10,17:19)]

# convert to long format again
df1 = data[data$method=="not subsampled",]
df2 = data[data$method=="bbtools",]
df3 = data[data$method=="picard",]

df = merge(df1, df2, by=c("reference_position","input","replicate"), all=TRUE)
df = merge(df, df3, by=c("reference_position","input","replicate"), all=TRUE)

df[is.na(df)] <- 0

# calculate regression between not subsampled and bbtools
r2.lm = lm(frequency.x~frequency.y, data=df)
r2.lm$residuals #get residuals
summary(r2.lm)

r2.lm = lm(frequency.x~frequency, data=df)
r2.lm$residuals #get residuals
summary(r2.lm)

df$input <- as.character(as.numeric(df$input))
df$input <- gsub("10000","1e+04",df$input)
df$input <- gsub("1000","1e+03",df$input)
df$input <- gsub("100","1e+02",df$input)
df$input <- factor(df$input, levels=c("1e+02","1e+03","1e+04","1e+05","1e+06","1e+07"))

# compare subsampled vs. not as a scatter plot (correlation)
ggplot(df, aes(x=frequency.y,y=frequency,color=input, alpha=0.8))+geom_point(size=2)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=13))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title=element_text(size=13, vjust=0.2))+theme(axis.text=element_text(size=11, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=11))+theme(legend.title=element_text(size=13, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="SNP frequency (%) after bbtools",y="SNP frequency (%) after picard")+scale_x_continuous(breaks=seq(0,100,10), limits=c(0,100))+scale_y_continuous(breaks=seq(0,100,10), limits=c(0,100))+scale_colour_manual(values=c(color1,color2,color3,color4,color5,color6))

# only SNPs 1-10%
df2 = df[df$frequency <= 10,]
df2 = df[df$frequency.x <= 10,]
df2 = df[df$frequency.y <= 10,]

# compare subsampled vs. not as a scatter plot (correlation)
ggplot(df2, aes(x=frequency.y,y=frequency,color=input, alpha=0.8))+geom_point(size=2)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=13))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title=element_text(size=13, vjust=0.2))+theme(axis.text=element_text(size=11, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=11))+theme(legend.title=element_text(size=13, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="SNP frequency (%) after bbtools",y="SNP frequency (%) after bbtools")+scale_x_continuous(breaks=seq(0,10,1), limits=c(0,10))+scale_y_continuous(breaks=seq(0,10,1), limits=c(0,10))+scale_colour_manual(values=c(color1,color2,color3,color4,color5,color6))
