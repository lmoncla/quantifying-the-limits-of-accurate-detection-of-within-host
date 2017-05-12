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
# Figure 5: Normalization does not impact SNP detection
###############################################################################
data=read.table("~/Documents/CA04_SNP_detection/subsampling_with_bbnorm/1000x_variant_calls/subsample_comparison.txt", header=T, stringsAsFactors=F, na.string=c(NA,"NA","na"), sep="\t", col.names=c("site","replicate","group","subsampled_freq","frequency"))

# remove indels (which is only site 11):
data <- data[data$site != 11,]

data$group <- as.character(as.numeric(data$group))
data$group <- gsub("10000","1e+04",data$group)
data$group <- gsub("1000","1e+03",data$group)
data$group <- gsub("100","1e+02",data$group)
data$group <- factor(data$group, levels=c("1e+02","1e+03","1e+04","1e+05","1e+06","1e+07"))

# find SNPs that were only detected in subsampled or not
data_2 <- data[data$frequency == 0 | data$subsampled_freq == 0,]
count(data_2, "frequency")
count(data_2, "subsampled_freq")
write.table(data_2,"~/Documents/CA04_SNP_detection/subsampling_with_bbnorm/SNPs_lost_and_gained_during_subsampling.txt", sep="\t")

# calculate regression between group and standard deviation
r2.lm = lm(frequency~subsampled_freq, data=data)
r2.lm$residuals #get residualsa
summary(r2.lm)

# compare subsampled vs. not as a scatter plot (correlation)
ggplot(data, aes(x=frequency,y=subsampled_freq,color=group, alpha=0.8))+geom_point(size=2)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+theme(plot.title=element_text(size=13))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title=element_text(size=13, vjust=0.2))+theme(axis.text=element_text(size=11, colour="black"))+theme(axis.text.x=element_text(hjust=0.5))+theme(legend.text=element_text(size=11))+theme(legend.title=element_text(size=13, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+labs(x="SNP frequency (%)",y="SNP frequency (%) after subsampling")+scale_x_continuous(breaks=seq(0,100,10), limits=c(0,100))+scale_y_continuous(breaks=seq(0,100,10), limits=c(0,100))+scale_colour_manual(values=c(color1,color2,color3,color4,color5,color6))


# calculate the average difference between subsampled vs. not SNP frequencies

data$difference <- abs(data$subsampled_freq - data$frequency)
mean(data$difference)

df_1e2 = data[data$group == "1e+02",]
df_1e3 = data[data$group == "1e+03",]
df_1e4 = data[data$group == "1e+04",]
df_1e5 = data[data$group == 1e5,]
df_1e6 = data[data$group == 1e6,]
df_1e7 = data[data$group == 1e7,]

# perform paired t-tests comparing frequencies detectedd before and after coverage normalization
t.test(data$frequency, data$subsampled_freq, paired=TRUE)
t.test(df_1e7$frequency,df_1e7$subsampled_freq,paired=TRUE)
t.test(df_1e6$frequency,df_1e6$subsampled_freq,paired=TRUE)
t.test(df_1e5$frequency,df_1e5$subsampled_freq,paired=TRUE)
t.test(df_1e4$frequency,df_1e4$subsampled_freq,paired=TRUE)
t.test(df_1e3$frequency,df_1e3$subsampled_freq,paired=TRUE)
t.test(df_1e2$frequency,df_1e2$subsampled_freq,paired=TRUE)
