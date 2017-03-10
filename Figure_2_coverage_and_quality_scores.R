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
data=read.table("~/Documents/CA04_SNP_detection/coverage depth/original_assemblies/combined2.pileup", header=F, sep="\t", col.names=c("dilution","reference","site","ref_base","coverage"))

data$dilution = gsub(".pileup","",data$dilution)
data["subsampled"] = "no"

data_1000x=read.table("~/Documents/CA04_SNP_detection/coverage depth/subsampled sam files/combined2.pileup", header=F, sep="\t", col.names=c("dilution","reference","site","ref_base","coverage"))

data_1000x$dilution = gsub(".1000x.pileup","", data_1000x$dilution)
data_1000x["subsampled"] = "1000x"

df = rbind(data, data_1000x)

ggplot(df, aes(x=site, y=coverage, color=subsampled))+geom_line(size=1.5)+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+labs(title="Q-score distribution", x="nucleotide site",y="Q-score")+theme(plot.title=element_text(size=13))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title=element_text(size=13, vjust=0.2))+theme(axis.text=element_text(size=11, colour="black"))+theme(axis.text.x=element_text(angle=45, hjust=1))+theme(legend.text=element_text(size=11))+theme(legend.title=element_text(size=13, face="plain"))+theme(panel.margin=unit(0.9, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+scale_y_log10(limits=c(100,1000000), breaks=c(100,1000,10000,100000,1000000))+scale_colour_manual(values=c("grey65","black"))+scale_x_continuous(breaks=seq(0,1800,300), limits=c(0,1800))+facet_wrap(~dilution, scales="free", ncol=3)

mean(data$coverage)
sd(data$coverage)

mean(data_1000x$coverage)
sd(data_1000x$coverage)


## to plot just the mean and standard deviation all on 1 plot

site_list = unique(data$site)
count = 0

for (s in site_list)
{
  df <- data[data$site == s,]
  mean_coverage = mean(df$coverage)
  sd_coverage = sd(df$coverage)
  df2 = data.frame(s,mean_coverage,sd_coverage)
  
  if (count == 0){
    coverage_df = df2
    count = count + 1
  } else {
    coverage_df = rbind(coverage_df,df2)
  }
}

coverage_df1 = coverage_df
coverage_df1["subsampled"] = "no"

coverage_df2 = coverage_df
coverage_df2["subsampled"] = "1000x"

coverage_df = rbind(coverage_df1, coverage_df2)
coverage_df["lower"] = coverage_df$mean_coverage - coverage_df$sd_coverage
coverage_df["upper"] = coverage_df$mean_coverage + coverage_df$sd_coverage

ggplot(coverage_df, aes(x=s, y=mean_coverage, color=subsampled))+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+labs(x="nucleotide site",y="coverage depth")+theme(plot.title=element_text(size=13))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title=element_text(size=13, vjust=0.2))+theme(axis.text=element_text(size=11, colour="black"))+theme(axis.text.x=element_text(angle=45, hjust=1))+theme(legend.text=element_text(size=11))+theme(legend.title=element_text(size=13, face="plain"))+theme(panel.margin=unit(0.9, "lines"))+theme(legend.key.size=unit(0.5, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+scale_y_log10(limits=c(100,1000000), breaks=c(100,1000,10000,100000,1000000))+scale_colour_manual(values=c("grey40","black"))+scale_x_continuous(breaks=seq(0,1800,300), limits=c(0,1800))+geom_ribbon(aes(x=s, ymin=lower, ymax=upper), fill="grey80", linetype=0, alpha=0.6)+geom_line(size=1)


####### FIGURE 2B: QSCORES PLOT ##############################################

data = read.table("~/Documents/CA04_SNP_detection/coverage depth and Qscores/quality_calculations/combined.aqhist2.txt", header=T, sep="\t", col.names=c("id","quality","count","fraction","read"))

total_reads = sum(data$count)
data$total_fraction = (data$count/total_reads)*100

# cast the data into long format
data = dcast(data, quality~id, value.var="count", fun.aggregate=sum)

data$total = rowSums(data)
data$corrected = data$total - data$quality
data$corrected_fraction = (data$corrected/total_reads)*100

# make a bar plot
ggplot(data, aes(x=quality, y=corrected_fraction))+geom_bar(position="dodge", stat = "identity")+theme(panel.grid.major=element_line(colour=NA,size=NA))+theme(panel.grid.minor=element_line(colour=NA,size=NA))+labs(x="average q-score",y="fraction of reads")+theme(plot.title=element_text(size=16))+theme(strip.background = element_rect(colour=NA, fill=NA))+theme(axis.line.x=element_line(colour="black"))+theme(axis.line.y=element_line(colour="black"))+theme(strip.text.x=element_text(size=11))+theme(axis.title.y=element_text(size=16, vjust=0.5))+theme(axis.title.x=element_text(size=16, vjust=0.5))+theme(axis.text=element_text(size=16, colour="black"))+theme(axis.text.x=element_text(hjust=1, angle=45))+theme(legend.text=element_text(size=16))+theme(legend.title=element_text(size=16, face="plain"))+theme(panel.margin=unit(1, "lines"))+theme(plot.margin=unit(c(1,1,1,1),"cm"))+theme(legend.key.size=unit(0.7, "cm"))+theme(panel.background=element_rect(fill=NA))+theme(legend.key=element_rect(fill=NA))+scale_x_continuous(limits=c(25,40), breaks=seq(25,40,1))+scale_y_continuous(breaks=seq(0,30,5))

