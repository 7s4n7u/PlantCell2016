#@author Luca Santuari
#@mail luca.santuari@wur.nl
#@group Plant Developmental Biology
#@university Wageningen University

#Copyright (C) 2016  Luca Santuari

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.


#Source functions file
source("../Rcode/PLT_regulated.genes_analysis_functions.R")

iterations <- 1e3
pval_cutoff <- 0.05

#Maximum distance to TSS
distance.tss <- 4e3

# Load annotation of PLT2 peaks by Homer
annotation.all.peaks <- read.csv("../Rdata/ChIP/PLT2/Homer/IP1_IP2_pooled_1FDR_chrNum_HomerAnnotation_4Kb_fromTSS.csv")

#Load PLT2 targets. Only genes that are interrogated both by the ath1 and the aragene11st arrays
PLT2_targets <- get.PLT2.direct.regulated.genes()
lapply(PLT2_targets, length)

for(regulation in c("activated", "repressed"))
{
  #regulation <- "repressed"
  print(regulation)
  
  regulatory.peaks.idx <- which( substr(annotation.all.peaks$Nearest.PromoterID,1,9)%in%PLT2_targets[[regulation]])
  annotation.regulatory.peaks <- annotation.all.peaks[regulatory.peaks.idx,]
  annotation.peaks.distance <- with(annotation.regulatory.peaks,
                                    table(cut(Distance.to.TSS[which(abs(Distance.to.TSS)<=distance.tss)], 
                                              seq(-distance.tss,distance.tss,by=500)))
  )
  
  output.shuffled.dir <- paste("../Rdata/peaks_to_TSS_analysis/shuffled/",regulation,"/",sep="")
  dir.create(output.shuffled.dir,
             recursive=T, showWarnings = F)
  
  for(peak.width in c(3e2,5e2))
  {
    
    peak.fileName <- paste("../Rdata/ChIP/PLT2/BED/",
                           regulation,"/",peak.width,"bp/",regulation,
                           "_",peak.width,"bp_peakCenter.bed",sep="")
    
      write.peaks.to.geneList(geneList = PLT2_targets[[regulation]],
                              regulation = regulation,
                              peakLength = peak.width,
                              output.directory = paste("../Rdata/ChIP/PLT2/BED/",regulation,"/",sep="")
      )
      
  }
  
  #Run sh shufflePeaks_[activated,repressed,total_peaks].sh from the folder Scripts/bash
  
  anno.iter.df <- data.frame()
  
  for(iter in 1:iterations)
  {
    #iter <- 1
    if(iter%%10==0){print(iter)}
    anno.iter <- read.delim(paste(output.shuffled.dir, "shuffledPeaks_",iter,".txt",sep=""))
    anno.iter.df <- rbind(anno.iter.df,
                          table(cut(anno.iter$Distance.to.TSS[which(abs(anno.iter$Distance.to.TSS)<=distance.tss)], seq(-distance.tss,distance.tss,by=500)))
    )
  }
  interval <- seq(-distance.tss,distance.tss,by=500)
  
  col.names.anno.iter.df <- paste(interval[1:(length(interval)-1)]/1e3, "", sep="")
  col.names.anno.iter.df[-grep("-", col.names.anno.iter.df)] <- paste("+",
                                                                      col.names.anno.iter.df[-grep("-", col.names.anno.iter.df)], sep="")
  col.names.anno.iter.df[col.names.anno.iter.df=="+0"] <- "TSS"
  
  names(anno.iter.df) <- col.names.anno.iter.df
  names(annotation.peaks.distance) <- col.names.anno.iter.df
  
  #calculate empirical pvalue
  pval <- sapply(1:length(annotation.peaks.distance), 
                 function(i){
                   length(which(annotation.peaks.distance[i] <= anno.iter.df[,i]))/iterations
                 } )
  #apply Bonferroni correction
  pval_cutoff.Bonferroni <- pval_cutoff/length(pval)
  
  #melt data
  longdata <- melt(anno.iter.df)
  longdata.anno <- melt(annotation.peaks.distance)
  names(longdata) <- names(longdata.anno) <- c("indices", "value")
  color.bars <- rep("black",length(annotation.peaks.distance))
  color.bars[which(pval<pval_cutoff.Bonferroni)] <- "red"
  
  bp1 <- ggplot(longdata, aes(x=indices, y=value)) + geom_boxplot(colour="black", size=0.5)
  bp1 <- bp1 + geom_bar(data=longdata.anno, aes(x=indices, y=value), colour=color.bars, fill=color.bars, stat="identity", alpha = 0.5, size=1)
  bp1 <- bp1 + labs(x = "position", y="peaks", title="") + 
    theme(plot.title = element_text(lineheight=1, face="bold"))
  bp1 <- bp1 + theme(axis.text=element_text(size=40, family="Helvetica"),
                     axis.title=element_text(size=40, family="Helvetica"), 
                     plot.title=element_text(size=40,face="bold", family="Helvetica"),
                     axis.text.x = element_text(size=40, family="Helvetica"),
                     axis.text.y = element_text(size=40, family="Helvetica")
  )
  
  xlabels <- c("-4","","-3","","-2","","-1",rep("",2),"+1","","+2","","+3","","+4")
  bp1 <- bp1 + scale_x_discrete(labels=xlabels)
  
  #Save plot
  fileName <- paste("../Figures/peaks_to_TSS_analysis/",
                    regulation,"_Peaks2TSS.eps",sep="")
  cairo_ps(filename = fileName,
           family="Helvetica",
           height = 12, width = 10)
  print(bp1)
  dev.off()
  
  fileName <- paste("../Figures/peaks_to_TSS_analysis/",
                    regulation,"_Peaks2TSS.svg",sep="")
  svg(filename = fileName,
      family="Helvetica",
      height = 12, width = 10)
  print(bp1)
  dev.off()
  
}

#########################################################################################################

###total_peaks

regulation <- "total_peaks"
output.shuffled.dir <- paste("../Rdata/peaks_to_TSS_analysis/shuffled/",regulation,"/",sep="")

annotation.peaks.distance <- with(annotation.all.peaks,
                                  table(cut(Distance.to.TSS[which(abs(Distance.to.TSS)<=distance.tss)], 
                                            seq(-distance.tss,distance.tss,by=500)))
)

anno.iter.df <- data.frame()

for(iter in 1:iterations)
{
  #iter <- 1
  if(iter%%10==0){print(iter)}
  anno.iter <- read.delim(paste(output.shuffled.dir, "shuffledPeaks_",iter,".txt",sep=""))
  anno.iter.df <- rbind(anno.iter.df,
                        table(cut(anno.iter$Distance.to.TSS[which(abs(anno.iter$Distance.to.TSS)<=distance.tss)], seq(-distance.tss,distance.tss,by=500)))
  )
}
interval <- seq(-distance.tss,distance.tss,by=500)

col.names.anno.iter.df <- paste(interval[1:(length(interval)-1)]/1e3, "", sep="")
col.names.anno.iter.df[-grep("-", col.names.anno.iter.df)] <- paste("+",
                                                                    col.names.anno.iter.df[-grep("-", col.names.anno.iter.df)], sep="")
col.names.anno.iter.df[col.names.anno.iter.df=="+0"] <- "TSS"

names(anno.iter.df) <- col.names.anno.iter.df
names(annotation.peaks.distance) <- col.names.anno.iter.df

#calculate empirical pvalue
pval <- sapply(1:length(annotation.peaks.distance), 
               function(i){
                 (1+length(which(annotation.peaks.distance[i] <= anno.iter.df[,i])))/(1+iterations)
               } )

#apply Bonferroni correction
pval_cutoff.Bonferroni <- pval_cutoff/length(pval)

#melt data
longdata <- melt(anno.iter.df)
longdata.anno <- melt(annotation.peaks.distance)
names(longdata) <- names(longdata.anno) <- c("indices", "value")
color.bars <- rep("black",length(annotation.peaks.distance))
color.bars[which(pval<pval_cutoff.Bonferroni)] <- "red"

bp1 <- ggplot(longdata, aes(x=indices, y=value/nrow(annotation.all.peaks)*100)) + geom_boxplot(colour="black", size=0.5)
bp1 <- bp1 + geom_bar(data=longdata.anno, aes(x=indices, y=value/nrow(annotation.all.peaks)*100), 
                      colour=color.bars, fill=color.bars, stat="identity", alpha = 0.5, size=1)
bp1 <- bp1 + labs(x = "position", y="peaks (%)", title="") + 
  theme(plot.title = element_text(lineheight=1, face="bold"))
bp1 <- bp1 + theme(axis.text=element_text(size=40, family="Helvetica"),
                   axis.title=element_text(size=40, family="Helvetica"), 
                   plot.title=element_text(size=40, family="Helvetica",face="bold"))
xlabels <- c("-4","","-3","","-2","","-1",rep("",2),"+1","","+2","","+3","","+4")
bp1 <- bp1 + scale_x_discrete(labels=xlabels)

#Save plot
fileName <- paste("../Figures/peaks_to_TSS_analysis/",
                  regulation,"_Peaks2TSS.eps",sep="")
cairo_ps(filename = fileName,
         family="Helvetica",
         height = 12, width = 10)
print(bp1)
dev.off()

fileName <- paste("../Figures/peaks_to_TSS_analysis/",
                  regulation,"_Peaks2TSS.svg",sep="")
svg(filename = fileName,
    family="Helvetica",
    height = 12, width = 10)
print(bp1)
dev.off()

#Proportion of PLT2 peaks in the 500 bp region upstream the TSS
paste(signif(
  nrow(annotation.all.peaks[
    annotation.all.peaks$Distance.to.TSS >= (-5e2) &
      annotation.all.peaks$Distance.to.TSS <= 0,
    ])/nrow(annotation.all.peaks)*100, digits = 2 ),
  "% of peaks in the 500 bp region upstream the TSS",sep="")

#Proportion of PLT2 peaks in the 3 Kb region upstream the TSS
paste(signif(
  nrow(annotation.all.peaks[
    annotation.all.peaks$Distance.to.TSS >= (-3e3) &
      annotation.all.peaks$Distance.to.TSS <= 0,
    ])/nrow(annotation.all.peaks)*100, digits = 2 ),
  "% of peaks in the 3 Kb region upstream the TSS",sep="")

#Write PLT2 peaks for PLT2 regulated (both activated and repressed) genes

regulation <- "activated_and_repressed"
print(regulation)

regulated.genes <- c(PLT2_targets[["activated"]], PLT2_targets[["repressed"]])
regulated.genes <- regulated.genes[order(regulated.genes)]
regulatory.peaks.idx <- which( substr(annotation.all.peaks$Nearest.PromoterID,1,9)%in%regulated.genes)
annotation.regulatory.peaks <- annotation.all.peaks[regulatory.peaks.idx,]
annotation.peaks.distance <- with(annotation.regulatory.peaks,
                                  table(cut(Distance.to.TSS[which(abs(Distance.to.TSS)<=distance.tss)], 
                                            seq(-distance.tss,distance.tss,by=500)))
)

output.shuffled.dir <- paste("../Rdata/peaks_to_TSS_analysis/shuffled/",regulation,"/",sep="")
dir.create(output.shuffled.dir,
           recursive=T, showWarnings = F)

for(peak.width in c(3e2,5e2))
{
  
  peak.fileName <- paste("../Rdata/ChIP/PLT2/BED/",
                         regulation,"/",peak.width,"bp/",regulation,
                         "_",peak.width,"bp_peakCenter.bed",sep="")
  
    write.peaks.to.geneList(geneList = regulated.genes,
                            regulation = regulation,
                            peakLength = peak.width,
                            output.directory = paste("../Rdata/ChIP/PLT2/BED/",regulation,"/",sep="")
    )
  
}
