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

library(ggplot2)
require(reshape2)

data <- read.table("../Rdata/Motif_analysis/Motif_occurrence_in_peaks/motif_distribution.txt",h=F)

#Use the motif center as reference position for plotting
motif.mid <- 16/2
#Region width
win <- 500
#Bin size
step.interval <- 50
interval <- seq(0,win,by=step.interval)
peak.pos <- table(cut(data[which(data[,1]==0),2]+motif.mid, interval, include.lowest=T))

peak.iter.df <- data.frame()

for(iter in 1:1000)
{
  #iter <- 1
  if(iter%%1e2==0){print(iter)}
  peak.iter <- data[which(data[,1]==iter),2]
  peak.iter.df <- rbind(peak.iter.df,
                        table(cut(peak.iter+motif.mid, interval, include.lowest=T))
  )
}

col.names.peak.iter.df <- as.character(seq(-win/2, win/2-50, by=50))
col.names.peak.iter.df[-grep("-", col.names.peak.iter.df)] <- paste("+",
                                                                    col.names.peak.iter.df[-grep("-", col.names.peak.iter.df)], sep="")
col.names.peak.iter.df[col.names.peak.iter.df=="+0"] <- "center"
names(peak.iter.df) <- col.names.peak.iter.df
names(peak.pos) <- col.names.peak.iter.df

pval_cutoff <- 5e-2
iterations <- 1e3

#calculate simulated p-values
peak.enrichment.matrix.pval <- sapply(1:length(peak.iter.df), function(i){
  (length(which(peak.iter.df[,i] >= peak.pos[i])))/(iterations)
} )
#apply multiple testing correction (Bonferroni, FWER)
pval_cutoff.Bonferroni <- pval_cutoff/length(peak.enrichment.matrix.pval)

#Generate the plot

#melt data
longdata <- melt(peak.iter.df)
longdata.pos <- melt(peak.pos)
longdata.pos$indices <- col.names.peak.iter.df
color.bars <- rep("black",length(peak.pos))

#Bonferroni
color.bars[which(peak.enrichment.matrix.pval<pval_cutoff.Bonferroni)] <- "red"

bp1 <- ggplot(longdata, aes(x=variable, y=value)) + geom_boxplot(colour="black", size=0.5)
bp1 <- bp1 + labs(x = "", y="", title="PLT2 motif occurrences in regulatory peaks") + 
  theme(plot.title = element_text(lineheight=1, face="bold"))
bp1 <- bp1 + theme(axis.text=element_text(size=30),
                   axis.title=element_text(size=30), plot.title=element_text(size=30,face="bold"))
bp1 <- bp1 + geom_bar(data=longdata.pos, aes(x=indices, y=value), colour=color.bars, fill=color.bars, stat="identity", alpha = 0.5, size=1)

#Save plot
fileName <- paste("../Figures/Motif_analysis/",
                  "PLT2_motif_occurrence_in_peaks.eps",sep="")
cairo_ps(filename = fileName,
         family="Helvetica",
         height = 12, width = 12)
bp1
dev.off()

svg(filename = paste("../Figures/Motif_analysis/",
                     "PLT2_motif_occurrence_in_peaks.svg", sep=""),
    family="Helvetica",
    height = 12, width = 12)
bp1
dev.off()

#Percentage of motif occurrences located within 50 bp from the peak summit
#sum(peak.pos[c("-50","center")])/sum(peak.pos)*100
