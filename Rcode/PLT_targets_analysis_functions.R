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


#Load libraries
library(Gviz)
options(ucscChromosomeNames=FALSE)

library(rtracklayer)
library(GenomicRanges)
library(biomaRt)
require(ggplot2)
require(reshape2)

library(BSgenome.Athaliana.TAIR.TAIR9)

#Get PLT2 regulated and bound genes in whole seedlings and QC
get.PLT2.regulated.genes.geneList <- function()
{
  PLT2.regulated.genes <- list()
  
  #Load PLT2 genes regulated in whole seedlings, only genes that are present on the ATH1 array
  load("../Rdata/PLT2_OE/limma_results/PLT2_regulated.genes.ath1.FC1.55_pval2e-2.rda")
  PLT2.regulated.genes[["Whole_seedling_activated"]] <- PLT2_regulated.genes.ath1[["DEXup"]]
  PLT2.regulated.genes[["Whole_seedling_repressed"]] <- PLT2_regulated.genes.ath1[["DEXdown"]]
  
  #Load PLT2 genes regulated in the QC
  load("../Rdata/pPLT2_QC/limma_results/PLT2_regulated.genes_QC.rda")
  PLT2_regulated.genes_QC_pooled <- list()
  PLT2_regulated.genes_QC_pooled[["activated"]] <- unique(c(
    PLT2_regulated.genes_QC[["1h.up"]],
    PLT2_regulated.genes_QC[["4h.up"]]
  ))
  PLT2_regulated.genes_QC_pooled[["repressed"]] <- unique(c(
    PLT2_regulated.genes_QC[["1h.down"]],
    PLT2_regulated.genes_QC[["4h.down"]]
  ))
  save(file="../Rdata/pPLT2_QC/limma_results/PLT2_regulated.genes_QC_pooled.rda", PLT2_regulated.genes_QC_pooled, compress = T)
  PLT2.regulated.genes[["QC_activated"]] <- PLT2_regulated.genes_QC_pooled[["activated"]]
  PLT2.regulated.genes[["QC_repressed"]] <- PLT2_regulated.genes_QC_pooled[["repressed"]]
  
  #Load PLT2 bound genes
  load("../Rdata/ChIP/PLT2/Homer/PLT2.bound.4Kb.rda")
  #Load ATH1 gene IDs
  genome.AGI <- as.vector(read.csv("../Rdata/BRAINARRAY/ath1121501/ath1.AGI.BRAINARRAY.csv")[,2]) 
  PLT2.regulated.genes[["bound"]] <- PLT2.bound.4Kb[which(PLT2.bound.4Kb%in%genome.AGI)]
  
  PLT2.regulated.genes
}

#Get PLT2 direct regulated.genes
get.PLT2.direct.regulated.genes <- function()
{
  PLT2.direct.regulated.genes <- list()
  
  #Load PLT2 genes regulated in whole seedlings, only genes that are present on the aragene11st array
  load("../Rdata/PLT2_OE/limma_results/PLT2_regulated.genes.ath1.FC1.55_pval2e-2.rda")
  
  #Load PLT2 genes regulated in QC
  load("../Rdata/pPLT2_QC/limma_results/PLT2_regulated.genes_QC_pooled.rda")
  #load PLT2 bound genes
  load("../Rdata/ChIP/PLT2/Homer/PLT2.bound.4Kb.rda")
  
  activated_genes <- unique(c(PLT2_regulated.genes.ath1[["DEXup"]],
                              PLT2_regulated.genes_QC_pooled[["activated"]]))
  PLT2.direct.regulated.genes[["activated"]] <- activated_genes[which(activated_genes%in%PLT2.bound.4Kb)]
  
  repressed_genes <- unique(c(PLT2_regulated.genes.ath1[["DEXdown"]],
                              PLT2_regulated.genes_QC_pooled[["repressed"]]))
  PLT2.direct.regulated.genes[["repressed"]] <- repressed_genes[which(repressed_genes%in%PLT2.bound.4Kb)]
  PLT2.direct.regulated.genes
}

#Write fasta file for region and BED files for regions and peak centers
write.Fasta.Bed <- function(dirName, fileName, myRange){
  
  seqs <- getSeq(Athaliana, myRange)
  names(seqs) <- paste(seqnames(myRange), start(myRange), end(myRange), sep="_")
  writeXStringSet(seqs, file=paste(dirName, "/", fileName, ".fasta", sep=""))
  
  myRange.num <- myRange
  seqlevels(myRange.num) <- substr(seqlevels(myRange.num), 4, 4)
  export.bed( myRange.num, paste(dirName, "/", fileName, ".bed", sep="") )
  
  myRangeCenter <- myRange.num
  start(myRangeCenter) <- myRangeCenter$peak
  width(myRangeCenter) <- 1
  export.bed( myRangeCenter, paste(dirName, "/", fileName, "_peakCenter.bed", sep="") )
  
}

#Write PLT2 peaks associated to a list of genes
write.peaks.to.geneList <- function(geneList, regulation, peakLength, output.directory){
  
  genome.len <- read.table("../Rdata/TAIR10/TAIR10.genome")
  seqinfo.TAIR10 <- Seqinfo(seqnames=c("Chr1", "Chr2", "Chr3", "Chr4",
                                       "Chr5", "ChrM", "ChrC"),
                            seqlengths=genome.len[,2],
                            isCircular=c(rep(F,5),rep(T,2)),
                            genome="TAIR9")
  ###Load annotation of PLT2 peaks
  annotation.peaks <- read.csv("../Rdata/ChIP/PLT2/Homer/IP1_IP2_pooled_1FDR_chrNum_HomerAnnotation_4Kb_fromTSS.csv")
  geneID <- substr(annotation.peaks$Nearest.PromoterID,1,9)
  names(annotation.peaks)[2] <- "peakID"
  annotation.peaks$Chr <- paste("Chr", annotation.peaks$Chr,sep="")
  annotation.peaks <- cbind(AGI=geneID, annotation.peaks)
  annotation.peaks <- annotation.peaks[which(geneID%in%geneList),]
  peaks.anno <- with(annotation.peaks,
                     GRanges(Chr,IRanges(Start-(peakLength/2)+1,Start+(peakLength/2), 
                                         names=make.names(paste(AGI,Chr,Start,sep="_"))),
                     score=Peak.Score, peak=Start,
                     seqinfo=seqinfo.TAIR10))
  peaks.anno <- unique(peaks.anno[order(peaks.anno),])
  length(peaks.anno)
  
  outDir <- paste(output.directory, peakLength,"bp",sep="")
  dir.create(outDir, recursive = T, showWarnings = F)
  write.Fasta.Bed(dirName = outDir, 
                  fileName = paste(regulation, "_", peakLength,"bp",sep=""), 
                  myRange = peaks.anno
  )
  
}


#Plot genomic view of ChIP-seq peaks
plot.Genomic.Tracks.from.wig <- function(geneList.IDs, dataTrack, num.rows, num.columns)
{
  #Load gene annotation
  load("../Rdata/TAIR10/TAIR10.GFF.genes.rda")
  
  #Prepare dataTrack object
  dataTrack.fileName <- "../Rdata/Gviz/dataTrack.rda"
  if(!file.exists(dataTrack.fileName))
  {
    # now create data track using new import function:
    wigFile = "../Rdata/Gviz/PLT2_IP1_IP2_pooled.bed"
    gr1 <- import(wigFile, asRangedData=FALSE)
    dataTrack = DataTrack(range=gr1, genome="tair10",
                          stream=TRUE, name="signal",
                          legend=TRUE, col=c("cornflowerblue")
    )
    save(file=dataTrack.fileName, dataTrack, compress=T)
  }else{
    load(dataTrack.fileName)
  }
  
  #use the Mart database
  #listMarts(host="plants.ensembl.org")
  mart = useMart(
    host="plants.ensembl.org",
    dataset="athaliana_eg_gene", biomart="plants_mart")
  
  #Get gene coordinates
  geneList <- TAIR10.GFF.genes[which(names(TAIR10.GFF.genes)%in%geneList.IDs)]
  
  #The plotted region is win*2
  win <- 6e3
  
  mgTrack <- GenomeAxisTrack(scale = 2e3, labelPos = "beside",
                             exponent = 3, col="black", lwd=8)
  displayPars(mgTrack) <- list(fill = "black", cex=2, cex.id=2, fontsize=12)
  
  output.dir <- "../Figures/Gviz/"
  dir.create(output.dir, showWarnings = F)
  
  #Generate plots
  cairo_ps(filename = paste(output.dir, "genomicPlot.eps", sep=""),
           width = 15, height = 15,
           family="Helvetica"
  )
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(num.rows, num.columns)))
  
  for (i in seq_along(geneList)) {
    #i <- 1
    myChr <- substr(as.vector(seqnames(geneList))[i],4,4)
    myStart <- start(geneList)[i] - win
    myEnd <- start(geneList)[i] + win
    
    #Gene model track
    biomTrack = BiomartGeneRegionTrack(genome="TAIR10", biomart=mart, 
                                       chromosome=myChr, start=myStart, end=myEnd,
                                       showId=T, geneSymbols=T, 
                                       rotate.title=TRUE, col.line=NULL, 
                                       col="black", fill="black",
                                       filters=list(biotype="protein_coding"),
                                       collapseTranscripts="longest", 
                                       name="",
                                       fontcolor.group="black", lwd=2, col.symbol="black",
                                       fontsize=14,
                                       utr5="grey", utr3="grey",
                                       protein_coding="grey",
                                       shape="arrow"
    )
    
    displayPars(biomTrack) <- list(fontcolor.title = "black", cex=2, fontcolor="black", cex.group=2)
    
    chromosome(dataTrack) <- myChr
    displayPars(dataTrack) <- list(fontcolor.title = "black", col.axis = "black",
                                   cex=2,
                                   fontsize=18,
                                   ylim=c(2,8)
    )
    
    pushViewport(
      viewport(layout.pos.col = ((i - 1)%%num.columns) + 1, 
               layout.pos.row = (((i) - 1)%/%num.columns) + 1)
    )
    
    plotTracks(list(mgTrack, dataTrack, biomTrack),
               type="hist", col.histogram="black", cex.title=2,
               cex.axis=1.5, title.width=2,
               col="black", fill="grey",
               filters=list(biotype="protein_coding"),
               collapseTranscripts="gene", 
               name="GENES",    
               chromosome = myChr, 
               from = myStart,
               to = myEnd,
               add = TRUE)
    
    popViewport(1)
  }
  dev.off()
  
}

