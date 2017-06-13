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
require(reshape2)
require(ggplot2)
require(seriation)
require(RColorBrewer)
require(grid)

#Source functions file
source("../Rcode/topGO_analysis_functions.R")

#ggplot2 color-blind-friendly colors
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#return gene expression matrix for gene list
get.gene.exp.dataset <- function(geneList, datasetID)
{
  switch(datasetID,
         rootTissueAtlas={
           load("../Rdata/Arex/RadialData.rda")
           widedata <- rootTissueAtlas[which(row.names(rootTissueAtlas)%in%geneList),]
         },
         rootLongAtlas={
           load("../Rdata/Arex/LongitudinalData.rda")
           widedata <- rootLongAtlas[which(row.names(rootLongAtlas)%in%geneList),]
         },
         stop("No dataset with such name")
  )
  widedata
}

#Plot heatmap for RootAtlas
getHeatMap.ggplot2 <- function(geneList,
                               datasetID,
                               main.title
)
{
  # Define palette
  myPalette <- colorRampPalette(rev(brewer.pal(11, "PuOr")))
  
  test.legend.range <- FALSE
  if(test.legend.range)
  {
    
    #Set legend limits for all heatmaps showing PLT regulated gene expression
    #Load PLT regulated.genes
    PLT_regulated.genes <- load.PLT.regulated.genes()
    PLT_regulated.genes.geneList <- as.vector(unlist(PLT_regulated.genes))
    
    #Load dataset
    #datasetID <- "rootLongAtlas"
    widedata.PLT_regulated.genes <- get.gene.exp.dataset(PLT_regulated.genes.geneList, datasetID)
    myData.PLT_regulated.genes <- as.matrix(widedata.PLT_regulated.genes)
    #get z-scores
    myData.PLT_regulated.genes <- t(apply(myData.PLT_regulated.genes,1,scale))
    colnames(myData.PLT_regulated.genes) <- colnames(widedata.PLT_regulated.genes)
    longData.PLT_regulated.genes <- melt(myData.PLT_regulated.genes)
    names(longData.PLT_regulated.genes) <- c("Var1", "Var2", "value")
    my.legend.min <- min(longData.PLT_regulated.genes$value)
    my.legend.max <- max(longData.PLT_regulated.genes$value)
    print(paste("Z-score range is from",
                signif(my.legend.min,digits = 2), "to",
                signif(my.legend.max,digits = 2)))
    
  }
  #Legend limits
  #if(datasetID=="rootLongAtlas")
  #{
  #  my.legend.min <- -2.7
  #  my.legend.max <- 3.5
  #}else if(datasetID=="rootTissueAtlas"){
  #  my.legend.min <- -2.6
  #  my.legend.max <- 4.2
  #}
  
  my.legend.min <- -2.7
  my.legend.max <- 4.2
  
  #Load dataset
  widedata <- get.gene.exp.dataset(geneList, datasetID)
  myData <- as.matrix(widedata)
  #get z-scores
  myData <- t(apply(myData,1,scale))
  colnames(myData) <- colnames(widedata)
  
  #Rename Columella as Col in Longitudinal RootAtlas dataset
  if(datasetID=="rootLongAtlas"){colnames(myData)[1] <- "Col"}
  
  # For melt() to work seamlessly, myData has to be a matrix.
  longData <- melt(myData)
  names(longData) <- c("Var1", "Var2", "value")
  longData$Var2 <- factor(longData$Var2, colnames(myData))
  
  min(longData$value)
  
  zp1 <- ggplot(longData,
                aes(x = Var1, y = Var2, fill = value))
  zp1 <- zp1 + geom_tile()
  
  breaks <- seq(-2,4,1)
  zp1 <- zp1 + scale_fill_gradientn(colours = myPalette(1000), name="",
                                    breaks = breaks, labels = format(breaks),
                                    limits=c(my.legend.min, my.legend.max))
  
  zp1 <- zp1 + scale_x_discrete(expand = c(0, 0))
  zp1 <- zp1 + scale_y_discrete(expand = c(0, 0))
  zp1 <- zp1 + ggtitle(main.title)
  
  zp1 <- zp1 + theme(
    axis.text.y = element_text(hjust = 1, size=40),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    #Only with genelist on x-axis:
    #axis.text.x = element_text(size=30, angle = 30, hjust = 1, vjust = 1),
    legend.text = element_text(size=30),
    legend.key.size = unit(3, "lines"), 
    plot.title = element_text(lineheight=.8, face="bold", size=40)
  )
  
  zp1 <- zp1 + labs(list(x = "", y = ""))
  #print(zp1)  # Here, the axes have their original order
  
  # "Optimally" reorder both the rows and columns
  optimalSeriation <- seriate(myData, method = "PCA_angle")
  if(datasetID=="rootTissueAtlas" & length(geneList)<5e2)
  {
    #For rootTissueAtlas and smaller list of genes
    #Anti-Robinson seriation by simulated annealing (Brusco et al 2007).
    optimalSeriation <- seriate(dist(myData), method = "ARSA")
  }
  
  longData$Var1 <- factor(longData$Var1, names(unlist(optimalSeriation[[1]][])))
  
  # The same plot, but with axes reordered according to optimal seriation
  zp1 <- zp1 %+% longData
  zp1
}

#Plot heatmap for both Longitudinal (Long) and Tissue RootAtlas, with gene order
#in the plot of the Tissue RootAtlas that is the same of the plot for 
#the Longitudinal RootAtlas
getHeatMap.ggplot2.Long.Tissue <- function(geneList,
                                           main.title,
                                           output.dir
)
{
  # Define palette
  myPalette <- colorRampPalette(rev(brewer.pal(11, "PuOr")))
  
  my.legend.min <- -2.7
  my.legend.max <- 4.2
  
  #Load dataset
  widedata <- get.gene.exp.dataset(geneList, "rootLongAtlas")
  myData <- as.matrix(widedata)
  #get z-scores
  myData <- t(apply(myData,1,scale))
  
  colnames(myData) <- colnames(widedata)
  row.names(myData) <- row.names(widedata)
  
  #Rename Columella as Col in Longitudinal RootAtlas dataset
  colnames(myData)[1] <- "Col"
  
  # For melt() to work seamlessly, myData has to be a matrix.
  longData <- melt(myData)
  names(longData) <- c("Var1", "Var2", "value")
  longData$Var2 <- factor(longData$Var2, colnames(myData))
  
  zp1 <- ggplot(longData,
                aes(x = Var1, y = Var2, fill = value))
  zp1 <- zp1 + geom_tile()
  
  breaks <- seq(-2,4,1)
  
  zp1 <- zp1 + scale_fill_gradientn(colours = myPalette(1000), name="",
                                    breaks = breaks, labels = format(breaks),
                                    limits=c(my.legend.min, my.legend.max)
  )
  
  zp1 <- zp1 + scale_x_discrete(expand = c(0, 0))
  zp1 <- zp1 + scale_y_discrete(expand = c(0, 0))
  zp1 <- zp1 + ggtitle(main.title)
  
  zp1 <- zp1 + theme(
    axis.text.y = element_text(hjust = 1, size=40),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    #Only with genelist on x-axis:
    #axis.text.x = element_text(size=30, angle = 30, hjust = 1, vjust = 1),
    legend.text = element_text(size=30),
    legend.key.size = unit(3, "lines"), 
    plot.title = element_text(lineheight=.8, face="bold", size=40)
  )
  
  zp1 <- zp1 + labs(list(x = "", y = ""))
  #print(zp1)  # Here, the axes have their original order
  
  # "Optimally" reorder both the rows and columns
  optimalSeriation <- seriate(myData, method = "PCA_angle")
  
  longData$Var1 <- factor(longData$Var1, names(unlist(optimalSeriation[[1]][])))
  
  # The same plot, but with axes reordered according to optimal seriation
  zp1 <- zp1 %+% longData
  
  ggsave(
    plot = zp1,
    filename=paste(main.title,"_Longitudinal_RootAtlas.eps", sep=""),
    width=8, height=10,
    path=output.dir
  )
  
  #Load dataset
  widedata.tissue <- get.gene.exp.dataset(geneList, "rootTissueAtlas")
  myData.tissue <- as.matrix(widedata.tissue)
  #get z-scores
  myData.tissue <- t(apply(myData.tissue,1,scale))
  colnames(myData.tissue) <- colnames(widedata.tissue)
  
  # For melt() to work seamlessly, myData has to be a matrix.
  longData.tissue <- melt(myData.tissue)
  names(longData.tissue) <- c("Var1", "Var2", "value")
  longData.tissue$Var2 <- factor(longData.tissue$Var2, colnames(myData.tissue))
  
  #min(longData$value)
  
  zp2 <- ggplot(longData.tissue,
                aes(x = Var1, y = Var2, fill = value))
  zp2 <- zp2 + geom_tile()
  
  breaks <- seq(-2,4,1)
  zp2 <- zp2 + scale_fill_gradientn(colours = myPalette(1000), name="",
                                    breaks = breaks, labels = format(breaks),
                                    limits=c(my.legend.min, my.legend.max))
  
  zp2 <- zp2 + scale_x_discrete(expand = c(0, 0))
  zp2 <- zp2 + scale_y_discrete(expand = c(0, 0))
  zp2 <- zp2 + ggtitle(main.title)
  
  zp2 <- zp2 + theme(
    axis.text.y = element_text(hjust = 1, size=40),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    #Only with genelist on x-axis:
    #axis.text.x = element_text(size=30, angle = 30, hjust = 1, vjust = 1),
    legend.text = element_text(size=30),
    legend.key.size = unit(3, "lines"), 
    plot.title = element_text(lineheight=.8, face="bold", size=40)
  )
  
  zp2 <- zp2 + labs(list(x = "", y = ""))
  #print(zp2)  # Here, the axes have their original order
  
  #Use same order as for rootLongAtlas
  longData.tissue$Var1 <- factor(longData.tissue$Var1, levels=levels(longData$Var1))
  
  # The same plot, but with axes reordered according to optimal seriation
  zp2 <- zp2 %+% longData.tissue
  
  ggsave(
    plot = zp2,
    filename=paste(main.title,"_Tissue_RootAtlas.eps", sep=""),
    width=16, height=16,
    path=output.dir
  )
  
}

#Create matrix with overlaps among PLT regulated.genes
create.overlap.matrix <- function(geneSetList)
{
  colNames <- names(geneSetList)
  geneList <- unique(as.vector(unlist(geneSetList)))
  olap.mat <- matrix(0, ncol=length(geneSetList), nrow=length(geneList))
  colnames(olap.mat) <- colNames
  row.names(olap.mat) <- geneList
  for(i in 1:length(geneSetList))
  {
    olap.mat[,i] <- as.integer(geneList%in%geneSetList[[i]])
  }
  olap.mat
}

#Get proportion of PLT regulated.genes shared by at least two PLT members and
#shared by at least three PLT members
get.overlap.proportions <- function()
{
  
  #Load PLT regulated.genes
  load("../Rdata/PLT_OE/ALLPLT_regulated.genes.rda")
  
  activated.overlap <- create.overlap.matrix(ALLPLTsetlistUP_DEX)
  repressed.overlap <- create.overlap.matrix(ALLPLTsetlistDOWN_DEX)
  
  activated.overlap.table <- table(apply(activated.overlap,1,sum))
  shared.by.table <- data.frame(regulation="activated",
                                type=c("unique",paste("shared by",2:6, "PLTs")),
                                counts=as.vector(activated.overlap.table))
  
  repressed.overlap.table <- table(apply(repressed.overlap,1,sum))
  shared.by.table <- rbind(shared.by.table,
                           data.frame(regulation="repressed",
                                      type=c("unique",paste("shared by",2:6, "PLTs")),
                                      counts=as.vector(repressed.overlap.table))
  )
  write.csv(file="../Rdata/PLT_OE/Shared_by_N_PLTs_table.csv", shared.by.table)
  
  overlap.mat <- rbind(
    activated.overlap,
    repressed.overlap
  )
  
  #Calculate how many PLT regulated.genes are shared by at least two PLT members and
  #how many PLT regulated.genes are shared by at least three PLT members
  olap <- apply(overlap.mat, 1, sum)
  olap.tab <- table(olap)
  #shared by at least two PLT members:
  print(
    paste(round(sum(olap.tab[names(olap.tab)!=1])/sum(olap.tab)*100),"% shared by at least two PLT members, ",
          #shared by at least three PLT members:
          round(sum(olap.tab[names(olap.tab)>2])/sum(olap.tab)*100),"% shared by at least three PLT members",sep="")
  )
  
}

#Plot percentage overlap with bars
plot.overlap.with.bars <- function(title.reg, legend.position)
{
  
  #Load PLT regulated.genes
  load("../Rdata/PLT_OE/ALLPLT_regulated.genes.rda")
  
  if(title.reg=="Activated")
  {
    reg.mat <- create.overlap.matrix(ALLPLTsetlistUP_DEX)
  }else if(title.reg=="Repressed")
  {
    reg.mat <- create.overlap.matrix(ALLPLTsetlistDOWN_DEX)
  }
  
  PLT.vec <- c(paste("PLT", c(1:5,7),sep=""),"PLTs")
  
  percent.mat <- matrix(nrow=length(PLT.vec), ncol=length(c(1:6)) )
  row.names(percent.mat) <- PLT.vec
  colnames(percent.mat) <- c("unique", paste("shared by ",c(2:6), " PLTs", sep="") )
  for(plt.ID in PLT.vec)
  {
    for(i in 1:6)
    {
      if(plt.ID=="PLTs")
      {
        percent.mat[plt.ID,] <- table(apply(reg.mat,1,sum))/nrow(reg.mat)*100
      }else{
        idx <- which(reg.mat[,which(colnames(reg.mat)==plt.ID)]==1)
        percent.mat[plt.ID,] <- table(apply(reg.mat[idx,],1,sum))/length(idx)*100
      }
    }
  }
  
  longdata <- melt(percent.mat)
  
  c <- ggplot(longdata, aes(x=Var1, y=value, fill=factor(Var2)) )
  c <- c + labs(x = "", y="regulated.genes (%)", title=title.reg) + 
    theme(plot.title = element_text(lineheight=1, face="bold"))
  c <- c + theme(
    axis.text.y=element_text(size=17, hjust=0),
    axis.text.x=element_text(size=15, hjust=1, angle = 30, face="bold"),
    axis.title=element_text(size=20), 
    plot.title=element_text(size=20,face="bold"),
    legend.title=element_blank(),
    legend.text = element_text(colour="black", size = 16, face = "bold"),
    legend.position=legend.position
  )
  
  c <- c + geom_bar(stat="identity")
  c <- c + guides(fill = guide_legend(reverse=TRUE))
  c <- c + scale_fill_manual(values=cbPalette)
  
  longdata$Var1 <- factor(longdata$Var1, levels(longdata$Var1)[order(levels(longdata$Var1),decreasing = F)])
  
  c %+% longdata
  
}

plot.regulated.genes.counts.with.bars <- function()
{
  
  #Load ALL PLT regulated.genes
  load(file="../Rdata/PLT_OE/ALLPLT_regulated.genes.rda")
  lapply(ALLPLTsetlistUP_DEX, length)
  lapply(ALLPLTsetlistDOWN_DEX, length)
  
  title.reg <- "PLT regulated.genes"
  legend.position <- "right"
  
  PLT.vec <- c(paste("PLT", c(1:5,7),sep=""),"PLTs")
  
  res <- data.frame()
  count.regulated.genes <- lapply(ALLPLTsetlistUP_DEX, length)
  res <- rbind(res, data.frame(PLT=names(count.regulated.genes),
                               regulated.genes=as.vector(unlist(count.regulated.genes)),
                               regulation="activated"))
  count.regulated.genes <- lapply(ALLPLTsetlistDOWN_DEX, length)
  res <- rbind(res, data.frame(PLT=names(count.regulated.genes),
                               regulated.genes=as.vector(unlist(count.regulated.genes)),
                               regulation="repressed"))
  
  require(reshape2)
  require(ggplot2)
  longdata <- melt(res)
  
  c <- ggplot(longdata, aes(x=PLT, y=value, fill=factor(regulation)) )
  c <- c + labs(x = "", y="", title=title.reg) + 
    theme(plot.title = element_text(lineheight=1, face="bold"))
  c <- c + theme(
    axis.text.y=element_text(size=30, hjust=1),
    axis.text.x=element_text(size=30, hjust=1, angle = 30, face="bold"),
    axis.title=element_text(size=30), 
    plot.title=element_text(size=40,face="bold"),
    legend.title=element_blank(),
    legend.text = element_text(colour="black", size = 25, face = "bold"),
    legend.position=legend.position
  )
  c <- c + geom_bar(stat="identity", position="dodge") #+ coord_flip()
  c <- c + guides(fill = guide_legend(reverse=FALSE))
  
  c <- c + scale_fill_manual(values=c(
    cbPalette[2],
    cbPalette[3]
  ))
  
  longdata$Var1 <- factor(longdata$PLT, levels(longdata$PLT)[order(levels(longdata$PLT),decreasing = T)])
  c %+% longdata
  
}

#Load PLT2 direct regulated.genes
load.PLT2.direct.regulated.genes <- function()
{
  PLT2.direct.regulated.genes <- list()
  #Load PLT2 regulated.genes
  load(file="../Rdata/pPLT2_QC/limma_results/PLT2_regulated.genes.seedlings.QC.rda")
  
  PLT2.bound.4Kb.fileName <- "../Rdata/ChIP/PLT2/Homer/PLT2.bound.4Kb.rda"
  if(!file.exists(PLT2.bound.4Kb.fileName))
  {
    #Get PLT2 bound genes
    PLT2.bound.4Kb <- unique(as.vector(substr(read.csv("../Rdata/ChIP/PLT2/Homer/IP1_IP2_pooled_1FDR_chrNum_HomerAnnotation_4Kb_fromTSS.csv")$Nearest.PromoterID,1,9)))
    save(file="../Rdata/ChIP/PLT2/Homer/PLT2.bound.4Kb.rda", PLT2.bound.4Kb, compress=T)
  }
  #Load PLT2 bound
  load(file="../Rdata/ChIP/PLT2/Homer/PLT2.bound.4Kb.rda")
  
  PLT2.direct.regulated.genes[["activated"]] <- PLT2_regulated.genes.seedlings.QC[["activated"]][
    which(PLT2_regulated.genes.seedlings.QC[["activated"]]%in%PLT2.bound.4Kb)
    ]
  PLT2.direct.regulated.genes[["repressed"]] <- PLT2_regulated.genes.seedlings.QC[["repressed"]][
    which(PLT2_regulated.genes.seedlings.QC[["repressed"]]%in%PLT2.bound.4Kb)
    ]
  
  #Print proportions
  print(paste(
    signif(length(PLT2.direct.regulated.genes[["activated"]])/length(PLT2_regulated.genes.seedlings.QC[["activated"]])*100, digits = 2),
    "% of PLT2 genes activated either in whole seedlings or QC are direct"
  ))
  print(paste(
    signif(length(PLT2.direct.regulated.genes[["repressed"]])/length(PLT2_regulated.genes.seedlings.QC[["repressed"]])*100, digits = 2),
    "% of PLT2 genes repressed either in whole seedlings or QC are direct"
  ))
  
  PLT2.direct.regulated.genes
}

#Load PLT2 direct regulated.genes, seedlings only
load.PLT2.direct.regulated.genes.seedlings <- function()
{
  PLT2.direct.regulated.genes <- list()
  #Load PLT2 regulated.genes
  load(file="../Rdata/PLT2_OE/limma_results/PLT2_regulated.genes.ath1.FC1.55_pval2e-2.rda")
  
  PLT2.bound.4Kb.fileName <- "../Rdata/ChIP/PLT2/Homer/PLT2.bound.4Kb.rda"
  if(!file.exists(PLT2.bound.4Kb.fileName))
  {
    #Get PLT2 bound genes
    PLT2.bound.4Kb <- unique(as.vector(substr(read.csv("../Rdata/ChIP/PLT2/Homer/IP1_IP2_pooled_1FDR_chrNum_HomerAnnotation_4Kb_fromTSS.csv")$Nearest.PromoterID,1,9)))
    save(file="../Rdata/ChIP/PLT2/Homer/PLT2.bound.4Kb.rda", PLT2.bound.4Kb, compress=T)
  }
  #Load PLT2 bound
  load(file="../Rdata/ChIP/PLT2/Homer/PLT2.bound.4Kb.rda")
  
  PLT2.direct.regulated.genes[["activated"]] <- PLT2_regulated.genes.ath1[["DEXup"]][
    which(PLT2_regulated.genes.ath1[["DEXup"]]%in%PLT2.bound.4Kb)
    ]
  PLT2.direct.regulated.genes[["repressed"]] <- PLT2_regulated.genes.ath1[["DEXdown"]][
    which(PLT2_regulated.genes.ath1[["DEXdown"]]%in%PLT2.bound.4Kb)
    ]
  
  #Print proportions
  print(paste(
    signif(length(PLT2.direct.regulated.genes[["activated"]])/length(PLT2_regulated.genes.ath1[["DEXup"]])*100, digits = 2),
    "% of PLT2 genes activated either in whole seedlings are direct"
  ))
  print(paste(
    signif(length(PLT2.direct.regulated.genes[["repressed"]])/length(PLT2_regulated.genes.ath1[["DEXdown"]])*100, digits = 2),
    "% of PLT2 genes repressed either in whole seedlings are direct"
  ))
  
  PLT2.direct.regulated.genes
}

#Load PLT2 direct regulated.genes, QC only
load.PLT2.direct.regulated.genes.QC <- function()
{
  PLT2.direct.regulated.genes <- list()
  #Load PLT2 regulated.genes
  load(file="../Rdata/pPLT2_QC/limma_results/PLT2_regulated.genes_QC_pooled.rda")
  
  PLT2.bound.4Kb.fileName <- "../Rdata/ChIP/PLT2/Homer/PLT2.bound.4Kb.rda"
  if(!file.exists(PLT2.bound.4Kb.fileName))
  {
    #Get PLT2 bound genes
    PLT2.bound.4Kb <- unique(as.vector(substr(read.csv("../Rdata/ChIP/PLT2/Homer/IP1_IP2_pooled_1FDR_chrNum_HomerAnnotation_4Kb_fromTSS.csv")$Nearest.PromoterID,1,9)))
    save(file="../Rdata/ChIP/PLT2/Homer/PLT2.bound.4Kb.rda", PLT2.bound.4Kb, compress=T)
  }
  #Load PLT2 bound
  load(file="../Rdata/ChIP/PLT2/Homer/PLT2.bound.4Kb.rda")
  
  PLT2.direct.regulated.genes[["activated"]] <- PLT2_regulated.genes_QC_pooled[["activated"]][
    which(PLT2_regulated.genes_QC_pooled[["activated"]]%in%PLT2.bound.4Kb)
    ]
  PLT2.direct.regulated.genes[["repressed"]] <- PLT2_regulated.genes_QC_pooled[["repressed"]][
    which(PLT2_regulated.genes_QC_pooled[["repressed"]]%in%PLT2.bound.4Kb)
    ]
  
  #Print proportions
  print(paste(
    signif(length(PLT2.direct.regulated.genes[["activated"]])/length(PLT2_regulated.genes_QC_pooled[["activated"]])*100, digits = 2),
    "% of PLT2 genes activated either in QC are direct"
  ))
  print(paste(
    signif(length(PLT2.direct.regulated.genes[["repressed"]])/length(PLT2_regulated.genes_QC_pooled[["repressed"]])*100, digits = 2),
    "% of PLT2 genes repressed either in QC are direct"
  ))
  
  PLT2.direct.regulated.genes
}
