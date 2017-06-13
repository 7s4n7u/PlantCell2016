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
require(VennDiagram)
require(GeneOverlap)

source("../Rcode/PLT_regulated.genes_analysis_functions.R")
source("../Rcode/PLT_targets_analysis_functions.R")

#Fisher's exact test
runFisher <- function(geneSet1, geneSet2, genome)
{
  
  if(genome=="aragene11st")
  {
    #Load aragene11st gene IDs
    genome.AGI <- as.vector(read.csv("../Rdata/aragene11st/aragene11st.AGI.csv")[,2])
  }else if(genome=="ath1")
  {
    #Load ATH1 gene IDs
    genome.AGI <- as.vector(read.csv("../Rdata/BRAINARRAY/ath1121501/ath1.AGI.BRAINARRAY.csv")[,2]) 
  }
  
  geneSet1.genome <- geneSet1[which(geneSet1%in%genome.AGI)]
  geneSet2.genome <- geneSet2[which(geneSet2%in%genome.AGI)]
  
  shared <- unique(geneSet2.genome[which(geneSet2.genome%in%geneSet1.genome)])
  not.geneSet1.not.geneSet2 <- genome.AGI[-which(genome.AGI%in%unique(c(geneSet1.genome, geneSet2.genome)))]
  counts = (matrix(data = c(length(shared), 
                            length(geneSet2.genome)-length(shared), 
                            length(geneSet1.genome)-length(shared), 
                            length(not.geneSet1.not.geneSet2)), nrow = 2,
                   dimnames= list(
                     Set1=c("set1", "not.set1"),
                     Set2=c("set2", "not.set2")
                   )))
  
  #Sanity check
  if(!(
    (length(geneSet2.genome)-length(shared)+length(not.geneSet1.not.geneSet2)+length(geneSet1.genome))==(length(genome.AGI))    
    &
    (length(geneSet1.genome)-length(shared)+length(not.geneSet1.not.geneSet2))==(length(genome.AGI)-length(geneSet2.genome))
  ))stop("Wrong conditions!")
  
  counts
  fisher.test(counts, alternative="greater")
  
}

#run Fisher's exact test on a list of gene sets
runFisher.geneSetList <- function(geneSetList, genome)
{
  stopifnot(genome == "ath1" | genome == "aragene11st")
  results.Fisher <- data.frame()
  #Avoid to repeat the same test again
  test.vec <- c()
  for(i in 1:length(geneSetList))
  {
    for(j in 1:length(geneSetList))
    {
      if(i!=j & !paste(i,j)%in%test.vec)
      {
        print(paste("Run test for", names(geneSetList)[i],names(geneSetList)[j]))
        my.fisher.result <- runFisher(geneSetList[[i]], 
                                      geneSetList[[j]],
                                      genome)
        results.Fisher <- rbind(results.Fisher,
                                data.frame(
                                  set1=names(geneSetList)[i],
                                  set2=names(geneSetList)[j],
                                  p.value=my.fisher.result$p.value,
                                  odds.ratio=as.vector(my.fisher.result$estimate)
                                )
        )
        
        test.vec <- c(test.vec, paste(i,j), paste(j,i))
      }
    }
  }
  #adjust p-values with BH
  results.Fisher <- cbind(results.Fisher,
                          adj.Pval=p.adjust(results.Fisher$p.value, method = "BH"))
  
  #Keep only two significant digits
  results.Fisher$p.value <- signif(results.Fisher$p.value, digits = 2)
  results.Fisher$adj.Pval <- signif(results.Fisher$adj.Pval, digits = 2)
  results.Fisher$odds.ratio <- signif(results.Fisher$odds.ratio, digits = 2)
  results.Fisher
}


#run Fisher's test to test the overlap between  PLT regulated.genes and
#clusters of temporal profiles from Lewis at al., 2013
plot.GeneOverlap.PLT.regulated.genes.Lewis.clusters <- function()
{
  #Load PLT regulated.genes
  PLT_regulated.genes <- load.PLT.regulated.genes()
  
  #Load ATH1 gene IDs
  genome.AGI <- as.vector(read.csv("../Rdata/BRAINARRAY/ath1121501/ath1.AGI.BRAINARRAY.csv")[,2]) 
  
  load("../Rdata/Lewis2013/temporal_profiles/lewisClusters_BRAINARRAY.rda")
  lewisClusters <- lewisList[paste("LewisCluster_",1:10,sep="")]
  names(lewisClusters) <- paste("cluster_",1:10,sep="")
  lapply(lewisClusters, length)
  
  gom.obj <- newGOM(lewisClusters, PLT_regulated.genes,
                    genome.size=length(genome.AGI))
  fileName <- "../Figures/PLT_auxin_response/auxin_PLTs_Lewis_clusters_Fishers_Test.eps"
  cairo_ps(filename = fileName,
           family="Helvetica")
  drawHeatmap(gom.obj, adj.p=T, cutoff = 5e-2,
              grid.col="Blues",
              note.col="orange")
  dev.off()
}

#Plot overlap with Lewis clusters
plot.lewisCluster.overlap.with.bars <- function(geneListType)
{
  #Load Lewis clusters
  load("../Rdata/Lewis2013/temporal_profiles/lewisClusters_BRAINARRAY.rda")
  lewisClusters <- lewisList[paste("LewisCluster_",1:10,sep="")]
  lapply(lewisClusters, length)
  
  if(geneListType=="PLT")
  {
    
    #Load PLT regulated.genes
    PLT.geneList <- load.PLT.regulated.genes()
    
  }else if(geneListType=="PLT2"){
    
    #Load PLT regulated.genes
    load("../Rdata/PLT_OE/ALLPLT_regulated.genes.rda")
    PLT.geneList$activated <- ALLPLTsetlistUP_DEX[["PLT2"]]
    PLT.geneList$repressed <- ALLPLTsetlistDOWN_DEX[["PLT2"]]
    
  }else if(geneListType=="PLT2_direct"){
    
    PLT.geneList <- get.PLT2.direct.regulated.genes()
    
  }
  
  #Plot overlap between PLT activated genes and Lewis cluster 1 along the root longitudinal axis
  gene.overlap <- PLT.geneList$activated[PLT.geneList$activated%in%lewisClusters$LewisCluster_1]
  if(length(gene.overlap)>5)
  {
  my.plot <- getHeatMap.ggplot2(geneList = gene.overlap, 
                                datasetID = "rootLongAtlas",
                                main.title = "Overlap with Lewis cluster 1")
  ggsave(
    plot = my.plot,
    filename=paste(geneListType, "_activated_genes_overlap_with_Lewis_cluster_1.eps", sep=""),
    width=8, height=10,
    path="../Figures/PLT_auxin_response"
  )
  }
  
  title.reg <- "Gene counts"
  legend.position <- "right"
  
  clust.vec <- paste("LewisCluster_",1:10,sep="")
  clust.names <- paste("cluster_",1:10,sep="")
  
  res <- data.frame()
  count.regulated.genes.activated <- as.vector(sapply(clust.vec, function(x){length(which(lewisClusters[[x]]%in%PLT.geneList$activated))}))
  res <- rbind(res, data.frame(Lewis=clust.names,
                               regulated.genes=as.vector(unlist(count.regulated.genes.activated)),
                               regulation="PLT activated"))
  count.regulated.genes.repressed <- as.vector(sapply(clust.vec, function(x){length(which(lewisClusters[[x]]%in%PLT.geneList$repressed))}))
  res <- rbind(res, data.frame(Lewis=clust.names,
                               regulated.genes=as.vector(unlist(count.regulated.genes.repressed)),
                               regulation="PLT repressed"))
  count.regulated.genes.non.reg <- as.vector(sapply(clust.vec, function(x){length(lewisClusters[[x]])}))-(count.regulated.genes.activated+count.regulated.genes.repressed)
  res <- rbind(res, data.frame(Lewis=clust.names,
                               regulated.genes=as.vector(unlist(count.regulated.genes.non.reg)),
                               regulation="non regulated"))
  
  require(reshape2)
  require(ggplot2)
  longdata <- melt(res)
  
  longdata$Lewis <- factor(longdata$Lewis, levels=levels(longdata$Lewis)[c(1,3:10,2)])
  
  c <- ggplot(longdata, aes(x=Lewis, y=value, fill=regulation) )
  c <- c + labs(x = "", y="", title=title.reg) + 
    theme(plot.title = element_text(lineheight=1, face="bold"))
  c <- c + theme(
    axis.text.y=element_text(size=30, hjust=1),
    axis.text.x=element_text(size=20, hjust=1, angle = 30, face="bold"),
    axis.title=element_text(size=30), 
    plot.title=element_text(size=30,face="bold"),
    legend.title=element_blank(),
    legend.text = element_text(colour="black", size = 20, face = "bold"),
    legend.position=legend.position
  )
  c <- c + geom_bar(stat="identity") #+ coord_flip()
  c <- c + guides(fill = guide_legend(reverse=FALSE))
  c <- c + scale_fill_manual(values=c(
    cbPalette[2],
    cbPalette[3],
    cbPalette[1]
  ))
  c
}


