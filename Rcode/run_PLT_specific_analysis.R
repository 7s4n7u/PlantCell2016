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
source("../Rcode/topGO_analysis_functions.R")
source("../Rcode/PLT_auxin_response_analysis_functions.R")

#Get a vector with the sections with maximum relative expression per gene
get.max.rel.exp <- function(geneList, datasetID)
{
  #Load dataset
  widedata <- get.gene.exp.dataset(geneList, datasetID)
  myData <- as.matrix(widedata)
  #get z-scores
  myData <- t(apply(myData,1,scale))
  colnames(myData) <- colnames(widedata)
  
  #Sections with maximum relative expression per gene
  max.rel.exp <- apply(myData, 1, function(gene.row){colnames(myData)[which(gene.row==max(gene.row))]})
  max.rel.exp
}

#get PLT specific genes divided by PLT member
get.PLT.specific.genes <- function()
{
  #Load PLT regulated.genes
  load("../Rdata/PLT_OE/ALLPLT_regulated.genes.rda")
  
  #List containing PLT specific genes, as genes regulated only by one PLT member and not by
  #any other PLT member.
  PLT_specific <- list()
  
  for(regulation in c("activated", "repressed"))
  {
    PLT_specific[[regulation]] <- list()
    
    if(regulation=="activated")
    {
      gene.overlap <- create.overlap.matrix(ALLPLTsetlistUP_DEX)
    }else if(regulation=="repressed")
    {
      gene.overlap <- create.overlap.matrix(ALLPLTsetlistDOWN_DEX)
    }
    
    gene.count <- data.frame()
    for(protein.name in colnames(gene.overlap))
    {
      PLT_specific[[regulation]][[protein.name]] <- row.names(gene.overlap)[
        which(gene.overlap[,which(colnames(gene.overlap)==protein.name)]==1 & 
                apply(gene.overlap,1,sum)==1)]
    }
  }
  PLT_specific
}


#Load ATH1 gene IDs
genome.AGI <- as.vector(read.csv("../Rdata/BRAINARRAY/ath1121501/ath1.AGI.BRAINARRAY.csv")[,2]) 

PLT_specific <- get.PLT.specific.genes()

for(regulation in c("activated", "repressed"))
{
  
  out.dir <- paste("../Rdata/topGO_specific/", regulation, sep="")
  dir.create(path = out.dir, 
             showWarnings = F, recursive = T)
  
  gene.count <- data.frame()
  for(protein.name in names(PLT_specific[[regulation]]))
  {
    gene.count <- rbind(gene.count, data.frame(protein=protein.name, gene.count=length(PLT_specific[[regulation]][[protein.name]])))
    
    topGO.results <- run.topGO.ath1(PLT_specific[[regulation]][[protein.name]])
    write.csv(file = paste(out.dir, "/", protein.name, "_", regulation,"_topGO_BP_BY_5FDR.csv", sep=""),
              topGO.results)
    
    my.plot <- getHeatMap.ggplot2(geneList = PLT_specific[[regulation]][[protein.name]], 
                                  datasetID = "rootLongAtlas",
                                  main.title = paste(protein.name, "specific genes"))
    ggsave(
      plot = my.plot,
      filename=paste(protein.name,"specific_LongitudinalRootAtlas.eps",sep="_"),
      width=8, height=10,
      path=out.dir
    )
    
    my.plot <- getHeatMap.ggplot2(geneList = PLT_specific[[regulation]][[protein.name]], 
                                  datasetID = "rootTissueAtlas",
                                  main.title = paste(protein.name, "specific genes"))
    ggsave(
      plot = my.plot,
      filename=paste(protein.name,"specific_TissueRootAtlas.eps",sep="_"),
      width=16, height=16,
      path=out.dir
    )
    
  }
  write.csv(file = paste(out.dir, "/Gene_count_", regulation,"_specific.csv", sep=""),
            gene.count)
  
  #pooled set
  topGO.geneList <- unique(unlist(PLT_specific[[regulation]]))
  length(topGO.geneList)
  topGO.results <- run.topGO.ath1(topGO.geneList)
  write.csv(file = paste(out.dir, "/ALLPLT_", regulation,"_topGO_BP_BY_5FDR.csv", sep=""),
            topGO.results)
  
}


#Create a matrix with the overlap among significant GO terms of the set of pooled PLT regulated genes (ALLPLT),
#and the set of PLT specific genes for each PLT member

for(regulation in c("activated", "repressed"))
{
  results <- data.frame()
  go.terms.allPLT <- read.csv(paste("../Rdata/topGO/",regulation,
                                    "/ALLPLT_",regulation,"_topGO_BP_BY_5FDR.csv",sep=""))
  results <- data.frame(PLT="ALLPLT", regulation=regulation,
                        GO.ID=go.terms.allPLT$GO.ID, 
                        Term=go.terms.allPLT$Term)
  
  for(proteinID in paste("PLT",c(1:5,7),sep=""))
  {
    #print(proteinID)
    go.terms.singlePLT.specific <- read.csv(paste("../Rdata/topGO_specific/",regulation,
                                                  "/", proteinID,
                                                  "_",regulation,"_topGO_BP_BY_5FDR.csv",sep=""))
    if(nrow(go.terms.singlePLT.specific)>0)
    {
      results <- rbind(results,
                       data.frame(PLT=proteinID, regulation=regulation,
                                  GO.ID=go.terms.singlePLT.specific$GO.ID, 
                                  Term=go.terms.singlePLT.specific$Term)
      )
    }
  }
  
  results.mat <- data.frame(unique(results[,3:4]))
  results.mat <- cbind(results.mat, matrix(0, nrow=nrow(results.mat), ncol=length(unique(results$PLT))))
  names(results.mat) <- c("GO.ID", "Term", as.vector(unique(results$PLT)) )
  
  for(protein.set in names(results.mat)[-c(1:2)])
  {
    results.mat[,protein.set] <- as.integer(results.mat$GO.ID%in%results$GO.ID[results$PLT==protein.set])
  }
  results.mat <- results.mat[order(results.mat$Term),]
  write.csv(file = paste("../Rdata/topGO_specific/",regulation,
                         "/PLT_",regulation,"_GO_overlap.csv",sep=""), results.mat, row.names = F)
}


#Analysis on PLT activated genes that have maximum expression in the sections from L7 to L12

#Load PLT regulated.genes
PLT_regulated.genes <- load.PLT.regulated.genes()

#Genes in L7-L12
regulation <- "activated"

out.dir <- paste("../Rdata/topGO_zones/", regulation, sep="")
dir.create(path = out.dir, 
           showWarnings = F, recursive = T)

max.rel.exp.activated <- get.max.rel.exp(PLT_regulated.genes[[regulation]], "rootLongAtlas")

#For activated genes
genes.L7_to_L12 <- names(max.rel.exp.activated)[max.rel.exp.activated%in%paste("L",7:12,sep="")]
length(genes.L7_to_L12)

topGO.results <- run.topGO.ath1(genes.L7_to_L12)
write.csv(file = paste(out.dir, "/L7_to_L12_", regulation,"_topGO_BP_BY_5FDR.csv", sep=""),
          topGO.results)

#Genes in columella section
genes.Col <- names(max.rel.exp.activated)[max.rel.exp.activated=="Col"]
length(genes.Col)

topGO.results <- run.topGO.ath1(genes.Col)
write.csv(file = paste(out.dir, "/Col_", regulation,"_topGO_BP_BY_5FDR.csv", sep=""),
          topGO.results)

#For repressed genes
regulation <- "repressed"

out.dir <- paste("../Rdata/topGO_zones/", regulation, sep="")
dir.create(path = out.dir, 
           showWarnings = F, recursive = T)

max.rel.exp.repressed <- get.max.rel.exp(PLT_regulated.genes[[regulation]], "rootLongAtlas")

genes.L1_to_L6 <- names(max.rel.exp.repressed)[max.rel.exp.repressed%in%paste("L",1:6,sep="")]
length(genes.L1_to_L6)

topGO.results <- run.topGO.ath1(genes.L1_to_L6)
write.csv(file = paste(out.dir, "/L1_to_L6_", regulation,"_topGO_BP_BY_5FDR.csv", sep=""),
          topGO.results)

#Genes in columella section
genes.Col <- names(max.rel.exp.repressed)[max.rel.exp.repressed=="Col"]
length(genes.Col)

topGO.results <- run.topGO.ath1(genes.Col)
write.csv(file = paste(out.dir, "/Col_", regulation,"_topGO_BP_BY_5FDR.csv", sep=""),
          topGO.results)
