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


#Set working directory
code.dir <- "~/Documents/ScheresLab/Manuscripts/PLT_Manuscript_code/Rcode/"
setwd(code.dir)

#Load libraries
require(topGO)
#citation("topGO")
require(ggplot2)

gene_ontology.dir <- "../Rdata/gene_ontology/Arabidopsis_thaliana/"

#Load PLT regulated.genes
load.PLT.regulated.genes <- function()
{
  load(file="../Rdata/PLT_OE/ALLPLT_regulated.genes.rda")
  PLT_regulated.genes <- list()
  PLT_regulated.genes[["activated"]] <- unique(as.vector(unlist(ALLPLTsetlistUP_DEX)))
  PLT_regulated.genes[["repressed"]] <- unique(as.vector(unlist(ALLPLTsetlistDOWN_DEX)))
  PLT_regulated.genes
}

#Generate the necessary data files from the GO annotation
create.GO.datafiles <- function()
{
  
  #Annotation file for Arabidopsis thaliana from http://geneontology.org
  #Annotation: http://geneontology.org/gene-associations/gene_association.tair.gz
  #README: http://geneontology.org/gene-associations/readme/tair.README
  
  #Create file gene_association.rda
  gene_association <- read.csv(file=paste(gene_ontology.dir,"gene_association.csv",sep=""), h=F)
  save(file=paste(gene_ontology.dir,"gene_association.rda",sep=""), gene_association, compress=T)
  
  #Create a mapping from locusID to TAIR10 AGI
  #From TAIR: https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_TAIRlocusaccessionID_AGI_mapping.txt
  locusIDtoAGI <- read.table("../Rdata/TAIR10/TAIR10_TAIRlocusaccessionID_AGI_mapping.txt", sep="\t",h=T)
  locusIDtoAGI$locus_tair_object_id <- paste("locus:", locusIDtoAGI$locus_tair_object_id, sep="")
  
  load(paste(gene_ontology.dir,"gene_association.rda",sep=""))
  gene_association_GO <- gene_association[,c(2,5)]
  names(gene_association_GO) <- c("locus_tair_object_id", "GO_term")
  gene_association_GO <- merge(locusIDtoAGI, gene_association_GO, by="locus_tair_object_id")
  dim(gene_association_GO)
  
  #Create a mapping from geneID to GO term
  geneID2GO <- aggregate(gene_association_GO$GO_term, by=list(gene_association_GO$locus_name), as.vector)
  #head(geneID2GO)
  my.geneID2GO <- geneID2GO[,2]
  names(my.geneID2GO) <- geneID2GO[,1]
  geneID2GO <- my.geneID2GO
  save(file=paste(gene_ontology.dir,"geneID2GO.geneontology.06.30.2015.rda",sep=""), geneID2GO, compress=T)
  
  #Create a mapping from GO term to geneID
  GO2geneID <- aggregate(gene_association_GO$locus_name, by=list(gene_association_GO$GO_term), as.vector)
  #head(GO2geneID)
  my.GO2geneID <- GO2geneID[,2]
  names(my.GO2geneID) <- GO2geneID[,1]
  GO2geneID <- my.GO2geneID
  save(file=paste(gene_ontology.dir,"GO2geneID.geneontology.06.30.2015.rda",sep=""), GO2geneID, compress=T)
  
}

############################
# Run enrichment analysis  #
############################

#run topGO on ath1 universe
run.topGO.ath1 <- function(geneList)
{
  geneID2GO.file <- paste(gene_ontology.dir,"geneID2GO.geneontology.06.30.2015.rda",sep="")
  #If GO files are not present, create them
  if(!file.exists(geneID2GO.file)){create.GO.datafiles()}
  load(geneID2GO.file)
  
  # background genes ( gene universe )
  all.genes <- as.vector(read.csv("../Rdata/BRAINARRAY/ath1121501/ath1.AGI.BRAINARRAY.csv")[,2])
  
  relevant.genes <- factor(as.integer(all.genes %in% geneList))
  names(relevant.genes) <- all.genes
  
  GOdata.BP <- new("topGOdata", ontology='BP', allGenes = relevant.genes, 
                   annotationFun = annFUN.gene2GO, gene2GO = geneID2GO)
  
  ############################
  # Run enrichment analysis  #
  ############################
  results <- runTest(GOdata.BP, algorithm = 'classic', statistic = 'fisher')
  #######################
  # Analysis of results #
  #######################
  ## Biological Processes
  # Basic information on input data can be accessed using the geneData function. The number of annotated
  # genes, the number of significant genes (if it is the case), the minimal size of a GO category as well as the
  # number of GO categories which have at least one signficant gene annotated are listed:
  
  # generate a summary of the enrichment analysis
  results.table <- GenTable(GOdata.BP, results, topNodes = length(results@score), numChar=100)
  
  # How many GO terms were tested?
  #dim(results.table)[1]
  # reduce results to GO terms passing Benjamini-Yekutieli multiple hypothesis corrected pval < 0.05, FDR < 5%
  #The pvalue results are in column 6, named "result1"
  names(results.table)[6] <- "p.value"
  
  #set "< 1e-30" values to 1e-30
  #to avoid generating NAs when applying multiple testing correction (Benjamini-Yekutieli method)
  idx <- grep("<", results.table$p.value)
  min.pvalue <- as.numeric(strsplit(results.table$p.value[idx[1]],"< ")[[1]][2])
  results.table$p.value[idx] <- min.pvalue
  
  results.table <- cbind(results.table, adj.P.val = p.adjust(results.table[,"p.value"],method="BY"))
  
  #Return results for a pvalue cutoff of 0.05
  results.table[which(results.table$adj.P.val<0.05),]
  
}

get.significant.genes.by.GO <- function(geneList, GO.term)
{
  #Load mapping from gene ID to GO term
  load("../Rdata/gene_ontology/Arabidopsis_thaliana/geneID2GO.geneontology.06.30.2015.rda")
  # background genes ( gene universe )
  all.genes <- as.vector(read.csv("../Rdata/BRAINARRAY/ath1121501/ath1.AGI.BRAINARRAY.csv")[,2])
  relevant.genes <- factor(as.integer(all.genes %in% geneList))
  names(relevant.genes) <- all.genes
  GOdata.BP <- new("topGOdata", ontology='BP', allGenes = relevant.genes, 
                   annotationFun = annFUN.gene2GO, gene2GO = geneID2GO)
  allGO = genesInTerm(GOdata.BP)
  my.annotation = lapply(allGO,function(x) x[x %in% geneList] ) 
  my.annotation[[GO.term]]
}

################################################################################################
