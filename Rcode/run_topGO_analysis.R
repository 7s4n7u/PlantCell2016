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

# ColorBrewer palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Load functions for GO analysis
source("../Rcode/topGO_analysis_functions.R")
source("../Rcode/PLT_regulated.genes_analysis_functions.R")
source("../Rcode/PLT_targets_analysis_functions.R")

################################################################################################

#Run topGO analysis for PLT regulated.genes

#Load PLT regulated.genes
PLT_regulated.genes <- load.PLT.regulated.genes()

for(regulation in c("activated", "repressed"))
{
  result.topGO <- run.topGO.ath1(PLT_regulated.genes[[regulation]])
  output.dir <- paste("../Rdata/topGO/",regulation,"/",sep="")
  dir.create(output.dir, showWarnings = F)
  write.csv(file=paste(output.dir,"ALLPLT_",regulation,"_topGO_BP_BY_5FDR.csv",sep=""), result.topGO,
            row.names = F)
  write.csv(file=paste(output.dir,"ALLPLT_",regulation,"_topGO_BP_BY_5FDR_REVIGO_input.csv",sep=""), 
            result.topGO[,c("GO.ID", "adj.P.val")],
            row.names = F)
}

#Run REVIGO from web: http://revigo.irb.hr/
#with the content of the file "*REVIGO_input*" for activated and repressed regulated.genes
#Add "%" to the header (first line with "GO.ID	adj.P.val") for REVIGO to accept it as a header

#REVIGO parameters:
#Allowed similarity:
#Small (0.5)
#If provided, the numbers associated to GO categories are...
#p-values
#Select a database with GO term sizes:
#Arabidopsis thaliana
#Select a semantic similarity measure to use:
#SimRel

#Information of the GO annotation used by REVIGO:
#The version of the Gene Ontology used is monthly release of Oct 2014 "go_201410-termdb.obo-xml.gz".
#The UniProt-to-GO mapping file "gene_association.goa_uniprot.gz" is dated 30th Sep 2014, downloaded from the EBI GOA project.
#Citation for REVIGO:
#Supek F, Bošnjak M, Škunca N, Šmuc T.
#"REVIGO summarizes and visualizes long lists of Gene Ontology terms"
#PLoS ONE 2011. doi:10.1371/journal.pone.0021800 

#Save REVIGO output as:
# ALLPLT_(regulation)_topGO_BP_BY_5FDR_REVIGO_output.csv where (regulation) is either
# "activated" or "repressed", accordingly

################################################################################################

#From the output of REVIGO, selected representative GO terms are plotted

regulation <- "activated"
topGO.results.activated <- read.csv(paste("../Rdata/topGO/",regulation,"/","ALLPLT_",regulation,"_topGO_BP_BY_5FDR.csv", sep=""))
topGO.results.activated <- cbind(topGO.results.activated, reg="activated")

#Select representative GO terms
selected.GO.terms <- c("GO:0007049", #cell cycle
                       "GO:0006259", #DNA metabolic process
                       "GO:0006260", #DNA replication
                       "GO:0000003"	 #reproduction
)
topGO.results.activated <- topGO.results.activated[topGO.results.activated$GO.ID%in%selected.GO.terms,]

regulation <- "repressed"
topGO.results.repressed <- read.csv(paste("../Rdata/topGO/",regulation,"/","ALLPLT_",regulation,"_topGO_BP_BY_5FDR.csv", sep=""))
topGO.results.repressed <- cbind(topGO.results.repressed, reg="repressed")

#Select representative GO terms
selected.GO.terms <- c("GO:0071554", #cell wall organization or biogenesis
                       "GO:0009719", #response to endogenous stimulus
                       "GO:0009605", #response to external stimulus
                       "GO:0019748", #secondary metabolic process
                       "GO:0010054"  #trichoblast differentiation
)
topGO.results.repressed <- topGO.results.repressed[topGO.results.repressed$GO.ID%in%selected.GO.terms,]

topGO.results.all <- rbind(topGO.results.activated, topGO.results.repressed)
#keep only the first two significant digits of the adjusted p-values
topGO.results.all$adj.P.val <- signif(topGO.results.all$adj.P.val, digits = 2)
#Calculated fold enrichment
topGO.results.all$Fold_enrichment <- with(topGO.results.all, Significant/Expected)

my.size <- 15
myplot <- ggplot(topGO.results.all, aes(Term, Fold_enrichment)) + 
  geom_bar(stat = "identity", aes(fill = reg)) + 
  geom_text(colour ='black', 
            aes(label = Term), hjust = 1, size=my.size, fontface="bold") +
  coord_flip() + labs(x = "", y = "") +
  scale_x_discrete(breaks = NULL) + theme_bw() +
  theme(
    axis.text.y=element_blank(),
    axis.text.x=element_text(size=my.size+30, hjust=0.5, family="Helvetica"),
    axis.title=element_text(size=my.size,  family="Helvetica"), 
    plot.title=element_text(size=my.size, face="bold", family="Helvetica"),
    legend.title=element_blank(),
    legend.text = element_text(colour="black", size = my.size+30,
                               face = "bold", family="Helvetica")
  ) + labs(x = "", y="", title="")

topGO.results.all$Term <- factor(
  topGO.results.all$Term, 
  levels = c(
    as.vector(with(topGO.results.all[topGO.results.all$reg=="repressed",],
                   Term[order(Fold_enrichment)])),
    as.vector(with(topGO.results.all[topGO.results.all$reg=="activated",],
                   Term[order(Fold_enrichment)]))
  )
)

myplot <- myplot + scale_fill_manual(values=c(
  cbPalette[2],
  cbPalette[3]
))
myplot%+%topGO.results.all

ggsave(file="ALLPLT_GO_bars.eps",
       width=30, height=20,
       path="../Figures/topGO"
)

#Run topGO analysis on PLT2 regulated genes in QC

PLT2_regulated.genes <- get.PLT2.regulated.genes.geneList()
names(PLT2_regulated.genes)

for(regulation in c("QC_activated", "QC_repressed"))
{
  result.topGO <- run.topGO.ath1(PLT2_regulated.genes[[regulation]])
  output.dir <- paste("../Rdata/topGO/",regulation,"/",sep="")
  dir.create(output.dir, showWarnings = F)
  write.csv(file=paste(output.dir,"PLT2_",regulation,"_topGO_BP_BY_5FDR.csv",sep=""), result.topGO,
            row.names = F)
  write.csv(file=paste(output.dir,"PLT2_",regulation,"_topGO_BP_BY_5FDR_REVIGO_input.csv",sep=""), 
            result.topGO[,c("GO.ID", "adj.P.val")],
            row.names = F)
}

#Get auxin biosynthesis genes regulated in the QC

activated.genes.QC <- PLT2_regulated.genes[["QC_activated"]]
#GO:0009851	auxin biosynthetic process
geneList.auxin.biosynthesis <- get.significant.genes.by.GO(geneList = activated.genes.QC,
                                                           GO.term = 'GO:0009851')

my.plot <- getHeatMap.ggplot2(geneList = geneList.auxin.biosynthesis, 
                              datasetID = "rootTissueAtlas",
                              main.title = "Auxin biosynthesis genes\nregulated by PLT2 in the QC")
ggsave(
  plot = my.plot,
  filename="PLT2_activated_QC_auxin_biosynthesis_Tissue_RootAtlas.eps",
  width=16, height=16,
  path="../Figures/PLT_regulated.genes_overlap/"
)

#Plot significant GO terms for genes regulated by PLT2 in the QC

#From the output of REVIGO, selected representative GO terms are plotted

regulation <- "QC_activated"
topGO.results.activated <- read.csv(paste("../Rdata/topGO/",regulation,"/","PLT2_",regulation,"_topGO_BP_BY_5FDR.csv", sep=""))
topGO.results.activated.REVIGO <- read.csv(paste("../Rdata/topGO/",regulation,"/","PLT2_",regulation,"_topGO_BP_BY_5FDR_REVIGO_output.csv", sep=""))
topGO.results.activated.REVIGO <- topGO.results.activated.REVIGO[topGO.results.activated.REVIGO$eliminated==0,]
topGO.results.activated <- topGO.results.activated[topGO.results.activated$GO.ID%in%topGO.results.activated.REVIGO$term_ID,]
topGO.results.activated <- cbind(topGO.results.activated, reg="activated")

regulation <- "QC_repressed"
topGO.results.repressed <- read.csv(paste("../Rdata/topGO/",regulation,"/","PLT2_",regulation,"_topGO_BP_BY_5FDR.csv", sep=""))
topGO.results.repressed.REVIGO <- read.csv(paste("../Rdata/topGO/",regulation,"/","PLT2_",regulation,"_topGO_BP_BY_5FDR_REVIGO_output.csv", sep=""))
topGO.results.repressed.REVIGO <- topGO.results.repressed.REVIGO[topGO.results.repressed.REVIGO$eliminated==0,]
topGO.results.repressed <- topGO.results.repressed[topGO.results.repressed$GO.ID%in%topGO.results.repressed.REVIGO$term_ID,]
topGO.results.repressed <- cbind(topGO.results.repressed, reg="repressed")

topGO.results.all <- rbind(topGO.results.activated, topGO.results.repressed)
#keep only the first two significant digits of the adjusted p-values
topGO.results.all$adj.P.val <- signif(topGO.results.all$adj.P.val, digits = 2)
#Calculated fold enrichment
topGO.results.all$Fold_enrichment <- with(topGO.results.all, Significant/Expected)

term.nchar <- sapply(as.vector(topGO.results.all$Term), nchar)
my.size <- 15
myplot <- ggplot(topGO.results.all, aes(Term, Fold_enrichment)) + 
  geom_bar(stat = "identity", aes(fill = reg)) + 
  geom_text(colour ='black', 
            aes(label = Term), y = 0, hjust=0,
            size=my.size, fontface="bold") +
  coord_flip() + labs(x = "", y = "") +
  scale_x_discrete(breaks = NULL) + theme_bw() +
  theme(
    axis.text.y=element_blank(),
    axis.text.x=element_text(size=my.size+30, hjust=0.5, family="Helvetica"),
    axis.title=element_text(size=my.size,  family="Helvetica"), 
    plot.title=element_text(size=my.size, face="bold", family="Helvetica"),
    legend.title=element_blank(),
    legend.text = element_text(colour="black", size = my.size+30,
                               face = "bold", family="Helvetica")
  ) + labs(x = "", y="", title="")

topGO.results.all$Term <- factor(
  topGO.results.all$Term, 
  levels = c(
    as.vector(with(topGO.results.all[topGO.results.all$reg=="repressed",],
                   Term[order(Fold_enrichment)])),
    as.vector(with(topGO.results.all[topGO.results.all$reg=="activated",],
                   Term[order(Fold_enrichment)]))
  )
)

myplot <- myplot + scale_fill_manual(values=c(
  cbPalette[2],
  cbPalette[3]
))
myplot%+%topGO.results.all

ggsave(file="PLT2_QC_GO_bars.eps",
       width=20, height=20,
       path="../Figures/topGO",
       limitsize=FALSE
)

