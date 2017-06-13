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


#ggplot2 color-blind-friendly colors
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Source functions file
source("../Rcode/topGO_analysis_functions.R")
source("../Rcode/PLT_auxin_response_analysis_functions.R")

#Load PLT regulated.genes
PLT_regulated.genes <- load.PLT.regulated.genes()

#Load auxin response genes
#Load Lewis time series dataset
load("../Rdata/Lewis2013/limma_results/lewisTime.rda")
#Load Bargmann intact root dataset
load("../Rdata/Bargmann2013/limma_results/Bargmann_intactRoot.rda")

#Create list of gene sets
geneSetList <- list(
  PLT_activated = PLT_regulated.genes[["activated"]],
  PLT_repressed = PLT_regulated.genes[["repressed"]],
  Auxin_root_induced = bargmannRoot[["auxinUp"]],
  Auxin_root_repressed = bargmannRoot[["auxinDown"]],
  Auxin_root_TS_induced = lewisTime[["auxinUp"]],
  Auxin_root_TS_repressed = lewisTime[["auxinDown"]]
)

results.Fisher <- runFisher.geneSetList(geneSetList, "ath1")
write.csv(file="../Rdata/PLT_auxin_response/PLT_auxin_response_Fishers_Test_results.csv",
          results.Fisher)

#Plot Venn diagrams
figure.dir <- "../Figures/PLT_auxin_response/"

my.venn <- venn.diagram(geneSetList[c("PLT_activated",
                           "Auxin_root_induced",
                           "Auxin_root_TS_induced")],
             imagetype="tiff",
             filename = NULL,
             category.names = c(
               expression( bold('PLT activated') ),
               expression( bold('Auxin induced') ),
               expression( bold('Auxin induced TS') )
             ),
             lwd=3,
             resolution=600,
             cex = 2.5,
             margin=0.01, 
             cat.cex = 2,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-20, 20, 180),
             cat.dist = c(0.07, 0.07, 0.05)
) 
fileName <- paste(figure.dir,"PLT_activated_auxin_induced","_3Venn.eps",sep="")
cairo_ps(filename = fileName,
         family="Helvetica")
grid.draw(my.venn)
dev.off()

my.venn <- venn.diagram(geneSetList[c("PLT_activated",
                           "Auxin_root_repressed",
                           "Auxin_root_TS_repressed")],
             imagetype="tiff",
             filename = NULL,
             category.names = c(
               expression( bold('PLT activated') ),
               expression( bold('Auxin repressed') ),
               expression( bold('Auxin repressed TS') )
             ),
             lwd=3,
             resolution=600,
             cex = 2.5,
             margin=0.01, 
             cat.cex = 2,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-20, 20, 180),
             cat.dist = c(0.07, 0.07, 0.05)
) 
fileName <- paste(figure.dir,"PLT_activated_auxin_repressed","_3Venn.eps",sep="")
cairo_ps(filename = fileName,
         family="Helvetica")
grid.draw(my.venn)
dev.off()

my.venn <- venn.diagram(geneSetList[c("PLT_repressed",
                           "Auxin_root_repressed",
                           "Auxin_root_TS_repressed")],
             imagetype="tiff",
             filename = NULL,
             category.names = c(
               expression( bold('PLT repressed') ),
               expression( bold('Auxin repressed') ),
               expression( bold('Auxin repressed TS') )
             ),
             lwd=3,
             resolution=600,
             cex = 2.5,
             margin=0.01, 
             cat.cex = 2,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-20, 20, 180),
             cat.dist = c(0.07, 0.07, 0.05)
)
fileName <- paste(figure.dir,"PLT_repressed_auxin_repressed","_3Venn.eps",sep="")
cairo_ps(filename = fileName,
         family="Helvetica")
grid.draw(my.venn)
dev.off()

my.venn <- venn.diagram(geneSetList[c("PLT_repressed",
                           "Auxin_root_induced",
                           "Auxin_root_TS_induced")],
             imagetype="tiff",
             filename = NULL,
             category.names = c(
               expression( bold('PLT repressed') ),
               expression( bold('Auxin induced') ),
               expression( bold('Auxin induced TS') )
             ),
             lwd=3,
             resolution=600,
             cex = 2.5,
             margin=0.01, 
             cat.cex = 2,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-20, 20, 180),
             cat.dist = c(0.07, 0.07, 0.05)
)
fileName <- paste(figure.dir,"PLT_repressed_auxin_induced","_3Venn.eps",sep="")
cairo_ps(filename = fileName,
         family="Helvetica")
grid.draw(my.venn)
dev.off()

#Overlap between PLT regulated.genes and temporal profiles from Lewis et al., 2013
#Draw heatmap with Fisher's exact test results
plot.GeneOverlap.PLT.regulated.genes.Lewis.clusters()

plot.lewisCluster.overlap.with.bars(geneListType="PLT")
ggsave(filename="Overlap_PLTs_Lewis_clusters_Fisher.eps",
       path="../Figures/PLT_auxin_response",
       height = 8, width = 12
)
