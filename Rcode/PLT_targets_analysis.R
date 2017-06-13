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
source("../Rcode/PLT_targets_analysis_functions.R")

PLT2_targets <- get.PLT2.regulated.genes.geneList()
names(PLT2_targets)
#Create list of gene sets
geneSetList <- list(
  Whole_seedling_activated = PLT2_targets[["Whole_seedling_activated"]],
  Whole_seedling_repressed = PLT2_targets[["Whole_seedling_repressed"]],
  QC_activated = PLT2_targets[["QC_activated"]],
  QC_repressed = PLT2_targets[["QC_repressed"]],
  Bound = PLT2_targets[["bound"]]
)

results.Fisher <- runFisher.geneSetList(geneSetList, "ath1")
dir.create("../Rdata/PLT2_direct_targets_analysis", showWarnings = F)
write.csv(file="../Rdata/PLT2_direct_targets_analysis/PLT2_direct_targets_Fishers_Test_results.csv",
          results.Fisher)

#Plot Venn diagrams
figure.dir <- "../Figures/PLT2_direct_targets_analysis/"
dir.create(figure.dir, showWarnings = F)

my.venn <- venn.diagram(geneSetList[c("Whole_seedling_activated",
                                      "QC_activated",
                                      "Bound")],
                        imagetype="tiff",
                        filename = NULL,
                        category.names = c(
                          expression( bold('Seedlings') ),
                          expression( bold('QC') ),
                          expression( bold('Bound') )
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
fileName <- paste(figure.dir,"PLT2_activated_QC_bound","_3Venn.eps",sep="")
cairo_ps(filename = fileName,
         family="Helvetica")
grid.draw(my.venn)
dev.off()

my.venn <- venn.diagram(geneSetList[c("Whole_seedling_repressed",
                                      "QC_repressed",
                                      "Bound")],
                        imagetype="tiff",
                        filename = NULL,
                        category.names = c(
                          expression( bold('Seedlings') ),
                          expression( bold('QC') ),
                          expression( bold('Bound') )
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
fileName <- paste(figure.dir,"PLT2_repressed_QC_bound","_3Venn.eps",sep="")
cairo_ps(filename = fileName,
         family="Helvetica")
grid.draw(my.venn)
dev.off()

my.venn <- venn.diagram(geneSetList[c(
  "Whole_seedling_activated",
  "QC_activated",
  "Whole_seedling_repressed",
  "QC_repressed")],
  imagetype="tiff",
  filename = NULL,
  category.names = c(
    expression( bold('Seedlings activated') ),
    expression( bold('QC activated') ),
    expression( bold('Seedlings repressed') ),
    expression( bold('QC repressed') )
  ),
  lwd=3,
  resolution=600,
  cex = 2.5,
  margin=0.01, 
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0,0,0,0),
  cat.dist = c(0.3,0.3,0.16,0.16)
)
fileName <- paste(figure.dir,"PLT2_regulated_QC","_4Venn.eps",sep="")
cairo_ps(filename = fileName,
         family="Helvetica")
grid.draw(my.venn)
dev.off()

#Load PLT2 genes regulated in the QC
load("../Rdata/pPLT2_QC/limma_results/PLT2_targets_QC.rda")

my.venn <- venn.diagram(PLT2_targets_QC,
  imagetype="tiff",
  filename = NULL,
  category.names = c(
    expression( bold('QC activated 1h') ),
    expression( bold('QC repressed 1h') ),
    expression( bold('QC activated 4h') ),
    expression( bold('QC repressed 4h') )
  ),
  lwd=3,
  resolution=600,
  cex = 2.5,
  margin=0.01, 
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0,0,0,0),
  cat.dist = c(0.3,0.3,0.16,0.16)
)
fileName <- paste(figure.dir,"PLT2_QC","_4Venn.eps",sep="")
cairo_ps(filename = fileName,
         family="Helvetica")
grid.draw(my.venn)
dev.off()

#Generate genomic plot
geneList.IDs <- c(
  "AT1G16070", #TLP8
  "AT3G24810", #ICK3
  "AT3G50870", #HAN
  "AT1G19850" #ARF5/MP
)

plot.Genomic.Tracks.from.wig(geneList.IDs, dataTrack, 
                             num.rows = 2, num.columns = 2)
