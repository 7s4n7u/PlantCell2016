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

output.dir <- "../Figures/PLT_regulated.genes_overlap"
dir.create(output.dir)

#Get proportion of shared PLT regulated.genes
get.overlap.proportions()

##Plot proportion of shared regulated.genes with bars
plot1 <- plot.overlap.with.bars(title.reg = "Activated", legend.position="left")
plot2 <- plot.overlap.with.bars(title.reg = "Repressed", legend.position="right")

plot1
ggsave(filename="Overlap_bars_Activated.eps",
       width=6, height=8,
       path=output.dir
)

plot2
ggsave(filename="Overlap_bars_Repressed.eps", 
       width=6, height=8,
       path=output.dir
)

#Plot bar with counts of PLT regulated.genes
plot.regulated.genes.counts.with.bars()
ggsave(filename="PLT_regulated.genes_counts.eps",
       path=output.dir
)

#Plot heatmaps with RootAtlas data

#Load PLT regulated.genes
PLT_regulated.genes <- load.PLT.regulated.genes()

my.plot <- getHeatMap.ggplot2(geneList = PLT_regulated.genes[["activated"]], 
                              datasetID = "rootLongAtlas",
                              main.title = "Activated regulated.genes")
ggsave(
  plot = my.plot,
  filename="PLT_activated_Longitudinal_RootAtlas.eps", 
  width=8, height=16,
  path=output.dir
)

my.plot <- getHeatMap.ggplot2(geneList = PLT_regulated.genes[["activated"]], 
                              datasetID = "rootTissueAtlas",
                              main.title = "Activated regulated.genes")
ggsave(
  plot = my.plot,
  filename="PLT_activated_Tissue_RootAtlas.eps",
  width=12, height=16,
  path=output.dir
)

widedata <- get.gene.exp.dataset(geneList = PLT_regulated.genes[["activated"]], 
                                 datasetID = "rootLongAtlas")

#Number of genes with highest expression in the L2
length(
  row.names(widedata)[
    which(apply(widedata,1,function(x){max(x)==x["L2"]}))
    ]
)
dim(widedata)

geneList.L2 <- row.names(widedata)[
  which(apply(widedata,1,function(x){max(x)==x["L2"]}))
  ]
geneList.L2 <- geneList.L2[order(geneList.L2)]
write.csv(file="../Rdata/Tables/PLT_activated_L2.csv",   
          geneList.L2)

PLT_regulated.genes[["activated"]][-which(PLT_regulated.genes[["activated"]]%in%row.names(widedata))]

#How many PLT activated genes peak in expression in Col, L1-6?
max.section <- apply(widedata,1,function(x){ names(x)[which(x==max(x))] })
length(which(max.section%in%c(paste("L", 1:6, sep=""))))
length(which(max.section%in%c("Col")))

#Plot PLT activated genes with highest expression in the L2
my.plot <- getHeatMap.ggplot2(geneList = row.names(widedata)[
  which(apply(widedata,1,function(x){max(x)==x["L2"]}))
  ], 
  datasetID = "rootLongAtlas",
  main.title = "Activated genes in L2")
ggsave(
  plot = my.plot,
  filename="PLT_L2_Longitudinal_RootAtlas.eps", 
  width=8, height=16,
  path=output.dir
)

my.plot <- getHeatMap.ggplot2(geneList = row.names(widedata)[
  which(apply(widedata,1,function(x){max(x)==x["L2"]}))
  ], 
  datasetID = "rootTissueAtlas",
  main.title = "Activated genes in L2")
ggsave(
  plot = my.plot,
  filename="PLT_L2_Tissue_RootAtlas.eps", 
  width=12, height=16,
  path=output.dir
)

my.plot <- getHeatMap.ggplot2(geneList = PLT_regulated.genes[["repressed"]], 
                              datasetID = "rootLongAtlas",
                              main.title = "Repressed regulated.genes")
ggsave(
  plot = my.plot,
  filename="PLT_repressed_Longitudinal_RootAtlas.eps", 
  width=8, height=16,
  path=output.dir
)

my.plot <- getHeatMap.ggplot2(geneList = PLT_regulated.genes[["repressed"]], 
                              datasetID = "rootTissueAtlas",
                              main.title = "Repressed regulated.genes")
ggsave(
  plot = my.plot,
  filename="PLT_repressed_Tissue_RootAtlas.eps",
  width=12, height=16,
  path=output.dir
)

widedata <- get.gene.exp.dataset(geneList = PLT_regulated.genes[["repressed"]], 
                                 datasetID = "rootLongAtlas")
dim(widedata)

PLT_regulated.genes[["repressed"]][-which(PLT_regulated.genes[["repressed"]]%in%row.names(widedata))]

#How many PLT repressed genes peak in expression in sections L7-12?
max.section <- apply(widedata,1,function(x){ names(x)[which(x==max(x))] })
length(which(max.section%in%c(paste("L", 7:12, sep=""))))
#How many PLT repressed genes peak in expression in Col section?
length(which(max.section%in%c("Col")))

getHeatMap.ggplot2.Long.Tissue(geneList = PLT_regulated.genes[["activated"]], 
                               main.title = "Activated.regulated.genes",
                               output.dir = output.dir)

getHeatMap.ggplot2.Long.Tissue(geneList = PLT_regulated.genes[["repressed"]], 
                               main.title = "Repressed.regulated.genes",
                               output.dir = output.dir)


#Load PLT2 direct regulated.genes
PLT2.direct.regulated.genes <- load.PLT2.direct.regulated.genes()

widedata <- get.gene.exp.dataset(geneList = PLT2.direct.regulated.genes[["activated"]], 
                                 datasetID = "rootLongAtlas")

#Number of genes with expression higher in the L2
length(
  row.names(widedata)[
    which(apply(widedata,1,function(x){max(x)==x["L2"]}))
    ]
)
dim(widedata)

my.plot <- getHeatMap.ggplot2(geneList = PLT2.direct.regulated.genes[["activated"]], 
                              datasetID = "rootLongAtlas",
                              main.title = "Directly activated regulated.genes")
ggsave(
  plot = my.plot,
  filename="PLT2_directly_activated_Longitudinal_RootAtlas.eps",
  width=8, height=10,
  path=output.dir
)

my.plot <- getHeatMap.ggplot2(geneList = PLT2.direct.regulated.genes[["activated"]], 
                              datasetID = "rootTissueAtlas",
                              main.title = "Directly activated regulated.genes")
ggsave(
  plot = my.plot,
  filename="PLT2_directly_activated_Tissue_RootAtlas.eps",
  width=16, height=16,
  path=output.dir
)

my.plot <- getHeatMap.ggplot2(geneList = PLT2.direct.regulated.genes[["repressed"]], 
                              datasetID = "rootLongAtlas",
                              main.title = "Directly repressed regulated.genes")
ggsave(
  plot = my.plot,
  filename="PLT2_directly_repressed_Longitudinal_RootAtlas.eps",
  width=8, height=10,
  path=output.dir
)

my.plot <- getHeatMap.ggplot2(geneList = PLT2.direct.regulated.genes[["repressed"]], 
                              datasetID = "rootTissueAtlas",
                              main.title = "Directly repressed regulated.genes")
ggsave(
  plot = my.plot,
  filename="PLT2_directly_repressed_Tissue_RootAtlas.eps",
  width=16, height=16,
  path=output.dir
)

#Load PLT2 direct regulated.genes, seedlings only
PLT2.direct.regulated.genes.seedlings <- load.PLT2.direct.regulated.genes.seedlings()

my.plot <- getHeatMap.ggplot2(geneList = PLT2.direct.regulated.genes.seedlings[["repressed"]], 
                              datasetID = "rootLongAtlas",
                              main.title = "Directly repressed regulated.genes, seedlings")
ggsave(
  plot = my.plot,
  filename="PLT2_directly_repressed_seedlings_Longitudinal_RootAtlas.eps",
  width=8, height=10,
  path=output.dir
)

#Load PLT2 direct regulated.genes, QC only
PLT2.direct.regulated.genes.QC <- load.PLT2.direct.regulated.genes.QC()

my.plot <- getHeatMap.ggplot2(geneList = PLT2.direct.regulated.genes.QC[["repressed"]], 
                              datasetID = "rootLongAtlas",
                              main.title = "Directly repressed regulated.genes, QC")
ggsave(
  plot = my.plot,
  filename="PLT2_directly_repressed_QC_Longitudinal_RootAtlas.eps",
  width=8, height=10,
  path=output.dir
)

#Plot PLTs expression

PLT.geneList <- c(
  PLT1="AT3G20840",
  PLT2="AT1G51190",
  PLT3="AT5G10510",
  PLT4="AT5G17430",
  PLT5="AT5G57390",
  PLT7="AT5G65510"
)

getHeatMap.ggplot2.Long.Tissue(geneList = PLT.geneList, 
                               main.title = "PLTs",
                               output.dir = output.dir)

getHeatMap.ggplot2.Long.Tissue(geneList = PLT2.direct.regulated.genes[["activated"]], 
                               main.title = "Directly activated regulated.genes",
                               output.dir = output.dir)
getHeatMap.ggplot2.Long.Tissue(geneList = PLT2.direct.regulated.genes[["repressed"]], 
                               main.title = "Directly repressed regulated.genes",
                               output.dir = output.dir)

dir.create(path = "../Rdata/geneLists", showWarnings = F)
write.csv(file="../Rdata/geneLists/PLT2.direct.regulated.genes.activated.csv", PLT2.direct.regulated.genes[["activated"]])
write.csv(file="../Rdata/geneLists/PLT2.direct.regulated.genes.repressed.csv", PLT2.direct.regulated.genes[["repressed"]])
