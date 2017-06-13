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
library(affy)
library(limma)
library(arrayQualityMetrics)

#Load BRAINARRAY array libraries for ath1121501 arrays
library(ath1121501attairg.db)
library(ath1121501attairgcdf)
library(ath1121501attairgprobe)

#load the auxiliary functions
source("diffExp_analysis_functions.R")

#Load experiment design information
expDetails <- read.csv("../Raw/arrays/pPLT2_QC/experiment_design.csv")
#Create labels
labels <- with(expDetails,
               make.names(paste(treatment, time, replicate, sep="_"))
)
#Create conditions
samples <- with(expDetails,
                make.names(paste(treatment, time, sep="_"))
)

#Load arrays
fileNames <- as.vector(expDetails$Sample)
data.affy=ReadAffy(filenames=paste("../Raw/arrays/pPLT2_QC/CEL/", fileNames, ".gz", sep=""), 
                   cdfname="ath1121501attairgcdf")

#perform RMA
data.rma=affy::rma(data.affy)

#perform array quality control
arrayQualityMetrics(data.rma, outdir="../Raw/arrays/pPLT2_QC/QC", force=TRUE)

#extract expression values
expr.rma=exprs(data.rma)

#write expression matrix in csv format
write.csv(file="../Rdata/pPLT2_QC/limma_results/expr_matrix.csv", expr.rma)

# convert into factors
samples<- as.factor(samples)
# check factors have been assigned
# set up the experimental design
design = model.matrix(~0+samples)
colnames(design)<-levels(samples)

#Make contrast
contrast.matrix <- makeContrasts(
  dex_1h = dex_1h - mock_1h,
  dex_4h = dex_4h - mock_1h,
  levels=design 
)

# check the contrast matrix
contrast.matrix

#Only select probes that target genes for the DE analysis
#Remove the 64 control probes
expr.rma <- expr.rma[grep("^AT.+\\_at", row.names(expr.rma)),]

#fit the linear model per gene 
fit <- lmFit(expr.rma, design)

# Now the contrast matrix is combined with the per-probeset linear model fit.
ath1_fits <- contrasts.fit(fit, contrast.matrix)
ath1_ebFit <- eBayes(ath1_fits)

#cutoffs for p.value and for logFC
logFC_cutoff <- log2(1.75)
pval_cutoff <- 5e-2

#Check results with the cutoffs used in the manuscript
results <- decideTests(ath1_ebFit, 
                       p.value=pval_cutoff,
                       lfc=logFC_cutoff)
summary(results)

#Load annotation data.frame generated from the BRAINARRAY file "ath1121501attairg.db_19.0.0.tar"
load("../Rdata/BRAINARRAY/ath1121501/annotation.df.rda")

#Create new directory for results
limma.results.dir <- "../Rdata/pPLT2_QC/limma_results/"
dir.create(limma.results.dir, showWarnings = F)

#Write annotated results
for(ctsID in colnames(ath1_ebFit))
{
  #ctsID <- colnames(ath1_ebFit)[1]
  print(paste("Writing results for",ctsID))
  res <- topTable(ath1_ebFit,coef=ctsID, adjust.method="BH", 
                  sort.by="p", number=Inf)
  res <- cbind(PROBEID=row.names(res), res)
  res <- merge(res, annotation.df,by="PROBEID")
  res <- res[order(res$adj.P.Val),]
  write.csv(file=paste(limma.results.dir, ctsID, ".csv", sep=""), res)
}


#Extract significant genes

PLT2_regulated.genes_QC <- list()
#QC 1h
data.1h <- read.csv(paste(limma.results.dir,"dex_1h.csv",sep=""))
data.1h <- with(data.1h, data.frame(AGI=TAIR, logFC=logFC, adj.P.Val=adj.P.Val, B=B))
#Remove genes covered by multiple probe sets
data.1h <- remove.dups(data.1h)

PLT2_regulated.genes_QC[["1h.up"]] <- with(data.1h, unique(as.vector(AGI[logFC >= logFC_cutoff &
                                                                   adj.P.Val<pval_cutoff])))
PLT2_regulated.genes_QC[["1h.down"]] <- with(data.1h, unique(as.vector(AGI[logFC <= (-logFC_cutoff) &
                                                                     adj.P.Val<pval_cutoff])))

#QC 4h
data.4h <- read.csv(paste(limma.results.dir,"dex_4h.csv",sep=""))
data.4h <- with(data.4h, data.frame(AGI=TAIR, logFC=logFC, adj.P.Val=adj.P.Val, B=B))
#Remove genes covered by multiple probe sets
data.4h <- remove.dups(data.4h)
PLT2_regulated.genes_QC[["4h.up"]] <- with(data.4h, unique(as.vector(AGI[logFC >= logFC_cutoff & 
                                                                   adj.P.Val<pval_cutoff])))
PLT2_regulated.genes_QC[["4h.down"]] <- with(data.4h, unique(as.vector(AGI[logFC <= (-logFC_cutoff) & 
                                                                     adj.P.Val<pval_cutoff])))

#Save results
save(file=paste(limma.results.dir,"PLT2_regulated.genes_QC.rda",sep=""), PLT2_regulated.genes_QC, compress=T)

#Save all PLT2 regulated genes
#Load PLT2 regulated.genes from the ath1 arrays
load("../Rdata/PLT2_OE/limma_results/PLT2_regulated.genes.ath1.FC1.55_pval2e-2.rda")
#Load PLT2 QC regulated genes
load(paste(limma.results.dir,"PLT2_regulated.genes_QC.rda",sep=""))

PLT2_regulated.genes.seedlings.QC <- list()
PLT2_regulated.genes.seedlings.QC[["activated"]] <- unique(c(
  PLT2_regulated.genes.ath1[["DEXup"]],
  PLT2_regulated.genes_QC[["1h.up"]],
  PLT2_regulated.genes_QC[["4h.up"]]
))
PLT2_regulated.genes.seedlings.QC[["activated"]] <- PLT2_regulated.genes.seedlings.QC[["activated"]][
  order(PLT2_regulated.genes.seedlings.QC[["activated"]])]

PLT2_regulated.genes.seedlings.QC[["repressed"]] <- unique(c(
  PLT2_regulated.genes.ath1[["DEXdown"]],
  PLT2_regulated.genes_QC[["1h.down"]],
  PLT2_regulated.genes_QC[["4h.down"]]
))
PLT2_regulated.genes.seedlings.QC[["repressed"]] <- PLT2_regulated.genes.seedlings.QC[["repressed"]][
  order(PLT2_regulated.genes.seedlings.QC[["repressed"]])]
lapply(PLT2_regulated.genes.seedlings.QC, length)

#Save results
save(file=paste(limma.results.dir,"PLT2_regulated.genes.seedlings.QC.rda",sep=""), PLT2_regulated.genes.seedlings.QC, compress=T)

