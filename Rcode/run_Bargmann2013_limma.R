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

#Set working directory
setwd("~/Documents/ScheresLab/Manuscripts/PLT_Manuscript_code/Rcode/")

#Load experiment design information
expDetails <- read.csv("../Raw/arrays/Bargmann2013/experiment_design.csv")
#Create labels
labels <- with(expDetails,
               make.names(paste(tissue, treatment, replicate, sep="_"))
)
#Create conditions
samples <- with(expDetails,
                make.names(paste(tissue, treatment, sep="_"))
)

#Load arrays
fileNames <- as.vector(expDetails$Sample)
data.affy=ReadAffy(filenames=paste("../Raw/arrays/Bargmann2013/CEL/", fileNames, sep=""), 
                   cdfname="ath1121501attairgcdf")

#perform RMA
data.rma=affy::rma(data.affy)

#perform array quality control
arrayQualityMetrics(data.rma, outdir="../Raw/arrays/Bargmann2013/QC", force=TRUE)

#extract expression values
expr.rma=exprs(data.rma)

# convert into factors
samples<- as.factor(samples)
# check factors have been assigned
# set up the experimental design
design = model.matrix(~0+samples)
colnames(design)<-levels(samples)

#Make contrast
contrast.matrix <- makeContrasts(
  columella = columella_auxin - columella_control,
  epidermis = epidermis_auxin - epidermis_control,
  intactRoot = intactRoot_auxin - intactRoot_control,
  pericycle = pericycle_auxin - pericycle_control,
  stele = stele_auxin - stele_control,
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
logFC_cutoff <- log2(1.5)
pval_cutoff <- 1e-2

#Check results with the cutoffs used in the manuscript
results <- decideTests(ath1_ebFit, 
                       p.value=pval_cutoff,
                       lfc=logFC_cutoff)
summary(results)

#Load annotation data.frame generated from the BRAINARRAY file "ath1121501attairg.db_19.0.0.tar"
load("../Rdata/BRAINARRAY/ath1121501/annotation.df.rda")

#Create new directory for results
limma.results.dir <- "../Rdata/Bargmann2013/limma_results/"
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

bargmannRoot <- vector(mode="list", length=2)
names(bargmannRoot) <- c("auxinUp","auxinDown")

data <- read.csv(paste(limma.results.dir, "intactRoot.csv",sep=""))
data <- with(data, data.frame(AGI=TAIR, logFC=logFC, adj.P.Val=adj.P.Val, B=B))
data <- remove.dups(data)

data <- data[which(data$adj.P.Val<pval_cutoff),]

bargmannRoot[["auxinUp"]] <- unique(c(
  bargmannRoot[["auxinUp"]],
  as.vector(data$AGI[data$logFC >= logFC_cutoff]))
)
bargmannRoot[["auxinDown"]] <- unique(c(
  bargmannRoot[["auxinDown"]],
  as.vector(data$AGI[data$logFC <= (-logFC_cutoff)]))
)
lapply(bargmannRoot, length)

#Save results
save(file="../Rdata/Bargmann2013/limma_results/Bargmann_intactRoot.rda", bargmannRoot, compress=T)

#Session info
devtools::session_info()
