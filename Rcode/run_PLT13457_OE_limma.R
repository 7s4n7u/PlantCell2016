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
expDetails <- read.csv("../Raw/arrays/PLT13457_OE/experiment_design.csv")
#Create labels
labels <- with(expDetails,
               make.names(paste(genotype, treatment, replicate, sep="_"))
)
#Create conditions
samples <- with(expDetails,
                   make.names(paste(genotype, treatment, sep="_"))
)
#List of PLT proteins
proteinList <- paste("PLT",c(1,3:5,7),sep="")

#Load arrays
fileNames <- paste(as.vector(expDetails$Sample),".gz", sep="")
data.affy=ReadAffy(filenames=paste("../Raw/arrays/PLT13457_OE/CEL/", fileNames, sep=""), 
                   cdfname="ath1121501attairgcdf")

#perform RMA
data.rma=affy::rma(data.affy)

#perform array quality control
arrayQualityMetrics(data.rma, outdir="../Raw/arrays/PLT13457_OE/QC", force=TRUE)

#extract expression values
expr.rma=exprs(data.rma)

#write expression matrix in csv format
write.csv(file="../Rdata/PLT13457_OE/limma_results/expr_matrix.csv", expr.rma)

# convert into factors
samples<- as.factor(samples)
# check factors have been assigned
# set up the experimental design
design = model.matrix(~0+samples)
colnames(design)<-levels(samples)

#Make contrast
contrast.matrix <- makeContrasts(
  PLT1_DEX = X35S.PLT1.GR_dex - X35S.PLT4.GR_mock,
  PLT3_DEX = X35S.PLT3.GR_dex - X35S.PLT4.GR_mock,
  PLT4_DEX = X35S.PLT4.GR_dex - X35S.PLT4.GR_mock,
  PLT5_DEX = X35S.PLT5.GR_dex - X35S.PLT4.GR_mock,
  PLT7_DEX = X35S.PLT7.GR_dex - X35S.PLT4.GR_mock,
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
pval_cutoff <- 2e-3

#Check results with the cutoffs used in the manuscript
results <- decideTests(ath1_ebFit, 
                       p.value=pval_cutoff,
                       lfc=logFC_cutoff)
summary(results)

#Load annotation data.frame generated from the BRAINARRAY file "ath1121501attairg.db_19.0.0.tar"
load("../Rdata/BRAINARRAY/ath1121501/annotation.df.rda")

#Create new directory for results
limma.results.dir <- "../Rdata/PLT13457_OE/limma_results/"
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
ALLPLTsetlistUP_DEX <- ALLPLTsetlistDOWN_DEX <- vector(mode="list", length=length(proteinList))
names(ALLPLTsetlistUP_DEX) <- names(ALLPLTsetlistDOWN_DEX) <- proteinList

for(proteinID in proteinList)
{
  
  data.DEX <- read.csv(paste(limma.results.dir,proteinID,"_DEX.csv",sep=""))
  data.DEX <- with(data.DEX, data.frame(AGI=TAIR, logFC=logFC, P.Value=P.Value, adj.P.Val=adj.P.Val, B=B))
  data.DEX <- remove.dups(data.DEX)
  
  data.DEX <- data.DEX[which(data.DEX$adj.P.Val < pval_cutoff),]
  ALLPLTsetlistUP_DEX[[proteinID]] <- unique(as.vector(data.DEX$AGI[which(data.DEX$logFC >= logFC_cutoff)]))
  ALLPLTsetlistDOWN_DEX[[proteinID]] <- unique(as.vector(data.DEX$AGI[which(data.DEX$logFC <= (-logFC_cutoff))]))
  
}

lapply(ALLPLTsetlistUP_DEX, length)
lapply(ALLPLTsetlistDOWN_DEX, length)

#Save results
save(file=paste(limma.results.dir,"PLT_13457_regulated.genes.ath1.FC1.75_pval2e-3.rda",sep=""), 
     ALLPLTsetlistUP_DEX, ALLPLTsetlistDOWN_DEX, compress=T)
