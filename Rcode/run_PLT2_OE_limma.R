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


#Requires Platform Design package for aragene11st from Bioconductor
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("pd.aragene.1.1.st")
#Carvalho B. pd.aragene.1.1.st: Platform Design Info for Affymetrix AraGene-1_1-st. R package version 3.12.0. 

library(oligo)
library(limma)
library(arrayQualityMetrics)

#Array package from http://nmg-r.bioinformatics.nl/NuGO_R.html WUR
library(aragene11st.db)

#load the auxiliary functions
source("diffExp_analysis_functions.R")

desc <- read.csv("../Raw/arrays/PLT2_OE/experiment_design.csv")
desc$filename <- paste(desc$filename, ".CEL", sep="")
rownames(desc) <- desc$filename

fileNames <- paste("../Raw/arrays/PLT2_OE/CEL/", desc$filename, ".gz", sep="")

affyGeneFS <- oligo::read.celfiles(fileNames)

#Perform RMA normalization
geneCore <- oligo::rma(affyGeneFS, target = "core")
featureData(geneCore) <- getNetAffx(geneCore, "transcript")

#perform array quality control
arrayQualityMetrics(geneCore, outdir="../Raw/arrays/PLT2_OE/QC", force=TRUE)

#extract expression values
eset <- exprs(geneCore)

#write expression matrix in csv format
write.csv(file="../Rdata/PLT2_OE/limma_results/expr_matrix.csv", eset)

colnames(eset) <- desc$replicates

#Create design matrix
design <- model.matrix(~0+factor(condition), desc)
colnames(design) <- as.vector(levels(factor(desc$condition)))

#Fit linear model per gene
fit <- lmFit(eset, design)
#Create contrast matrix
contrast.matrix <- makeContrasts(PLT2_DEX=DEX-mock, levels=design)

aragene11.fit <- contrasts.fit(fit, contrast.matrix)
aragene11.ebFit <- eBayes(aragene11.fit)

results <- decideTests(aragene11.fit, p.value=2e-2,lfc=log2(1.55))
summary(results)

#Load ACCNUMDF for aragene11st
load("../Rdata/aragene11st/ACCNUMDF.rda")

#Create new directory for results
limma.results.dir <- "../Rdata/PLT2_OE/limma_results/"
dir.create(limma.results.dir, showWarnings = F)

results <- topTable(aragene11.ebFit, coef=1, adjust="BH", sort.by="B", number=39000)
results <- cbind(AffyID=row.names(results), results)
results <- merge(results, ACCNUMDF, by="AffyID")
results <- results[order(results$adj.P.Val,decreasing = F),]
write.csv(results, file = paste(limma.results.dir, "PLT2_DEX.csv",sep="") )


#Extract relevant genes

#cutoffs for aragene11st
logFC_cutoff <- log2(1.55)
pval_cutoff <- 2e-2

dataset <- read.csv(paste(limma.results.dir,"PLT2_DEX.csv",sep=""))
dataset <- with(dataset, data.frame(AGI=as.vector(AGI),
                                    logFC=logFC,
                                    P.Value=P.Value,
                                    adj.P.Val=adj.P.Val,
                                    B=B
))
#Remove genes covered by multiple probe sets
dataset <- remove.dups(dataset)

DEXup <- with(dataset,
              as.vector(AGI[which(adj.P.Val < pval_cutoff &
                                    logFC >= logFC_cutoff ) ])
)
DEXdown <- with(dataset,
                as.vector(AGI[which(adj.P.Val < pval_cutoff &
                                      logFC <= (-logFC_cutoff) ) ])
)


PLT2_regulated.genes.aragene11st <- list(DEXup=unique(DEXup),
                                 DEXdown=unique(DEXdown))
lapply(PLT2_regulated.genes.aragene11st, length)

save(file=paste(limma.results.dir,"PLT2_regulated.genes.aragene11st.FC1.55_pval2e-2.rda",sep=""), PLT2_regulated.genes.aragene11st)

write.csv(file="../Rdata/aragene11st/aragene11st.AGI.csv", unique(as.vector(dataset$AGI[order(dataset$AGI)])))

#Load ATH1 gene AGI
ath1.AGI <- as.vector(read.csv("../Rdata/BRAINARRAY/ath1121501/ath1.AGI.BRAINARRAY.csv")[,2])

PLT2_regulated.genes.ath1 <- list(DEXup=unique(PLT2_regulated.genes.aragene11st$DEXup[which(PLT2_regulated.genes.aragene11st$DEXup%in%ath1.AGI)]),
                          DEXdown=unique(PLT2_regulated.genes.aragene11st$DEXdown[which(PLT2_regulated.genes.aragene11st$DEXdown%in%ath1.AGI)])
)
lapply(PLT2_regulated.genes.ath1, length)

save(file=paste(limma.results.dir,"PLT2_regulated.genes.ath1.FC1.55_pval2e-2.rda",sep=""), PLT2_regulated.genes.ath1)

###PLT2 bound aragene11st genes
anno.4Kb <- read.csv("../Rdata/ChIP/PLT2/Homer/IP1_IP2_pooled_1FDR_chrNum_HomerAnnotation_4Kb_fromTSS.csv")
#Load aragene11st gene AGI
aragene11st.AGI <- as.vector(read.csv("../Rdata/aragene11st/aragene11st.AGI.csv")[,2])
bound <- unique(as.vector(substr(anno.4Kb$Nearest.PromoterID,1,9)))
bound <- bound[bound%in%aragene11st.AGI]

#Proportions
#Bound and activated
print(paste(
  signif(length(PLT2_regulated.genes.aragene11st$DEXup[which(PLT2_regulated.genes.aragene11st$DEXup%in%bound)])/length(PLT2_regulated.genes.aragene11st$DEXup)*100,
         digits = 2),
  "% of PLT2 genes activated in whole seedling are directly regulated.genes", sep=""
))
#Bound and repressed
print(paste(
  signif(length(PLT2_regulated.genes.aragene11st$DEXdown[which(PLT2_regulated.genes.aragene11st$DEXdown%in%bound)])/length(PLT2_regulated.genes.aragene11st$DEXdown)*100,
         digits = 2),
  "% of PLT2 genes repressed in whole seedling are direct regulated.genes", sep=""
))

#Create ALLPLT target lists
load("../Rdata/PLT13457_OE/limma_results/PLT_13457_regulated.genes.ath1.FC1.75_pval2e-3.rda")
#Load PLT2 ath1 genes
load("../Rdata/PLT2_OE/limma_results/PLT2_regulated.genes.ath1.FC1.55_pval2e-2.rda")

ALLPLTsetlistUP <- list(PLT1=ALLPLTsetlistUP_DEX$PLT1,
                        PLT2=PLT2_regulated.genes.ath1$DEXup,
                        PLT3=ALLPLTsetlistUP_DEX$PLT3,
                        PLT4=ALLPLTsetlistUP_DEX$PLT4,
                        PLT5=ALLPLTsetlistUP_DEX$PLT5,
                        PLT7=ALLPLTsetlistUP_DEX$PLT7
)

ALLPLTsetlistDOWN <- list(PLT1=ALLPLTsetlistDOWN_DEX$PLT1,
                          PLT2=PLT2_regulated.genes.ath1$DEXdown,
                          PLT3=ALLPLTsetlistDOWN_DEX$PLT3,
                          PLT4=ALLPLTsetlistDOWN_DEX$PLT4,
                          PLT5=ALLPLTsetlistDOWN_DEX$PLT5,
                          PLT7=ALLPLTsetlistDOWN_DEX$PLT7
)

ALLPLTsetlistUP_DEX <- ALLPLTsetlistUP
ALLPLTsetlistDOWN_DEX <- ALLPLTsetlistDOWN

lapply(ALLPLTsetlistUP_DEX, length)
lapply(ALLPLTsetlistDOWN_DEX, length)

save(file="../Rdata/PLT_OE/ALLPLT_regulated.genes.rda", ALLPLTsetlistUP_DEX, ALLPLTsetlistDOWN_DEX)

#
