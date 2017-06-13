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

#Create the compendium table

#remove genes covered by multiple probesets
remove.dups <- function(data)
{
  #remove dups
  AGI.tab <- table(data$AGI)
  AGI.dup <- names(AGI.tab[AGI.tab>1])
  
  idx <- which(data$AGI%in%AGI.dup)
  if(length(idx)>0){
    data <- data[-idx,]
  }
  data
}

#merge expression table with existing compendium
merge.exp.data <- function(fileName, id, data)
{
  data.tmp <- read.csv(fileName)
  data.tmp <- with(data.tmp, data.frame(AGI=TAIR, logFC=logFC, adj.P.Val=adj.P.Val, B=B))
  data.tmp <- remove.dups(data.tmp)
  names(data.tmp)[-1] <- paste(id,names(data.tmp)[-1],sep="")
  data <- merge(data,data.tmp, by="AGI")
  data
}

#merge expression table with existing compendium
merge.exp.data.aragene <- function(fileName, id, data)
{
  data.tmp <- read.csv(fileName)
  data.tmp <- with(data.tmp, data.frame(AGI=AGI, logFC=logFC, adj.P.Val=adj.P.Val, B=B))
  data.tmp <- remove.dups(data.tmp)
  names(data.tmp)[-1] <- paste(id,names(data.tmp)[-1],sep="")
  data <- merge(data,data.tmp, by="AGI",all.x=T)
  data
}

#if(FALSE)
#{

working.dir <- "../Rdata/PLT13457_OE/limma_results/"

res <- data.frame()

regulation <- "DEX"
protein.ID <- "PLT1"
data.start <- read.csv(paste(working.dir, protein.ID,"_",regulation,".csv",sep=""))
data.start <- with(data.start, data.frame(AGI=TAIR, logFC=logFC, adj.P.Val=adj.P.Val, B=B))
data.start <- remove.dups(data.start)
names(data.start)[-1] <- paste(protein.ID,regulation,names(data.start)[-1],sep="_")
res <- data.start

for(regulation in c("DEX"))
{
  for(protein.ID in paste("PLT",c(1,3:5,7),sep=""))
  {
    if(!(regulation=="DEX"&protein.ID=="PLT1"))
    {
      fileName <- paste(working.dir, protein.ID,"_",regulation,".csv",sep="")
      id <- paste(protein.ID,"_",regulation,"_",sep="")
      res <- merge.exp.data(fileName,id,res)
    }
  }
}

working.dir <- "../Rdata/PLT2_OE/limma_results/"

fileName <- paste(working.dir,"PLT2_DEX.csv",sep="")
id <- paste("PLT2","_","DEX","_",sep="")
res <- merge.exp.data.aragene(fileName,id,res)

working.dir <-"../Rdata/pPLT2_QC/limma_results/"

fileName <- paste(working.dir,"dex_1h.csv",sep="")
id <- paste("pPLT2_QC","_","DEX_1h","_",sep="")
res <- merge.exp.data(fileName,id,res)

fileName <- paste(working.dir,"dex_4h.csv",sep="")
id <- paste("pPLT2_QC","_","DEX_4h","_",sep="")
res <- merge.exp.data(fileName,id,res)

#Reordering columns by logFC and adjusted p-values
logFC_names <- names(res)[grep("logFC", names(res))]
logFC_names <- logFC_names[order(logFC_names)]
pval_names <- names(res)[grep("adj.P.Val", names(res))]
pval_names <- pval_names[order(pval_names)]
B_names <- names(res)[grep("_B", names(res))]
B_names <- B_names[order(B_names)]

res.ordered <- res[,c(
  "AGI",
  logFC_names,
  pval_names
)]
names(res.ordered)
head(res.ordered)
dim(res.ordered)

load("../Rdata/BRAINARRAY/ath1121501/annotation.df.rda")
names(annotation.df)[1] <- "AGI"
res.anno <- merge(annotation.df,res.ordered,by="AGI")
compendium <- res.anno
#write.csv(file="../Rdata/Tables/compendium.csv", compendium)
dim(compendium)

compendium.sig <- cbind(compendium,
                    PLT1.sig=0,
                    PLT2.sig=0,
                    PLT3.sig=0,
                    PLT4.sig=0,
                    PLT5.sig=0,
                    PLT7.sig=0,
                    pPLT2_QC_1h.sig=0,
                    pPLT2_QC_4h.sig=0
)

pval_cutoff.ath1 <- 2e-3
logFC_cutoff.ath1 <- log2(1.75)

pval_cutoff.aragene11st <- 2e-2
logFC_cutoff.aragene11st <- log2(1.55)

pval_cutoff.QC <- 5e-2
logFC_cutoff.QC <- log2(1.75)

compendium.sig$PLT1.sig[
  compendium.sig$PLT1_DEX_logFC >= logFC_cutoff.ath1 &
  compendium.sig$PLT1_DEX_adj.P.Val < pval_cutoff.ath1
    ] <- 1
compendium.sig$PLT1.sig[
  compendium.sig$PLT1_DEX_logFC <= (-logFC_cutoff.ath1) &
    compendium.sig$PLT1_DEX_adj.P.Val < pval_cutoff.ath1
  ] <- (-1)

compendium.sig$PLT2.sig[
  compendium.sig$PLT2_DEX_logFC >= logFC_cutoff.aragene11st &
    compendium.sig$PLT2_DEX_adj.P.Val < pval_cutoff.aragene11st
  ] <- 1
compendium.sig$PLT2.sig[
  compendium.sig$PLT2_DEX_logFC <= (-logFC_cutoff.aragene11st) &
    compendium.sig$PLT2_DEX_adj.P.Val < pval_cutoff.aragene11st
  ] <- (-1)
compendium.sig$PLT2.sig[is.na(compendium.sig$PLT2_DEX_logFC)] <- "NA"

compendium.sig$PLT3.sig[
  compendium.sig$PLT3_DEX_logFC >= logFC_cutoff.ath1 &
    compendium.sig$PLT3_DEX_adj.P.Val < pval_cutoff.ath1
  ] <- 1
compendium.sig$PLT3.sig[
  compendium.sig$PLT3_DEX_logFC <= (-logFC_cutoff.ath1) &
    compendium.sig$PLT3_DEX_adj.P.Val < pval_cutoff.ath1
  ] <- (-1)

compendium.sig$PLT4.sig[
  compendium.sig$PLT4_DEX_logFC >= logFC_cutoff.ath1 &
    compendium.sig$PLT4_DEX_adj.P.Val < pval_cutoff.ath1
  ] <- 1
compendium.sig$PLT4.sig[
  compendium.sig$PLT4_DEX_logFC <= (-logFC_cutoff.ath1) &
    compendium.sig$PLT4_DEX_adj.P.Val < pval_cutoff.ath1
  ] <- (-1)

compendium.sig$PLT5.sig[
  compendium.sig$PLT5_DEX_logFC >= logFC_cutoff.ath1 &
    compendium.sig$PLT5_DEX_adj.P.Val < pval_cutoff.ath1
  ] <- 1
compendium.sig$PLT5.sig[
  compendium.sig$PLT5_DEX_logFC <= (-logFC_cutoff.ath1) &
    compendium.sig$PLT5_DEX_adj.P.Val < pval_cutoff.ath1
  ] <- (-1)

compendium.sig$PLT7.sig[
  compendium.sig$PLT7_DEX_logFC >= logFC_cutoff.ath1 &
    compendium.sig$PLT7_DEX_adj.P.Val < pval_cutoff.ath1
  ] <- 1
compendium.sig$PLT7.sig[
  compendium.sig$PLT7_DEX_logFC <= (-logFC_cutoff.ath1) &
    compendium.sig$PLT7_DEX_adj.P.Val < pval_cutoff.ath1
  ] <- (-1)

compendium.sig$pPLT2_QC_1h.sig[
  compendium.sig$pPLT2_QC_DEX_1h_logFC >= logFC_cutoff.QC &
    compendium.sig$pPLT2_QC_DEX_1h_adj.P.Val < pval_cutoff.QC
  ] <- 1
compendium.sig$pPLT2_QC_1h.sig[
  compendium.sig$pPLT2_QC_DEX_1h_logFC <= (-logFC_cutoff.QC) &
    compendium.sig$pPLT2_QC_DEX_1h_adj.P.Val < pval_cutoff.QC
  ] <- (-1)

compendium.sig$pPLT2_QC_4h.sig[
  compendium.sig$pPLT2_QC_DEX_4h_logFC >= logFC_cutoff.QC &
    compendium.sig$pPLT2_QC_DEX_4h_adj.P.Val < pval_cutoff.QC
  ] <- 1
compendium.sig$pPLT2_QC_4h.sig[
  compendium.sig$pPLT2_QC_DEX_4h_logFC <= (-logFC_cutoff.QC) &
    compendium.sig$pPLT2_QC_DEX_4h_adj.P.Val < pval_cutoff.QC
  ] <- (-1)

compendium.sig$PLT1_DEX_logFC <- signif(compendium.sig$PLT1_DEX_logFC, digits = 2)
compendium.sig$PLT2_DEX_logFC <- signif(compendium.sig$PLT2_DEX_logFC, digits = 2)
compendium.sig$PLT3_DEX_logFC <- signif(compendium.sig$PLT3_DEX_logFC, digits = 2)
compendium.sig$PLT4_DEX_logFC <- signif(compendium.sig$PLT4_DEX_logFC, digits = 2)
compendium.sig$PLT5_DEX_logFC <- signif(compendium.sig$PLT5_DEX_logFC, digits = 2)
compendium.sig$PLT7_DEX_logFC <- signif(compendium.sig$PLT7_DEX_logFC, digits = 2)
compendium.sig$pPLT2_QC_DEX_1h_logFC <- signif(compendium.sig$pPLT2_QC_DEX_1h_logFC, digits = 2)
compendium.sig$pPLT2_QC_DEX_4h_logFC <- signif(compendium.sig$pPLT2_QC_DEX_4h_logFC, digits = 2)

write.csv(file="../Rdata/Tables/compendium_sig.csv", compendium.sig)
