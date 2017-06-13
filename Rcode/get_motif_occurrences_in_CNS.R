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


get.stats <- function(fileName, type, cns.type)
{
  data <- as.vector(read.table(fileName,h=F)[,1])
  data.frame(
    type=type,
    cns.type=cns.type,
    observed=data[1],
    expected.median=median(data[-1]),
    expected.MAD=mad(data[-1])
  )
}

get.pvalue <- function(fileName, type, cns.type)
{
  data <- as.vector(read.table(fileName,h=F)[,1])
  iterations <- length(data[-1])
  #calculate empirical pvalue
  pval <- length(which(data[1] <= data[-1]))/iterations
}

res <- data.frame()
res <- rbind(res,
             get.stats("../Rdata/CNS_analysis/Haudry_2013/total_peaks/results.txt", "total_peaks", "Haudry2013"),
             get.stats("../Rdata/CNS_analysis/Haudry_2013/activated_and_repressed/results.txt", "peaks_assigned_to_activated_and_repressed_genes", "Haudry2013"),
             get.stats("../Rdata/CNS_analysis/Van_de_Velde_2014/total_peaks/results.txt", "total_peaks", "Van_de_Velde2014"),
             get.stats("../Rdata/CNS_analysis/Van_de_Velde_2014/activated_and_repressed/results.txt", "peaks_assigned_to_activated_and_repressed_genes", "Van_de_Velde2014")
)

#Calculate empirical p-values
pval <- c(
  get.pvalue("../Rdata/CNS_analysis/Haudry_2013/total_peaks/results.txt", "total_peaks", "Haudry2013"),
  get.pvalue("../Rdata/CNS_analysis/Haudry_2013/activated_and_repressed/results.txt", "peaks_assigned_to_activated_and_repressed_genes", "Haudry2013"),
  get.pvalue("../Rdata/CNS_analysis/Van_de_Velde_2014/total_peaks/results.txt", "total_peaks", "Van_de_Velde2014"),
  get.pvalue("../Rdata/CNS_analysis/Van_de_Velde_2014/activated_and_repressed/results.txt", "peaks_assigned_to_activated_and_repressed_genes", "Van_de_Velde2014")
)

pval_cutoff <- 5e-2
pval_cutoff.Bonferroni <- pval_cutoff/nrow(res)
res <- cbind(res,
             enriched=(res$observed/res$expected.median),
             significant=as.character(pval<pval_cutoff.Bonferroni)
)

write.csv(file="../Rdata/CNS_analysis/motif_occurrences_in_CNS_stats.csv", res)
