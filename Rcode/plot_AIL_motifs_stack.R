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
library(motifStack)

motif.dir <- "../Rdata/Motif_analysis/AIL_motifs/"
output.dir <- "../Figures/Motif_analysis/"

motifs <- list()
for(fileName in list.files(path = motif.dir, pattern = "txt"))
{
  
  motifName <- strsplit(fileName,"\\.")[[1]][1]
  motif <- t(read.table(paste(motif.dir,fileName,sep="")))
  rownames(motif) <- c("A","C","G","T")
  motif.pfm <- new("pfm", mat=motif, name=paste(motifName,"motif",sep="_"))
  if(motifName=="PLT2_SELEX"|motifName=="PLT5_SELEX"){
    #reverse complement
    motif.pfm <- matrixReverseComplement(motif.pfm)
  }
  
  motifs[[motifName]] <- motif.pfm

}

## plot stacks
cairo_ps(filename = paste(output.dir, "AIL_motifStack.eps", sep=""))
motifStack(rev(motifs), layout="stack", ncex=1.0)
dev.off()

svg(filename = paste(output.dir, "AIL_motifStack.svg", sep=""),
    family="Helvetica")
    #height = 12, width = 10)
motifStack(rev(motifs), layout="stack", ncex=1.0)
dev.off()

