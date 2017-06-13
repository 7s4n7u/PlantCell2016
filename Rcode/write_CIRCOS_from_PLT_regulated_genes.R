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

require(systemPipeR)

#Load ALLPLT lists
load("../Rdata/PLT_OE/ALLPLT_regulated.genes.rda")

for(regulation in c("activated", "repressed"))
{
  
  #regulation <- "activated"
  
  if(regulation=="activated")
  {
    geneSetList <- ALLPLTsetlistUP_DEX
  }else{
    geneSetList <- ALLPLTsetlistDOWN_DEX
  }
  
  circos.output.dir <- paste("../Figures/Circos/",regulation,"/",sep="")
  
  lapply(geneSetList, length)
  OLlist <- overLapper(setlist=geneSetList, sep="_", type="vennsets")
  
  PLT.vec <- paste("PLT", c(1:5,7),sep="")
  start.vec <- vector(mode="integer", length=6)
  names(start.vec) <- PLT.vec
  segdup.vec <- vector(mode="integer", length=length(2:6))
  names(segdup.vec) <- 2:6
  
  #remove files
  for(len in 2:6)
  {
    #len <- 2
    fileName <- paste(circos.output.dir,"segdup",len,"_chx.txt",sep="")
    if(file.exists(fileName)){file.remove(fileName)}
  }
  
  #write cyto file
  fileName <- paste(circos.output.dir,"xx_chx.txt",sep="")
  if(file.exists(fileName)){file.remove(fileName)}
  for(plt.ID in PLT.vec)
  {
    write(
      paste(
        "chr - ", plt.ID, " ", plt.ID, " ", 0 , " ", length(geneSetList[[plt.ID]]), " grey"
      ), append = T, sep = "\n",
      file=fileName
    )
  }
  
  for(setName in rev(names(vennlist(OLlist))) )
  {
    #setName <- rev(names(OLlist$Venn_List))[2]
    PLT.vec.set <- strsplit(setName, "_")[[1]]
    set.len <- length(vennlist(OLlist)[[setName]])
    len <- length(PLT.vec.set)
    if(len>1)
    {
      #print(len)
      #print(PLT.vec.set)
      
      fileName <- paste(circos.output.dir, "segdup",len,"_chx.txt",sep="")
      
      for(i in 1:(length(PLT.vec.set)-1))
      {
        #i <- 1
        start.val.i <- as.vector(start.vec[which(names(start.vec)==PLT.vec.set[i])])
        for(j in (i+1):length(PLT.vec.set))
        {
          #j <- (i+1)
          start.val.j <- as.vector(start.vec[which(names(start.vec)==PLT.vec.set[j])])
          
          write(
            paste(
              paste("segdup", segdup.vec[len-1], sep="_"), "\t",
              PLT.vec.set[i], "\t", 
              start.val.i, "\t",
              start.val.i+set.len
            ), append = T, sep = "\n",
            file=fileName
          )
          write(
            paste(
              paste("segdup", segdup.vec[len-1], sep="_"), "\t",
              PLT.vec.set[j], "\t", 
              start.val.j, "\t",
              start.val.j+set.len
            ), append = T, sep = "\n",
            file=fileName
          )
          segdup.vec[len-1] <- segdup.vec[len-1]+1
        }
      }
      for(plt.ID in PLT.vec.set)
      {
        start.vec[which(names(start.vec)==plt.ID)] <- start.vec[which(names(start.vec)==plt.ID)]+set.len
      }
    }
  }
  
  ###############################################################################################################
  
  #To run it:
  setwd(circos.output.dir)
  cmd.line <- "/usr/local/ActivePerl-5.16/bin/perl /Applications/Bioinformatics/Circos/circos-0.67-7/bin/circos -conf circos.conf"
  try(system(cmd.line, intern = TRUE, ignore.stderr = TRUE))
  #cmd.line <- "/usr/local/ActivePerl-5.16/bin/perl /Applications/Bioinformatics/Circos/circos-0.67-7/bin/circos -conf circos_2way.conf"
  #try(system(cmd.line, intern = TRUE, ignore.stderr = TRUE))
  setwd(code.dir)
  # 
  
}
