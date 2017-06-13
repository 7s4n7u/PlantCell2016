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

#Load PLT regulated.genes
PLT_regulated.genes <- load.PLT.regulated.genes()

#datasetID <- "rootLongAtlas"
#regulation <- "activated"

for(datasetID in c("rootLongAtlas", "rootTissueAtlas"))
{
  for(regulation in c("activated", "repressed"))
  {
    
    #Load dataset
    widedata <- get.gene.exp.dataset(PLT_regulated.genes[[regulation]], datasetID)
    myData <- as.matrix(widedata)
    #get z-scores
    myData <- t(apply(myData,1,scale))
    
    colnames(myData) <- colnames(widedata)
    row.names(myData) <- row.names(widedata)
    
    if(datasetID=="rootLongAtlas")
      #Rename Columella as Col in Longitudinal RootAtlas dataset
      colnames(myData)[1] <- "Col"
    
    myData <- as.data.frame(myData)
    myTable <- cbind(myData, 
                     highest.z.score=as.vector(apply(myData,1,function(x){names(x)[which(x==max(x))]}))
    )
    write.csv(file=paste("../Rdata/Tables/z_score_",regulation,"_",datasetID,".csv",sep=""), myTable)
    
  }
}
