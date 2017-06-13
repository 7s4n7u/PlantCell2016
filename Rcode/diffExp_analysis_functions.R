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


#remove genes covered by multiple probesets
#data is a data.frame with the output from limma and a column with TAIR10 AGI gene identifiers
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
