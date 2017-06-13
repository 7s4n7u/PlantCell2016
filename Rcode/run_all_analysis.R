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

#Set working directory
code.dir <- "~/Documents/ScheresLab/Manuscripts/PLT_Manuscript_code/Rcode/"
setwd(code.dir)

#Run gene expression analysis for PLT1, PLT3, PLT4, PLT5 and PLT7
source("../Rcode/run_PLT13457_OE_limma.R")

#Run gene expression analysis for PLT2
source("../Rcode/run_PLT2_OE_limma.R")

#Run gene expression analysis for the induction of PLT2 in QC
source("../Rcode/run_pPLT2_QC_limma.R")

#Run analysis on PLT regulated genes
source("../Rcode/run_PLT_regulated.genes_analysis.R")

#Generate Circos plots
source("../Rcode/write_CIRCOS_from_PLT_regulated_genes.R")

#Run GO enrichment analysis
source("../Rcode/run_topGO_analysis.R")

#Run gene expression analysis for the dataset Lewis et al., 2013
source("../Rcode/run_Lewis2013_limma.R")

#Run gene expression analysis for the dataset Bargmann et al., 2013
source("../Rcode/run_Bargmann2013_limma.R")

#Run analysis on the overlap between PLT response and auxin response
source("../Rcode/PLT_auxin_response_analysis.R")

#Run analysis on the distances of PLT2 peaks to regulated genes
source("../Rcode/run_peaks_to_TSS_analysis.R")

#Plot PLT2 motif distribution in PLT2 peaks
source("../Rcode/plot_motif_occurrence_in_peaks.R")

#Plot stack of ANT-like (AIL) motifs
source("../Rcode/plot_AIL_motifs_stack.R")

#Generate table with PLT2 binding sites in CNS regions
source("../Rcode/get_motif_occurrences_in_CNS.R")

#Create compendium table
source("../Rcode/make_compendium.R")

#Session info
devtools::session_info()
