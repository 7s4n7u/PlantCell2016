#!/bin/bash

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

#Dependencies:
#shuffledBed from bedtools v2.17.0

#Calculate the motif distribution in peak regions

ACTIVATED_PEAKS="../../Rdata/ChIP/PLT2/BED/activated/500bp/activated_500bp.bed"
REPRESSED_PEAKS="../../Rdata/ChIP/PLT2/BED/repressed/500bp/repressed_500bp.bed"
PEAK_REGIONS="../../Rdata/ChIP/PLT2/BED/activated_and_repressed/500bp/activated_and_repressed_500bp.bed"
OUTPUT_FILE="../../Rdata/Motif_analysis/Motif_occurrence_in_peaks/motif_distribution.txt"

rm $OUTPUT_FILE
cat $ACTIVATED_PEAKS $REPRESSED_PEAKS | sortBed -i - > $PEAK_REGIONS

SELEX_SITES="../../Rdata/Motif_analysis/Motif_matches/PLT2_SELEX_motif_WG_FIMO_pval_1e-3.score.bed"
GENOME="../../Rdata/TAIR10/TAIR10.genome"

#At least half the PLT2 motif (8 bases) should match for a match to be counted. This applies in particular for the bins at the edges of the considered regions.
#The output contains two columns. The first column shows the iteration number. Iteration 0 records the results for the PLT2 peak regions, while iterations from 1 to iterN record the results for the random selection of genomic regions. Column two contains the distances from the start of the region for the PLT2 motif matches.

intersectBed -a $PEAK_REGIONS -b $SELEX_SITES -wao | awk '$13 >= 8' | uniq | awk -v ITER=0 '{print ITER "\t" $8-$2}' > $OUTPUT_FILE

iterN=1000

for iter in $(seq 1 $iterN);
do
    shuffleBed -i  $PEAK_REGIONS  -g $GENOME -chrom | intersectBed -a - -b $SELEX_SITES -wao | awk '$13 >= 8' | uniq | awk -v ITER=$iter '{print ITER "\t" $8-$2}' >> $OUTPUT_FILE
done
