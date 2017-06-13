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
#annotatePeaks.pl from Homer v3.12

REGULATED_PEAKS="../../Rdata/ChIP/PLT2/BED/repressed/500bp/repressed_500bp_peakCenter.bed"
GENOME="../../Rdata/TAIR10/TAIR10.genome"
SHUFFLED_BED=`which shuffleBed`
OUTPUT_DIR="../../Rdata/peaks_to_TSS_analysis/shuffled/repressed"
iterN=1000

for iter in $(seq 1 $iterN);
do
    echo "Shuffling iter:" $iter
    $SHUFFLED_BED -i $REGULATED_PEAKS -g $GENOME -chrom > $OUTPUT_DIR/shuffledPeaks_$iter\.bed
    annotatePeaks.pl $OUTPUT_DIR/shuffledPeaks_$iter\.bed tair10 > $OUTPUT_DIR/shuffledPeaks_$iter\.txt 2> $OUTPUT_DIR/shuffledPeaks_$iter\.log   
done
