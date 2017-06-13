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

#Dependencies
#gshuf
#bedtools v2.17.0

SELEX_SITES="../../Rdata/Motif_analysis/Motif_matches/PLT2_SELEX_motif_WG_FIMO_pval_1e-3.score.bed"
MIN_OVERLAP=`head -n 1 $SELEX_SITES | awk '{print $3-$2}'`


CNS_REGIONS_HAUDRY="../../Rdata/CNS_datasets/Haudry_2013/AT_CNSv5_1.bed"
CNS_REGIONS_VANDEVELDE="../../Rdata/CNS_datasets/Van_de_Velde_2014/AllFootPrintsFDR0.10_scores.bed"

iterN=1000

#Create the dataset of PLT2 peaks assigned to PLT2 regulated genes
cat "../../Rdata/ChIP/PLT2/BED/activated/300bp/activated_300bp.bed" "../../Rdata/ChIP/PLT2/BED/repressed/300bp/repressed_300bp.bed" |  sortBed -i >  "../../Rdata/ChIP/PLT2/BED/activated_and_repressed/300bp/activated_and_repressed_300bp.bed"

for peaks in "total_peaks" "activated_and_repressed"; do

    echo "Running CNS analysis for " $peaks

    if [ $peaks == "total_peaks" ]; then
        PEAK_DATASET="../../Rdata/ChIP/PLT2/BED/"$peaks"/IP1_IP2_pooled_1FDR_300bp.bed"
    else 
        PEAK_DATASET="../../Rdata/ChIP/PLT2/BED/"$peaks"/300bp/"$peaks"_300bp.bed"
    fi

    BEDINPUT_LINES=`intersectBed -a $SELEX_SITES -b $PEAK_DATASET -wao | awk -v overlap=$MIN_OVERLAP '$13 >= overlap' | uniq | cut -f1-3 | wc -l`
    echo "Motif occurrences in PLT2 peaks: " $BEDINPUT_LINES

    OVERLAP_FILE="../../Rdata/CNS_analysis/Motifs_in_"$peaks".bed"

    intersectBed -a $SELEX_SITES \
        -b $PEAK_DATASET \
        -wao | awk -v overlap=$MIN_OVERLAP '$13 >= overlap' | uniq | cut -f1-3 > $OVERLAP_FILE

    #analysis with CNS from Haudry et al., 2013

    OUTPUT_FILE="../../Rdata/CNS_analysis/Haudry_2013/"$peaks"/results.txt"
    #If the result file exists, delete it
    if [ -a $OUTPUT_FILE ]; then
        rm $OUTPUT_FILE
    fi

    intersectBed -a $OVERLAP_FILE -b $CNS_REGIONS_HAUDRY -wao | awk -v overlap=$MIN_OVERLAP '$12 >= overlap' | uniq | wc -l > $OUTPUT_FILE


    for iter in $(seq 1 $iterN); do

        gshuf  -n $BEDINPUT_LINES $SELEX_SITES | sortBed -i | intersectBed -a - -b $CNS_REGIONS_HAUDRY -wao | awk  -v overlap=$MIN_OVERLAP '$15 >= overlap' | uniq | wc -l >> $OUTPUT_FILE

    done

    #analysis with CNS from Van de Velde et al., 2014

    OUTPUT_FILE="../../Rdata/CNS_analysis/Van_de_Velde_2014/"$peaks"/results.txt"
    #If the result file exists, delete it
    if [ -a $OUTPUT_FILE ]; then
        rm $OUTPUT_FILE
    fi

    intersectBed -a $OVERLAP_FILE -b $CNS_REGIONS_VANDEVELDE -wao | awk -v overlap=$MIN_OVERLAP '$9 >= overlap' | uniq | wc -l > $OUTPUT_FILE

    for iter in $(seq 1 $iterN); do
    
        gshuf  -n $BEDINPUT_LINES $SELEX_SITES | sortBed -i | intersectBed -a - -b $CNS_REGIONS_VANDEVELDE -wao | awk  -v overlap=$MIN_OVERLAP '$12 >= overlap' | uniq | wc -l >> $OUTPUT_FILE

    done

done
