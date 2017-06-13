#!/bin/bash
#
#$ -wd PLT2/XSQ/Input
#$ -j n
#$ -S /bin/bash
#$ -m abe
#$ -M luca.santuari@wur.nl
#$ -e PLT2_Input.log
#$ -o PLT2_Input.out
#$ -N BWA_Input
#$ -pe orte 6

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

bwa samse -r "@RG\tID:PLT2_Input\tSM:Input\tPU:Input\tLB:Input\tPL:SOLID" /home/santu001/genome/BWA_index_CS/TAIR10_chr_all.fas aln.sai p1.BasRutjensInput_T_F3.fq.gz | samtools view -bSh - > aln_unsorted.bam

samtools sort -o aln_unsorted.bam deleteme > PLT2_Input.bam
samtools index PLT2_Input.bam
samtools flagstat PLT2_Input.bam > PLT2_Input.flagstat.txt

samtools view -b -F 1548 -q 30 PLT2_Input.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c >  PLT2_Input.tagAlign.gz
samtools view -b -F 1548 -q 30 PLT2_Input.bam | bamToBed -i stdin | awk '{print 1 "\t" $3-$2 "\t" $6 "\t" $1 "\t" $2+1}' | gzip -c > PLT2_Input.csar.gz

samtools rmdup -s PLT2_Input.bam PLT2_Input.rmdup.bam
samtools index PLT2_Input.rmdup.bam

samtools view -b -F 1548 -q 30 PLT2_Input.rmdup.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c >  PLT2_Input.rmdup.tagAlign.gz
samtools view -b -F 1548 -q 30 PLT2_Input.rmdup.bam | bamToBed -i stdin | awk '{print 1 "\t" $3-$2 "\t" $6 "\t" $1 "\t" $2+1}' | gzip -c > PLT2_Input.rmdup.csar.gz

