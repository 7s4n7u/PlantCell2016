#!/bin/bash
#
#$ -wd PLT2/XSQ/IP1
#$ -j n
#$ -S /bin/bash
#$ -m abe
#$ -M luca.santuari@wur.nl
#$ -e PLT2_IP1.log
#$ -o PLT2_IP1.out
#$ -N BWA_IP1
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

bwa samse -r "@RG\tID:PLT2_IP1\tSM:IP1\tPU:IP1\tLB:IP1\tPL:SOLID" genome/BWA_index_CS/TAIR10_chr_all.fas aln.sai p1.BasRutjensIP1_T_F3.fq.gz | samtools view -bSh - > aln_unsorted.bam

samtools sort -o aln_unsorted.bam deleteme > PLT2_IP1.bam
samtools index PLT2_IP1.bam
samtools flagstat PLT2_IP1.bam > PLT2_IP1.flagstat.txt

samtools rmdup -s PLT2_IP1.bam PLT2_IP1.rmdup.bam
samtools index PLT2_IP1.rmdup.bam

samtools view -b -F 1548 -q 30 PLT2_IP1.rmdup.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c >  PLT2_IP1.rmdup.tagAlign.gz
samtools view -b -F 1548 -q 30 PLT2_IP1.rmdup.bam | bamToBed -i stdin | awk '{print 1 "\t" $3-$2 "\t" $6 "\t" $1 "\t" $2+1}' | gzip -c > PLT2_IP1.rmdup.csar.gz

