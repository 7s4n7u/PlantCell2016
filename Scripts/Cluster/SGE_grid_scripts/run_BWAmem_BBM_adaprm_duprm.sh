#!/bin/bash
#
#$ -wd BBM
#$ -j n
#$ -S /bin/bash
#$ -m abe
#$ -M luca.santuari@wur.nl
#$ -o BBM.out
#$ -e BBM.log
#$ -N BBM
#$ -pe orte 10

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
#fastq-mcf Version: 1.04.676
#bwa Version: 0.7.5a-r405
#samtools Version: 0.1.19-96b5f2294a
#bedtools v2.17.0

#The file illumina-adapters.fasta contains Illumina sequencing paired-end and mate pair adapters

fastq-mcf illumina-adapters.fasta ../FASTQ/pBBM-BBM-YFP_NoIndex_L001_R1_001.fastq.gz 2> BBM_adaprm.log | gzip -c > ../FASTQ/BBM_adaprm.fastq.gz
bwa mem -t 10 -R '@RG\tID:BBM\tSM:BBM' ../BWA_index/TAIR10_chr_all.fas ../FASTQ/BBM_adaprm.fastq.gz | samtools view -bSh - > aln_unsorted.bam
samtools sort -o aln_unsorted.bam deleteme > BBM_adaprm.bam

samtools index BBM_adaprm.bam
samtools flagstat BBM_adaprm.bam > BBM_adaprm.flagstat

samtools rmdup -s BBM_adaprm.bam BBM_adaprm.rmdup.bam
samtools index BBM_adaprm.rmdup.bam
samtools flagstat BBM_adaprm.rmdup.bam

samtools view -b -F 1548 -q 30 BBM_adaprm.rmdup.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c >  BBM_adaprm.rmdup.tagAlign.gz
samtools view -b -F 1548 -q 30 BBM_adaprm.rmdup.bam | bamToBed -i stdin | awk '{print 1 "\t" $3-$2 "\t" $6 "\t" $1 "\t" $2+1}' | gzip -c > BBM_adaprm.rmdup.csar.gz

