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

require(CSAR)

Input <- read.table("NLS_adaprm.rmdup.csar.gz", h=F)
IP <- read.table("BBM_adaprm.rmdup.csar.gz", h=F)

names(Input) <- names(IP) <- c("Nhits", "lengthRead", "strand","chr", "pos")

library(BSgenome.Athaliana.TAIR.TAIR9)
chr.names <- c(as.character(c(1:7)))
seq.length <- seqlengths(Athaliana)
names(seq.length) <- chr.names

Input$chr <- as.vector(Input$chr)
IP$chr <- as.vector(IP$chr)

save(file="CSAR_BED_dataset.rda", Input, IP, compress=T)

nhitsIP <- mappedReads2Nhits(IP, file = "IP", chr = chr.names, chrL = seq.length, w = 200)
nhitsInput <- mappedReads2Nhits(Input, file = "Input", chr = chr.names, chrL = seq.length, w = 200)

BBM <- ChIPseqScore(control = nhitsInput, sample = nhitsIP, file = "BBM", times = 100000)
winIP <- sigWin(BBM)
score2wig(BBM, file = "BBM.wig", name="BBM", t=2, times = 2e6)

#dIP <- distance2Genes(win = winIP, gff = TAIR10_genes)
#genesIP <- genesWithPeaks(dIP)

save(file="CSAR_results.rda", BBM, winIP, compress=T)
#load("PLT2_CSAR_results.rda")
write.csv(file="BBM.csv", as.data.frame(winIP))

iterations <- 10
for(i in 1:10)
{
permutatedWinScores(nn = i, sample = IP,
                    control = Input, fileOutput = "BBM",
                    chr = c(1:7), chrL = seq.length)
}

nulldist <- getPermutatedWinScores(file = "BBM",
                                      nn = 1:iterations)

res <- data.frame()
for(cutoff in c(.05,.01,.005,.001,.0005,.0001))
{
  thrs <- getThreshold(winscores = values(winIP)$score, permutatedScores = nulldist,
                       FDR = cutoff)
  res <- rbind(res, data.frame(thrs, NPeaks=length(winIP[values(winIP)$score>=thrs$threshold])))
}
write.csv(file="BBM_FDR_thresholds.csv", res)

