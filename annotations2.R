install.packages("tidyverse")
library(tidyverse)

ENSCINP = read.delim("ENSCINP to T.txt", header = T)
KH = read.delim("KH-ENS-conversion2013.txt", header = T)
Names = read.delim("KH2013-UniqueNAME.txt", header = T)
probe = read.delim("probeset-KH.txt", header = T)
probes = read.delim("probes_to_kh.txt", header = F)
Mouse = read.delim("~/Documents/Ciona conversions/Ciona_mouse_orthologs_final.txt", header = T)

colnames(KH)[colnames(KH)=="ENS"] <- "Transcript.stable.ID"
colnames(KH)[colnames(KH)=="KHID"] <- "KH.complete"
colnames(probe)[colnames(probe)=="probeannotation11.kh"] <- "KH.short.short"
colnames(probes)[colnames(probes)=="V1"] <- "probeset"
colnames(probes)[colnames(probes)=="V2"] <- "KH.complete"
probes$V3<- NULL

KH <- KH %>% distinct(Transcript.stable.ID, KH.complete, .keep_all = TRUE)
probes$probeset<-gsub(".*>","",probes$probeset)
probes$KH.short  <- as.vector( sub( "(^[^.]+[.][^.]+[.][^.]+)(.+$)", "\\1", as.character(probes$KH.complete)) )
Names$KH.short.short<-gsub(".*:","",Names$GeneID)
probes$GeneName  <- Names$UniqueNAME[match(probes$KH.short, Names$KH.short.short)]
write.csv(probes, file="annotations_PROBE.csv")

YAY<-merge(probes,KH, all=TRUE)
YAY<- merge(KH,ENSCINP, all=TRUE)
YAY$KH.complete<- NULL

YAY$KH.short  <- as.vector( sub( "(^[^.]+[.][^.]+[.][^.]+)(.+$)", "\\1", as.character(YAY$KH.ID)) )
YAY$GeneName  <- Names$UniqueNAME[match(YAY$KH.short, Names$GeneID)]
YAY$KH.short.short<-gsub(".*:","",YAY$KH.short)

YAY2<- merge(YAY,probe, all=TRUE)
write.table(YAY2, file="ALL_annotations_PROBE.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

probe.short<-probe[!(is.na(probe$KH.short.short) | df$start_pc==""), ]

##KHID TO ENS##
setwd("~/Documents/Ciona conversions")
KH = read.delim("KH-ENS-conversion2013.txt", header = T)
KH_unique <- unique(KH)

write.csv(KH_unique, file="KH-ENS-conversion2013_unique.txt", 
            append = FALSE, sep = "\t", dec = ".",
            row.names = TRUE, col.names = TRUE)

#simplify KH
KH_unique$KH.short<- as.vector(sub("(^[^.]+[.][^.]+[.][^.]+)(.+$)", "\\1", as.character(KH_unique$KHID)) )
TOPTAGS<- read.csv("~/Documents/RNF149 RNASEQ results/toptags_rnfko_KH_UniqieName.csv")
#convert
TOPTAGS$ENS  <- KH_unique$ENS[match(TOPTAGS$X, KH_unique$KH.short)]
write.csv(TOPTAGS, file="TOPTAGS_ENS.csv")

TOPTAGS_inclusive<-merge(TOPTAGS, KH_unique, by.x = "X", by.y = "KH.short")
write.csv(TOPTAGS_inclusive, file="TOPTAGS_ENS_inclusive.csv")

##counts convert

geneCounts <- read.csv("geneCounts_KH.csv", header = T, stringsAsFactors = FALSE)
