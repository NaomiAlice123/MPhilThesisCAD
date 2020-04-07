remove(list = ls())

library(stringr)
library(tidyverse)

#import data 
CAD_Loci <- read.csv("./1.Loci /Data/CAD_Loci.csv")[, 1:9]
  CAD_Loci$rsID <- as.character(CAD_Loci$rsID)
  CAD_Loci$rsID <- str_replace_all(CAD_Loci$rsID, pattern=" ", repl="")
  CAD_Loci$Chromosome <- as.numeric(CAD_Loci$Chromosome)
  CAD_Loci$Position_B37 <- as.numeric((CAD_Loci$Position_B37))
  
  attach(CAD_Loci)
  CAD_Loci <- CAD_Loci[order(Chromosome, Position_B37), ]
  detach(CAD_Loci)

  CAD_Loci$rsID_Upper <- CAD_Loci$Position_B37 + 250000
  CAD_Loci$rsID_Lower <- CAD_Loci$Position_B37 - 250000

  Loci <- data.frame(matrix(ncol=44, nrow=500))
  x <- c("Loci", "Chromosome", "Start", "End", str_c("SNP", c(1:40)))
  colnames(Loci) <- x

  CAD_Loci_remove <- CAD_Loci
  
  i=1
  while (nrow(CAD_Loci_remove)>0){
    nrow=nrow(CAD_Loci_remove)
    for (j in 1:nrow){
      CAD_Loci_remove$Count[j] <- length(which(CAD_Loci_remove$Chromosome==CAD_Loci_remove$Chromosome[j] & CAD_Loci_remove$Position_B37>CAD_Loci_remove$rsID_Lower[j] & CAD_Loci_remove$Position_B37<CAD_Loci_remove$rsID_Upper[j]))
    }
    max <-which.max(CAD_Loci_remove$Count)
    Loci$Chromosome[i] <- CAD_Loci_remove$Chromosome[max]
    Loci$Start[i] <- CAD_Loci_remove$rsID_Lower[max]
    Loci$End[i] <- CAD_Loci_remove$rsID_Upper[max]
    Take <- which(CAD_Loci_remove$Chromosome==CAD_Loci_remove$Chromosome[max] & CAD_Loci_remove$Position_B37>CAD_Loci_remove$rsID_Lower[max] & CAD_Loci_remove$Position_B37<CAD_Loci_remove$rsID_Upper[max])
    k <- length(Take)
    Loci[i, 5:(k+4)] <- CAD_Loci_remove[Take, 1]
    Loci$Loci[i]=i
    i=i+1
    CAD_Loci_remove <- CAD_Loci_remove[-Take, ]
  }
  remove(CAD_Loci_remove)
  
  Loci <- Loci[which(rowSums(is.na(Loci)) != 44),]
  
  # check no duplicates across multiple columns:
  sum(apply(Loci[, 5:44], 2, function(x) length(which(!is.na(x))))) #508 :) 
  
  # Are there any overlapping loci? (need to merge 72 to 36)
  Loci <- Loci[order(Loci$Chromosome, Loci$Start), ]
  
  Merge <- data.frame(Loci1=c(), Loci2=c())
  for (i in 1:22){
    Data <- Loci[which(Loci$Chromosome==i),]
    Data$Okay <- ifelse(Data$Start < Data$End[c(NA,1:(nrow(Data)-1))], 1, 0)
    Loci1 <- Data[which(Data$Okay==1),1]
    Loci2 <- Data[(which(Data$Okay==1))-1, 1]
    Merge2 <- cbind(Loci1, Loci2)
    Merge <- rbind(Merge, Merge2)
  }
  
  # Combine overlapping loci:
  for (i in 1:nrow(Merge)){
    Overlap <- subset(Loci, Loci == Merge$Loci1[i] | Loci == Merge$Loci2[i])
    newEnd <- Overlap[which.max(Overlap$Start),4 ]
    rsIDs <- Overlap[which.max(Overlap$Start),5:44 ]
    rsIDs <- rsIDs[!is.na(rsIDs)]
    Loci[Loci$Loci==Overlap[which.min(Overlap$Start), 1 ], 4] <- newEnd
    Options <- which(is.na(Loci[Loci$Loci==Overlap[which.min(Overlap$Start), 1 ], ]))
    Fill <- Options[1:length(rsIDs)]
    Loci[Loci$Loci==Overlap[which.min(Overlap$Start), 1 ], Fill] <- rsIDs
    Loci[Loci$Loci==Overlap[which.max(Overlap$Start), 1 ], 1:44]  <-0
  }
  
  # delete rows with zeros:
  Loci <- Loci[-which(Loci$End==0),]
  
  # check no overlapping now (same code as above):
  Loci <- Loci[order(Loci$Chromosome, Loci$Start), ]
  
  Merge <- data.frame(Loci1=c(), Loci2=c())
  for (i in 1:22){
    Data <- Loci[which(Loci$Chromosome==i),]
    Data$Okay <- ifelse(Data$Start < Data$End[c(NA,1:(nrow(Data)-1))], 1, 0)
    Loci1 <- Data[which(Data$Okay==1),1]
    Loci2 <- Data[(which(Data$Okay==1))-1, 1]
    Merge2 <- cbind(Loci1, Loci2)
    Merge <- rbind(Merge, Merge2)
  }
  # Merge is empty 
  
  #re-number Loci numbers:
  Loci$Loci <- seq(1:176)
  
  #remove all NA columns:
  Loci <- Loci[, !apply(is.na(Loci), 2, all)]
  
  write.csv(Loci, "./1.Loci /Data/Loci.csv")
  
  
  
 
