#!/usr/bin/env Rscript

#Plots coverage in BED, positions in pileup and genome
#Run with "Rscript Coverage_bed_pileup_genome.R [infile.txt] [SampleName1] [xcoverage] [out_directory]" from bash
#©2011 Henrik Stranneheim

#1st argument infile
#2nd argument sample name 1
#3rd argument x scale for coverage
#4th argument out directory

args <- commandArgs(TRUE) #Collects arguments

infile <-args[1] #Collect infile

sampleName1 <- args[2] #Collects sample name 1

xcov <- as.numeric(args[3]) #How much to plot on the x and how many points on the y axis

od <- args[4] #Collect location of out directory

#Read infile

cov <-read.table(paste(infile), sep="\t", header=T, row.names=NULL)
#cov <-read.table("10-7050_Cov_BED_Genome.txt", header=T)

rownames( cov ) <- cov$cvr

cov <- cov[ , -1 ]

colNamecov<-names(cov)

plotcov <- cov[0:xcov+1,] #c(0:xcov) yields x steps but cov[0:xcov] only yields x-1

#set working directory
setwd( paste(od, sep="") )

pdf(paste( paste( sampleName1, "Cov_BED_pileup_Genome_BaseVScov", sep="_"),"pdf", sep=".") )

plot(c(0:xcov), seq(0,1,(1/xcov)) , col="white", xlab="Coverage (x)", ylab="Fraction bases ") 
title(paste("Fraction of Bases Vs Coverage ", "(", sampleName1, ")", sep="") )

points(c(0:xcov), plotcov[, 2], col=1, pch=21)
points(c(0:xcov), plotcov[, 4], col=2, pch=22)
points(c(0:xcov), plotcov[, 5], col=3, pch=23)
legend(x="topright", c( paste(colNamecov[2]), paste(colNamecov[4]), paste(colNamecov[5]) ), cex=0.8, 
col=1:3, pch=21:23 )

dev.off()