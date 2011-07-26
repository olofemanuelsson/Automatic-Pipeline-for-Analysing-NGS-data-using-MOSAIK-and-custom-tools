#!/usr/bin/env Rscript

#Plots coverage across chromosomes
#Run with "Rscript Coverage_hist_by_Chr.R [infile.txt] [SampleName1] [xcoverage] [out_directory]" from bash
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

covChr <-read.table(paste(infile), header=T)

rownames( covChr ) <- covChr$cvr

covChr <- covChr[ , -1 ]

colNamecovChr <-names(covChr)

plotcov <- covChr[0:xcov+1,] #c(0:xcov) yields x steps but covChr[0:xcov] only yields x-1

#set working directory
setwd( paste(od, sep="") )

pdf(paste( paste( sampleName1, "Chr_all_perBaseVScov", sep="_"),"pdf", sep=".") )

k <- 0
ktrack <-0

plot(c(0:xcov), seq(0,1,(1/xcov) ) , col="white", xlab="Coverage (x)", ylab="Fraction bases ") 
title(paste("Fraction of Bases Vs Coverage ", "(", sampleName1, ")", sep="") )
legend(x="bottomleft", c("chrX", "chrY", "chrM", "Chr1-22"), cex=0.8, 
col=c("red","blue","black", "black"), pch=c(22,22,22,21) )


for(i in 1: length(covChr)  )  {

ktrack <-0

	if ( colNamecovChr[i] == "chrX" ) {
	
	ktrack <-1
	k <- k+1
	points(c(0:xcov), plotcov[, i], col="red", pch=22)
	}
	if ( colNamecovChr[i] == "chrY" ) {
	
	ktrack <-1
	k <- k+1
	points(c(0:xcov), plotcov[, i], col="blue", pch=22)
	}
	if ( colNamecovChr[i] == "chrM" ) {
	
	ktrack <-1
	k <- k+1
	points(c(0:xcov), plotcov[, i], col="black", pch=22)
	}
	if (ktrack < 1) {
	
	k <- k+1
	points(c(0:xcov), plotcov[, i], col=k)
 	}

}

dev.off()
