#!/usr/bin/env Rscript

#Plots average coverage in across all chromosomes
#Run with "Rscript Average_cov_chr.R [infile.txt] [SampleName1] [out_directory]" from bash
#©2011 Henrik Stranneheim

#1st argument infile
#2nd argument sample name 1
#3rd argument out directory

args <- commandArgs(TRUE) #Collects arguments

infile <-args[1] #Collect infile

sampleName1 <- args[2] #Collects sample name 1

od <- args[3] #Collect location of out directory

#Read infile

cov <-read.table(paste(infile), sep="\t", header=T, row.names=NULL)
#cov <-read.table(paste(infile), header=T)

#set working directory
setwd( paste(od, sep="") )

sum1<-summary(cov[1:24, 5])

pdf(paste( paste( sampleName1, "Average_Cov_per_chr", sep="_"),"pdf", sep=".") )

barplot(cov[1:24, 5], xlab="Chromosomes", ylab="Average Coverage (x)", 
main=paste("Average Coverage per Chromosome ", "(", sampleName1, ")", sep=""), names.arg=as.character( cov[1:24,1] ), cex.names=0.8 ) 
#title(paste("Percentage of Bases Vs Coverage ", "(", sampleName1, ")", sep="") )

dev.off()