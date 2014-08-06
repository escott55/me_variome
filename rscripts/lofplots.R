#!/usr/bin/env Rscript

#library(gplots)
library(ggplot2)
library(reshape)
library(grid)
library(gridExtra)

setwd("..")
source("rscripts/themes.R")

simpleCap <- function(x) {
    s <- strsplit(x, "[. ]")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
}

#args <- commandArgs( trailingOnly = TRUE )
#prefix <- ifelse(length(args)>=1,args[1],"default")
##prefix <- "merged.chimp.regions.filt.samp.samp_plink"
#ibcdatafile <- paste("results/",prefix,"_continent_ibcrates.txt",sep="")
##countryfile <- paste("results/",prefix,"_country_ibcrates.txt",sep="")
#
#print(paste("Prefix:",prefix))
#print(paste("IBC data file:",ibcdatafile))
##print(paste("Country file:",countryfile))
#
#basename <- basename(gsub( "(.*)\\.(.*)", "\\1", ibcdatafile))


prefix <- "merged"
countsfile <- "./results/merged.clean_regioncounts.txt"
regioncounts <- read.table(countsfile,header=TRUE,sep="\t")

regioncounts <- regioncounts[regioncounts$region != "Oceania",]

outimage <- paste("results/figures/",prefix,"_violin_lof.png",sep="")
print(outimage)
png(outimage)
m2 <- ggplot(regioncounts, aes(x=region, y=total))
m2 <- m2 + geom_violin()
m2 <- m2 + theme(axis.text.x=element_text(angle = 45, hjust = 1))
m2 <- m2 + ng1
m2 <- m2 + labs(title="Loss of Function counts")
m2 <- m2 + xlab("Region")+ ylab("Total")
print(m2)
dev.off()

outimage <- paste("results/figures/",prefix,"_violin_lof_het.png",sep="")
print(outimage)
png(outimage)
m2 <- ggplot(regioncounts, aes(x=region, y=het))
m2 <- m2 + geom_violin()
m2 <- m2 + theme(axis.text.x=element_text(angle = 45, hjust = 1))
m2 <- m2 + ng1
m2 <- m2 + labs(title="Het Loss of Function counts")
m2 <- m2 + xlab("Region")+ ylab("Hets")
print(m2)
dev.off()

outimage <- paste("results/figures/",prefix,"_violin_lof_hom.png",sep="")
print(outimage)
png(outimage)
m2 <- ggplot(regioncounts, aes(x=region, y=hom))
m2 <- m2 + geom_violin()
m2 <- m2 + theme(axis.text.x=element_text(angle = 45, hjust = 1))
m2 <- m2 + ng1
m2 <- m2 + labs(title="Hom Loss of Function counts")
m2 <- m2 + xlab("Region")+ ylab("Homs")
print(m2)
dev.off()


