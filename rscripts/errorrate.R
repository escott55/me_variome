#!/usr/bin/env Rscript

library(ggplot2)

######################################################################
# Main
######################################################################

setwd( "/home/escott/workspace/variome/rscripts/" )
print(paste("Working dir:",getwd()))

targetdir <- "errorrateplot"
infile <- "logs/merged_snps.plink.mod.filt.uniq_error.log"
# Parse Basename
args <- commandArgs( trailingOnly = TRUE )
if( length(args) == 0 ){ q() }
#gsub("(.*)\\.(.*)", "\\1\\2", N)
infile <- args[1]
basename <- sub( "_error\\.log.*","", basename(ifelse(length(args)>=1,args[1],basename)))
targetdir <- ifelse(length(args)>=2,args[2],targetdir)
print(paste("Basename",basename))
print(paste("Targetdir",targetdir))

data <- read.csv(infile,header=TRUE,sep="\t") #row.names=NULL
print(dim(data))
print(rownames(data))
print(head(data))
colnames(data) <- c("K","cve")
outfile <- paste(targetdir,"/",basename,"_error.png",sep="")
print(outfile)
png(outfile)
#plot(data, main="K Selection")
#abline(data)
p <- ggplot(data, aes(x=K, y=cve ) )
print(p + geom_line() +
    geom_point() +
    ylab("Cross-validation Error") +
    ggtitle("K Selection for Admixture"))
dev.off()

