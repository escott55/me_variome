#!/usr/bin/env Rscript


library(RColorBrewer)
palette(brewer.pal(8, "Dark2"))
#library(gplots, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(reshape, quietly=TRUE)

setwd( "/home/escott/workspace/variome/rscripts/" )
source("themes.R")
source("myutils.R")

#start <- grep("^[[:digit:]]+$", x)

testfile <- "/home/escott/workspace/inbreed/rawdata/turks/turks.recode.chimp_stats.txt"
args <- commandArgs( trailingOnly = TRUE )
targetfile <- ifelse( length(args) > 0, args[1], testfile )

finfo <- getBasename( targetfile )
outputdir <- "./refstatsplot"

df <- readMultiTable( targetfile )

#plotAF( df$

# Allele Frequency
af <- t(df$`Allele Frequency`)
colnames(af) <- af[1,]
af <- as.data.frame(af[2:dim(af)[1],])
af$chimp <- as.numeric(as.character(af$chimp))
af$human <- as.numeric(as.character(af$human))
af <- melt(af, id=c("Rtype"))

outimage <- paste(outputdir,"/",finfo[2],"_af.png",sep="")
print(paste( "Writing to file:",outimage ))
png(outimage)
m <- ggplot(af, aes(x=Rtype, y=value ))
m <- m + geom_bar(stat="identity")
m <- m + ylab("Allele Frequency Bins")+ xlab("Counts")
m <- m + labs(title = "Number of variants by Allele Frequency")
m <- m + facet_grid( variable~. )
m <- m + ng1
print(m)
dev.off()

# Missingness
miss <- t(df$`Missingness`)
colnames(miss) <- miss[1,]
miss <- as.data.frame(miss[2:dim(miss)[1],])
miss$chimp <- as.numeric(as.character(miss$chimp))
miss$human <- as.numeric(as.character(miss$human))
miss <- melt(miss, id=c("Rtype"))

outimage <- paste(outputdir,"/",finfo[2],"_missing.png",sep="")
print(paste( "Writing to file:",outimage ))
png(outimage)
m <- ggplot(miss, aes(x=Rtype, y=value ))
m <- m + geom_bar(stat="identity")
m <- m + ylab("Allele Frequency Bins")+ xlab("Counts")
m <- m + labs(title = "Number of variants by Missingness")
m <- m + facet_grid( variable~. )
m <- m + ng1
print(m)
dev.off()

# DBsnp
dbsnp <- t(df$`DBsnp`)
colnames(dbsnp) <- dbsnp[1,]
dbsnp <- as.data.frame(dbsnp[2:dim(dbsnp)[1],])
dbsnp$chimp <- as.numeric(as.character(dbsnp$chimp))
dbsnp$human <- as.numeric(as.character(dbsnp$human))
dbsnp <- melt(dbsnp, id=c("Rtype"))
dbsnp <- dbsnp[dbsnp$Rtype == "Ratio",]
dbsnp$Rtype <- factor(dbsnp$Rtype)

outimage <- paste(outputdir,"/",finfo[2],"_dbsnp.png",sep="")
print(paste( "Writing to file:",outimage ))
png(outimage)
m <- ggplot(dbsnp, aes(x=variable, y=value ))
m <- m + geom_bar(stat="identity")
m <- m + ylab("Ratio (dbsnp/total)")+ xlab("Organism Ref")
m <- m + labs(title = "Ratio of variants overlaping with dbsnp")
#m <- m + facet_grid( variable~. )
m <- m + ng1
print(m)
dev.off()




