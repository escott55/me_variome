#!/usr/bin/env Rscript

#source("http://bioconductor.org/biocLite.R")
#biocLite("quantsmooth")
library(quantsmooth)
library(scales)
if (!require("RColorBrewer")) {
    install.packages("RColorBrewer")
    library(RColorBrewer)
}

ideogram <- function (chrompos, cols = "grey50", bleach = 0.3, topspace = 1, 
                sexChromosomes = FALSE, showaxis=TRUE,point.vals=NULL,
                paintCytobands=FALSE,units = c("bases", "cM", "ISCN"), 
                ylab="Chromosomes", palette="Blues", ...) {
    cytobandWidth <- 0.2
    #par(mar = c(1, 4, 2, 3) + 0.1)
    #units <- "bases"
    units <- match.arg(units)
    chrom.n <- 22
    chrs2 <- factor(numericCHR(chrompos[, "CHR"]), levels = c(1:chrom.n, 
                    if (sexChromosomes) c(98, 99) else NULL))

    lens <- lengthChromosome(levels(chrs2), units = units)

    names(lens) <- characterCHR(names(lens))
    cols <- rep(cols, length.out = length(lens))
    names(cols) <- names(lens)
    dheight <- NULL
    for (i in 1:chrom.n) dheight[i] <- lens[i] + lens[chrom.n + 1 - i]
    maxdheight <- max(dheight) * 1.0001
    plot(  c(1, 1.0 + length(dheight)/2+1), c(0, maxdheight),
         type = "n", ylab="", xlab = "", axes = FALSE, 
         las = 2) # ...)cex.lab=2, 
    title(ylab=ylab,cex.lab=2,mgp=c(1.5,1,0))
    if( showaxis )
        axis(3, c(1:length(dheight))/2+.5, characterCHR(names(cols)), las = 0)

    if( paintCytobands ){
        for (i in 1:length(dheight)) {
            paintCytobands(i, c(.5+i/2 + cytobandWidth/2,lens[i]), "bases", 
                width = cytobandWidth, length.out = lens[i], orientation="v",
                legend = FALSE, bleach = bleach)
        }
    }
    chrompos$fixed <- 0
    for( chr in names(chroms) ){
        chrlen <- as.numeric(chroms[chr])
        pos <- chrompos[chrompos$CHR == as.character(chr),c("POS1")]
        chrompos[chrompos$CHR == as.character(chr),c("fixed")] <- chrlen - pos
    }

    # on top of - cex =5, chrompos = .5
    # next to - cex= anything, chrompos += .72
    if( length(point.vals) > 0 ) {
        if( sum(point.vals <0) > 0 ) {
            labels <- paste(c('>-30',as.character(-3:3*10),'>30'),"%",sep="")
            breaks <- c(-1,-3:-1/10,-.01,.01,1:3/10,1)
            print(labels)
            print(breaks)
            binned.x=cut(point.vals, breaks=breaks, labels=labels)
            hcols <- brewer.pal(length(levels(binned.x)),"RdBu")
            pcolors <- hcols[as.numeric(binned.x)]
        }else {
            labels <- paste(c(as.character(0:5*10),'>50'),"%",sep="")
            breaks <- c(0:6/10,1)
            binned.x=cut(point.vals, breaks=breaks, labels=labels)
            #hcols <- adjustcolor(brewer.pal(length(levels(binned.x)),"OrRd"),.2)
            hcols <- brewer.pal(length(levels(binned.x)),palette) # OrRd
            pcolors <- hcols[as.numeric(binned.x)]
        }
        pointloc <- 
        points(chrompos[,1]/2+.5,chrompos[,"fixed"],pch="-",col=pcolors,cex=4.8)
        legend(length(dheight)/2-.5, maxdheight, labels, col=hcols, text.col="black", lwd=2, 
               merge=TRUE, bg="gray90") #lty = c(2, -1, 1)
    }else{
        points(chrompos[,1]/2+.5,chrompos[,"fixed"],pch="-",col="red",cex=4.8)
    }
    if( ! paintCytobands ){
        for (i in 1:length(dheight)) {
            rect( i/2+.5-cytobandWidth/2, 0, i/2+.5+cytobandWidth/2, lens[i]+1000 )
        }
    }
}


#chrs2 <- factor(numericCHR(chrompos[, "CHR"]), levels = c(1:chrom.n,
                   #if (sexChromosomes) c(98, 99) else NULL))

chroms <- lengthChromosome(1:22,units="bases")

#positions$fixed <- 0
#for( chr in names(chroms) ){
    #chrlen <- as.numeric(chroms[chr])
    #pos <- positions[positions$CHR == as.character(chr),c("POS1")]
    #positions[positions$CHR == as.character(chr),c("fixed")] <- chrlen - pos
#}
#print(head(positions))

#"cM","Score",
header = c("Marker","dbSNP","CHR","POS1","End","GeneticDist","U_MR49.2","A_MR49.3","A_MR49.4","A_MR49.5")
data <- read.table("chipdata.txt",sep="",header=F, na.strings="",skip=1,col.names=header)

tokeep <- with(data, (U_MR49.2 == "R") & (A_MR49.3 == "L") 
               & (A_MR49.4 == "L") & (A_MR49.5 =="L") & CHR != "X")
positions <- data[tokeep,c("CHR","POS1")]
positions$CHR <- as.numeric(as.character(positions$CHR))
#within(data, (U_MR42.2 == "R") & (A_MR49.3 == "L"))

png( "ideograms/chip_ideogram_cyto.png",  width=12, height=5, units="in", res=200 )
ideogram(positions, ylab="MR49 Homozygosity", units="bases",paintCytobands=F)
dev.off()

