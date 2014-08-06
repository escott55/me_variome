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
         type = "n", ylab = "", cex.lab=2, xlab = "", axes = FALSE, 
         las = 2) # ...)
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
    for( chr in names(lens) ){
        chrlen <- as.numeric(lens[chr])
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

# make a single ideogram
singleIdeogramPlot <- function( bedfile, popsize, prefix ){
    data <- read.table(bedfile,sep="",header=F,
                       na.strings="",stringsAsFactors=F,
                       col.names=c("CHR","POS1","POS2","count"))

    data$percent <- data$count / popsize

    positions <- data[,c("CHR","POS1","count")]
    positions <- positions[positions$CHR != 23,]

    plotfile <- paste("ideograms/",prefix,"_ideogram.png",sep="")
    print( paste("Writing file:",plotfile))
    png( plotfile,  width=12, height=5, units="in", res=200 )
    ideogram(positions, ylab=prefix, units="bases",point.vals=data$percent)
    dev.off()
}

differenceIdeogramPlot <- function( bedfile1,bedfile2,popsize1,popsize2,prefix ){
    print( "Running differenceIdeogramPlot" )
    print(bedfile1)
    print(bedfile2)
    data1 <- read.table(bedfile1,sep="",header=F,
                       na.strings="",stringsAsFactors=F,
                       col.names=c("CHR","POS1","POS2","count"))
    data1$percent <- data1$count / popsize1

    data2 <- read.table(bedfile1,sep="",header=F,
                       na.strings="",stringsAsFactors=F,
                       col.names=c("CHR","POS1","POS2","count"))
    data2$percent <- data2$count / popsize2

    both <- merge(data1,data2, by=c("CHR","POS1"), all = TRUE)
    both <- both[both$CHR != 23,]
    both[is.na(both)] <- 0
    both$newpercent <- both$percent.x - both$percent.y
    #print( head(both) )
    #print( head(both[order(-both$newpercent),]) )

    plotfile <- paste("ideograms/",prefix,"_ideogram.png",sep="")
    print( paste("Writing file:",plotfile))
    png( plotfile,  width=12, height=5, units="in", res=200 )
    ideogram(both, ylab=prefix, units="bases",point.vals=both$newpercent)
    dev.off()
}

runTests <- function (){
    singleIdeogramPlot( "../results/rohcov/testing.bed",14,"test")

    singleIdeogramPlot( "../results/rohcov/test_daly.bed",1488,"daly")

    singleIdeogramPlot( "../results/rohcov/test_variome.bed",2008,"variome")
    testbed1 <- "../results/rohcov/test_daly.bed"
    testbed2 <- "../results/rohcov/test_variome.bed"
    differenceIdeogramPlot( testbed1,testbed2,1488,2008,"CEU-ME" )
}


setwd("~/workspace/variome/rscripts")
# Parse Basename
args <- commandArgs( trailingOnly = TRUE )
#basename <- sub( "\\..*","", ifelse(length(args)>=1,args[1],basename))
#basename <- sub( "_rel\\.kin.*","", basename(ifelse(length(args)>=1,args[1],basename)))
#targetdir <- dirname(ifelse(length(args)>=1,args[1],targetdir))
#targetdir <- ifelse(length(args)>=1,dirname(args[1]),targetdir)
#print(paste("Basename",basename))
#print(paste("Targetdir",targetdir))

targetcommand <- args[1]
if( targetcommand == "single" ){
    bedfile <- args[2]
    popsize <- as.numeric(args[3])
    prefix  <- args[4]
    singleIdeogramPlot( bedfile, popsize, prefix )
}else if(targetcommand == "difference"){
    bedfile1 <- args[2]
    bedfile2 <- args[3]
    popsize1 <- as.numeric(args[4])
    popsize2 <- as.numeric(args[5])
    prefix   <- args[6]
    differenceIdeogramPlot(bedfile1,bedfile2,popsize1,popsize2,prefix)
}






















#data <- read.table("../rawdata/onekg/old2/roh/onekg.clean.plink.filt.hom", sep = "" , header=T,
                                                         #na.strings ="", stringsAsFactors= F)
#cat variome_JGG-157*exomeroh.bed | awk '$1 ~ /[0-9]/ {print $1,$2,$3}' | sort -V -k1,1 -k2,2n | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$1}' > testing.bed
#nfiles <- 1488
#cat daily_*exomeroh.bed | awk '$1 ~ /[0-9]/ {print $1,$2,$3}' | sort -V -k1,1 -k2,2n | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$1}' > test_daly.bed
#nfiles <- 2008
#cat variome_*exomeroh.bed | awk '$1 ~ /[0-9]/ {print $1,$2,$3}' | sort -V -k1,1 -k2,2n | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$1}' > test_variome.bed



# construct genomic positions
#CHR<-sample(19,40,replace=TRUE)  # Chromosomes
#MapInfo<-lengthChromosome(CHR,"bases")*runif(length(CHR)) # position on chromosome
#chrompos <- data.frame(CHR,MapInfo)
#chrompos<-prepareGenomePlot(data.frame(CHR,MapInfo),paintCytobands = TRUE, organism="hsa")

#par(mfrow=c(2,1))
#ideogram(data.frame(CHR,MapInfo), ylab="Europe", units="bases")
#ideogram(data.frame(CHR,MapInfo), ylab="Asia", units="bases", showaxis=FALSE)

# Chrompos returns a matrix with the positions of the elements on the plot
# You can use all kinds of base graphics functions to annotate the chromosomes
#points(chrompos[,2],chrompos[,1]+0.1,pch="x",col="red")
# Show connection between 3rd and 4th element
#segments(chrompos[3,2],chrompos[3,1],chrompos[4,2],chrompos[4,1],col="blue",lwd=2)
#cols="grey50"
#sexChromosomes=FALSE
#units="bases"





    #dchrompos <- matrix(0, nrow = nrow(chrompos), ncol = 1, 
        #dimnames = list(rownames(chrompos), c("CHR", "MapInfo")))
    #for (i in 1:length(dheight)) if (dheight[i] != "") {
        ##probes <- characterCHR(chrompos[, "CHR"]) == as.character(i)
        #dchrompos[, 2] <- chrompos[probes, "MapInfo"] + maxdwidth - lens[dheight[i]]
        #dchrompos[, 1] <- i
    #}

    #for (i in 1:length(leftrow)) {
        #probes <- characterCHR(chrompos[, "CHR"]) == leftrow[i]
        #dchrompos[probes, 2] <- chrompos[probes, "MapInfo"]
        #dchrompos[probes, 1] <- i
    #}

    #leftrow <- c(if (sexChromosomes) "X" else NULL, ((chrom.n + 1)%/%2):1)
    #rightrow <- c(if (sexChromosomes) "Y" else NULL, if (chrom.n%%2 == 
        #1) "" else NULL, ((chrom.n + 1)%/%2 + 1):chrom.n)
    #plot(c(0.5, 0.5 + length(dheight)), c(0, maxdheight), 
        #type = "n", ylab = "Chromosome", xlab = "", axes = FALSE, 
        #las = 2) # ...)
    #axis(1, c(1:length(dheight)), characterCHR(names(cols)), las = 0)
    #axis(4, c(1:length(dheight)), characterCHR(names(cols)), las = 2)
    #for (i in 1:length(dheight)) {
        #paintCytobands(i, c(i + cytobandWidth/2,0),orientation="v", 
            #"bases", width = cytobandWidth, length.out = lens[i], 
            #legend = FALSE, bleach = bleach)
    #}

    # x <- chromosomes, y <- height
   
    #for (i in 1:(chrom.n%/%2)) dwidth[i] <- lens[i] + lens[chrom.n + 1 - i]
    #if (chrom.n%%2 == 1) 
        #dwidth <- c(dwidth, lens[chrom.n%/%2 + 1])
    #if (sexChromosomes) 
        #dwidth <- c(dwidth, lens["X"] + lens["Y"])
        #if (lens[leftrow[i]] > 0) 
        #if (rightrow[i] != "" && lens[rightrow[i]] > 
          #0) 
          #paintCytobands(rightrow[i], c(maxdwidth - lens[rightrow[i]], 
            #i + cytobandWidth/2), "bases", width = cytobandWidth, 
            #length.out = lens[rightrow[i]], legend = FALSE, 
            #bleach = bleach)

#ideogram <- function (chrompos, cols = "grey50", paintCytobands = FALSE, 
    #bleach = 0, topspace = 1, organism, sexChromosomes = FALSE, 
    #units = c("bases", "cM", "ISCN"), ...) {
#
    #units <- match.arg(units)
    #organism <- match.arg(organism, c("hsa", "mmu", "rno"))
    #chrom.n <- switch(organism, hsa = 22, mmu = 19, rno = 20)
    #chrs2 <- factor(numericCHR(chrompos[, "CHR"]), levels = c(1:chrom.n, 
        #if (sexChromosomes) c(98, 99) else NULL))
    #if (organism == "hsa") 
        #lens <- lengthChromosome(levels(chrs2), units = units)
    #else lens <- sapply(split(chrompos[, "MapInfo"], chrs2), 
        #function(x) max(c(0, x)))
    #names(lens) <- characterCHR(names(lens))
    #cols <- rep(cols, length.out = length(lens))
    #names(cols) <- names(lens)
    #dwidth <- NULL
    #for (i in 1:(chrom.n%/%2)) dwidth[i] <- lens[i] + lens[chrom.n + 
        #1 - i]
    #if (chrom.n%%2 == 1) 
        #dwidth <- c(dwidth, lens[chrom.n%/%2 + 1])
    #if (sexChromosomes) 
        #dwidth <- c(dwidth, lens["X"] + lens["Y"])
    #maxdwidth <- max(dwidth) * 1.05
    #leftrow <- c(if (sexChromosomes) "X" else NULL, ((chrom.n + 
        #1)%/%2):1)
    #rightrow <- c(if (sexChromosomes) "Y" else NULL, if (chrom.n%%2 == 
        #1) "" else NULL, ((chrom.n + 1)%/%2 + 1):chrom.n)
    #plot(c(0, maxdwidth), c(0.5, 0.5 + length(dwidth) + topspace), 
        #type = "n", ylab = "Chromosome", xlab = "", axes = FALSE, 
        #las = 2, ...)
    #axis(2, c(1:length(dwidth)), characterCHR(leftrow), las = 2)
    #axis(4, c(1:length(dwidth)), characterCHR(rightrow), 
        #las = 2)
    #if (paintCytobands && organism == "hsa") {
        #for (i in 1:length(dwidth)) {
            #if (lens[leftrow[i]] > 0) 
              #paintCytobands(leftrow[i], c(0, i + cytobandWidth/2), 
                #"bases", width = cytobandWidth, length.out = lens[leftrow[i]], 
                #legend = FALSE, bleach = bleach)
            #if (rightrow[i] != "" && lens[rightrow[i]] > 
              #0) 
              #paintCytobands(rightrow[i], c(maxdwidth - lens[rightrow[i]], 
                #i + cytobandWidth/2), "bases", width = cytobandWidth, 
                #length.out = lens[rightrow[i]], legend = FALSE, 
                #bleach = bleach)
        #}
    #}
#}
#
