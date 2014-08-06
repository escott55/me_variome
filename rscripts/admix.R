#!/usr/bin/env Rscript

library(RColorBrewer)
palette(brewer.pal(8, "Dark2"))
#library(gplots)
library(ggplot2)
library(ape)
library(reshape)
library(grid)
library(gridExtra)
source("/home/escott/workspace/variome/rscripts/themes.R")
#tbl=read.table("ALL.genome.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.vcf.plink_unrel.5.Q")

#####################################################################
# printHelp
#####################################################################
printHelp <- function(args){
    print( "Hello, you have entred the help text")
    print( args )
    print("length")
    print( length(args) )
    print( "./admix.R <basename> <annotfile.ped>" )
    print( "Goodbye")
    #q()
} # End printHelp

######################################################################
# readAnnotation
######################################################################
readAnnotation <- function(file, filterfile=NULL)
{
    print("readAnnotation")
    #print(file)
    allannotdata <- read.table(file, header=TRUE, sep="\t", comment.char="")
    rownames(allannotdata) <- allannotdata[,2]
    normethnicity <- read.table("resources/HGDP_ethnicgroups.txt", header=TRUE, sep="\t", comment.char="")
    #print("done reading...")
    allannotdata <- merge( allannotdata, normethnicity, by.x="ethnicity",by.y="Ethnicity",all.x=TRUE)
    #print(paste("Filterfile",filterfile))
    if( !is.null(filterfile) ){
        #namelist <- read.table(filterfile,header=FALSE,sep=" ")
        #print(head(namelist))
        #print(head(namelist[,2]))
        #matchingnames <-allannotdata$Individual.ID %in% intersect(namelist[,2], allannotdata$Individual.ID)
        #print(head(matchingnames))
        #newallannotdata <- allannotdata[matchingnames,]
        #allannotdata <- newallannotdata[match(namelist[,2],newallannotdata$Individual.ID),]
        namelist <- read.table(filterfile,header=FALSE,sep=" ",comment.char="",col.names=c("V1","Individual.ID","V3","V4","V5","V6"))
        allannotdata <- merge( namelist, allannotdata, by="Individual.ID", all.x=TRUE )
        allannotdata <- allannotdata[namelist$Individual.ID,]
        #allannotdata <- allannotdata[with(allannotdata,
        print(head(namelist$Individual.ID))
        print(head(allannotdata$Individual.ID))
        #print("Missing samples")
        #print(allannotdata[is.null(allannotdata["Family.ID"]),])
    }
    #print(tail(allannotdata))
    return(allannotdata)
} # End readAnnotation

######################################################################
# clusterAdmix
######################################################################
clusterAdmix <- function(mydata,currK,outfile)
{
    print("clusterAdmix")
    ## Ward Hierarchical Clustering
    d <- dist(mydata, method = "euclidean") # distance matrix
    fit <- hclust(d, method="ward") 
    #png( outfile )
    #plot(fit) # display dendogram
    groups <- cutree(fit, k=currK) # cut tree into 5 clusters
    # draw dendogram with red borders around the 5 clusters 
    #rect.hclust(fit, k=currK, border="red")
    #dev.off()
    hc <- hclust(dist(mydata))
    png( outfile )
    plot( as.phylo(hc), type="radial")
    dev.off()
    return(groups)
}

######################################################################
# chooseColors
######################################################################
chooseColors <- function( classlist ){
    print(head(classlist))
    classes <- data.frame(unique(classlist))
    #print(classes)
    #print(rainbow(dim(classes)[1]))
    mypalette <- cbind(classes, rainbow( dim(classes)[1] ))
    colnames(mypalette) <- c("classes","colors")
    print(mypalette)

    classlist <- as.data.frame(classlist)
    colnames(classlist) <- c("classes")
    print("classlist")
    #print(head(classlist))
    eachcolors <- merge(mypalette,classlist,by="classes")
    #print(eachcolors)
    return(eachcolors)
}


######################################################################
# makeHeatmap
######################################################################
makeHeatmap <- function( mydata, outfile, classes=NULL ){
    print("Running makeHeatmap")
    print(paste("Writing file:",outfile))
    hcol<-colorRampPalette(c("Red","Green"))(32)
    if( !is.null(classes) ){
        #colorPalette <- chooseColors( newannot$Population )
        colorPalette <- chooseColors( classes )
        scol <- as.character(colorPalette[,2])
        #print(scol)
        #print(dim(mydata))
        #print(length(scol))
        png(outfile)
        heatmap.2(as.matrix(mydata), col=redgreen(75),  key=TRUE, symkey=TRUE, RowSideColors=scol, density.info="none", trace="none", cexRow=0.5)
        dev.off()
        #scale="row",
    }else{
        png(outfile)
        heatmap.2(as.matrix(mydata), col=redgreen(75), scale="row", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
        dev.off()
    }
    #ColSideColors=scol, 
}

######################################################################
# makeBarplot
######################################################################
makeBarplot <- function( mydata, outfile, continents, populations ){
    print("Running makeBarplot")
    conts <- unique(continents)
    cont_index <- mat.or.vec(length(conts), 1)
    for(i in 1:length(conts)) {
        cont_index[i] <- mean(which(continents == conts[i]));
    }
    #cont_index <- 1360/nrow(mydata)*cont_index -55;

    regs <- unique(populations);
    reg_index <- mat.or.vec(length(regs), 1);
    for(i in 1:length(regs)) {
        reg_index[i] <- mean(which(populations == regs[i]));
    }
    #reg_index <- 500/nrow(mydata)*reg_index - 55;
    #coef <- 1980/nrow(mydata)
    if( any(grep("daily", outfile)) ){
        print( "Daily specific coefficient" )
        coef <- 1460/nrow(mydata)
    }else if( any(grep("everything", outfile)) ){
        print("found everything")
        coef <- 1058/nrow(mydata)
    }else if( any(grep("merged", outfile)) ){
        print("found merged")
        coef <- 2900/nrow(mydata)
    }else if( any(grep("onekg", outfile)) ){
        print("found onekg")
        coef <- 1460/nrow(mydata)
    }else{
        coef <- 710/nrow(mydata)
    }
    print( paste("Coef:",coef) )
    #print( paste("Coef:",coef) )
    #print( paste("atts:",reg_index, length(mydata), mat.or.vec(length(regs),1) ))
    #print(paste("reg_index",reg_index))
    #print(paste("cont_index",cont_index))
    #print( length(mydata[,1]) )
    #newtbl <- mydata
    #newtbl$sample <- rownames(mydata)
    #newtbl <- melt( newtbl, id="sample" ) 
    #m <- ggplot(newtbl, aes(x=sample, y=value, fill=variable ))
    #m <- m + geom_bar(stat="identity")
    #m <- m + theme(legend.position = "none")
    #m <- m + theme(axis.text.x=element_blank())
    #m <- m + scale_x_continuous( breaks=(reg_index-reg_index[1])*coef,labels=regs)
    #m <- m + ng1
    #barwidth = 7/nrow(mydata)
    print( paste("Writing file:",outfile) )
    png( outfile, width = 7, height = 3, units="in", res=400) 
    barplot(t(as.matrix(mydata)), col=rainbow( currK ),
           ylab="Ancestry", border=NA, xaxt="n" )
        #xlab="Individual", main="Population Admxiture"
    #axis(1, at = (reg_index-reg_index[1])*coef, labels = regs, hadj=2.5)
    axis(1, at = (cont_index-cont_index[1])*coef, labels=conts)#, padj = 2.5)#,tick=FALSE)
    #text((cont_index-cont_index[1])*coef, labels=conts, xpd=TRUE, srt=45)
    dev.off()
} # makeBarplot

######################################################################
# Main
######################################################################

setwd( "/home/escott/workspace/variome/" )

annotfile <- "./resources/annotation/patientannotation.ped"
figuredir <- "./rscripts/admixplot"

#targetdir <- "../data/merged"
#basename <- "ciliopathies.unfilt.recode.plink.mod.filt.uniq"

#targetdir <- "./rawdata/daily/admixture"
#basename <- "daily.chimp.regions.filt.samp.plink.mod.filt.uniq.fam"

#targetdir <- "./rawdata/test/admixture"
#basename <- "everything_set1.chr1.snp.chimp.regions.filt.samp.plink.mod.filt.uniq.fam"

#targetdir <- "./rawdata/merged/admixture"
#basename <- "admixture/merged.chimp.regions.filt.samp.plink.filt.fam"

targetdir <- "./rawdata/onekg/admixture"
#basename <- "onekg.chimp.regions.filt.samp.plink.mod.filt.uniq.fam"
basename <- "onekg.chimp.regions.filt.samp.plink.filt.bed"

#targetdir <- "../data/daily1"
#basename <- "daily.recode.plink.mod.filt.uniq"

# Parse Basename
args <- commandArgs( trailingOnly = TRUE )
if( length(args) == 0 ){ printHelp(args) }
#gsub("(.*)\\.(.*)", "\\1\\2", N)
basename <- basename(gsub( "(.*)\\.(.*)", "\\1", ifelse(length(args)>=1,args[1],basename)))
targetdir <- ifelse(length(args)>=1,dirname(args[1]),targetdir)
annotfile <- ifelse(length(args)>=2,args[2],annotfile)
maxK <- ifelse(length(args)>=3,args[3],6)
print(paste("Max K",maxK))


print(paste("Basename",basename))
print(paste("Targetdir",targetdir))
print(paste("Annotation",annotfile))


filterfile <- paste(targetdir,"/",basename,".fam",sep="")
print(paste("Filterfile:",filterfile))
annotdata <- readAnnotation( annotfile, filterfile )

print(paste("Using",dim(annotdata)[1],"Samples"))
print(paste("Unique Samples:",length(unique(annotdata$Individual.ID))))
print(paste("Duplicated Samples:",annotdata[duplicated(annotdata$Individual.ID),]))
print(head(annotdata))

suffix <- ".Q"

for( currK in 4:maxK ){
    currfile <- paste(targetdir,"/",basename,".",currK,suffix,sep="")
    print("read table")
    tbl=read.table( currfile, header=FALSE, sep=" ", comment.char="#" )
    print(paste("Table Dimensions Rows:",dim(tbl)[1],"Cols:",dim(tbl)[2]))

    rownames(tbl) <- annotdata$Individual.ID

    print("Sort Data!")
    #newtbl <- tbl[order(annotdata$Continent, annotdata$Country, annotdata$ethnicity, annotdata$Individual.ID),]
    newtbl <- tbl[order(annotdata$Continent, annotdata$ethnicity),]
    print(paste("New table", dim(newtbl)))
    #newannot <- annotdata[order(annotdata$Continent, annotdata$Country, annotdata$ethnicity, annotdata$Individual.ID),]
    newannot <- annotdata[order(annotdata$Continent, annotdata$ethnicity),]
    print(paste("New Annot", dim(newannot)))
    outfile <- paste(targetdir,"/",basename,"_",currK,".sorted.txt",sep="")
    print(paste("Writing file:",outfile))
    write.table( cbind(newannot,newtbl), outfile, quote=FALSE,row.names=FALSE, col.names=TRUE, sep="\t")

    mydata <- newtbl[!(is.na(newannot$ethnicity)),]
    newannot <- newannot[!(is.na(newannot$ethnicity)),]
    continents <- newannot$Continent
    populations <- newannot$ethnicity

    #m <- m + theme(axis.text.x=element_blank())
    #m <- m + scale_x_continuous( breaks=(reg_index-reg_index[1])*coef,labels=regs)
 
    outfile <- paste(figuredir,"/",basename,currK,".png",sep="")
    makeBarplot( mydata, outfile, newannot$Continent, newannot$ethnicity )

    outfile <- paste(figuredir,"/",basename,".",currK,"_clust.png",sep="")
    groups <- clusterAdmix(mydata, currK, outfile) 

    print("making groupdata")
    groupdata <- cbind(newannot,groups)
    outfile <- paste(targetdir,"/",basename,".",currK,".groups.txt",sep="")
    print(paste("Writing file:",outfile))
    write.table(groupdata,file=outfile, quote=FALSE,row.names=FALSE, col.names=TRUE, sep="\t")

    #outfile <- paste(figuredir,"/",basename,currK,"_hmap.png",sep="")
    #makeHeatmap( newtbl, outfile, newannot$ethnicity ) #newannot$Continent
}

#filterfile <- "all_plates_allchr.plink.mod.filt_unrelated.txt"
#annotfile <- "phase1_samples_integrated_20101123_clean1.ped"
#filterfile <- "ALL.genome.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.plink.mod.filt_unrelated.txt"
#filterfile <- "ALL.genome.phase1.plink.mod.filt.uniq.fam"
#filterfile <- "ALL.genome.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.plink.mod.filt.uniq.fam"

#q()
#basename <- "ALL.genome.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.vcf.plink_unrel."
##colorPalette <- chooseColors( newannot$Population )
##makeHeatmap( newtbl, outfile, newannot$Population )
## K = 5
#tbl=read.table("plateII_snps.vcf.plink_unrel.2.Q")
#tbl.sort <- tbl[with(tbl,order(-V1,-V2)),]
#barplot(t(as.matrix(tbl.sort)), col=rainbow(2),
      #xlab="Individual #", ylab="Ancestry", border=NA)

