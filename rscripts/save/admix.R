#!/usr/bin/env Rscript

library(RColorBrewer, quietly=TRUE)
palette(brewer.pal(8, "Dark2"))
library(gplots, quietly=TRUE)
library(ggplot2, quietly=TRUE)
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
    #print(paste("readAnnotation",file))
    allannotdata <- read.table(file, header=TRUE, sep="\t",comment.char="")
    #print(head(allannotdata))
    #print(allannotdata[allannotdata["Individual.ID"] == "<NA>",])
    rownames(allannotdata) <- allannotdata[,2]
    #print(paste("Filterfile",filterfile))
    if( !is.null(filterfile) ){
        namelist <- read.table(filterfile,header=FALSE,sep=" ",comment.char="",col.names=c("V1","Individual.ID","V3","V4","V5","V6"))
        allannotdata <- merge( namelist, allannotdata, by="Individual.ID", all.x=TRUE )
        #print(tail(allannotdata[,colnames(allannotdata) %in% c("Family.ID","Individual.ID","Paternal.ID","Maternal.ID","Gender","Phenotype")]))
        print("Missing samples")
        print(allannotdata[is.null(allannotdata["Family.ID"]),])
    }

    return(allannotdata)
} # End readAnnotation
#matchingnames <-allannotdata$Individual.ID %in% intersect(namelist[,2],allannotdata$Individual.ID)
#newallannotdata <- allannotdata[matchingnames,]
#allannotdata <- newallannotdata[match(namelist[,2],newallannotdata$Individual.ID),]

######################################################################
# clusterAdmix
######################################################################
clusterAdmix <- function(mydata,currK,outfile)
{
    # Ward Hierarchical Clustering
    d <- dist(mydata, method = "euclidean") # distance matrix
    fit <- hclust(d, method="ward") 
    png( outfile )
    plot(fit) # display dendogram
    groups <- cutree(fit, k=currK) # cut tree into K clusters
    print(head(groups))
    print(paste("K:",currK))
    # draw dendogram with red borders around the K clusters 
    rect.hclust(fit, k=as.numeric(currK), border="red")
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
    hcol<-colorRampPalette(c("Red","Green"))(32)
    if( !is.null(classes) ){
        #colorPalette <- chooseColors( newannot$Population )
        colorPalette <- chooseColors( classes )
        scol <- as.character(colorPalette[,2])
        #print(scol)
        #print(dim(mydata))
        #print(length(scol))
        png(outfile)
        heatmap.2(as.matrix(mydata), col=redgreen(75), scale="row", key=TRUE, symkey=FALSE, RowSideColors=scol, density.info="none", trace="none", cexRow=0.5)
        dev.off()
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
    print(outfile)
    conts <- unique(continents);
    cont_index <- mat.or.vec(length(conts), 1);
    for(i in 1:length(conts))
    {
            cont_index[i] <- mean(which(continents == conts[i]));
    }
    #cont_index <- 1360/length(mydata[,1])*cont_index -55;

    regs <- unique(populations);
    reg_index <- mat.or.vec(length(regs), 1);
    for(i in 1:length(regs))
    {
            reg_index[i] <- mean(which(populations == regs[i]));
    }
    #reg_index <- 1360/length(mydata[,1])*reg_index - 55;
    coef <- 1980/length(mydata[,1])
    print(paste("reg_index",reg_index))
    print(paste("cont_index",cont_index))
    print( length(mydata[,1]) )
    png( outfile, width = 11, height = 5, units="in", res=1200) 
    barplot(t(as.matrix(newtbl)), col=rainbow( currK ),
           ylab="Ancestry", border=NA, xaxt = "n", main="Population Admixture")
        #xlab="Individual",
    axis(1, at = (reg_index-1)*coef, labels = regs)
    axis(1, at = (cont_index-1)*coef, labels = conts, padj = 2.5,tick=FALSE)
    dev.off()
} # END makeBarplot
    #axis(1, at = reg_index, labels = regs)
    #axis(1, at = cont_index, labels = conts, padj = 2.5)

######################################################################
# Main
######################################################################

setwd( "/home/escott/workspace/variome/" )
print(paste("Working dir:",getwd()))

figuredir <- "./rscripts/admixplot"
targetdir <- "./rawdata/daily"
#annotfile <- paste(targetdir,"/","inhouse_samples_annotation.ped",sep="")
annotfile <- "./resources/annotation/patientannotation.ped"
#basename <- "ALL.genome.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.plink.mod.filt.uniq."
basename <- "daily.regions.filt.samp.plink.mod.filt.uniq."
currK <- 10

# Parse Basename
args <- commandArgs( trailingOnly = TRUE )
if( length(args) == 0 ){ printHelp(args) }
#gsub("(.*)\\.(.*)", "\\1\\2", N)
basename <- basename(gsub( "(.*)\\.(.*)", "\\1", ifelse(length(args)>=1,args[1],basename)))
targetdir <- ifelse(length(args)>=1,dirname(args[1]),targetdir)
annotfile <- ifelse(length(args)>=2,args[2],annotfile)
currK <- ifelse(length(args)>=3,args[3],currK)

print(paste("Basename",basename))
print(paste("Targetdir",targetdir))
print(paste("Annotation",annotfile))
print(paste("Curr K",currK))


filterfile <- paste(targetdir,"/",basename,".fam",sep="")
annotdata <- readAnnotation( annotfile, filterfile )
print(head(annotdata))
suffix <- ".Q"

currfile <- paste(targetdir,"/",basename,".",currK,suffix,sep="")
tbl=read.table( currfile ,comment.char="")
#print(head(tbl))

print( "IDs to rownames" )
#annotdata$Individual.ID <- factor( annotdata$Individual.ID, levels=sort(unique(annotdata$Individual.ID)))
#s <- data.frame(table(annotdata$Individual.ID))
#print( head(s[s$Freq > 1,]) )
print( dim(tbl) )
print( dim(annotdata) )
rownames(tbl) <- annotdata$Individual.ID
#print(head(tbl))
print("Sort Data!")
newtbl <- tbl[order(annotdata$continent, annotdata$ethnicity, annotdata$Individual.ID),]
newannot <- annotdata[order(annotdata$continent, annotdata$ethnicity, annotdata$Individual.ID),]
#print(head(newtbl))
#write.table( cbind(newtbl,newannot), "sorteddata.txt", quote=FALSE,row.names=FALSE, col.names=FALSE, sep="\t")
#print(head(annotdata))
#print(head(newannot))

outfile <- paste(figuredir,"/",basename,currK,".png",sep="")
print(paste("Making image:",outfile))
makeBarplot( newtbl, outfile, newannot$continent, newannot$ethnicity )

outfile <- paste(figuredir,"/",basename,".",currK,"_clust.png",sep="")
print(paste("Making image:",outfile))
groups <- clusterAdmix(newtbl, currK, outfile) 

groupdata <- cbind(newannot,groups)
outfile <- paste(targetdir,"/",basename,".",currK,".groups.txt",sep="")
write.table(groupdata,file=outfile, quote=FALSE,row.names=FALSE, col.names=TRUE, sep="\t")

outfile <- paste(figuredir,"/",basename,currK,"_hmap.png",sep="")
print(paste("Making image:",outfile))
makeHeatmap( newtbl, outfile, newannot$ethnicity ) #newannot$continent
#for( currK in 2:maxK ){
#}


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

