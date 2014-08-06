#!/usr/bin/env Rscript

library(RColorBrewer)
#palette(brewer.pal(8, "Dark2"))
#library(gplots)
library(ggplot2)
library(ape)
library(reshape)
library(grid)
library(gridExtra)
library(plyr) 
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
        namelist <- read.table(filterfile,header=FALSE,sep=" ",comment.char="",col.names=c("V1","Individual.ID","V3","V4","V5","V6"))
        allannotdata <- merge( namelist, allannotdata, by="Individual.ID", all.x=TRUE )
        allannotdata <- allannotdata[namelist$Individual.ID,]
        #allannotdata <- allannotdata[with(allannotdata,
        print(head(namelist$Individual.ID))
        print(head(allannotdata$Individual.ID))
        #print(allannotdata[is.null(allannotdata["Family.ID"]),])
    }
    return(allannotdata)
} # End readAnnotation


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
# Main
###################################################################### 
setwd( "/home/escott/workspace/variome/" )

#annotfile <- "./resources/annotation/patientannotation.ped"
#figuredir <- "./rscripts/relatedness"
#targetdir <- "./rawdata/merged/admixture"
#basename <- "admixture/merged.chimp.regions.filt.samp.plink.filt.fam"


#relatedfile <- "/home/escott/workspace/variome/rawdata/merge1kg/qc/me_1000G_imp.clean.relatedness"
relatedfile <- "/media/data/workspace/variome/rawdata/mevariome/qc/variome.clean.relatedness"
unrelatedfile <- "/media/data/workspace/variome/rawdata/mevariome/qc/variome.clean.relatedness"

relness <- read.table( relatedfile, header=TRUE, sep="\t", comment.char="" )
unrelated <- read.table( unrelatedfile, header=F, sep="\t", comment.char="" )

subdata = relness[(relness$INDV1 %in% unrelated$V1) & (relness$INDV2 %in% unrelated$V1),]
                                                      #c("ID1","ID2","IBS0")]
#subdata2 <- subdata
#colnames(subdata2) <- c("ID2","ID1","IBS0")
#subdata2 <- subdata2[,c("ID1","ID2","IBS0")]
#combined <- rbind( subdata, subdata2 )
#head(combined)
#subdata[subdata$RELATEDNESS_AJK > .25,c("RELATEDNESS_AJK")] <- 0.0
casted = cast(subdata, INDV1 ~ INDV2 )
submat <- data.matrix(casted[,2:dim(casted)[2]])
#submat[is.na(submat)] <- 0.0
#submat[submat)] <- 0.0
is.numeric( submat )

pdf("test2.pdf")
heatmap.2(submat, col=redgreen(75), scale="none", key=TRUE, 
                    symkey=FALSE, density.info="none", trace="none", cexRow=0.5, symm=T)
dev.off()


# KING 
relatedfile = "/home/escott/workspace/variome/rawdata/mevariome/king/variome_rel.kin0"
unrelatedfile = "/home/escott/workspace/variome/rawdata/mevariome/king/variome_unrelated.txt"

relness <- read.table( relatedfile, header=TRUE, sep="\t", comment.char="" )
unrelated <- read.table( unrelatedfile, header=F, sep="\t", comment.char="" )


subdata = relness[(relness$ID1 %in% unrelated$V2) 
                  & (relness$ID2 %in% unrelated$V1),
                  c("ID1","ID2","IBS0")]

subdata2 <- subdata
colnames(subdata2) <- c("ID2","ID1","IBS0")
subdata2 <- subdata2[,c("ID1","ID2","IBS0")]
combined <- rbind( subdata, subdata2 )
head(combined)
casted = cast(combined, ID1 ~ ID2 )
submat <- data.matrix(casted[,2:dim(casted)[2]])
submat[is.na(submat)] <- 0.0
is.numeric( submat )
#heatmap.2(submat)

#boxplot(apply(submat,1,mean))

#s.means <- submat.apply(mean)

#pdf("test.pdf")
heatmap.2(submat, col=redgreen(75), scale="none", key=TRUE, 
          symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
#dev.off()
#
#symm=TRUE, 

outs <- boxplot(apply(submat,1,mean))$out

q()
print(head(relness))
#relness$INDV1 <- sub( ".variant2?", "", relness$INDV1 )
#relness$INDV2 <- sub( ".variant2?", "", relness$INDV2 )

#mymat <- cast( relness, INDV1 ~ INDV2, value="RELATEDNESS_AJK" )

#heatmap.2(as.matrix(mymat), col=redgreen(75), key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)

relness_sub <- relness[relness$RELATEDNESS_AJK > 0. & relness$RELATEDNESS_AJK < .5,]
(m <- ggplot(relness_sub, aes(x=RELATEDNESS_AJK) ) +
            geom_histogram() )#+ scale_y_log10() )

(p <- ggplot(relness, aes(INDV1, INDV2)) + 
                geom_tile(aes(fill = RELATEDNESS_AJK),colour = "white") + 
                scale_fill_gradient(low = "white", high = "steelblue"))


tocut <- c(-1, .12,.13, .25, .5, 1.)
cutlabels <- c("unrelated", "error?","first cousin","sibling","identical")
relness["rtype"] <- cut( relness$RELATEDNESS_AJK, tocut, labels=cutlabels)
table(relness$rtype)

rel_sub <- relness[relness$rtype != "unrelated",]

#scale="row", 
# Parse Basename
#args <- commandArgs( trailingOnly = TRUE )
#if( length(args) == 0 ){ printHelp(args) }

#basename <- basename(gsub( "(.*)\\.(.*)", "\\1", ifelse(length(args)>=1,args[1],basename)))
#targetdir <- ifelse(length(args)>=1,dirname(args[1]),targetdir)
#annotfile <- ifelse(length(args)>=2,args[2],annotfile)
#maxK <- ifelse(length(args)>=3,args[3],6)
#print(paste("Max K",maxK))

#print(paste("Basename",basename))
#print(paste("Targetdir",targetdir))
#print(paste("Annotation",annotfile))




