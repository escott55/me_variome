#!/usr/bin/env Rscript

library(RColorBrewer, quietly=TRUE)
palette(brewer.pal(8, "Dark2"))
#library(gplots, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(reshape, quietly=TRUE)
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
# Main
######################################################################
#annotfile <- "patientannotation.tsv"

setwd( "/home/escott/workspace/variome/" )
source("rscripts/themes.R")

print(paste("Working dir:",getwd()))

annotfile <- "resources/annotation/patientannotation.ped"
outdir <- "rscripts/pseqfigs"
#datadir <- "rscripts/istatsfiles"
#istatsfile <- paste(datadir,"/ciliopathies.unfilt_istats.txt",sep="")
#istatsfile <- paste(datadir,"/eichler_istats.txt",sep="")
#istatsfile <- "/home/escott/workspace/inbreed/rawdata/onekg/ALL.genome.phase1.names.chimp.regions.recode_istats.txt"
#istatsfile <- "/home/escott/workspace/variome/rawdata/test/everything_set1.chr1.snp_istats.txt"
#istatsfile <- "/home/escott/workspace/variome/rawdata/variome/variome_snps_istats.txt"
istatsfile <- "/home/escott/workspace/variome/rawdata/merged/merged_istats.txt"
istatsfile <- "/home/escott/workspace/variome/rawdata/onekg/onekg.chimp.regions.filt.samp.samp_istats.txt"
#istatsfile <- "/home/escott/workspace/variome/rawdata/onekg/"

# Arguments
args <- commandArgs( trailingOnly = TRUE )
#if( length(args) == 0 ){ printHelp(args) }
istatsfile <- ifelse(length(args)>=1,args[1],istatsfile)
outdir <- ifelse(length(args)>=2,args[2],outdir)
annotfile <- ifelse(length(args)>=3,args[3],annotfile)

print(paste("Istats File",istatsfile))
print(paste("Annotation",annotfile))
print(paste("Outdir",outdir))

# Parse Basename
basename <- basename(gsub( "(.*)\\.(.*)", "\\1", istatsfile))
filepath <- dirname( istatsfile )
print(paste("Basename:",basename))

# Read Data
istats <- read.table( istatsfile, header=TRUE, sep="\t",comment.char="")
rownames(istats) <- istats[,1]
#istats <- istats[order(istats$TITV),]

lvls <- istats$ID[order(istats$TITV)]
istats$ID <- factor( istats$ID, levels=lvls)

# Histogram of Transiton Translation Ratio
outimage <- paste(outdir,"/",basename,"_bar.png",sep="")
print(paste("Writing:",outimage))
png(outimage)
m <- ggplot(istats,aes(x=ID, y=TITV)) + geom_bar(stat="identity")
m <- m + ylab("transition/transversion ratio")+ xlab("Individual")
m <- m + labs(title = "Transition Transversion Ratio per Sample")
print(m)
dev.off()

# Histogram of Transiton Translation Ratio
outimage <- paste(outdir,"/",basename,"_hist.png",sep="")
print(paste("Writing:",outimage))
png(outimage)
m <- ggplot(istats,aes(x=TITV)) + geom_histogram(colour = "darkgreen",fill = "white")
m <- m + ylab("transition/transversion ratio")+ xlab("Individual")
m <- m + labs(title = "Transition Transversion Ratio per Sample")
print(m)
dev.off()

# Get Patient Annotation
#patannot <- read.table(annotfile,header=TRUE,sep="\t",comment.char="")
fullannotation <- read.table(annotfile,header=TRUE,sep="\t",comment.char="")
patannot <- fullannotation[c("Individual.ID","Other.info","continent","ethnicity")]
colnames(patannot) <- c("Individual.ID","Other.info","continent","ethnicity")
istats.merged <- merge( istats, patannot, by.x="ID", by.y="Individual.ID" )

if( dim(istats.merged)[1] == 0 ){
    patannot <- fullannotation[c("Normalized","Other.info","continent","ethnicity")]
    colnames(patannot) <- c("Individual.ID","Other.info","continent","ethnicity")
    istats.merged <- merge( istats, patannot, by.x="ID", by.y="Individual.ID" )
}

istats.merged$continent <- factor( istats.merged$continent, levels=sort(unique(istats.merged$continent) ))

# Transition / Translation Ratio per Ethnicity
outimage <- paste(outdir,"/",basename,"_grpdensity.png",sep="")
print(paste("Writing:",outimage))
print( head(istats.merged$ethnicity) )
png(outimage)
m <- ggplot(istats.merged, aes(TITV, fill = ethnicity))
m <- m + geom_density(alpha = 0.2)
m <- m + ylab("transition/transversion ratio")+ xlab("Individual")
m <- m + labs(title = "Transition Transversion Ratio per Ethnicity")
print(m)
dev.off()

# Transition / Translation Ratio per Plate
outimage <- paste(outdir,"/",basename,"_srcdensity.png",sep="")
print(paste("Writing:",outimage))
png(outimage)
m <- ggplot(istats.merged, aes(TITV, fill = continent))
m <- m + geom_density(alpha = 0.2)
m <- m + ylab("transition/transversion ratio")+ xlab("Individual")
m <- m + labs(title = "Transition Transversion Ratio per Plate")
print(m)
dev.off()


# Number of Singletons
outimage <- paste(outdir,"/",basename,"_singdensity.png",sep="")
print(paste("Writing:",outimage))
png(outimage)
m <- ggplot(istats.merged, aes(SING, fill = ethnicity))
m <- m + geom_density(alpha = 0.2)
m <- m + ylab("Number of Singletons")+ xlab("Individual")
m <- m + labs(title = "Transition Transvertion Ratio per Ethnicity")
print(m)
dev.off()

outliers <- list()
istats.melt <- melt( istats.merged, id=c("ID","Other.info","continent","ethnicity") )
outimage <- paste(outdir,"/",basename,"_attsbox.png",sep="")
print(paste("Writing:",outimage))
png(outimage)
m <- ggplot(istats.melt, aes(x=variable, y=value,fill=variable))
m <- m + geom_boxplot()
m <- m + ylab("Frequency")+ xlab("Attribute")
m <- m + labs(title = "Boxplots for every Attribute")
print(m)
dev.off()

alloutliers <- NULL
#colnames(alloutliers) <- c( "ID", "Value", "Attribute" )
#istats.merged[,-which(names(istats.merged) %in% c("ID"))]
for( att in names(istats.merged[,-which(names(istats.merged) %in% c("ID","Other.info","continent","ethnicity"))])) {
    print( att)
    outliers <- boxplot(istats.merged[[att]],plot=FALSE)$out
    if( length(outliers) == 0 ){
        print(paste("No outliers",att))
        next
    }
    couts <- istats.merged[(istats.merged[[att]] %in% outliers),][[att]]
    cids <- as.character(istats.merged[(istats.merged[[att]] %in% outliers),]$ID)
    new <- as.data.frame(cbind(cids,couts, att ))
    colnames(new) <- c( "ID", "Value", "Attribute" )
    rownames(new) <- NULL
    alloutliers <- rbind( alloutliers, new )
}

#colnames(alloutliers) <- c( "ID","Attribute","Value")
outlierstowrite <- alloutliers[!(alloutliers$Attribute %in% c("RATE","DP","NVAR")),]
outfile <- paste(filepath,"/",basename,"_outliers.txt",sep="")
print(paste("Writing:",outfile))
write.table( outlierstowrite,outfile,quote=FALSE,row.names=FALSE,sep="\t") 
#istats.merged <- istats.merged[istats.merged$continent != "BGI" & istats.merged$continent != "Unknown",]
istats.merged$continent <- factor(istats.merged$continent, levels=sort(unique(istats.merged$continent)))

plotlist <- list()
affplotset <- list()
i <- 1
for( att in names(istats.merged[,-which(names(istats.merged) %in% c("ID","Other.info","continent","ethnicity"))])) {
    outliers <- boxplot(istats.merged[[att]],plot=FALSE)$out
    if( length(outliers) == 0 ){
        print(paste("No outliers",att))
        next
    }

    qs <- quantile( istats.merged[[att]], na.rm=T )
    irq <- qs[4] - qs[2]
    upperb <- qs[4] + irq*1.5
    lowerb <- qs[2] - irq*1.5
    print( paste("IQR",irq,upperb,lowerb) )
    if( irq == upperb && irq == lowerb ){
        print("No variation!")
        next
    }

    failedatt <- NULL
    for (dsource in unique(istats.merged$continent)) {
        values <- istats.merged[istats.merged$continent == dsource,att]
        print(length(unique(values)))
        print(length(values))
        if (length(unique(values)) == 1 & length(values) > 1) {
            failedatt <- dsource
        }
    }
    if( !is.null(failedatt) ){ print(failedatt); next }

    print(head(istats.merged))
    outimage <- paste(outdir,"/",basename,"_",att,"_violin.png",sep="")
    print(outimage)
    png(outimage)
    m <- ggplot(istats.merged, aes_string(x="continent", y=att, fill="continent"))
    m <- m + geom_violin(show_guide=FALSE)
    m <- m + theme(legend.position = "none")
    m <- m + theme(axis.title.y=element_blank())
    m <- m + theme(axis.title.x=element_blank())
    m <- m + theme(axis.text.x=element_text(angle = 45, hjust = 1))
    m <- m + geom_hline(yintercept=upperb, linetype="dashed")
    m <- m + geom_hline(yintercept=lowerb, linetype="dashed")
    m <- m + labs(title = att)
    print("Got here")
    plotlist[[i]] <- m
    #m <- m + scale_y_log10()
    m <- m + ylab("Variant rate (log)")+ xlab("All rates")
    print(m)
    dev.off()


    #outimage <- paste(outdir,"/",basename,"_",att,"_affviolin.png",sep="")
    #print(outimage)
    #png(outimage)
    #p <- ggplot(istats.merged, aes_string(x="Other.info", y=att, fill="Other.info"))
    #p <- p + geom_violin(show_guide=FALSE)
    #p <- p + theme(legend.position = "none")
    #p <- p + theme(axis.title.y=element_blank())
    #p <- p + theme(axis.title.x=element_blank())
    #p <- p + theme(axis.text.x=element_text(angle = 45, hjust = 1))
    #p <- p + geom_hline(yintercept=upperb, linetype="dashed")
    #p <- p + geom_hline(yintercept=lowerb, linetype="dashed")
    #p <- p + labs(title = att)
    print("got here!")
    #affplotset[[i]] <- p
    #p <- p + ylab("Variant rate (log)")+ xlab("All rates")
    #print(p)
    #dev.off()
 
    print("Next att")
    i <- i + 1
}
if(length(plotlist) == 1){print("Only one plot!")}
if(length(plotlist) > 1){
    print(paste("More than one plot!",length(plotlist)))
    print("About to multiplot!")
    outimage <- paste(outdir,"/",basename,"_allatts_violin.png",sep="")
    print(paste("Writing:",outimage))
    png(outimage)
    multiplot(plotlist=plotlist, cols=3 )
    dev.off()
}
if(length(affplotset) > 1){
    print(paste("More than one plot!",length(affplotset)))
    outimage <- paste(outdir,"/",basename,"_affatts_violin.png",sep="")
    print(paste("Writing:",outimage))
    png(outimage)
    multiplot(plotlist=affplotset, cols=3 )
    dev.off()
}


# Numbers of outlier types
#outlierlist <- as.data.frame(table(alloutliers$Attribute))
alloutliers <- merge(alloutliers,patannot, by.x="ID", by.y="Individual.ID")
outimage <- paste(outdir,"/",basename,"_outliers_bar.png",sep="")
print(paste("Writing:",outimage))
m <- ggplot(alloutliers,aes(x=Attribute, fill=continent)) + geom_bar()
m <- m + ylab("Freq")+ xlab("QC Class")
m <- m + labs(title = "Number of outliers ") + ng1
m <- m + theme(axis.text.x=element_text(angle = 45, hjust = 1))
png(outimage)
print(m)
dev.off()

# Proportions per source
#m <- ggplot(alloutliers, aes(x = )) + 
     #geom_bar(aes(y = (..count..)/sum(..count..))) + 
     #scale_y_continuous(formatter = 'percent')
istats.merged <- merge( istats, patannot, by.x="ID", by.y="Individual.ID" )
istats.merged$continent <- factor( istats.merged$continent, levels=sort(unique(istats.merged$continent) ))

alloutliers$continent <- factor(alloutliers$continent, levels=sort(unique(alloutliers$continent)))
attbreakdown <- as.data.frame(table(alloutliers[,c("ID","continent","Attribute")]))
attbreakdown <- attbreakdown[attbreakdown$Freq !=0,]

s <- as.data.frame( table(attbreakdown[,c("continent","Attribute","Freq")]) )
#s <- s[s$continent != "BGI" & s$continent != "Unknown",]
colnames(s) <- c("continent","Attribute","Failed","Freq")
sourcecounts <- as.data.frame(table(istats.merged$continent))
colnames(sourcecounts) <- c("continent","Total")
attbreakdown <- merge(s, sourcecounts, on="continent", all.x=TRUE)
attbreakdown$Percent <- attbreakdown$Freq / attbreakdown$Total*100
attbreakdown <- attbreakdown[attbreakdown$Freq > 0,]
#attbreakdown <- attbreakdown[!(attbreakdown$Attribute %in% c("RATE","DP","NVAR")),]
#attbreakdown$Failed <- as.numeric(as.character(attbreakdown$Failed))
#attbreakdown$Failed[attbreakdown$Failed>3] <- '>3'
require(scales)
m <- ggplot(attbreakdown, aes(x = continent, y=Percent, fill=continent)) + geom_bar(stat="identity")#, fill=factor(Failed)
m <- m + facet_wrap(~Attribute)
m <- m + labs(title = "Distribution of Outliers") + ng1
m <- m + theme(axis.text.x=element_text(angle=45, hjust=1))
outimage <- paste(outdir,"/",basename,"_source_bar.png",sep="")
print(paste("Writing:",outimage))
png(outimage)
print(m)
dev.off()

#m <- m + scale_y_continuous(labels="percent")
#m <- m + scale_y_continuous("",formatter="percent")
