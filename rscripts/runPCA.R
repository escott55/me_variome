#!/usr/bin/env Rscript

#source("http://bioconductor.org/biocLite.R")
#biocLite("snpStats")

library(RColorBrewer)
palette(brewer.pal(8, "Dark2"))
#library(gplots)
library(ggplot2)
library(reshape)
#library(MultiPhen)
library(snpStats)

source("themes.R")


######################################################################
# LAZER data
######################################################################
plotLazerData <- function(outputdir ){
    pcs <- read.table("data/reference_pca.tsv",header=TRUE,sep="\t")
    png(paste(outputdir,"lazerpcaref.png",sep=""))
    m <- ggplot(pcs,aes(x=PC1, y=PC2)) + geom_point(aes(colour=popID))
    m <- m + labs(title = "Lazer PCA Reference colored by Ethnicity")
    print(m)
    dev.off()

    inhouse.pcs <- read.table("targetannotated.ped",header=TRUE,sep="\t")
    png(paste(outputdir,"lazerpcacontinent.png",sep=""))
    m <- ggplot(inhouse.pcs,aes(x=PC1, y=PC2)) + geom_point(aes(colour=Continent))
    m <- m + labs(title = "Lazer PCA colored by Continent")
    print(m)
    dev.off()
}

######################################################################
# hclusterPCA
######################################################################
hclusterPCA <- function(mydata,bigK,outputdir)
{
    # Ward Hierarchical Clustering
    d <- dist(mydata, method = "euclidean") # distance matrix
    fit <- hclust(d, method="ward")
    outfile <- paste(outputdir,"clustergram_K",bigK,".png",sep="")
    print(paste("Writing file:",outfile))
    png( outfile )
    plot(fit) # display dendogram
    groups <- cutree(fit, k=bigK) # cut tree into K clusters
    # draw dendogram with red borders around the K clusters 
    rect.hclust(fit, k=bigK, border="red")
    dev.off()

    groupdata <- cbind(mydata,groups)
    outfile <- paste(outputdir,"/pcagrouping.",bigK,"_groups.txt",sep="")
    write.table(groupdata,file=outfile, quote=FALSE,row.names=FALSE, col.names=TRUE, sep="\t")

    outfile <- paste(outputdir,"pca_clustK",bigK,"grpbreakdown.txt",sep="")
    write("Group Breakdown",file=outfile,append=FALSE)
    for( group in unique(groupdata$groups) ){
        sub <- groupdata[groupdata$groups == group,]$Population
        s <- as.data.frame(table(sub)) 
        write(paste("#Group",group,":",length(sub)),file=outfile,append=TRUE)
        write.table(s[s$Freq != 0,],file=outfile,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }

    return(groups)
} # End hclusterPCA

######################################################################
# runPCAnalysis
######################################################################
runPCAnalysis <- function( bedprefix, annotfile, outputdir, fileprefix="" ) {
    print("Run function runPCAnalysis")
    bedfile <- paste(bedprefix,".bed",sep="")
    bimfile <- paste(bedprefix,".bim",sep="")
    famfile <- paste(bedprefix,".fam",sep="")

    # select.snps works on a list
    plink.full <- read.plink(bedfile,bimfile,famfile)
    allpat.annot <- read.table(annotfile, header=TRUE,sep="\t",comment.char="")
    plink.annot.full <- merge( plink.full$fam, allpat.annot, by.x="member", by.y="Individual.ID", sort=FALSE )
    print(head(plink.annot.full))
    print(head(plink.full$fam$member))
    print(head(plink.annot.full$member))
    print(tail(plink.full$fam$member))
    print(tail(plink.annot.full$member))

    #if( fileprefix == "collective" ){
        #pat.sub <- plink.annot.full$member[plink.annot.full$Source %in% c("Plate5","PlateIII","PlateIV","PlateV","PlateVI")]
    #}else if( fileprefix == "daily" ) {
        #pat.sub <- plink.annot.full$member[plink.annot.full$Source %in% c("Daly")]
    #}
    geno.mat <- plink.full$genotypes[rownames(plink.full$genotypes) %in% plink.annot.full$member,]
    plink.annot <- plink.annot.full#[plink.annot.full$member %in% plink.annot.full$member,]

    s <- plink.full$fam$member
    t <- plink.annot$member
    s[!(s %in% t)]

    ###################################################
    ### code chunk number 2: xxt-matrix
    ###################################################
    xxmat <- xxt(geno.mat, correct.for.missing=FALSE)

    ###################################################
    ### code chunk number 3: eigen
    ###################################################
    evv <- eigen(xxmat, symmetric=TRUE)
    numvectors <- 260
    pcs <- evv$vectors[,1:numvectors]
    pcs.cols <- list()
    for( i in 1:numvectors ){ pcs.cols[i] = (paste("PC",i,sep="")) }
    #colnames(pcs) <- c("PC1","PC2","PC3","PC4","PC5")
    colnames(pcs) <- pcs.cols
    rownames(pcs) <- plink.annot$member
    evals <- evv$values[1:24]
    evals
    all.pcs <- cbind(plink.annot,pcs[,1:5])

    #all.pcs.filt <- all.pcs[all.pcs$Country != "Unknown",]
    all.pcs.filt <- all.pcs
    all.pcs.filt$Origin <- factor(all.pcs.filt$Origin, levels = sort(unique(all.pcs.filt$Origin)))
    head(all.pcs.filt$Origin)

    outfilename <- paste(outputdir,fileprefix,"_pcs.txt",sep="")
    print(paste("Writing file:",outfilename))
    write.table(all.pcs.filt,file=outfilename,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
    # Compute proportion of the total variance explained by the first eigenvalues
    cumprop <- cumsum(evv$values)/sum(evv$values)
    #pc.sds <- apply(pcs,2,sd)
    #variances <- data.frame(variances=pc.sds**2, pcomp=1:length(pc.sds))
    variances <- data.frame(variance=cumprop[1:25], pcomp=1:25)
    outfilename <- paste(outputdir,fileprefix,"pcavariances.png",sep="")
    print(paste("Writing file:",outfilename))
    png(outfilename)
    varPlot <- ggplot(variances, aes(pcomp, variance)) + geom_bar(stat="identity", 
            fill="grey") + geom_line()
    print(varPlot)
    dev.off()

    # Plot Principle Components
    outfilename <- paste(outputdir,fileprefix,"_pcacontinent.png",sep="")
    print(paste("Writing file:",outfilename))
    png(outfilename)
    m <- ggplot(all.pcs.filt,aes(x=PC1, y=PC2)) + geom_point(aes(colour=continent))
    m <- m + labs(title = "PCA colored by Continent") + ng1
    print(m)
    dev.off()

    outfilename <- paste(outputdir,fileprefix,"_pcacountry.png",sep="")
    print(paste("Writing file:",outfilename))
    png(outfilename)
    m <- ggplot(all.pcs.filt,aes(x=PC1, y=PC2)) + geom_point(aes(colour=Origin))
    m <- m + labs(title = "PCA colored by Origin") + ng1
    print(m)
    dev.off()

    outfilename <- paste(outputdir,fileprefix,"_pcaethnicity.png",sep="")
    print(paste("Writing file:",outfilename))
    png(outfilename)
    m <- ggplot(all.pcs.filt,aes(x=PC1, y=PC2)) + geom_point(aes(colour=ethnicity))
    m <- m + labs(title = "PCA colored by Ethnicity") + ng1
    print(m)
    dev.off()

    outfilename <- paste(outputdir,fileprefix,"_pcaplate.png",sep="")
    print(paste("Writing file:",outfilename))
    png(outfilename)
    m <- ggplot(all.pcs.filt,aes(x=PC1, y=PC2)) + geom_point(aes(colour=Plate))
    m <- m + labs(title = "PCA colored by Plate") + ng1
    print(m)
    dev.off()

    outfilename <- paste(outputdir,fileprefix,"_pcasource.png",sep="")
    print(paste("Writing file:",outfilename))
    png(outfilename)
    m <- ggplot(all.pcs.filt,aes(x=PC1, y=PC2)) + geom_point(aes(colour=Source))
    m <- m + labs(title = "PCA colored by Source") + ng1
    print(m)
    dev.off()

    #outfilename <- paste(outputdir,fileprefix,"_pcabatch.png",sep="")
    #print(paste("Writing file:",outfilename))
    #png(outfilename)
    #m <- ggplot(all.pcs.filt,aes(x=PC1, y=PC2)) + geom_point(aes(colour=Batch))
    #m <- m + labs(title = "PCA colored by Batch") + ng1
    #print(m)
    #dev.off()


    return( list(all.pcs.filt,pcs) )
} # End 

######################################################################
# Main
######################################################################

#bedprefix <- "plinkfiles/merged.1kg.percentcoverage.plink.mod.filt.uniq"
#bedprefix <- "plinkfiles/merged.eichler.plink.mod.filt.uniq"
#bedprefix <- "../data/merged/ciliopathies.unfilt.recode.plink.mod.filt.uniq"
#bedprefix <- "../rawdata/daily1/daily.recode.plink.recode12.mod.filt.uniq"
#bedprefix <- "../data/collective/ciliopathies.unfilt.recode.plink.recode12.mod.filt.uniq"
#fileprefix <- "collective"

setwd( "/home/escott/workspace/variome/rscripts/" )
outputdir <- "pcafigures/"

bedprefix <- "../rawdata/test/everything_set1.chr1.snp.chimp.regions.recode.recode.plink.mod.filt.uniq"
fileprefix <- "chr1test"
#bedprefix <- "../rawdata/variome/variome_snps.chimp.regions.recode.recode.plink.mod.filt.uniq"
bedprefix <- "../rawdata/variome1/variome.chimp.regions.filt.samp.plink.mod.filt.uniq"
fileprefix <- "variome"
annotfile <- "../resources/annotation/patientannotation.ped"

# parse args
args <- commandArgs( trailingOnly = TRUE )
bedprefix <- ifelse(length(args)>=1,dirname(args[1]),bedprefix)
fileprefix <- ifelse(length(args)>=2,dirname(args[2]),fileprefix)


s <- runPCAnalysis( bedprefix, annotfile, outputdir, fileprefix )
all.pcs.filt <- s[[1]]
pcs <- s[[2]]

q()
bigK <- 11
#basename <- paste(outputdir,"clustergram",sep="")
groups <- hclusterPCA( all.pcs.filt, bigK, outputdir )

print("Kmeans")
t <- kmeans( pcs, bigK, nstart=25 )

selectK <- list()
bestK <- 0
for( i in 1:50 ){
    t <- kmeans( pcs[,1:50], bigK, nstart=10 )
    #selectK[[i]] <- t$betweenss / t$totss *100
    #selectK[i] <- sum(t$withinss)
    selectK[i] <- sum(t$tot.withinss)
}
selectK.df <- data.frame(Within=as.numeric(as.character(selectK)),K=1:length(selectK))

outfilename <- paste(outputdir,fileprefix,"kmeansselection.png",sep="")
print(paste("Writing file:",outfilename))
png(outfilename)
varPlot <- ggplot(selectK.df, aes(K, Within)) + geom_point() + geom_line()
# + geom_bar()
print(varPlot)
dev.off()

outfile <- paste(outputdir,"pca_kmeans",bigK,"grpbreakdown.txt",sep="")
write("Group Breakdown",file=outfile,append=FALSE)
for( grp in sort(unique(t$cluster)) ){
    pnames <- names(t$cluster[t$cluster == grp])
    sub <- all.pcs.filt[all.pcs.filt$member %in% pnames, ]$Population

    s <- as.data.frame(table(sub)) 
    write(paste("#Group",grp,":",length(sub)),file=outfile,append=TRUE)
    write.table(s[s$Freq != 0,],file=outfile,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
}




# End Main
