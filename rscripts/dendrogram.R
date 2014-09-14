#!/usr/bin/env Rscript

library(RColorBrewer)
palette(brewer.pal(8, "Dark2"))
library(ggdendro)
library(ggplot2)
library(reshape)
library(plyr)

source("/home/escott/workspace/variome/rscripts/themes.R")

plotAdmixDendrogram <- function(tfile,K, clustmethod="ward") {
    print(tfile)

    fbase <- strsplit(basename(tfile),'\\.')[[1]][1]

    figname <- paste("dendrograms/",fbase,"_admix_",clustmethod,".pdf",sep="") 

    if( clustmethod == "ward" ){
        clustmethod <- "ward.D2"
    }

    admix <- read.table(tfile, header=TRUE, sep="\t", comment.char="")
    rownames(admix) <- admix$IID

    x <- t(as.matrix(admix[,paste("K",1:K,sep="")]))

    dd.row <- as.dendrogram(hclust(dist(t(x)), method="ward.D2"))
    ddata_x <- dendro_data(dd.row)

    grlabs <- label(ddata_x)
    grlabs$group <- admix[as.vector(grlabs$label),"GeographicRegions"]

    xsegments <- segment(ddata_x)
    xsegments$newy <- xsegments$y / max(xsegments$y)*10
    xsegments$newyend <- xsegments$yend / max(xsegments$yend)*10

    allcolors <- c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(10,"Set3"))
    p2 <- (ggplot(xsegments) +
           geom_segment(aes(x=x, y=newy, xend=xend, yend=newyend)) + 
           geom_tile(data=grlabs, aes(x=x, y=-.5, fill=group), size=.25 ) +
           ggtitle( paste(fbase,":",clustmethod) ) +
           scale_y_continuous("Distance") +
           scale_x_continuous("Sample") +
           #scale_fill_brewer(palette="Set3") +
           scale_fill_manual(values=allcolors[1:length(unique(grlabs$group))]) +
           hclusttheme )

    print(paste("Writing file:",figname))
    pdf( figname, width=10)
    plot(p2)
    dev.off()

} # plotAdmixDendrogram


plotPCAdendro <- function( pcafile, withdaisy=F, clustmethod="ward" ){
    require(cluster)

    fbase <- strsplit(basename(pcafile),'\\.')[[1]][1]

    pcatable <- read.table(pcafile, header=F, sep="", comment.char="", na.strings="", fill=T)
    rownames(pcatable) <- pcatable$V1
    pcatable <- pcatable[,2:(ncol(pcatable))]
    eigvals <- pcatable[1,1:(ncol(pcatable)-1)]
    pcatable <- pcatable[2:(nrow(pcatable)),]
    tannot <- pcatable[,ncol(pcatable)]
    contcompare <- c("LWK","CHB","CHS","JPT","YRI","CEU","GBR","IBS","TSI","FIN")
    cont <- c("Africa","East.Asia","Europe")
    pca <- pcatable[tannot %in% contcompare,1:(ncol(pcatable)-1)]

    figname <- paste("dendrograms/",fbase,"_",clustmethod,".pdf",sep="") 

    if( clustmethod == "ward" ){
        clustmethod <- "ward.D2"
    }

    if( withdaisy )
    {
        dd.row <- as.dendrogram(hclust(daisy(pca,metric="euclidean",weights=eigvals), method=clustmethod))
        ddata_x <- dendro_data(dd.row)

        p2 <- ggplot(segment(ddata_x)) +
              geom_segment(aes(x=x, y=y, xend=xend, yend=yend))

        grlabs <- label(ddata_x)
        grlabs$group <- pcatable[as.vector(grlabs$label),ncol(pcatable)]#"GeographicRegions"

        p2 <- p2 + geom_tile(data=grlabs, aes(x=x, y=-.5, fill=group), size=.25 )#label=label, 
    }
    else 
    {
        # ward.D, ward.D2, single, complete, average, mcquitty
        pcamat <- sweep( as.matrix(pca), MARGIN=2, as.numeric(as.vector(eigvals)), `*` )
        dd.row <- as.dendrogram(hclust(dist(pcamat), method=clustmethod))
        #pcamat <- as.matrix(pca)
        #dd.row <- as.dendrogram(hclust(daisy(pca,metric="manhattan",weights=eigvals), method=clustmethod))
        ddata_x <- dendro_data(dd.row)

        xsegments <- segment(ddata_x)
        xsegments$newy <- xsegments$y / max(xsegments$y)*10
        xsegments$newyend <- xsegments$yend / max(xsegments$yend)*10
        #p2 <- ggplot(segment(ddata_x)) +
              #geom_segment(aes(x=x, y=y, xend=xend, yend=yend))

        grlabs <- label(ddata_x)
        grlabs$group <- pcatable[as.vector(grlabs$label),ncol(pcatable)]#"GeographicRegions"

        allcolors <- c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(10,"Set3"))
        p2 <- (ggplot(xsegments) +
               geom_segment(aes(x=x, y=newy, xend=xend, yend=newyend)) +
               geom_tile(data=grlabs, aes(x=x, y=-.5, fill=group), size=.25 ) + 
               ggtitle( paste(fbase,":",clustmethod) ) +
               scale_y_continuous("Distance") +
               scale_x_continuous("Sample") +
               #scale_fill_brewer(palette="Set3") +
               scale_fill_manual(values=allcolors[1:length(unique(grlabs$group))]) +
               hclusttheme )
    }
    print(paste("Writing file:",figname))
    pdf( figname, width=10 )
    plot(p2)
    dev.off()
}

# plotIBSdendrogram
plotIBSdendrogram <- function( famfile, distfile, clustmethod="ward", 
                              nclust=7, targetdir="./dendrograms" ){

    fbase <- strsplit(basename(distfile),'\\.')[[1]][1]
    famtable <- read.table(famfile, header=T, sep="\t", comment.char="", na.strings="")
    rownames(famtable) <- famtable$IID

    ibsdist <- read.table(distfile, header=F, sep="", comment.char="", na.strings="", fill=T)
    rownames(ibsdist) <- famtable$IID
    colnames(ibsdist) <- famtable$IID
    print( paste("Any missing values:",any(is.na(as.matrix(ibsdist))) ))

    onekgpops <- c("ASW","CEU","CHB","CHS","CLM","FIN","GBR",
                   "IBS","JPT","LWK","MXL","PUR","TSI","YRI")
    tokeep <- !(is.na(famtable$GeographicRegions) | 
                      famtable$GeographicRegions %in% c("America","Oceania","")  |
                      famtable$GeographicRegions2 %in% c(c("East Asia"),onekgpops) ) 
    #"Africa","Europe","Unknown"
    tokeep <- !( is.na(famtable$GeographicRegions) ) 

    famtable <- famtable[tokeep,]
    ibsdist <- ibsdist[tokeep,tokeep]

    figname <- paste(targetdir,"/ibs_",fbase,"_",clustmethod,".pdf",sep="") 
    labfile <- paste(targetdir,"/ibs_",fbase,"_",clustmethod,".txt",sep="") 
    if( clustmethod == "ward" ){
        clustmethod <- "ward.D2"
    }

    # ward.D, ward.D2, single, complete, average, mcquitty
    hc <- hclust(as.dist(ibsdist), method=clustmethod)

    # Remove Africans!
    if( F ) {
        myclusters <- as.data.frame(cutree( hc, 2 ))
        colnames(myclusters) <- c("clust")
        myclusters$labels <- rownames(myclusters)
        myclusters$group <- famtable[as.vector(myclusters$label),"Continent2"]
        gcounts <- melt(table( myclusters[,c("group","clust")] ) )
        clustrename <- ddply(gcounts, "group", function(x) x[which.max(x$value),])
        rownames(clustrename) <- clustrename$group
        meclust <- clustrename["Middle East","clust"]
        tokeep <- as.vector(myclusters[myclusters$clust == meclust,"labels"])
        ibsdist_filt <- ibsdist[tokeep,tokeep]
        hc <- hclust(as.dist(ibsdist_filt), method=clustmethod)
    }

    myclusters <- cutree( hc, nclust )

    dd.row <- as.dendrogram(hc)
    ddata_x <- dendro_data(dd.row)

    xsegments <- segment(ddata_x)
    xsegments$newy <- xsegments$y / max(xsegments$y)*10
    xsegments$newyend <- xsegments$yend / max(xsegments$yend)*10

    grlabs <- label(ddata_x)
    grlabs$group <- famtable[as.vector(grlabs$label),"GeographicRegions2"]#"GeographicRegions"
    grlabs$eth <- famtable[as.vector(grlabs$label),"ethnicity"]
    grlabs$source <- famtable[as.vector(grlabs$label),"Source"]
    grlabs$origin <- famtable[as.vector(grlabs$label),"Origin"]
    grlabs$clust <- as.factor(myclusters[as.vector(grlabs$label)])

    print(table( grlabs[,c("group","clust")] ))
    gcounts <- melt(table( grlabs[,c("group","clust")] ) )
    #maxclust <- data.frame(tapply(gcounts[,3], gcounts[,2], max))
    clustrename <- ddply(gcounts, "clust", function(x) x[which.max(x$value),])
    colnames( clustrename ) <- c("closest","clust","nsamp")
    #dups <- ddply(clustrename,.(group),nrow)
    #clustrename <- merge( clustrename, dups, by="group" )
    grlabs <- merge( grlabs, clustrename, by="clust" )

    write.table( grlabs, file=labfile,quote=F, row.names=F, col.names=T, sep="\t" )

    allcolors <- c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(10,"Set3"))
    p2 <- (ggplot(xsegments) +
           geom_segment(aes(x=x, y=newy, xend=xend, yend=newyend)) +
           geom_tile(data=grlabs, aes(x=x, y=-.5, fill=group), size=.25 ) +
           geom_tile(data=grlabs, aes(x=x, y=-1.5, fill=clust), size=.25 ) +
           #geom_tile(data=grlabs, aes(x=x, y=-.5, fill=eth), size=.25 ) +
           scale_fill_manual("Clusters",values=allcolors) +
           ggtitle( paste(fbase,":",clustmethod) ) +
           scale_y_continuous("Distance") +
           scale_x_continuous("Sample") +
           hclusttheme )
           #scale_fill_manual(values=allcolors[1:length(unique(grlabs$group))]) +
    p2

    print(paste("Writing file:",figname))
    pdf( figname, width=10 )
    plot(p2)
    dev.off()
} # END plotIBSdendrogram

plotMDSdendrogram <- function( famfile, mdsfile, clustmethod="ward", targetdir="./dendrograms" ){
    require(cluster)

    fbase <- strsplit(basename(mdsfile),'\\.')[[1]][1]
    famtable <- read.table(famfile, header=T, sep="\t", comment.char="", na.strings="")
    rownames(famtable) <- famtable$IID

    mdstable <- read.table(mdsfile, header=T, sep="", comment.char="", na.strings="", fill=T)
    rownames(mdstable) <- mdstable$IID
    mds <- mdstable[,4:ncol(mdstable)]
    #contcompare <- c("LWK","CHB","CHS","JPT","YRI","CEU","GBR","IBS","TSI","FIN")
    #cont <- c("Africa","East.Asia","Europe")
    #mds <- mdstable[tannot %in% contcompare,1:(ncol(mdstable)-1)]

    tokeep <- !(is.na(famtable$GeographicRegions) | 
                      famtable$GeographicRegions %in% c("Unknown","America","Oceania","") |
                      famtable$GeographicRegions2 %in% c("Africa","Europe","East Asia") )

    famtable <- famtable[tokeep,]
    mds <- mds[tokeep,]
    figname <- paste(targetdir,"/mds_",fbase,"_",clustmethod,".pdf",sep="") 

    if( clustmethod == "ward" ){
        clustmethod <- "ward.D2"
    }

    # ward.D, ward.D2, single, complete, average, mcquitty
    dd.row <- as.dendrogram(hclust(dist(as.matrix(mds)), method=clustmethod))
    ddata_x <- dendro_data(dd.row)

    xsegments <- segment(ddata_x)
    xsegments$newy <- xsegments$y / max(xsegments$y)*10
    xsegments$newyend <- xsegments$yend / max(xsegments$yend)*10

    grlabs <- label(ddata_x)
    grlabs$group <- famtable[as.vector(grlabs$label),"GeographicRegions2"]#"GeographicRegions"

    allcolors <- c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(10,"Set3"))
    p2 <- (ggplot(xsegments) +
           geom_segment(aes(x=x, y=newy, xend=xend, yend=newyend)) +
           geom_tile(data=grlabs, aes(x=x, y=-.5, fill=group), size=.25 ) + 
           ggtitle( paste(fbase,":",clustmethod) ) +
           scale_y_continuous("Distance") +
           scale_x_continuous("Sample") +
           #scale_fill_brewer(palette="Set3") +
           scale_fill_manual(values=allcolors[1:length(unique(grlabs$group))]) +
           hclusttheme )
    
    print(paste("Writing file:",figname))
    pdf( figname, width=10 )
    plot(p2)
    dev.off()
} # plotMDSdendrogram


runTests <- function(){
    tfile <- "/media/data/workspace/variome/rawdata/merge1kg/main/admixture/Middle_East/me1000G.clean_Middle_East.5_annot.tsv"
    K <- 5
    tfile <- "/media/data/workspace/variome/rawdata/merge1kg/main/admixture/me1000G.clean.10_annot.tsv"
    K <- 10
    #tfile <- "/media/data/workspace/variome/rawdata/onekg/main/admixture/onekg.clean.6_annot.tsv"
    #K <- 6
    #plotAdmixDendrogram( tfile, K )

    #pcafile <- "/media/data/workspace/variome/rawdata/merge1kg/main/pca/pca/GeographicRegions2/me1000G.clean.recode12.exclude.evec"
    pcafile <- "/media/data/workspace/variome/rawdata/onekg/main/pca/pca/GeographicRegions2/onekg.clean.recode12.exclude.evec"
    #plotPCAdendro( pcafile )


    # ONEKG
    famfile <- "/media/data/workspace/variome/rawdata/onekg/main/clust/onekg.clean_annot.fam"
    distfile <- "/media/data/workspace/variome/rawdata/onekg/main/clust/onekg.clean.mdist"
    mdsfile <- "/media/data/workspace/variome/rawdata/onekg/main/clust/onekg.clean.mds"

    # Merge 1KG variome
    famfile <- "/media/data/workspace/variome/rawdata/merge1kg/main/clust/me1000G.clean_annot.fam"
    distfile <- "/media/data/workspace/variome/rawdata/merge1kg/main/clust/me1000G.clean.mdist"
    mdsfile <- "/media/data/workspace/variome/rawdata/merge1kg/main/clust/me1000G.clean.mds"

    famfile <- "/media/data/workspace/variome/rawdata/mevariome/main/clust/variome.clean_annot.fam"
    distfile <- "/media/data/workspace/variome/rawdata/mevariome/main/clust/variome.clean.mdist"
    mdsfile <- "/media/data/workspace/variome/rawdata/mevariome/main/clust/variome.clean.mds"

    famfile <- "/media/data/workspace/variome/rawdata/daily/main/clust/daily.clean_annot.fam"
    distfile <- "/media/data/workspace/variome/rawdata/daily/main/clust/daily.clean.mdist"
    mdsfile <- "/media/data/workspace/variome/rawdata/daily/main/clust/daily.clean.mds"

    plotIBSdendrogram( famfile, distfile, clustmethod="ward", nclust=10 )

    plotMDSdendrogram( famfile, mdsfile, clustmethod="ward" )

    allclustmethods <- c( "ward.D", "ward.D2", "complete", "average", "mcquitty" )
    for( clustm in allclustmethods ){
        print(paste("Running method",clustm))
        plotIBSdendrogram( famfile, distfile, clustm )
        plotMDSdendrogram( famfile, mdsfile, clustm )
    }
}
#genotypes <- read.table(pedfile, header=F, sep="", comment.char="", na.strings="", fill=T)

#labs$group <- with( labs, admix [ cbind(label, "GeographicRegions") ] )
#labs$group <- c(rep("Clust1", 5), rep("Clust2", 2), rep("Clust3", 4))
#labs
#p3<-ggplot(labs,aes(x,y,fill=factor(group)))+geom_tile()+
  #scale_y_continuous(expand=c(0,0))+
  #theme(axis.title=element_blank(),
        #axis.ticks=element_blank(),
        #axis.text=element_blank(),
        #legend.position="none")
#
#plot(p2)

