#!/usr/bin/env Rscript

library(RColorBrewer)
palette(brewer.pal(8, "Dark2"))
library(ggdendro)
library(ggplot2)

source("themes.R")

plotAdmixDendrogram <- function(tfile,K) {
    print(tfile)

    admix <- read.table(tfile, header=TRUE, sep="\t", comment.char="")
    rownames(admix) <- admix$IID

    x <- t(as.matrix(admix[,paste("K",1:K,sep="")]))

    dd.row <- as.dendrogram(hclust(dist(t(x)), method="average"))#, method="canberra"
    ddata_x <- dendro_data(dd.row)

    p2 <- ggplot(segment(ddata_x)) +
          geom_segment(aes(x=x, y=y, xend=xend, yend=yend))

    labs <- label(ddata_x)
    labs$group <- admix[as.vector(labs$label),"GeographicRegions"]

    p2 + geom_tile(data=labs, aes(x=x, y=-.5, fill=group), size=.25 )#label=label, 
}


plotPCAdendro <- function( pcafile, withdaisy=F, clustmethod="ward" ){
    require(cluster)

    basename <- strsplit(basename(pcafile),'\\.')[[1]][1]

    pcatable <- read.table(pcafile, header=F, sep="", comment.char="", na.strings="", fill=T)
    rownames(pcatable) <- pcatable$V1
    pcatable <- pcatable[,2:(ncol(pcatable))]
    eigvals <- pcatable[1,1:(ncol(pcatable)-1)]
    pcatable <- pcatable[2:(nrow(pcatable)),]
    tannot <- pcatable[,ncol(pcatable)]
    contcompare <- c("LWK","CHB","CHS","JPT","YRI","CEU","GBR","IBS","TSI","FIN")
    cont <- c("Africa","East.Asia","Europe")
    pca <- pcatable[tannot %in% contcompare,1:(ncol(pcatable)-1)]

    figname <- paste("dendrograms/",basename,"_",clustmethod,".pdf",sep="") 

    if( clustmethod == "ward" ){
        clustmethod <- "ward.D2"
    }

    if( withdaisy )
    {
        dd.row <- as.dendrogram(hclust(daisy(pca,metric="euclidean",weights=eigvals), method=clustmethod))
        ddata_x <- dendro_data(dd.row)

        p2 <- ggplot(segment(ddata_x)) +
              geom_segment(aes(x=x, y=y, xend=xend, yend=yend))

        labs <- label(ddata_x)
        labs$group <- pcatable[as.vector(labs$label),ncol(pcatable)]#"GeographicRegions"

        p2 <- p2 + geom_tile(data=labs, aes(x=x, y=-.5, fill=group), size=.25 )#label=label, 
    }
    else 
    {
        #pcamat <- as.matrix(pca)
        pcamat <- sweep( as.matrix(pca), MARGIN=2, as.numeric(as.vector(eigvals)), `*` )
        # ward.D, ward.D2, single, complete, average, mcquitty
        dd.row <- as.dendrogram(hclust(dist(pcamat), method=clustmethod))
        ddata_x <- dendro_data(dd.row)

        xsegments <- segment(ddata_x)
        xsegments$newy <- xsegments$y / max(xsegments$y)*10
        xsegments$newyend <- xsegments$yend / max(xsegments$yend)*10
        #p2 <- ggplot(segment(ddata_x)) +
              #geom_segment(aes(x=x, y=y, xend=xend, yend=yend))

        labs <- label(ddata_x)
        labs$group <- pcatable[as.vector(labs$label),ncol(pcatable)]#"GeographicRegions"

        allcolors <- c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(10,"Set3"))
        p2 <- (ggplot(xsegments) +
               geom_segment(aes(x=x, y=newy, xend=xend, yend=newyend)) +
               geom_tile(data=labs, aes(x=x, y=-.5, fill=group), size=.25 ) + 
               ggtitle( paste(basename,":",clustmethod) ) +
               scale_y_continuous("Distance") +
               scale_x_continuous("Sample") +
               #scale_fill_brewer(palette="Set3") +
               scale_fill_manual(values=allcolors[1:length(unique(labs$group))]) +
               hclusttheme )
    }
    print(paste("Writing file:",figname))
    pdf( figname )
    plot(p2)
    dev.off()
}

tfile <- "/media/data/workspace/variome/rawdata/merge1kg/main/admixture/Middle_East/me1000G.clean_Middle_East.5_annot.tsv"
K <- 5
tfile <- "/media/data/workspace/variome/rawdata/merge1kg/main/admixture/me1000G.clean.10_annot.tsv"
K <- 10
#plotAdmixDendrogram( tfile, K )

#pcafile <- "/media/data/workspace/variome/rawdata/merge1kg/main/pca/pca/GeographicRegions2/me1000G.clean.recode12.exclude.evec"
pcafile <- "/media/data/workspace/variome/rawdata/onekg/main/pca/pca/GeographicRegions2/onekg.clean.recode12.exclude.evec"
#plotPCAdendro( pcafile )


pedfile <- "/media/data/workspace/variome/rawdata/onekg/main/pca/onekg.clean.recode12.mod.ped"


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

