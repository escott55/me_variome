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


ancientPCA <- function( snpdata, outputdir, fileprefix="" ){

    geno.mat <- snpdata$genotypes[rownames(snpdata$genotypes) %in% snpdata$annot$member,]

    if( fileprefix == "hgdp" ){
        humanpops <- c("Adygei","Balochi","BantuKenya","BantuSouthAfrica","Basque","Bedouin",
                   "BiakaPygmy", "Brahui","Burusho","Cambodian","Colombian","Druze",
                  "French","Han","Han-NChina","Hazara","Hezhen","Href","Italian","Japanese",
                  "Kalash","Karitiana","Lahu","Makrani","Mandenka","Maya",
                  "MbutiPygmy","Melanesian","Miao","Mongola","Mozabite","Naxi","Orcadian",
                  "Oroqen","Palestinian","Papuan","Pathan","Pima","Russian","San",
                  "Sardinian","She","Sindhi","Surui","Tu","Tujia","Tuscan","Uygur",
                  "Xibo","Yakut","Yi","Yoruba")
    }else{
        print( "Error! havent done this yet" )
        q()
    }

    #"San_discover","Neander","Orang","Chimp","Denisova","Marmoset","Dai","Gorilla","Daur","Macaque",

    tmat <- as(geno.mat[row.names(geno.mat) %in% c("Denisova","Vindija","Clint"),],"numeric")
    tmat <- replace(tmat, is.na(tmat), 0)

    #test <- as(geno.mat[row.names(geno.mat) %in% c("Denisova","Vindija","Clint"),1:100],"numeric")
    humansamps <- as.vector(snpdata$annot$member[snpdata$annot$region %in% humanpops])

    testmat <- as(geno.mat[row.names(geno.mat) %in% humansamps,],"numeric")
    testmat <- replace(testmat, is.na(testmat), 0)

    #preform principal components analysis
    pca<-prcomp(na.omit(tmat)) 

    #project new data onto the PCA space
    x <- scale(testmat,pca$center,pca$scale)%*%pca$rotation 

    projected <- as.data.frame(rbind( x, pca$x ))
    projected$member <- rownames(projected)
    projected <- merge(projected, snpdata$annot, on="member")
    projected$region <- factor(projected$region, levels=sort(unique(projected$region)))

    primary <- projected[projected$member %in% c("Denisova","Vindija","Clint"),]
    cnames <- aggregate(cbind(PC1, PC2) ~ region, data=primary, 
                                    FUN=function(x)mean(range(x)))
    cnames$PC1 <- cnames$PC1 *.90
    cnames$PC2 <- cnames$PC2 *.90
    centroid <- aggregate( cbind(PC1, PC2) ~ region, data=primary, 
                                    FUN=function(x) mean(range(x)))
    centroid <- as.data.frame(t(apply( primary[,c("PC1","PC2","PC3")], 2, mean )))

    m <- (ggplot(projected,aes(x=PC1, y=PC2)) + 
          geom_point(aes(colour=region))+
          labs(title = "Projected PCA") +
          geom_text(data=cnames, aes(PC1, PC2, label = region), size=10) +
          #geom_point(data=centroid,aes(PC1,PC2), colour="black",size=15)+
          theme(legend.position="none") +
          ng1 )

    outfilename <- paste(outputdir,fileprefix,"_ancientfull.png",sep="")
    print(paste("Writing file:",outfilename))
    png(outfilename)
    plot(m)
    dev.off()

    projected <- as.data.frame(x)
    projected$member <- rownames(projected)
    projected <- merge(projected, snpdata$annot, on="member")
    projected$region <- factor(projected$region, levels=sort(unique(projected$region)))

    cnames <- aggregate(cbind(PC1, PC2) ~ region, data=projected, 
                                    FUN=function(x)mean(range(x)))

    m <- (ggplot(projected,aes(x=PC1, y=PC2)) + geom_point(aes(colour=region)) +
          geom_text(data=cnames, aes(PC1, PC2, label = region), size=3) +
          theme(legend.position="none") +
          labs(title = "Projected PCA") + ng1 )

    outfilename <- paste(outputdir,fileprefix,"_ancientzoom.png",sep="")
    print(paste("Writing file:",outfilename))
    png(outfilename)
    plot(m)
    dev.off()
}


######################################################################
# runPCAnalysis
######################################################################
runPCAnalysis <- function( snpdata, outputdir, fileprefix="" ) {
    print("Run function runPCAnalysis")

    geno.mat <- snpdata$genotypes[rownames(snpdata$genotypes) %in% snpdata$annot$member,]
    #plink.annot <- plink.annot.full#[plink.annot.full$member %in% plink.annot.full$member,]

    s <- plink.full$fam$member
    t <- snpdata$annot$member
    s[!(s %in% t)]

    ###################################################
    ### code chunk number 2: xxt-matrix
    ###################################################
    xxmat <- xxt(geno.mat, correct.for.missing=FALSE)

    ###################################################
    ### code chunk number 3: eigen
    ###################################################
    evv <- eigen(xxmat, symmetric=TRUE)
    numvectors <- 20
    pcs <- evv$vectors[,1:numvectors]
    pcs.cols <- list()
    #for( i in 1:numvectors ){ pcs.cols[i] = (paste("PC",i,sep="")) }
    #colnames(pcs) <- c("PC1","PC2","PC3","PC4","PC5")
    colnames(pcs) <- paste("PC",1:numvectors,sep="")
    rownames(pcs) <- snpdata$fam$member
    #evals <- evv$values[1:24]
    pcs <- cbind(snpdata$annot,as.data.frame(pcs[,1:5]))

    #all.pcs.filt <- all.pcs[all.pcs$Country != "Unknown",]
    #all.pcs.filt <- all.pcs
    #all.pcs.filt$Origin <- factor(all.pcs.filt$Origin, levels = sort(unique(all.pcs.filt$Origin)))
    #head(all.pcs.filt$Origin)

    outfilename <- paste(outputdir,fileprefix,"_pcs.txt",sep="")
    print(paste("Writing file:",outfilename))
    write.table(pcs,file=outfilename,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
    # Compute proportion of the total variance explained by the first eigenvalues
    #cumprop <- cumsum(evv$values)/sum(evv$values)
    #pc.sds <- apply(pcs,2,sd)
    #variances <- data.frame(variances=pc.sds**2, pcomp=1:length(pc.sds))
    #variances <- data.frame(variance=cumprop[1:25], pcomp=1:25)
    #outfilename <- paste(outputdir,fileprefix,"pcavariances.png",sep="")
    #print(paste("Writing file:",outfilename))
    #png(outfilename)
    #varPlot <- ggplot(variances, aes(pcomp, variance)) + geom_bar(stat="identity", 
            #fill="grey") + geom_line()
    #print(varPlot)
    #dev.off()

    # Plot Principle Components
    outfilename <- paste(outputdir,fileprefix,"_pcacontinent.png",sep="")
    print(paste("Writing file:",outfilename))
    png(outfilename)
    m <- ggplot(pcs,aes(x=PC1, y=PC2)) + geom_point(aes(colour=region))
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

readHgdpData <- function(){
    filepath <- "/media/data/workspace/variome/resources/hgdppanel4/"
    plinkprefix <- "panel4"
    annotfile <- "/media/data/workspace/variome/resources/hgdppanel4/sample_panel4.txt"

    bedfile <- paste(filepath,plinkprefix,".bed",sep="")
    bimfile <- paste(filepath,plinkprefix,".bim",sep="")
    famfile <- paste(filepath,plinkprefix,".fam",sep="")

    annot <- read.table(annotfile,sep="\t",header=F, 
                        col.names=c("member","gender","region"))
    row.names(annot) <- annot$member

    hgdpdata <- read.plink(bedfile,bimfile,famfile)
    hgdpdata$annot <- annot
    return(hgdpdata)
}
######################################################################
# Main
######################################################################

#bedprefix <- "plinkfiles/merged.1kg.percentcoverage.plink.mod.filt.uniq"
#bedprefix <- "plinkfiles/merged.eichler.plink.mod.filt.uniq"
#bedprefix <- "../data/merged/ciliopathies.unfilt.recode.plink.mod.filt.uniq"
#bedprefix <- "../rawdata/daily1/daily.recode.plink.recode12.mod.filt.uniq"
#bedprefix <- "../data/collective/ciliopathies.unfilt.recode.plink.recode12.mod.filt.uniq"
#fileprefix <- "collective"
# parse args
#args <- commandArgs( trailingOnly = TRUE )
#bedprefix <- ifelse(length(args)>=1,dirname(args[1]),bedprefix)
#fileprefix <- ifelse(length(args)>=2,dirname(args[2]),fileprefix)


#setwd( "/home/escott/workspace/variome/rscripts/" )
outputdir <- "/home/escott/workspace/variome/rscripts/pcafigures/"

bedprefix <- "../rawdata/test/everything_set1.chr1.snp.chimp.regions.recode.recode.plink.mod.filt.uniq"
fileprefix <- "chr1test"
#bedprefix <- "../rawdata/variome/variome_snps.chimp.regions.recode.recode.plink.mod.filt.uniq"
bedprefix <- "../rawdata/variome1/variome.chimp.regions.filt.samp.plink.mod.filt.uniq"
fileprefix <- "variome"
annotfile <- "../resources/annotation/patientannotation.ped"

hgdpdata <- readHgdpData()

ancientPCA( hgdpdata, outputdir, fileprefix="hgdp" )

#s <- runPCAnalysis( bedprefix, annotfile, outputdir, fileprefix )
#all.pcs.filt <- s[[1]]
#pcs <- s[[2]]

q()
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
