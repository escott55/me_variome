#!/usr/bin/env Rscript

library(Geneland)
library(ggplot2)
library(reshape)

setwd( "/home/escott/workspace/variome/rscripts/" )
print(paste("Working dir:",getwd()))

#genotypes <- read.table("../rawdata/variome1/ibc/variome.regions.filt.samp.samp.Middle_East.recode12.mod.exclude.ped",sep=" ",header=FALSE)
#genotypes <- read.table("./rawdata/variome1/ibc/variome.regions.filt.samp.samp.plink.recode12.mod.exclude.ped",sep=" ",header=FALSE)
#tfile <- "../rawdata/merged/ibc/old1/merged.chimp.regions.filt.samp.samp.plink.recode12.ped"
#tfile <- "../rawdata/mevariome/main/ibc/variome.clean.Middle_East.recode12.ped"
tfile <- "../rawdata/mevariome/main/ibc/variome.clean.Middle_East.recode12.mod.exclude.ped"
genotypes <- read.table(tfile,sep=" ",header=FALSE)
annot <- read.table("../resources/annotation/patientannotation.ped",sep="\t",header=TRUE)

tsamp <- genotypes$V2
tfam <- genotypes$V1
rownames(genotypes) <- genotypes$V2

annot_filt <- annot[annot$Individual.ID %in% tsamp,]

imatch <- cbind(match(tsamp, annot$Individual.ID),match(tsamp, annot$Family.ID))
annot_filt <- annot[apply(imatch,1, function(x){ min(x,na.rm=TRUE) }),]

genotypes <- genotypes[,!(colnames(genotypes) %in% c("V1","V2","V3","V4","V5","V6"))]

print(paste("Genotypes Shape",ncol(genotypes),"x",nrow(genotypes)))
if( ncol(genotypes) > 5000 ){
    genotypes <- genotypes[,1:5000]
}
print(paste("Genotypes Shape",ncol(genotypes),"x",nrow(genotypes)))

#ethtargets <- c("Bedouin","Palestinian","Druze","Adygei","Mozabite")
#annot_filt <- annot_filt[annot_filt$ethnicity %in% ethtargets,]
#genotypes <- genotypes[annot_filt$ethnicity %in% ethtargets,]
#
#s <- FormatGenotypes(genotypes,2)
#npop <-  length(unique(annot_filt$ethnicity))
#mbrship <- factor(annot_filt$ethnicity,levels=unique(annot_filt$ethnicity))

regiontargets <- c("Bedouin","Palestinian","Druze","Adygei","Mozabite")
regiontargets <- c("Algeria","Egypt","Iran","Iraq","Jordan","Kuwait",
                   "Lebanon","Libya","Morocco","North Africa","Oman",
                   "Palestine","Qatar","Saudi Arabia","Syria","Tunisia","Turkey","UAE","Yemen")
annot_filt <- annot_filt[annot_filt$Origin %in% regiontargets,]
annot_filt$Origin <- factor( annot_filt$Origin, levels=unique(annot_filt$Origin))
genotypes <- genotypes[annot_filt$Origin %in% regiontargets,]

s <- FormatGenotypes(genotypes,2)
npop <-  length(unique(annot_filt$Origin))
mbrship <- factor(annot_filt$Origin,levels=unique(annot_filt$Origin))

x <- Fstat(s$genotypes,npop,as.integer(mbrship),2)

fst <- x$Fst
rownames(fst) <- levels(mbrship)
colnames(fst) <- levels(mbrship)

fst.m <- as.data.frame(melt(fst[rownames(fst) %in% regiontargets, colnames(fst) %in% regiontargets]))
(p <- ggplot(fst.m, aes(X1, X2)) + geom_tile(aes(fill = value),
       colour = "white") + scale_fill_gradient(low = "white", high = "steelblue"))

png( "fst/origin_fst_heat2.png", width=5, height=3, units="in", res=300)
p
dev.off()

#base_size <- 9
#p + theme_grey(base_size = base_size) + labs(x = "",
       #y = "") + scale_x_discrete(expand = c(0, 0)) +
       #scale_y_discrete(expand = c(0, 0)) + opts(legend.position = "none",
       #axis.ticks = theme_blank(), axis.text.x = theme_text(size = base_size *
       #0.8, angle = 330, hjust = 0, colour = "grey50"))


#library(adegenet)
#rawfile  <- "/home/escott/workspace/consang/data/collective/plink.raw"
#s <- read.PLINK( rawfile )
#
#pca1 <- glPca(s)
#png("figures/scatter.png")
#scatter(pca1, posi="bottomright")
#title("PCA of the US influenza data\n axes 1-2")
#dev.off()
#
#library(ape)
#tre <- nj(dist(as.matrix(s)))
#
#png("figures/njtree.png")
#plot(tre, typ="fan", cex=0.7)
#title("NJ tree of the US influenza data")
#dev.off()
#
#png("figures/njscatter.png")
#myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4)
#abline(h=0,v=0, col="grey")
#add.scatter.eig(pca1$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)
#dev.off()

#png("figures/njtree2.png")
#plot(tre, typ="fan", show.tip=FALSE)
#tiplabels(pch=20, col=myCol, cex=4)
#title("NJ tree of the US influenza data")
#dev.off()

#library(PopGenome)
#vcffile <- "rawdata/ciliopathies/ciliopathies.unfilt.vcf.gz"
#GENOME.class <- readData( "rawdata/ciliopathies", format="VCF" )
#
#s@region.names
#s <- neutrality.stats(GENOME.class,FAST=TRUE)
#get.sum.data( s )



