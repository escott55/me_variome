#!/usr/bin/env Rscript


setwd( "/home/escott/workspace/variome/rscripts/" )
print(paste("Working dir:",getwd()))

#source("http://bioconductor.org/biocLite.R")
#biocLite("snpStats")
#install.packages(c("gPCA","reshape","ggplot2"))

library(RColorBrewer, quietly=TRUE)
palette(brewer.pal(8, "Dark2"))
#library(gplots, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(reshape, quietly=TRUE)
#library(MultiPhen)
library(snpStats, quietly=TRUE)

library( gPCA , quietly=TRUE)

testgpca <- function( ){
    data(caseDat)
    batch<-caseDat$batch
    data<-caseDat$data
    out<-gPCA.batchdetect(x=data,batch=batch,center=FALSE,nperm=250)
    print( out$delta )
    print( out$p.val )

    gDist(out)
    png("batchfigs/test_cumvar.png")
    CumulativeVarPlot(out,ug="unguided",col="blue")
    dev.off()
    png("batchfigs/test_1v2.png")
    PCplot(out,ug="unguided",type="1v2")
    dev.off()
    png("batchfigs/test_comp.png")
    PCplot(out,ug="unguided",type="comp",npcs=4)
    dev.off()
}

#testgpca()

testsnpstats <- function( ){
    data(for.exercise)
    show(snps.10)
    data <- apply(snps.10@.Data,2,as.numeric)
    batch <- sample(1:5, dim(data)[1], replace=T)
    out<-gPCA.batchdetect(x=data,batch=batch,center=FALSE,nperm=250)
    print(out$delta)
    print(out$p.val)

    gDist(out)
    png("batchfigs/test1_cumvar.png")
    CumulativeVarPlot(out,ug="unguided",col="blue")
    dev.off()
    png("batchfigs/test1_1v2.png")
    PCplot(out,ug="unguided",type="1v2")
    dev.off()
    png("batchfigs/test1_comp.png")
    PCplot(out,ug="unguided",type="comp",npcs=4)
    dev.off()
}

#testsnpstats()


#annotfile <- "../patientannot/annotation/patientannotation_laser.txt"
annotfile <- "../patientannot/annotation/patientannotation.ped"
allpat.annot <- read.table(annotfile, header=TRUE,sep="\t",comment.char="")

#bedprefix <- "../data/daily1/daily.recode.plink.recode12.mod.filt.uniq"
bedprefix = "/home/escott/workspace/inbreed/rawdata/turks/turks.names.chimp.regions.recode.recode.plink.mod.filt.uniq"

bedfile <- paste(bedprefix,".bed",sep="")
bimfile <- paste(bedprefix,".bim",sep="")
famfile <- paste(bedprefix,".fam",sep="")

plink.full <- read.plink(bedfile,bimfile,famfile)
plink.annot.full <- merge( plink.full$fam, allpat.annot, by.x="member", by.y="Individual.ID", sort=FALSE )
plink.annot.full$Source <- factor(plink.annot.full$Source, levels=unique(plink.annot.full$Source))
plink.annot.full$ethnicity <- factor(plink.annot.full$ethnicity, levels=unique(plink.annot.full$ethnicity))

geno.mat <- plink.full$genotypes[rownames(plink.full$genotypes) %in% plink.annot.full$member,]

data <- apply(geno.mat@.Data,2,as.numeric)
#batch <- sample(1:5, dim(data)[1], replace=T)
batch <- plink.annot.full$Plate

plink.annot.full <- transform( plink.annot.full, 
          batch  = as.integer(factor(Source, unique(as.character(Source)))) - 1L)
#batch <- as.vector(plink.annot.full$Source)
#batch <- data.matrix(as.data.frame(plink.annot.full$Source))
batch <- plink.annot.full$batch
out<-gPCA.batchdetect(x=data,batch=batch,center=FALSE,nperm=250)
print(out$delta)
print(out$p.val)

gDist(out)
png("batchfigs/forreal_cumvar.png")
CumulativeVarPlot(out,ug="unguided",col="blue")
dev.off()
png("batchfigs/forreal_1v2.png")
PCplot(out,ug="unguided",type="1v2")
dev.off()
png("batchfigs/forreal_comp.png")
PCplot(out,ug="unguided",type="comp",npcs=3)
dev.off()


cplate <- "PlateI"
mydata <- as.data.frame(cbind( sample=plink.annot.full$member, ethnicity=plink.annot.full$ethnicity, as.data.frame(data[,1:10]) ))

Y <- as.data.frame( data[,1:1000] )
Ethnicity <- plink.annot.full$ethnicity
Source <- plink.annot.full$Source
Samples <- plink.annot.full$member

#fit <- aov( mydata$rs8997 ~ Source + Ethnicity )
#summary(fit)[[1]][["Pr(>F)"]][1]

fit <- aov( data[,1] ~ Source + Ethnicity )
summary(fit)[[1]][["Pr(>F)"]][1]

fit <- kruskal.test( data[,1] ~ Source + Ethnicity )
summary(fit)[[1]][["Pr(>F)"]][1]

test <- function( x ) {
    fit <- aov( x ~ Source + Ethnicity )
    return(summary(fit)[[1]][["Pr(>F)"]][1:2])
}

#varfit <- as.data.frame(t(apply( Y, 2, test )))
#varfit.fdr <- cbind( p.adjust(varfit[,0]), p.adjust(varfit[,1]) )
#colnames(varfit.fdr) <- c("Source","Ethnicity")

#snpstocheck <- varfit.fdr[varfit.fdr$Source < .05 & varfit.fdr$Ethnicity >= .05,]
#dim(snpstocheck)

test1 <- function( x ) {
    fit <- kruskal.test( x ~ Source + Ethnicity )
    return(fit$p.value)
}

varfit <- t(apply( Y, 2, test1 ))
varfit.fdr <- as.data.frame(colnames(varfit))
colnames(varfit.fdr) <- c("Var")
varfit.fdr$Pvalue <- p.adjust(varfit)
#varfit.fdr$Pvalue <- as.numeric(as.character(varfit.fdr$Pvalue))
#colnames(varfit.fdr) <- c("Pvalue")
#varfit.fdr$var <- rownames(varfit)

snpstocheck <- varfit.fdr[varfit.fdr$Pvalue < .05,]
length(snpstocheck)

mysnp <- "rs74980814"
mysnp <- "rs1801274"
mysnp <- "1.22156531"
mydata <- cbind( plink.annot.full[,colnames(plink.annot.full) %in% c("member","ethnicity","Source")], snp=as.data.frame(data)[[mysnp]] )
table(mydata$ethnicity,mydata$snp)
s <- table(mydata$Source,mydata$snp)


s <- s/apply(s,1,sum)

#inplate <- ifelse( plink.annot.full$Source == cplate, 0, 1)
#mydata <- as.data.frame(cbind( inplate, ethnicity=plink.annot.full$ethnicity, as.data.frame(data[,1:10]) ))
#x <- TukeyHSD( fit )

#fit <- manova( as.matrix(Y) ~ Ethnicity + Source )
#summary(fit)
#
#summary(mod <- Anova(lm(as.matrix(Y) ~ 1),
        #type=3, 
        #idesign=~Ethnicity+Source),
    #multivariate=FALSE)
#
#require(heplots)
#etasq(mod, anova = TRUE)
#
#fit1 <- lm( as.matrix(Y) ~ mydata$ethnicity + mydata$inplate )
#fit2 <- lm( as.matrix(Y) ~ mydata$inplate )
#
#summary( fit1) 
#summary( fit2) 
#
#anova(fit1,fit2,test="Wilks")
#
#require(car)
#fit3 <- Anova(lm(as.matrix(Y) ~ mydata$ethnicity + mydata$inplate))
#summary( fit3) 

# END
