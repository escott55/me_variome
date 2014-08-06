#!/usr/bin/env Rscript

library(RColorBrewer)
palette(brewer.pal(8, "Dark2"))
library(gplots)
library(ggplot2)
library(reshape)
library(VennDiagram)

######################################################################
# Theme
######################################################################
immigration_theme <- theme(
    panel.grid.major = element_line(colour = "grey80"),
    panel.grid.minor = element_blank(), 
    #panel.background = theme_blank(),
    panel.background = element_rect(fill="grey90",colour="white"), 
    plot.title = element_text(lineheight=.8, face="bold", size=20,colour="black"),
    axis.text.x = element_text(colour="black",size=13,angle=45),
    axis.text.y = element_text(colour="black",size=13),
    axis.title.x = element_text(colour="black",size=20),
    axis.title.y = element_text(colour="black",size=20,angle=90)
    #axis.ticks = element_blank()
    #legend.position = "none"
    )

ng1 = theme(panel.background = element_rect(fill = "white",colour = "white"),
    panel.grid.major = element_line(colour = "grey90"),
    axis.line = element_line(size = 1.2, colour="black"),
    axis.ticks=element_line(color="black"),
    axis.text=element_text(color="black",size=15),
    axis.title=element_text(color="black",size=20),
    plot.title = element_text(face="bold", size=20,colour="black"),
    panel.grid.minor = element_line(colour = NA),
    #legend.position = "top", 
    #legend.direction="horizontal", 
    legend.text = element_text(size=13),
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(size=15)
    )

ng2 = theme(panel.background = element_rect(fill = "grey90",colour = "white"),
    panel.grid.major = element_line(colour = "white"),
    axis.line = element_line(size = 1.2, colour="black"),
    axis.ticks=element_line(color="black"),
    axis.text=element_text(color="black",size=15),
    axis.text.x = element_text(colour="black",angle=45,hjust=1,vjust=1),
    axis.title=element_text(color="black",size=20),
    plot.title = element_text(face="bold", size=20,colour="black"),
    panel.grid.minor = element_line(colour = NA),
    #legend.position = "top", 
    #legend.direction="horizontal", 
    legend.text = element_text(size=13),
    legend.key = element_rect(fill = "white"),
    legend.title = element_blank(),
    strip.text.y = element_text(size=13,face="bold"),
    strip.text.x = element_text(size=13,face="bold")
    )


######################################################################
# removeSingletons
######################################################################
removeSingletons <- function( cdata ){
    annotationfile <- "../processeddata/annotation/patientannotation_laser.txt"
    patannot <- read.table(annotationfile, header=TRUE,sep="\t",comment.char="")
    newcdata <- NULL
    for( plate in unique(patannot$Source) ){
        patlist <- patannot[patannot$Source == plate,]$Individual.ID
        cdata.sub <- cdata[cdata$Patient %in% patlist,]
        print(paste("cdata.sub:",dim(cdata.sub)))
        if( length( cdata.sub$Patient ) == 0 ){ next; }
        print( paste("Num patients:",length(cdata.sub$Patient) ))
        minaf <- min(cdata.sub$AF)
        print(names(cdata.sub))
        print(head(cdata.sub))
        print(paste("Min AF:",plate,minaf))
        dim(cdata.sub)
        dim(cdata.sub[cdata.sub$AF > minaf,])
        newcdata <- rbind(newcdata,cdata.sub[cdata.sub$AF > minaf,])
    }
    return(newcdata)
} #END removeSingletons


######################################################################
# recalculateAF
######################################################################
recalculateAF <- function( cdata ){
    varid <- paste(cdata$chrom,cdata$pos,cdata$ref,sep=":")
    cdata$varid <- varid
    counts <- t(table(cdata[,colnames(cdata) %in% c("varid","Genotype")]))
    numpats <- length(unique(cdata$Patient))
    newAF <- as.data.frame(cbind(varid=rownames(counts),AF=((counts[,1]+2*counts[,2]) / numpats)))
    newAF$AF <- as.numeric(as.character(newAF$AF))
    print(paste("Min:",min(newAF$AF),"Max:",max(newAF$AF)))
    #colnames(newAF) <- c("varid","newAF")
    print("right before")
    cdata.new <- merge(cdata[,!(colnames(cdata) %in% c("AF","X.HomShares","X.HetShares","scorePhastCons","consScoreGERP"))], newAF,by="varid")
    print("right after")
    print(dim(cdata.new))
    cdata.new <- cdata.new[cdata.new$AF != 0,]
    print(dim(cdata.new))
    return(cdata.new)
} # END recalculateAF

######################################################################
# recalculateAFPops
######################################################################
recalculateAFPops <- function( cdata, patset1, patset2 ){
    print(" - Running recalculateAFPops")
    print(dim(cdata))
    print(head(cdata))
    cdata.sub <- cdata[cdata$Patient %in% patset1,]
    varid <- paste(cdata.sub$chrom,cdata.sub$pos,cdata.sub$ref,sep=":")
    cdata.sub$varid <- varid
    counts <- t(table(cdata.sub[,colnames(cdata.sub) %in% c("varid","Genotype")]))
    numpats <- length(unique(cdata.sub$Patient))
    newAF <- as.data.frame(cbind(varid=rownames(counts),AF=((counts[,1]+2*counts[,2]) / numpats)))
    newAF$AF <- as.numeric(as.character(newAF$AF))
    print(paste("Min:",min(newAF$AF),"Max:",max(newAF$AF)))
    #colnames(newAF) <- c("varid","newAF")
    print("right before")
    cdata.1 <- merge(cdata.sub[,!(colnames(cdata.sub) %in% c("AF","X.HomShares","X.HetShares","scorePhastCons","consScoreGERP"))], newAF,by="varid")
    print("right after")
    print(dim(cdata.1))
    cdata.1 <- cdata.1[cdata.1$AF != 0,]
    print(dim(cdata.1))

    cdata.sub <- cdata[cdata$Patient %in% patset2,]
    varid <- paste(cdata.sub$chrom,cdata.sub$pos,cdata.sub$ref,sep=":")
    cdata.sub$varid <- varid
    counts <- t(table(cdata.sub[,colnames(cdata.sub) %in% c("varid","Genotype")]))
    numpats <- length(unique(cdata.sub$Patient))
    newAF <- as.data.frame(cbind(varid=rownames(counts),AF=((counts[,1]+2*counts[,2]) / numpats)))
    newAF$AF <- as.numeric(as.character(newAF$AF))
    print(paste("Min:",min(newAF$AF),"Max:",max(newAF$AF)))
    #colnames(newAF) <- c("varid","newAF")
    print("right before")
    cdata.2 <- merge(cdata.sub[,!(colnames(cdata.sub) %in% c("AF","X.HomShares","X.HetShares","scorePhastCons","consScoreGERP"))], newAF,by="varid")
    print("right after")
    print(dim(cdata.2))
    cdata.2 <- cdata.2[cdata.2$AF != 0,]
    print(dim(cdata.2))

    cdata.new <- rbind( cdata.1, cdata.2 )

    return(cdata.new)
} # END recalculateAFPops

#makeVarCountsGraph <- function( tfile, figuredir ){
#} # End makeVarCountsGraph

 # readClassData
   #dev.off()

 # End inbredOutbredRatios

   #dev.off()

 # End groupplot
 # END affectedplot
 # END inbredOutbredBurden

######################################################################
# getPatlist
######################################################################
getPatlist <- function( patannot, psource ){
    patlist <- NULL
    print(head(patannot))
    annot <- patannot[patannot$Unrelated == "True",]
    #print(head(annot))
    print( paste("Psource:",psource) )
    if( psource == "internal" ){
        patlist <- as.vector(annot[annot$Lab == "Broad" & annot$continent == "Middle East" , ]$Individual.ID)
        #& annot$Phenotype == 0
    #}else if( psource == "eichler" ){
        #patlist <- as.vector(annot[annot$Lab == "Eichler" & annot$continent == "Europe",]$Individual.ID)
    #}else if(psource == "onekg" ){
        #patlist <- as.vector(annot[annot$Lab == "1KG",]$Individual.ID)
    #}else if(psource == "daily" ){
        #patlist <- as.vector(annot[annot$Lab == "Daly" & annot$continent == "Europe",]$Individual.ID)
    }
    print( paste("Number of patients:",length(patlist)) )
    return(patlist)
} # END getPatlist

 # End variantclass_analysis
 # END allVariantAnalysis 

######################################################################
# tabulateGeneBurden
######################################################################
tabulateGeneBurden <- function( tfile, patlisttype ){

    patset1 <- getPatlist( patannot, patlisttype )

    #tfile <- paste(targetdir,"gaiadata_class",varclass,"_all.tsv",sep="")
    print(paste("Reading file:",tfile))
    cdata <- read.table(tfile, header=TRUE, sep="\t", comment.char = "" )
    originalpatients <- cdata$Patient
    #print(dim(cdata))

    #keeppatients <- c(patset1, patset2)
    keeppatients <- patset1
    print(paste("Keep Patients:",patlisttype))
    print(length(keeppatients))
    if( !(is.null(keeppatients)) ){
        cdata <- cdata[cdata$Patient %in% keeppatients,]
        cdata$Patient <- factor( cdata$Patient, levels=unique(cdata$Patient) )
    }else{
        print( "Keep patients is null" )
    }
    print(dim(cdata))
    print("done")

    homs <- cdata[cdata$Genotype == "hom",]
    hets <- cdata[cdata$Genotype == "het",]

    hetburden <- as.data.frame(table(hets[,colnames(hets) %in% c("Patient","gene")]))
    #chets <- as.data.frame( table( hetburden$gene[hetburden$Freq > 1] ) )
    #kos <- as.data.frame( table( homs$gene ) )
    chets <-  table( hetburden$gene[hetburden$Freq > 1] )
    kos <- table( homs$gene )
    allkos <- as.data.frame( chets + kos )
    return( list( chets=as.data.frame(chets),kos=as.data.frame(kos),allkos=allkos ) )
} # END tabulateGeneBurden 

######################################################################
# knockoutAnalysis
######################################################################
knockoutAnalysis <- function( targetdir, figuredir, varclasstype, varclass ){
    #tfile <- paste(targetdir,varclasstype,"_class",varclass,"_all.tsv",sep="")
    tfile <- paste(targetdir,varclasstype,"_",varclass,".tsv",sep="")
   
    tmp <- tabulateGeneBurden( tfile, "internal" )
    set1.chets <- tmp$chets
    set1.kos <- tmp$kos
    set1.allkos <- tmp$allkos

    tmp <- tabulateGeneBurden( tfile, "daily" )
    set2.chets <- tmp$chets
    set2.kos <- tmp$kos
    set2.allkos <- tmp$allkos

    tmp <- tabulateGeneBurden( tfile, "eichler" )
    set3.chets <- tmp$chets
    set3.kos <- tmp$kos
    set3.allkos <- tmp$allkos

    venn.plot <- venn.diagram(list(ME = set1.chets$Var1[set1.chets$Freq > 0], 
                      CEU = set2.chets$Var1[set2.chets$Freq > 0]),
              fill = c("red", "green"),
              main=paste("Compound Het Class",varclass,"Variants"),
              sub="Middle Easterns Versus Europeans",
              main.cex=3,sub.cex=2.3,cat.cex=2.1,
              alpha = c(0.5, 0.5), cex = 2, cat.fontface = 4,lty =2, fontfamily =3, 
              filename=NULL)
              #filename = "figures/kos/chets.png");
    png(paste("figures/kos/chets_class",varclass,".png"));
    grid.draw(venn.plot);
    dev.off();

    venn.plot <- venn.diagram(list(ME = set1.kos$Var1[set1.kos$Freq > 0], 
                      CEU = set2.kos$Var1[set2.kos$Freq > 0]),
              fill = c("red", "green"),
              main=paste("Homozygous Class",varclass,"Variants"),
              sub="Middle Easterns Versus Europeans",
              main.cex=3,sub.cex=2.3,cat.cex=2.1,
              alpha = c(0.5, 0.5), cex = 2, cat.fontface = 4,lty =2, fontfamily =3, 
              filename=NULL)
              #filename = "figures/kos/kos.png");
    png(paste("figures/kos/kos_class",varclass,".png"));
    grid.draw(venn.plot);
    dev.off();

    venn.plot <- venn.diagram(list(ME = set1.allkos$Var1[set1.allkos$Freq > 0], 
                      CEU = set2.allkos$Var1[set2.allkos$Freq > 0]),
              fill = c("red", "green"),
              main=paste("All Class",varclass,"Variants"),
              sub="Middle Easterns Versus Europeans",
              main.cex=3,sub.cex=2.3,cat.cex=2.1,
              alpha = c(0.5, 0.5), cex = 2, cat.fontface = 4,lty =2, fontfamily =3, 
              filename=NULL)
              #filename = "figures/kos/allkos.png");
    png(paste("figures/kos/allkos_class",varclass,".png"));
    grid.draw(venn.plot);
    dev.off();

    venn.plot <- venn.diagram(list(JG = set1.allkos$Var1[set1.allkos$Freq > 0], 
                      EE = set2.allkos$Var1[set2.allkos$Freq > 0],
                      MD = set3.allkos$Var1[set3.allkos$Freq > 0]),
              fill = c("red", "green", "blue"),
              main=paste("All Class",varclass,"Variants"),
              sub="Lab Comparison",
              main.cex=3,sub.cex=2.3,cat.cex=2.1,
              alpha = c(0.5, 0.5, 0.5), cex = 2, cat.fontface = 4,lty =2, fontfamily =3, 
              filename=NULL)
              #filename = "figures/kos/allkos.png");
    png(paste("figures/kos/allkos_alllabs_class",varclass,".png"));
    grid.draw(venn.plot);
    dev.off();


} # END knockoutAnalysis

######################################################################
# Main
######################################################################
#varclasstype <- "gaiadata"
#varclasstype <- "rawdatatenp"
#varclasstype <- "eichler"
#varclasstype <- "daily"
#varclasstype <- "rawdata"
setwd("/home/escott/workspace/variome/")
varclasstype <- "variome_snps.chimp.regions.recode.annot_"
figuredir <- paste("results/affunaff/",varclasstype,"/",sep="")
dir.create(file.path(figuredir), showWarnings = FALSE)

targetdir <- "classdata/"
#removepatients <- read.table("toremove.txt",sep="\t")$V1

genelist <- read.table("resources/allgene_classes.tsv",header=FALSE,sep="\t",comment.char="")
#patannot <- read.table("patientannotation.tsv",header=TRUE,sep="\t",comment.char="")
annotationfile <- "../processeddata/annotation/patientannotation.ped"
patannot <- read.table(annotationfile,header=TRUE,sep="\t",comment.char="")
patannot$Affected <- ifelse( patannot$Phenotype == 1, "Affected", "Unaffected" )


varclass <- "HIGH"
knockoutAnalysis( targetdir, figuredir, varclasstype, varclass )

## Test readClassData
#varclass <- 1
#patset1 <- getPatlist( patannot, "eichler")
#patset2 <- getPatlist( patannot, "internal")
#keeppatients <- c(patset1, patset2)
#
#tmp <- readClassData( varclass, targetdir, figuredir, varclasstype, keeppatients )


# END Main


