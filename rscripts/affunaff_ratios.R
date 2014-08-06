#!/usr/bin/env Rscript

library(RColorBrewer)
palette(brewer.pal(8, "Dark2"))
#library(gplots)
library(ggplot2)
library(reshape)

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
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
######################################################################
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
} # End multiplot


######################################################################
# removeSingletons
######################################################################
removeSingletons <- function( cdata ){
    annotationfile <- "patientannot/annotation/patientannotation_laser.txt"
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

######################################################################
# readClassData
######################################################################
readClassData <- function( tfile, figuredir, patset1=NULL, patset2=NULL, 
                          genelist=NULL, geneclass=NULL ){
    #tfile <- paste(targetdir,varclasstype,"_class",varclass,"_all.tsv",sep="")

    print(paste("Reading file:",tfile))
    cdata <- read.table(tfile, header=TRUE, sep="\t", comment.char = "" )
    originalpatients <- cdata$Patient
    #print(dim(cdata))

    keeppatients <- c(patset1, patset2)
    print("Keep Patients:")
    print(keeppatients)
    if( !(is.null(keeppatients)) ){
        cdata <- cdata[cdata$Patient %in% keeppatients,]
        cdata$Patient <- factor( cdata$Patient, levels=unique(cdata$Patient) )
    }else{
        print( "Keep patients is null" )
    }
    print(dim(cdata))
    print("done")

    keepgenes <- list()
    if( !is.null( geneclass) ){
        keepgenes <- genelist[genelist$Class == geneclass,]$Gene
    }else if(!is.null(genelist)){
        keepgenes <- genelist$Gene
    }
    if( length(keepgenes) > 0 ){
        print(length(keepgenes))
        cdata <- cdata[cdata$gene %in% keepgenes,]
    }
    cdata.old <- cdata
    print(head(cdata))
    #cdata <- recalculateAF( cdata )
    cdata <- recalculateAFPops( cdata, patset1, patset2 )
    #cdata <- removeSingletons( cdata )
    #cdata <- cdata[cdata$AF > min(cdata$AF),]
    print(head(cdata))


    Het <- tapply(cdata$Genotype,cdata$Patient,function(x){ sum(x=="het")})
    Hom <- tapply(cdata$Genotype,cdata$Patient,function(x){ sum(x=="hom")})

    print(paste("Hets:",sum(Het),"Homs",sum(Hom)))
    x <- as.data.frame(table(cdata[,colnames(cdata) %in% c("Patient","Genotype")]))
    x$Freq <- as.numeric(x$Freq)
    #print(head(x))
    #print(dim(x))
    homhet <- cast(x , Patient~Genotype, value='Freq',fill=0)

    #homhet <- cbind(Patient=as.character(sort(unique(cdata$Patient))),
               #Het=tapply(cdata$Genotype,cdata$Patient,function(x){ sum(x=="het")}),
               #Hom=tapply(cdata$Genotype,cdata$Patient,function(x){ sum(x=="hom")}))
    print("up next")
    #homhet <- as.data.frame(homhet)
    homhet$hom <- as.numeric(homhet$hom)
    homhet$het <- as.numeric(homhet$het)
    #print(head(homhet))
    homhet[["homhet"]] <- homhet[["hom"]]/homhet[["het"]]
    print(head(homhet))

    print("TiTv ratio")
    transitions <- c("AG","GA","CT","TC")
    transversion <- c("AT","TA", "GC", "CG","GT","TG","AC","CA")
    cdata$combined <- apply(cdata,1,function(x){paste(x[["ref"]],x[["mut"]],sep="")})
    print(unique(cdata$combined))
    print(dim(cdata))
    #print(table(cdata$Patient))
    print(head(tapply( cdata$combined, cdata$Patient, function(x){ sum(x %in% transitions)})))
    titv <- cbind( Patient=as.character(sort(unique(cdata$Patient))),
           TI=tapply( cdata$combined, cdata$Patient, function(x){ sum(x %in% transitions)}),
           TV=tapply( cdata$combined, cdata$Patient, function(x){ sum(x %in% transversion)}))
    print("afterwards?")
    titv <- as.data.frame( titv )
    titv$TI <- as.numeric(titv$TI)
    titv$TV <- as.numeric(titv$TV)
    titv[["titv"]] <- titv[["TI"]]/titv[["TV"]]
    patientstats <- merge( homhet, titv, by="Patient" )
    print(head(patientstats))
    patientstats$Burden <- patientstats$het + 2*patientstats$hom
    print(head(patientstats))
    #patientstats$Burden <- patientstats$het + patientstats$hom
    if( !is.null(geneclass) ){
        patientstats$Type <- geneclass
    }
    #print(head(patientstats))

    afcuts <- NULL
    #cdata.hets <- cdata[cdata$Genotype == "het",]
    #cdata.homs <- cdata[cdata$Genotype == "hom",]
    print( "Calculate Cutoffs" )

    smallestaf <- min( 1/(2*length(patset1)), 1/(2*length(patset2)) )
    for( cutoff in c(1:200/2000) ){
        if( cutoff < smallestaf ){ next }
        tmp <- tapply( cdata$AF, cdata$Patient, function(x){ sum(x<cutoff) })
        # Multiply homozygous by 2
        #tmp <- apply(table(cdata[cdata$AF<cutoff, 1:2]),1,function(x){x[1]+2*x[2]})
        afcuts <- rbind(afcuts,cbind(Class=varclass,Cutoff=cutoff,Patient=as.character(sort(unique(cdata$Patient))),Count=tmp,Type=geneclass))
    }
    afcuts <- as.data.frame( afcuts )
    #print(head(afcuts))

    #colnames(afcuts) <- c("Class","Cutoff","Patient","Count")
    afcuts$Cutoff <- as.numeric(as.character(afcuts$Cutoff))
    afcuts$Count <- as.numeric(as.character(afcuts$Count))

    print("End readClassData")
    return(list(afcuts,patientstats))
} # readClassData
    #outfile <- ifelse( !is.null(geneclass) ,
        #paste(figuredir,"class",varclass,"_",geneclass,"_patvarcounts.png",sep=""),
        #paste(figuredir,"class",varclass,"_patvarcounts.png",sep="")
        #)
    #png(outfile)
    #p <- ggplot(afcuts, aes(x=Cutoff, y=Count, group=Patient))
    #p <- p + geom_line()
    #print(p)
    #dev.off()


######################################################################
# inbredOutbredRatios
######################################################################
inbredOutbredRatios <- function( patannot, patset1, patset2, classdata, figuredir, prefix=NULL ){
    print(" - Running inbredOutbredRatios")
    #inbredpats <- patannot[ patannot$Group %in% c(1,9,10) 
                           #& patannot$Phenotype==0, ]$Individual.ID
    #outbredpats <- patannot[ patannot$Group %in% c(4,5,6,11)
                           #& patannot$Phenotype==0, ]$Individual.ID
    #print( unique(patannot$Region) )
    #print( head(patannot) )
    #inbredpats <- patannot[patannot$Region == "Middle Eastern",]$Individual.ID
                           #& patannot$Phenotype==0, ]$Individual.ID
    #outbredpats <- patannot[ patannot$Source == "Eichler"
                           #& patannot$Phenotype==0, ]$Individual.ID
    #outbredpats <- patannot[ patannot$Source == "KG"
                           #& patannot$Phenotype==0, ]$Individual.ID
    numoutbred <- sum(unique(classdata$Patient) %in% patset2)
    numinbred <- sum(unique(classdata$Patient) %in% patset1)
    print(paste("Outbred #:",numoutbred))
    print(paste("Inbred #:",numinbred))
    #print(patset1)
    #print(unique(classdata$Patient))
    print(length(unique(classdata$Patient)))

    breedingdata <- list()
    i <- 1
    for( class in unique( classdata$Class ) ){
        for( cutoff in unique(classdata$Cutoff) ){
            inbred <- sum( classdata[classdata$Class == class & 
                     classdata$Cutoff == cutoff & !is.na(classdata$Count) &
                     classdata$Patient %in% patset1,]$Count ) / numinbred
            outbred <- sum( classdata[classdata$Class == class & 
                     classdata$Cutoff == cutoff & !is.na(classdata$Count) &
                     classdata$Patient %in% patset2,]$Count ) / numoutbred
            new <- c(class,cutoff,inbred,outbred)
            #print(new)
            breedingdata[[i]] <- new
            i <- i +1
        }
    }
    breedingdata <- as.data.frame(do.call(rbind,breedingdata))
    colnames(breedingdata) <- c("VariantClass", "Cutoff", "Inbred","Outbred")
    breedingdata$Cutoff <- as.numeric(as.character(breedingdata$Cutoff))
    breedingdata$Inbred <- as.numeric(as.character(breedingdata$Inbred))
    breedingdata$Outbred <- as.numeric(as.character(breedingdata$Outbred))
    #print(breedingdata)

    #print(unique(breedingdata$Cutoff))
    #outfile <- paste(figuredir,"breedingcompare_",prefix,"_ratio.png",sep="")
    outfile <- ifelse( !is.null(prefix), 
        paste(figuredir,"breedingcompare_",prefix,"_ratio.png",sep=""),
        paste(figuredir,"breedingcompare_ratio.png",sep=""))
    logfile <- ifelse( !is.null(prefix), 
        paste(figuredir,"breedingcompare_",prefix,"_ratio.txt",sep=""),
        paste(figuredir,"breedingcompare_ratio.txt",sep=""))
    print(outfile)
    print(head(breedingdata))
    print(head(classdata))


    print(paste("Writing file:",logfile))
    write.table(breedingdata, logfile, quote=FALSE,sep="\t",row.names=FALSE)

    png(outfile, width=520)
    p <- ggplot(breedingdata, aes(x=Cutoff, y=Inbred/Outbred, group=VariantClass))#, colour=VariantClass))
    p <- p + labs(title = "Burden Ratio across Variant Classes with\nDecreasing Allele Frequency Cutoffs")
    p <- p + ylab("Inbred/Outbred burden ratio")+ xlab("Allele Frequency Cutoff")
    #p <- p + geom_line()
    p <- p + geom_line(size=1,aes(colour=VariantClass)) #linetype=VariantClass,
    #p <- p + geom_point(size=3,aes(shape=VariantClass))
    p <- p + ng1
    print(p)
    dev.off()
    print( "END inbredOutbredRatios" )
} # End inbredOutbredRatios

    #breeddata.melt <- melt(breedingdata, id=c("Class","Cutoff"))
    #outfile <- ifelse( !is.null(prefix), 
        #paste(figuredir,"breedingcompare_",prefix,"_separate.png",sep=""),
        #paste(figuredir,"breedingcompare_separate.png",sep=""))
    #print(outfile)
    #png(outfile)
    #p <- ggplot(breeddata.melt, aes(x=Cutoff, y=value, group=variable, colour=variable))
    #p <- p + geom_line()
    #print(p)
    #dev.off()


######################################################################
# groupplot
######################################################################
groupplot <- function( pstats.annot, column, figuredir ){
    print( " - Running groupplot" )
    barplotstack <- list()
    i <- 1
    #for( class in sort(unique(pstats.annot$Class)) ){
        #pstats.sub <- pstats.annot[pstats.annot$Class == class,]
        #p <- ggplot(pstats.sub, aes_string(x="Group", y=column, group="Group"))
        #p <- p + geom_boxplot(aes(fill = factor(Type)))
        #p <- p + labs(title = paste("Class",class) )
        #outfile <- paste(figuredir,column,"group_class",class,".png",sep="")
        #png(outfile)
        #print(p)
        #dev.off()
        #barplotstack[[i]] <- p
        #i <- i+1
    #}
    #multiplot(plotlist=barplotstack, cols=2 )
    png( paste(figuredir,column,"group_byclass.png",sep="") )
    p <- ggplot(pstats.annot, aes_string(x="Group", y=column, group="Group"))
    p <- p + geom_boxplot(aes(fill = factor(Type)))
    p <- p + labs(title = paste("Distribution of",column,"Across Groups") )
    p <- p + facet_grid(Class ~ . ) + theme(legend.position="none")
    p <- p + ng2
    dev.off()
} # End groupplot

######################################################################
# affectedplot
######################################################################
affectedplot <- function( pstats.annot, column, figuredir ){
    print(" - running affectedplot")
    #barplotstack <- list()
    #i <- 1
    #for( class in sort(unique(pstats.annot$Class)) ){
        #pstats.sub <- pstats.annot[pstats.annot$Class == class,]
        #p <- ggplot(pstats.sub, aes_string(x="Affected", y=column ))
        #p <- p + geom_boxplot(aes(fill = factor(Type)))
        #p <- p + labs(title = paste("Class",class) )
        #outfile <- paste(figuredir,column,"affected_class",class,".png",sep="")
        #png(outfile)
        #print(p)
        #dev.off()
        #barplotstack[[i]] <- p
        #i <- i+1
    #}
    #multiplot(plotlist=barplotstack, cols=2 )
    outfile <- paste(figuredir,column,"affected_byclass.png",sep="")
    p <- ggplot(pstats.annot, aes_string(x="Affected", y=column ))
    p <- p + geom_boxplot(aes(fill = factor(Affected)))
    p <- p + labs(title = paste("Affectedness versus",column,"By Variant Class") )
    p <- p + facet_grid(. ~ Class ) + theme(legend.position="none")
    p <- p + ng2
    png( outfile )
    print(p)
    dev.off()
} # END affectedplot

######################################################################
# inbredOutbredBurden
######################################################################
inbredOutbredBurden <- function( pstats.annot, patset1, patset2, figuredir, prefix=NULL ){
    print(" - Running inbredOutbredBurden ")
    #inbredpats <- pstats.annot[ pstats.annot$Group %in% c(1,9,10) 
                           #& pstats.annot$Phenotype==0, ]$Patient
    #outbredpats <- pstats.annot[ pstats.annot$Group %in% c(4,5,6,11)
                           #& pstats.annot$Phenotype==0, ]$Patient
    #inbredpats <- pstats.annot[ pstats.annot$Region == "Middle Eastern", ]$Patient
                           #& pstats.annot$Phenotype==0, ]$Patient
    #outbredpats <- pstats.annot[ pstats.annot$Source == "Eichler"
                           #& pstats.annot$Phenotype==0, ]$Patient
    numoutbred <- sum(levels(pstats.annot$Patient) %in% patset2)
    numinbred <- sum(levels(pstats.annot$Patient) %in% patset1)
    print(paste("Outbred #:",numoutbred))
    print(paste("Inbred #:",numinbred))
    #print(patset1)
    #print(unique(classdata$Patient))
    #print(length(unique(classdata$Patient)))

    pstats.annot$Breeding <- ifelse( pstats.annot$Patient %in% patset1,"Inbred",
                            ifelse(pstats.annot$Patient %in% patset2, "Outbred", "other" ))

    #pstats.sub <- pstats.annot[pstats.annot$Group != "#N/A" &
    pstats.sub <- pstats.annot[ pstats.annot$Breeding != "other",]
    #pstats.sub$Group <- factor( pstats.sub$Group, levels=sort(unique(pstats.sub$Group)))
    pstats.sub$Patient <- factor( pstats.sub$Patient, levels=sort(unique(pstats.sub$Patient)))


    print(table(pstats.sub$Breeding))
    importantcolumns <- c("homhet", "Burden", "titv" )
    comparisons = NULL
    for( column in importantcolumns ){
        for( class in unique(pstats.sub$Class) ){
            s <- wilcox.test( x=pstats.sub[pstats.sub$Breeding=="Inbred"
                             & pstats.sub$Class==class,][[column]],
                         y=pstats.sub[pstats.sub$Breeding=="Outbred"
                             & pstats.sub$Class==class,][[column]])
                         #,alternative="greater" )
            #print(s)
            comparisons <- rbind( comparisons, c("InbredvOutbred",column,class,
                                         numoutbred,numinbred,s$statistic,s$p.value))
        }
        p <- ggplot(pstats.sub, aes_string(x="Breeding", y=column ))
        p <- p + geom_boxplot(aes(fill=Breeding))
        p <- p + ylab("Variant Count") + xlab("Breeding Class")
        p <- p + labs(title = paste("Breeding Type",column,"Comparison\nBetween Variant Classes") )
        p <- p + facet_grid(. ~ Class ) + theme(legend.position = "none")
        p <- p + ng2
        #outfile <- paste(figuredir,column,"_breeding.png",sep="")
        outfile <- ifelse( !is.null(prefix), 
            paste(figuredir,"breeding_",column,"_",prefix,".png",sep=""),
            paste(figuredir,"breeding_",column,".png",sep=""))
        print(outfile)
        png(outfile)
        print(p)
        dev.off()
    }


    psums <- tapply( as.numeric(pstats.sub$Burden), pstats.sub$Patient, function(x){ sum(x) } )
    psums <- tapply( as.numeric(pstats.sub$Burden), pstats.sub$Patient, function(x){ sum(x) } )
    totalcounts <- cbind( unique(pstats.sub[colnames(pstats.sub) %in% c("Patient","Affected","Breeding")]), Burden=psums )

    s <- wilcox.test( x=totalcounts[totalcounts$Breeding=="Inbred",]$Burden,
                 y=totalcounts[totalcounts$Breeding!="Outbred",]$Burden)
    comparisons <- rbind( comparisons, c("InbredvOutbred","Burden","Total",
                                         numoutbred,numinbred,s$statistic,s$p.value))
    p <- ggplot(totalcounts, aes(x=Breeding, y=Burden ))
    p <- p + geom_boxplot()
    p <- p + labs(title = paste("Breeding Type Burden Comparison") )
    p <- p + ng1
    #p <- p + facet_grid(Class ~ .) + theme(legend.position = "none")
    #outfile <- paste(figuredir,column,"_breeding.png",sep="")
    outfile <- ifelse( !is.null(prefix),
        paste(figuredir,"breeding_Burden_total_",prefix,".png",sep=""),
        paste(figuredir,"breeding_Burden_total.png",sep=""))
    print(outfile)
    png(outfile)
    print(p)
    dev.off()

    colnames(comparisons) <- c("Type","Factor","Class","Outbred","Inbred","Statistic", "Pvalue" )
    outfile <- ifelse( !is.null(prefix),
        paste(figuredir,"breeding_mannwhitneyu_",prefix,".txt",sep=""),
        paste(figuredir,"breeding_mannwhitneyu.txt",sep=""))
    write.table( comparisons, outfile, quote=FALSE,sep="\t",row.names=FALSE )
    print( "End inbredOutbredBurden" )
} # END inbredOutbredBurden

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
        patlist <- as.vector(annot[annot$Lab == "Broad" & annot$continent == "Middle East", ]$Individual.ID)
        #& annot$Phenotype == 0
    }else if( psource == "eichler" ){
        patlist <- as.vector(annot[annot$Lab == "Eichler" & annot$continent == "Europe",]$Individual.ID)
    }else if(psource == "onekg" ){
        patlist <- as.vector(annot[annot$Lab == "1KG",]$Individual.ID)
    }else if(psource == "daily" ){
        patlist <- as.vector(annot[annot$Lab == "Daly" & annot$continent == "Europe",]$Individual.ID)
    }
    print( paste("Number of patients:",length(patlist)) )
    return(patlist)
} # END getPatlist


######################################################################
# variantclass_analysis
######################################################################
variantclass_analysis <- function( genelist,figuredir,targetdir,removepatients,patannot,varclasstype, classlist){
    print( " - Running variantclass_analysis" )
    # Essential - e, complex - c, mendelian - m
    colnames(genelist) <- c( "Gene", "Entrez", "Class" )
    geneclasses <- c("Essential","Complex","Mendelian", "Childhood Recessive", "Fertility")

    geneclass <- "Essential"
    classdata <- NULL
    pstats <- NULL
    print(names(patannot))
    print(head(patannot))
    print(unique(patannot$Region))
    print(unique(patannot$Source))
    #print(sum(patannot$Region == "Middle Eastern"))
    #print(sum(patannot$Source == "Eichler"))
    #patset1 <- getPatlist( patannot, "eichler")
    patset1 <- getPatlist( patannot, "internal")
    #patset2 <- getPatlist( patannot, "daily")
    patset2 <- getPatlist( patannot, "eichler")
    #keeppatients <- patannot[(patannot$Region == "Middle Eastern" | patannot$Source == "Eichler"),]$Individual.ID
    print("patlist1")
    print(head(patset1))
    print(length(patset1))
    print("patset2")
    print(head(patset2))
    print(length(patset2))
    if( length(patset1) == 0 | length(patset2) == 0 ){ 
        print("Error: no patients found!")
        q()
    }
    for( geneclass in geneclasses ){
        for( varclass in classlist ){
            tfile <- paste(targetdir,varclasstype,"_class",varclass,"_all.tsv",sep="")
            tmp <- readClassData( tfile, figuredir, patset1,patset2, genelist, geneclass )
            classdata <- rbind( classdata, tmp[[1]] )
            pstats <- rbind( pstats, tmp[[2]] )
            pstats$Class <- varclass
        }
    }
    classdata <- classdata[!(classdata$Patient %in% removepatients),]
    pstats <- pstats[!(pstats$Patient %in% removepatients),]
    #print(head(pstats))
    #print(head(patannot))

    #patannot$Affected <- ifelse( patannot$Phenotype %in% c(1,9,10),"Affected", 
                   #ifelse(patannot$Phenotype %in% c(4,5,6,11),"Unaffected","other"))
    #patannot$Affected <- ifelse( patannot$Phenotype == 1, "Affected", "Unaffected" )

    # Ratio Analysis
    print(table(classdata$Type))
    for( geneclass in geneclasses ){
        classdata.sub <- classdata[classdata$Type == geneclass,]
        inbredOutbredRatios( patannot, patset1, patset2, classdata.sub, figuredir, geneclass )
    }

    # Patient stats
    pstats.annot <- merge(pstats,patannot[names(patannot) %in% c("Individual.ID", "Affected","Phenotype","Source","Region")], by.x="Patient", by.y="Individual.ID" )
    #"Inbreeding.Coeff", "Group",
    #pstats.annot$Inbreeding.Coeff <- as.numeric(as.character(pstats.annot$Inbreeding.Coeff))

    # Burden Analysis
    for( geneclass in geneclasses ){
        print(paste("Geneclass",geneclass))
        pstats.annot.sub <- pstats.annot[pstats.annot$Type == geneclass,]
        inbredOutbredBurden( pstats.annot.sub, patset1, patset2, figuredir, geneclass )
    }

    print("Make Group plots")
    # TI/TV ratio plot
    #groupplot( pstats.annot, "titv", figuredir )
    # Burden plot
    #groupplot( pstats.annot, "Burden", figuredir )
    # Hom Het plot
    #groupplot( pstats.annot, "homhet", figuredir )
    # Affected titv ratio
    affectedplot( pstats.annot, "titv", figuredir )
    # Burden Affected plot
    affectedplot( pstats.annot, "Burden", figuredir )
    # Hom Het Affected plot
    affectedplot( pstats.annot, "homhet", figuredir )
    print( "END variantclass_anaysis" )
} # End variantclass_analysis

######################################################################
# allVariantsAnalysis
######################################################################
allVariantAnalysis <- function( targetdir, figuredir, varclasstype, classlist ){
    print("All variant analysis")

    patset1 <- getPatlist( patannot, "internal")
    #patset2 <- getPatlist( patannot, "daily")
    patset2 <- getPatlist( patannot, "eichler")
    #keeppatients <- c(patset1, patset2)
    print("patlist1")
    print(head(patset1))
    print(length(patset1))
    print("patset2")
    print(head(patset2))
    print(length(patset2))
    if( length(patset1) == 0 | length(patset2) == 0 ){ 
        print("Error: no patients found!")
        q()
    }
    #keeppatients <- patannot[patannot$Region == "Middle Eastern" | patannot$Source == "Eichler",]$Individual.ID
    #print( keeppatients )

    classdata <- NULL
    pstats <- NULL
    for( varclass in classlist ){
        tfile <- paste(targetdir,varclasstype,"_",varclass,".tsv",sep="")
        tmp <- readClassData( tfile,figuredir,patset1,patset2 )
        print(length(tmp))
        classdata <- rbind( classdata, tmp[[1]] )
        print("classdata")
        print(dim(classdata))
        pstats <- rbind( pstats, tmp[[2]] )
        print("pstats")
        print(dim(pstats))
    }
# Patient stats
    classdata <- classdata[!(classdata$Patient %in% removepatients),]
    pstats <- pstats[!(pstats$Patient %in% removepatients),]

    pstats.annot <- merge(pstats,patannot[,colnames(patannot) %in% c("Individual.ID",  "Region","Source","Affected","Other.info","Phenotype")], by.x="Patient", by.y="Individual.ID" )
    #"Inbreeding.Coeff","Group",
    #pstats.annot$Inbreeding.Coeff <- as.numeric(as.character(pstats.annot$Inbreeding.Coeff))

    inbredOutbredRatios( patannot, patset1, patset2, classdata, figuredir )
    inbredOutbredBurden( pstats.annot, patset1, patset2, figuredir )

    # Hom Het Affected plot

    #print( "Class 6 analysis" )
    #if( maxclass == 6 ){
        #pstats.sub <- pstats.annot[pstats.annot$Class == 6,]
    #} else {
        #pstats.sub <- pstats.annot[pstats.annot$Class == 4,]
    #}
    pstats.sub <- pstats.annot[pstats.annot$Class == "HIGH",]
    print(head(pstats.sub))
    #tmp <- pstats.sub[pstats.sub$Region == "Middle Eastern",]
    tmp <- pstats.sub[pstats.sub$Patient %in% patset1,]
    tmp$Burden <- as.numeric(tmp$Burden)
    print(head(tmp))

    p <- ggplot(tmp, aes(x=Affected, y=Burden, fill=Affected ))
    p <- p + geom_boxplot()
    p <- p + labs(title = "Affectedness vs Burden of Deleterious Variants")
    p <- p + xlab("Affectedness")
    p <- p + ng1
    png(paste(figuredir,"affectednessandburden.png",sep=""))
    print(p)
    dev.off()

    # Burden
    print("Burdens")
    #print(paste("Max class:",maxclass))
    #pstats.sub <- pstats.annot[pstats.annot$Class %in% c(6),]
    #if( maxclass == 6 ){
        #pstats.tmp <- pstats.annot[pstats.annot$Class == 6 & pstats.annot$Burden < 1000,]
    #}else{
        #pstats.tmp <- pstats.annot[pstats.annot$Class == 4 & pstats.annot$Burden < 1000,]
    #}
    pstats.tmp <- pstats.annot[pstats.annot$Class == "HIGH" & pstats.annot$Burden < 1000,]
    #pstats.tmp <- pstats.tmp[pstats.tmp$Group %in% c(1,9,10),]
    pstats.tmp <- pstats.tmp[pstats.tmp$Patient %in% patset1,]

    pstats.tmp$Burden <- as.numeric(pstats.tmp$Burden)

    p <- ggplot(pstats.tmp, aes(x=Affected, y=Burden, fill=Affected ))
    p <- p + geom_boxplot()
    p <- p + labs(title = "Affectedness vs Homozygous Burden\nof Deleterious Variants")
    p <- p + xlab("Affectedness")
    p <- p + ng1
    figurefile <-paste(figuredir,"affectedness_allburden.png",sep="")
    print( paste("Writing file:",figurefile))
    png(figurefile)
    print(p)
    dev.off()

    s <- wilcox.test( x=pstats.tmp[pstats.tmp$Affected=="Affected",]$Burden, y=pstats.tmp[pstats.tmp$Affected!="Affected",]$Burden, alternative="greater" )
    print(s)

    print("Homozygous Burden")
    p <- ggplot(pstats.tmp, aes(x=Other.info, y=hom, fill=Other.info ))
    p <- p + geom_boxplot()
    p <- p + labs(title = "Affectedness vs Homozygous Burden\nof Deleterious Variants")
    p <- p + xlab("Affectedness")
    p <- p + ng1
    figurefile <- paste(figuredir,"otherinfo_homburden.png",sep="")
    print(paste("Writing file:",figurefile))
    png(figurefile)
    print(p)
    dev.off()

    s <- wilcox.test( x=pstats.tmp[pstats.tmp$Other.info=="Affected",]$Burden, y=pstats.tmp[pstats.tmp$Other.info!="Affected",]$Burden, alternative="greater" )
    print(s)
    print(paste("Num patients:",length(pstats.tmp[,1])))

    print("Total Burden")
    p <- ggplot(pstats.tmp, aes(x=Other.info, y=Burden, fill=Other.info ))
    p <- p + geom_boxplot()
    p <- p + labs(title = "Affectedness vs Mutational Burden\nof Deleterious Variants")
    p <- p + xlab("Affectedness")
    p <- p + ng1
    figurefile <- paste(figuredir,"otherinfo_allburden.png",sep="")
    print(paste("Writing file:",figurefile))
    png(figurefile)
    print(p)
    dev.off()

    print("Homhet")
    p <- ggplot(pstats.tmp, aes(x=Other.info, y=homhet, fill=Other.info ))
    p <- p + geom_boxplot()
    p <- p + labs(title = "Affectedness vs Hom/Het Ratio\nof Deleterious Variants")
    p <- p + xlab("Affectedness")
    p <- p + ng1
    figurefile <- paste(figuredir,"otherinfo_homhetburden.png",sep="")
    print(paste("Writing file:",figurefile))
    png(figurefile)
    print(p)
    dev.off()
    print("End allVariantAnalysis")
} # END allVariantAnalysis 

######################################################################
# Main
######################################################################
#varclasstype <- "gaiadata"
#varclasstype <- "rawdatatenp"
#varclasstype <- "eichler"
#varclasstype <- "daily"
#varclasstype <- "rawdata"
#varclasstype <- "alldata"

varclasstype <- "ciliopathies.chimp.annot"
figuredir <- paste("figures/affunaff/",varclasstype,"/",sep="")
dir.create(file.path(figuredir), showWarnings = FALSE)

targetdir <- "classes/"
#removepatients <- read.table("toremove.txt",sep="\t")$V1

genelist <- read.table("resources/allgene_classes.tsv",header=FALSE,sep="\t",comment.char="")
annotationfile <- "patientannot/annotation/patientannotation_laser.txt"
patannot <- read.table(annotationfile,header=TRUE,sep="\t",comment.char="")
patannot$Affected <- ifelse( patannot$Phenotype == 1, "Affected", "Unaffected" )

#maxclass <- 5
classlist <- c("LOW","MODIFIER","MODERATE","HIGH")
# Compar Essential Complex and Mendelian
#variantclass_analysis( genelist, figuredir, targetdir, removepatients, patannot,varclasstype )

# allVariantsAnalysis
allVariantAnalysis( targetdir, figuredir, varclasstype, classlist )

## Test readClassData
#varclass <- 1
#patset1 <- getPatlist( patannot, "eichler")
#patset2 <- getPatlist( patannot, "internal")
#keeppatients <- c(patset1, patset2)
#
#tfile <- paste(targetdir,varclasstype,"_class",varclass,"_all.tsv",sep="")
#tmp <- readClassData( varclass, targetdir, figuredir, varclasstype, keeppatients )


# END Main


