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
        #namelist <- read.table(filterfile,header=FALSE,sep=" ")
        #print(head(namelist))
        #print(head(namelist[,2]))
        #matchingnames <-allannotdata$Individual.ID %in% intersect(namelist[,2], allannotdata$Individual.ID)
        #print(head(matchingnames))
        #newallannotdata <- allannotdata[matchingnames,]
        #allannotdata <- newallannotdata[match(namelist[,2],newallannotdata$Individual.ID),]
        namelist <- read.table(filterfile,header=FALSE,sep=" ",comment.char="",col.names=c("V1","Individual.ID","V3","V4","V5","V6"))
        allannotdata <- merge( namelist, allannotdata, by="Individual.ID", all.x=TRUE )
        allannotdata <- allannotdata[namelist$Individual.ID,]
        #allannotdata <- allannotdata[with(allannotdata,
        print(head(namelist$Individual.ID))
        print(head(allannotdata$Individual.ID))
        #print("Missing samples")
        #print(allannotdata[is.null(allannotdata["Family.ID"]),])
    }
    #print(tail(allannotdata))
    return(allannotdata)
} # End readAnnotation

######################################################################
# clusterAdmix
######################################################################
clusterAdmix <- function(mydata,currK,outfile)
{
    print("clusterAdmix")
    ## Ward Hierarchical Clustering
    d <- dist(mydata, method = "euclidean") # distance matrix
    fit <- hclust(d, method="ward") 
    #png( outfile )
    #plot(fit) # display dendogram
    groups <- cutree(fit, k=currK) # cut tree into 5 clusters
    # draw dendogram with red borders around the 5 clusters 
    #rect.hclust(fit, k=currK, border="red")
    #dev.off()
    hc <- hclust(dist(mydata))
    png( outfile )
    plot( as.phylo(hc), type="radial")
    dev.off()
    return(groups)
}

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
# makeSummaryBarplot
######################################################################
makeSummaryBarplot <- function( figuredir, basename, currK, mydata, newannot ){
    print("###############################################################")
    print("Running makeSummaryBarplot")
    print(paste("K is",currK))

    continents <- newannot$Continent
    populations <- newannot$ethnicity

    meltdat <- mydata
    meltdat$Individual.ID <- rownames(mydata)
    meltdat <- melt( meltdat, id="Individual.ID" ) 
    mergedat <- merge(newannot[,c("Individual.ID","Ethnicity","Source","ethnicity","Continent","Origin")],meltdat,on=Individual.ID,ALL.y=TRUE)
    #dd[with(dd, order(-z, b)), ]
    mergedat <- mergedat[with(mergedat,order(variable,ethnicity,Individual.ID,value)),]
    #mergedat <- cbind(newannot[,c("Individual.ID","Ethnicity","Source","ethnicity","Continent")],newtbl)
    mergedat.filt <- mergedat[mergedat$Continent %in% c("Europe","Africa","East Asia","Middle East","South Asia"),]
    s <- cast(melt(mergedat.filt[,c("value","variable","Continent")]),Continent ~ variable, mean)
    new <- melt(s, id="Continent")
    t <- do.call("rbind", by(new, new$Continent, function(x) x[which.max(x$value), ]))
    print(t)
    print(head(new))
    newlevels <- unique(c(as.character(t$variable), as.character(levels(new$variable))))
    print("newlevels:")
    print(newlevels)
    #new$variable <- factor(new$variable, levels=newlevels)
    renamedlvls = paste("A",1:currK,sep="")
    new$variable <- mapvalues(new$variable, from = newlevels, to=renamedlvls)
    new$variable <- factor(new$variable, levels=sort(unique(new$variable)))
    new <- new[order(new$variable,new$Continent),]
    print( head(new))
    print( levels(new$variable) )
    mergedat$variable <- factor(mergedat$variable, levels=newlevels)
    mergedat <- mergedat[order(mergedat$variable,mergedat$Continent),]

    m1 <- ggplot(new, aes(x=Continent, y=value, fill=variable ))
    m1 <- m1 + geom_bar(stat="identity")
    m1 <- m1 + scale_fill_brewer(palette="Set1")
    m1 <- m1 + theme(legend.position = "none")
    m1 <- m1 + theme(axis.text.x=element_text(angle = 45, hjust = 1))
    m1 <- m1 + xlab("Continents") #ylab("Country")+ 
    m1 <- m1 + ng1
    outfile1 <- paste(figuredir,"/",basename,currK,"_summary.png",sep="")
    print( paste("Writing file:",outfile1) )
    png( outfile1 )
    print(m1)
    dev.off()

    mergedat.filt <- mergedat[mergedat$Continent == "Middle East",]
    s <- cast(melt(mergedat.filt[,c("value","variable","ethnicity")]),ethnicity ~ variable, mean)
    new <- melt(s, id="ethnicity")
    m2 <- ggplot(new, aes(x=ethnicity, y=value, fill=variable ))
    m2 <- m2 + geom_bar(stat="identity")
    m2 <- m2 + theme(legend.position = "none")
    m2 <- m2 + theme(axis.text.x=element_text(angle = 45, hjust = 1))
    m2 <- m2 + scale_fill_brewer(palette="Set1")
    m2 <- m2 + ng1
    #m2 <- m2 + labs(title="Middle East IBC")
    m2 <- m2 + xlab("Middle East Ethnicities") #ylab("Country")+ 
    #print(m2)
    outfile2 <- paste(figuredir,"/",basename,currK,"_summary_ethnicity.png",sep="")
    print( paste("Writing file:",outfile2) )
    png( outfile2 )
    print(m2)
    dev.off()

    #grid.newpage()
    #pushViewport(viewport(layout = grid.layout(1, 2)))
    #print(m1, vp = viewport(1, 1))
    #print(m2, vp = viewport(1, 2))
    #multiplot( m1, m2, cols=2 )

    mergedat.filt <- mergedat[mergedat$Continent == "South Asia",]
    s <- cast(melt(mergedat.filt[,c("value","variable","ethnicity")]),ethnicity ~ variable, mean)
    new <- melt(s, id="ethnicity")
    m3 <- ggplot(new, aes(x=ethnicity, y=value, fill=variable ))
    m3 <- m3 + geom_bar(stat="identity")
    m3 <- m3 + theme(legend.position = "none")
    m3 <- m3 + theme(axis.text.x=element_text(angle = 45, hjust = 1))
    m3 <- m3 + scale_fill_brewer(palette="Set1")
    m3 <- m3 + ng1
    m3 <- m3 + labs(title="South Asia IBC")
    m3 <- m3 + xlab("South Asia") #ylab("Country")+ 
    outfile3 <- paste(figuredir,"/",basename,currK,"_summary_SA.png",sep="")
    print( paste("Writing file:",outfile3) )
    png( outfile3 )
    print(m3)
    dev.off()

    #mecountries <- c("Morocco","Algeria","Tunisia","Libya","North Africa","Egypt","Saudi Arabia","Oman","Qatar","UAE","Yemen","Jordan","Palestine","Lebanon","Syria","Kuwait","Iraq","Turkey")
    meethnicities <- c("Mozabite","Bedouin","Palestinian","Druze","Adygei")

    mergedat.filt <- mergedat[(mergedat$ethnicity %in% meethnicities),] #(mergedat$Continent == "Middle East") & 
    mergedat.filt$ethnicity <- factor(mergedat.filt$ethnicity, levels=meethnicities)

    mergedat.filt <- mergedat.filt[with(mergedat.filt,order(variable,value)),]
    d <- dist(mydata, method = "euclidean") # distance matrix
    fit <- hclust(d, method="ward") 

    #row.names(mydata[with(mydata,order("V1","V2","V3","V4")),])
    #mergedat.filt$Individual.ID <- factor(mergedat.filt$Individual.ID, levels=row.names(mydata[with(mydata,order(V2,V3,V4,V1)),]))

    m4 <- ggplot(mergedat.filt, aes(x=Individual.ID, y=value, fill=variable ))
    m4 <- m4 + geom_bar(stat="identity")
    m4 <- m4 + theme(legend.position = "none")
    m4 <- m4 + ng1 + facet_grid( .~ethnicity, scales="free_x" )
    m4 <- m4 + scale_fill_brewer(palette="Set1")
    m4 <- m4 + theme(axis.text.x=element_blank())
    m4 <- m4 + theme(axis.title.y=element_blank()) + theme(axis.text.y=element_blank())
    outfile4 <- paste(figuredir,"/",basename,currK,"_summary_admix.png",sep="")
    print( paste("Writing file:",outfile4) )
    png( outfile4 )
    print(m4)
    dev.off()
    finalfig <- paste("results/figures/",basename,"_",currK,"_ethnicity_final.png",sep="")
    print( paste("Writing file:",finalfig) )
    png(finalfig, width=7, height=4, units="in", res=400)
    grid.arrange(m1,m4, widths=c(1/3,2/3), ncol=2) 
    dev.off()
    
    #mergedat.filt <- mergedat[(mergedat$Continent == "Middle East"),]
    mecountries <- c("Morocco","Algeria","Tunisia","Libya","North Africa","Egypt","Saudi Arabia","Oman","Qatar","UAE","Yemen","Jordan","Palestine","Lebanon","Syria","Kuwait","Iraq","Turkey")

    mergedat.filt <- mergedat[(mergedat$Continent == "Middle East") & (mergedat$Origin %in% mecountries),]
    mergedat.filt$Origin <- factor(mergedat.filt$Origin, levels=mecountries)
    s <- cast(melt(mergedat.filt[,c("value","variable","Origin")]),Origin ~ variable, mean)
    new <- melt(s, id="Origin")
    new$variable <- mapvalues(new$variable, from = newlevels, to=renamedlvls)
    new$variable <- factor(new$variable, levels=sort(unique(new$variable)))
    m5 <- ggplot(new, aes(x=Origin, y=value, fill=variable ))
    m5 <- m5 + geom_bar(stat="identity")
    m5 <- m5 + theme(legend.position = "none")
    m5 <- m5 + theme(axis.text.x=element_text(angle = 45, hjust = 1))
    m5 <- m5 + scale_fill_brewer(palette="Set1")
    m5 <- m5 + ng1
    m5 <- m5 + theme(axis.title.y=element_blank()) + theme(axis.text.y=element_blank())
    m5 <- m5 + xlab("Middle East Countries") #ylab("Country")+ 
    print(m5)
    outfile5 <- paste(figuredir,"/",basename,currK,"_summary_ME.png",sep="")
    print( paste("Writing file:",outfile5) )
    png( outfile5 )
    print(m5)
    dev.off()

    finalfig <- paste("results/figures/",basename,"_",currK,"_country_final.png",sep="")
    print( paste("Writing file:",finalfig) )
    png(finalfig, width=7, height=4, units="in", res=400)
    grid.arrange(m1,m5, widths=c(1/3,2/3), ncol=2) 
    dev.off()

    return(new)
} # makeSummaryBarplot
    #m <- m + scale_x_continuous( breaks=(reg_index-reg_index[1])*coef,labels=regs)

    #print( paste("Coef:",coef) )
    #print( paste("atts:",reg_index, length(mydata), mat.or.vec(length(regs),1) ))
    #print(paste("reg_index",reg_index))
    #print(paste("cont_index",cont_index))
    #print( length(mydata[,1]) )
    #newtbl <- mydata
    #newtbl$sample <- rownames(mydata)
    #newtbl <- melt( newtbl, id="sample" ) 
    #m <- ggplot(newtbl, aes(x=sample, y=value, fill=variable ))
    #m <- m + geom_bar(stat="identity")
    #m <- m + theme(legend.position = "none")
    #m <- m + theme(axis.text.x=element_blank())
    #m <- m + scale_x_continuous( breaks=(reg_index-reg_index[1])*coef,labels=regs)
    #m <- m + ng1
    #barwidth = 7/nrow(mydata)
    #print( paste("Writing file:",outfile) )
    #png( outfile, width = 7, height = 3, units="in", res=400) 
    #barplot(t(as.matrix(mydata)), col=rainbow( currK ),
           #ylab="Ancestry", border=NA, xaxt="n" )
        ##xlab="Individual", main="Population Admxiture"
    ##axis(1, at = (reg_index-reg_index[1])*coef, labels = regs, hadj=2.5)
    #axis(1, at = (cont_index-cont_index[1])*coef, labels=conts)#, padj = 2.5)#,tick=FALSE)
    #dev.off()


######################################################################
# makeBarplot
######################################################################
makeBarplot <- function( mydata, outfile, continents, populations ){
    print("Running makeBarplot")
    conts <- unique(continents)
    cont_index <- mat.or.vec(length(conts), 1)
    for(i in 1:length(conts)) {
        cont_index[i] <- mean(which(continents == conts[i]));
    }
    #cont_index <- 1360/nrow(mydata)*cont_index -55;

    regs <- unique(populations);
    reg_index <- mat.or.vec(length(regs), 1);
    for(i in 1:length(regs)) {
        reg_index[i] <- mean(which(populations == regs[i]));
    }
    #reg_index <- 500/nrow(mydata)*reg_index - 55;
    #coef <- 1980/nrow(mydata)
    if( any(grep("daily", outfile)) ){
        print( "Daily specific coefficient" )
        coef <- 1460/nrow(mydata)
    }else if( any(grep("everything", outfile)) ){
        print("found everything")
        coef <- 1058/nrow(mydata)
    }else if( any(grep("merged", outfile)) ){
        print("found merged")
        coef <- 2900/nrow(mydata)
    }else if( any(grep("onekg", outfile)) ){
        print("found onekg")
        coef <- 1460/nrow(mydata)
    }else{
        coef <- 710/nrow(mydata)
    }
    print( paste("Coef:",coef) )
    #print( paste("Coef:",coef) )
    #print( paste("atts:",reg_index, length(mydata), mat.or.vec(length(regs),1) ))
    #print(paste("reg_index",reg_index))
    #print(paste("cont_index",cont_index))
    #print( length(mydata[,1]) )
    #newtbl <- mydata
    #newtbl$sample <- rownames(mydata)
    #newtbl <- melt( newtbl, id="sample" ) 
    #m <- ggplot(newtbl, aes(x=sample, y=value, fill=variable ))
    #m <- m + geom_bar(stat="identity")
    #m <- m + theme(legend.position = "none")
    #m <- m + theme(axis.text.x=element_blank())
    #m <- m + scale_x_continuous( breaks=(reg_index-reg_index[1])*coef,labels=regs)
    #m <- m + ng1
    #barwidth = 7/nrow(mydata)
    print( paste("Writing file:",outfile) )
    png( outfile, width = 7, height = 3, units="in", res=400) 
    barplot(t(as.matrix(mydata)), col=rainbow( currK ),
           ylab="Ancestry", border=NA, xaxt="n" )
        #xlab="Individual", main="Population Admxiture"
    #axis(1, at = (reg_index-reg_index[1])*coef, labels = regs, hadj=2.5)
    axis(1, at = (cont_index-cont_index[1])*coef, labels=conts)#, padj = 2.5)#,tick=FALSE)
    #text((cont_index-cont_index[1])*coef, labels=conts, xpd=TRUE, srt=45)
    dev.off()
} # makeBarplot

######################################################################
# Main
######################################################################

setwd( "/home/escott/workspace/variome/" )

annotfile <- "./resources/annotation/patientannotation.ped"
figuredir <- "./rscripts/admixplot"

#targetdir <- "../data/merged"
#basename <- "ciliopathies.unfilt.recode.plink.mod.filt.uniq"

#targetdir <- "./rawdata/daily/admixture"
#basename <- "daily.chimp.regions.filt.samp.plink.mod.filt.uniq.fam"

#targetdir <- "./rawdata/test/admixture"
#basename <- "everything_set1.chr1.snp.chimp.regions.filt.samp.plink.mod.filt.uniq.fam"

targetdir <- "./rawdata/merged/admixture"
basename <- "admixture/merged.chimp.regions.filt.samp.plink.filt.fam"

#targetdir <- "./rawdata/onekg/admixture"
#basename <- "onekg.chimp.regions.filt.samp.plink.mod.filt.uniq.fam"

#targetdir <- "../data/daily1"
#basename <- "daily.recode.plink.mod.filt.uniq"

# Parse Basename
args <- commandArgs( trailingOnly = TRUE )
if( length(args) == 0 ){ printHelp(args) }
#gsub("(.*)\\.(.*)", "\\1\\2", N)
basename <- basename(gsub( "(.*)\\.(.*)", "\\1", ifelse(length(args)>=1,args[1],basename)))
targetdir <- ifelse(length(args)>=1,dirname(args[1]),targetdir)
annotfile <- ifelse(length(args)>=2,args[2],annotfile)
maxK <- ifelse(length(args)>=3,args[3],6)
print(paste("Max K",maxK))


print(paste("Basename",basename))
print(paste("Targetdir",targetdir))
print(paste("Annotation",annotfile))


filterfile <- paste(targetdir,"/",basename,".fam",sep="")
print(paste("Filterfile:",filterfile))
annotdata <- readAnnotation( annotfile, filterfile )
#print(head(annotdata))
#print( "End annotdata" )

print(dim(annotdata))
#x <- as.data.frame( table( annotdata$Individual.ID ) )
#print(x[x$Freq > 1,])
print(length(unique(annotdata$Individual.ID)))
print(annotdata[duplicated(annotdata$Individual.ID),])

suffix <- ".Q"

allKdata <- NULL
for( currK in 2:maxK ){
    currfile <- paste(targetdir,"/",basename,".",currK,suffix,sep="")
    print("read table")
    tbl=read.table( currfile, header=FALSE, sep=" ", comment.char="#" )
    print(dim(tbl))

    rownames(tbl) <- annotdata$Individual.ID

    print("Sort Data!")
    #newtbl <- tbl[order(annotdata$Continent, annotdata$Country, annotdata$ethnicity, annotdata$Individual.ID),]
    newtbl <- tbl[order(annotdata$Continent, annotdata$ethnicity),]
    print(paste("New table", dim(newtbl)))
    #newannot <- annotdata[order(annotdata$Continent, annotdata$Country, annotdata$ethnicity, annotdata$Individual.ID),]
    newannot <- annotdata[order(annotdata$Continent, annotdata$ethnicity),]
    print(paste("New Annot", dim(newannot)))
    outfile <- paste(targetdir,"/",basename,"_",currK,".sorted.txt",sep="")
    print(paste("Writing file:",outfile))
    write.table( cbind(newannot,newtbl), outfile, quote=FALSE,row.names=FALSE, col.names=TRUE, sep="\t")

    mydata <- newtbl[!(is.na(newannot$ethnicity)),]
    newannot <- newannot[!(is.na(newannot$ethnicity)),]

    #m <- m + theme(axis.text.x=element_blank())
    #m <- m + scale_x_continuous( breaks=(reg_index-reg_index[1])*coef,labels=regs)
 
    fplotdata <- makeSummaryBarplot( figuredir, basename, currK, mydata, newannot )
    fplotdata$currK <- currK
    print(head(fplotdata))
    allKdata <- rbind(allKdata,fplotdata)
}
head(allKdata)

m6 <- ggplot(allKdata, aes(x=Origin, y=value, fill=variable ))
m6 <- m6 + geom_bar(stat="identity")
m6 <- m6 + theme(legend.position = "none")
m6 <- m6 + theme(axis.text.x=element_text(angle = 45, hjust = 1))
m6 <- m6 + scale_fill_brewer(palette="Set1")
m6 <- m6 + ng1
m6 <- m6 + theme(axis.title.y=element_blank()) + theme(axis.text.y=element_blank())
m6 <- m6 + xlab("Middle East Countries") #ylab("Country")+ 
m6 <- m6 + facet_grid( currK~. )
print(m6)

outimage <- paste(figuredir,"/",basename,"_alladmix_final.png",sep="")
print( paste("Writing file:",outimage) )
png( outimage )
print(m6)
dev.off()


#print(paste("Writing:",outimage))
#png(outimage)
#multiplot(plotlist=fplotlist, cols=1 )
#dev.off()

    #multiplot( m1, m2, cols=2 )
#filterfile <- "all_plates_allchr.plink.mod.filt_unrelated.txt"
#annotfile <- "phase1_samples_integrated_20101123_clean1.ped"
#filterfile <- "ALL.genome.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.plink.mod.filt_unrelated.txt"
#filterfile <- "ALL.genome.phase1.plink.mod.filt.uniq.fam"
#filterfile <- "ALL.genome.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.plink.mod.filt.uniq.fam"

#q()
#basename <- "ALL.genome.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.vcf.plink_unrel."
##colorPalette <- chooseColors( newannot$Population )
##makeHeatmap( newtbl, outfile, newannot$Population )
## K = 5
#tbl=read.table("plateII_snps.vcf.plink_unrel.2.Q")
#tbl.sort <- tbl[with(tbl,order(-V1,-V2)),]
#barplot(t(as.matrix(tbl.sort)), col=rainbow(2),
      #xlab="Individual #", ylab="Ancestry", border=NA)

