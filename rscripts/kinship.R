#!/usr/bin/env Rscript

library(ggplot2)
#tbl=read.table("ALL.genome.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.vcf.plink_unrel.5.Q")

######################################################################
# readAnnotation
######################################################################
readAnnotation <- function(file, filterfile=NULL)
{
    print("readAnnotation")
    print(file)
    allannotdata <- read.table(file, header=TRUE, sep="\t", comment.char="")
    print("done reading...")
    rownames(allannotdata) <- allannotdata[,2]
    print(paste("Filterfile",filterfile))
    if( !is.null(filterfile) ){
        namelist <- read.table(filterfile,header=FALSE,sep=" ", comment.char="")
        #print(head(namelist))
        #print(head(namelist[,2]))
        matchingnames <-allannotdata$Individual.ID %in% intersect(namelist[,2], allannotdata$Individual.ID)
        print(head(matchingnames))
        newallannotdata <- allannotdata[matchingnames,]
        allannotdata <- newallannotdata[match(namelist[,2],newallannotdata$Individual.ID),]
        #tmpannot <- allannotdata[order(namelist),]
        #write.table(allannotdata$Individual.ID,file="test.txt",quote=FALSE,row.names=FALSE, col.names=FALSE)
        print(tail(namelist))
    }
    print(tail(allannotdata))
    return(allannotdata)
} # End readAnnotationData

#####################################################################
# printHelp
#####################################################################
printHelp <- function(args){
    print( "Hello, you have entred the help text")
    print( args )
    print("length")
    print( length(args) )
    print( "./kinship.R <basename>" )
    print( "Goodbye")
    q()
} # End printHelp


######################################################################
# MAIN
######################################################################

setwd( "/home/escott/workspace/variome/rscripts/" )
print(paste("Working dir:",getwd()))

targetdir <- "../rawdata/daily/"
outdir <- "kinshipplots"
#basename <- "all_plates_allchr.plink.mod.filt_rel"
#basename <- "ALL.genome.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.plink.mod.filt_rel"
#basename <- "ALL.genome.phase1.plink.mod.filt"
#basename <- "ALL.genome.phase1.plink.mod.filt.uniq"
#basename <- "chr1.phase1.plink.mod.filt_rel"
#basename <- "merged_snps.plink.mod.filt"
basename <- "daily.chimp.plink.filt" #_rel.kin

# Parse Basename
args <- commandArgs( trailingOnly = TRUE )
#basename <- sub( "\\..*","", ifelse(length(args)>=1,args[1],basename))
basename <- sub( "_rel\\.kin.*","", basename(ifelse(length(args)>=1,args[1],basename)))
#targetdir <- dirname(ifelse(length(args)>=1,args[1],targetdir))
targetdir <- ifelse(length(args)>=1,dirname(args[1]),targetdir)
print(paste("Basename",basename))
print(paste("Targetdir",targetdir))

kinship0 <- paste(targetdir,"/",basename,"_rel.kin",sep="")
kinship1 <- paste(targetdir,"/",basename,"_rel.kin0",sep="")
print( paste( "Kin:",kinship0,"Kin1:",kinship1) )
tbl1 <- read.table(kinship0, header=TRUE, comment.char="")
tbl2 <- read.table(kinship1, header=TRUE, comment.char="")

#Kinship <- c( tbl1$Kinship, tbl2$Kinship )
#IBS0 <- c( tbl1$IBS0, tbl2$IBS0 )

outfile <- paste(outdir,"/",basename,"_rel.png",sep="")
print(paste("Writing file:", outfile))
png(outfile)
par(mfrow=c(1,2))
plot(tbl1$IBS0,tbl1$Kinship, ylim=c(0,.5), xlim=c(0,0.06), main="Family Relationships", xlab="Pr(IBS=0)", ylab="Kinship coefficient", col = ifelse(tbl1$Kinship < 0.0884,'black','red'))
abline(h=c(.354,.177,.0884,.0442),lty=2,col="blue")
plot(tbl2$IBS0,tbl2$Kinship, ylim=c(0,.5), xlim=c(0,0.06), main="Inter-family Relationships", xlab="Pr(IBS=0)", ylab="Kinship coefficient", col = ifelse(tbl2$Kinship < 0.0884,'black','red'))
abline(h=c(.354,.177,.0884,.0442),lty=2,col="blue")
dev.off()

#>0.354, [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884] 
# duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree
labels <- c("duplicates","1st degree","2nd degree","3rd degree","Possible")
start <-c(.354,.177,0.0884,.0442,0)
end <-c(1,.354,.177,0.0884,.0442)
ranges <- cbind( labels, start, end )

outfile <- paste(targetdir,"/",basename,"_interfam.txt",sep="")
append <- FALSE
for( i in 1:length(ranges[,1]) ){
    #print(ranges[i,])
    kin <- tbl1[tbl1$Kinship >= as.numeric(ranges[i,2]) & tbl1$Kinship < as.numeric(ranges[i,3]),]
    if( dim(kin)[1] == 0 ) {
        next
    }
    kin <- cbind(kin,relationship=ranges[i,1])
    #print(kin)
    write.table(kin,outfile,append=append, quote=FALSE,row.names=FALSE, col.names=(!append), sep="\t")
    append <- TRUE
}

outfile <- paste(targetdir,"/",basename,"_btwfam.txt",sep="")
append <- FALSE
for( i in 1:length(ranges[,1]) ){
    #print(ranges[i,])
    kin <- tbl2[tbl2$Kinship >= as.numeric(ranges[i,2]) & tbl2$Kinship < as.numeric(ranges[i,3]),]
    if( dim(kin)[1] == 0 ) {
        next
    }
    kin <- cbind(kin,relationship=ranges[i,1])
    #print(head(kin))
    write.table(kin,outfile,append=append, quote=FALSE,row.names=FALSE, col.names=(!append), sep="\t")
    append <- TRUE
}

#q()
#annotfile <- "inhouse_samples_annotation.ped"
#annotfile <- "phase1_samples_integrated_20101123_clean1.ped"
#filterfile <- paste(targetdir,"/",basename,".fam",sep="")
#annotdata <- readAnnotation( annotfile, filterfile )



#tmp <- merge(kin, annotdata, by.y="Individual.ID",by.x="ID1")
#kin.annot <- merge(tmp, annotdata,by="Individual.ID")
#print(head(kin.annot))
#kin.dup <- tbl1[tbl1$Kinship >= .354,]
#kin.1st <- tbl1[tbl1$Kinship >= .177 & tbl1$Kinship < .354,]
#kin.2nd <- tbl1[tbl1$Kinship >= 0.0884 & tbl1$Kinship < .177,]
#kin.3rd <- tbl1[tbl1$Kinship >= 0.0442 & tbl1$Kinship < .0884,]
#kin.pos <- tbl1[tbl1$Kinship >= 0.0 & tbl1$Kinship < .0442,]
#tmp <- merge(kin.3rd, annotdata, by.y="Individual.ID",by.x="ID1")
#kin.3rd.annot <- merge(tmp, annotdata,by.x="Individual.ID",by.y="Individual.ID")


#plot(IBS0,Kinship, ylim=c(0,.5), xlim=c(0,0.06), main="Kinship", xlab="Pr(IBS=0)", ylab="Kinship coefficient", col = ifelse(Kinship < 0.0884,'black','red'))
#

#abline(.354,0)
#abline(.177,0)
#abline(.0884,0)
#abline(.0442,0)
#estimated kinship coefficient range >0.354, [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884]

#tbl <- read.table("plateII_snps.vcf.plink_homo.kin0", header=TRUE)
#plot(tbl$IBD0,tbl$Kinship)

#tbl=read.table("ALL.genome.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.vcf.plink_unrel.1.Q")

