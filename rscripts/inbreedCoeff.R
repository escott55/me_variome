#!/usr/bin/env Rscript

####!/usr/bin/env Rscript
#install.packages("reshape")

#library(gplots)
library(ggplot2)
library(reshape)

#####################################################################
# Globals
#####################################################################
outputdir <- "./images/"

#####################################################################
# printHelp
#####################################################################
printHelp <- function(args){
    print( "Hello, you have entred the help text")
    print( args )
    print("length")
    print( length(args) )
    print( "./inbreedingCoeff.R <targetdir>" )
    print( "Goodbye")
    q()
} # End printHelp

#--------------------------------------------------------------------#
# Simulation Data
#--------------------------------------------------------------------#

#####################################################################
# resultsStatsMain
#####################################################################
resultsStatsMain <- function(exDir, prefix=NULL){
    print( "running - resultsStatsMain" )
    fpattern <- paste(prefix,".*het",sep="")
    data.files <- list.files(exDir, pattern=fpattern, full.names = TRUE)

    print(data.files)
    inbreeddata <- data.frame()
    for( cfile in data.files ){
        components <- strsplit(cfile,'\\.')
        group <- as.character(components[[1]][length(components[[1]])-1])
        group <- as.numeric(as.character(substr(group,4,nchar(group))))
        #print( group )

        alldata <- read.table( cfile, header=TRUE, comment.char="" )
        #colnames(alldata) <- c("popSize","AF","genMax","Type","Rep","Deaths","Purged")
        #print( head(alldata) )
        print( paste("Group",group,":",mean(alldata$F)))
        inbreeddata <- rbind( inbreeddata, cbind(alldata,group) )
    }

    ibcoeff.sort <- inbreeddata[order(inbreeddata$group),]
    write.table(ibcoeff.sort, file=paste(exDir,"/combined.tsv",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
    #head( ibcoeff.sort )
    #ibcoeff.sort <- within(ibcoeff.sort, 
          #group <- factor(group,
          #levels=names(sort(table(group), 
          #decreasing=FALSE))))
    outfile <- paste(outputdir,prefix,"_ibcbox.png",sep="")
    print(paste("Making Figure:",outfile))
    png(outfile)
    p <- ggplot(ibcoeff.sort, aes(factor(group), F))
    p <- p + geom_boxplot(aes(fill = factor(group)),notch = TRUE)
    p <- p + ylab("Fsnp")+ xlab("Groups")
    p <- p + labs(title = "Inbreeding Coefficient Ranges Across Groups")
    print(p)
    dev.off()
    #png(paste(outputdir,"deathhist_",breedingtype,".png",sep=""))
    #m <- ggplot(purgeddata, aes(Deaths))
    #print(m + geom_histogram(colour = "darkgreen", fill="white")) #+ scale_y_log10())
    #dev.off()
} # End resultsStatsMain


#####################################################################
# Main
#####################################################################


setwd( "/home/escott/workspace/variome/rscripts/" )
print(paste("Working dir:",getwd()))

#exDir <- "./grp11"
exDir <- "./merged/cohort"

# Parse arguments
prefix <- NULL
prefix <- "/home/escott/workspace/variome/rawdata/test/ibc everything_set1.chr1.snp.chimp.regions.filt.samp.samp.plink.recode12.mod"
args <- commandArgs( trailingOnly = TRUE )
if( length(args) == 0 ){ printHelp(args) }
exDir <- ifelse(length(args)>=1,args[1],exDir)
if(length(args)>=2){
    prefix <- args[2]
}

print(exDir)
print(prefix)
resultsStatsMain( exDir, prefix )

