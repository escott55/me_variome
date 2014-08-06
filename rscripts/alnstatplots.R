#!/usr/bin/env Rscript


#library(gplots)
library(ggplot2)
library(reshape)

source("themes.R")

simpleCap <- function(x) {
    s <- strsplit(x, "[. ]")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
}

alnstats <- read.table("../resources/alignmentstats.txt", sep="\t",header=TRUE)
print( head(alnstats) )

#gleesonannot <- read.table("/home/escott/nextjensystem/src/results/gleesonannot.ped",sep="\t",header=TRUE)
gleesonannot <- read.table("../resources/gleesonannot.ped",sep="\t",header=TRUE)
fowzanannot <- read.table("../resources/fowzanannot.ped",sep="\t",header=TRUE)
casanovaannot <- read.table("../resources/casanovaannot.ped",sep="\t",header=TRUE)

allannot <- rbind(fowzanannot,gleesonannot,casanovaannot)

stats.merged <- merge( alnstats, allannot, by.x="Patient", by.y="Individual.ID" )

head(alnstats[!(alnstats$Patient %in% allannot$Individual.ID),]["Patient"],50)

print( colnames(stats.merged) )
# Stats 
# Mapping
#  - How many reads mapped violin
#  - Paired reads mapped violin
# Coverage
#  - Box MeanCov
#  - Percent over 8x.... somehow
# Ethnicities
#  - Pie chart ethnicities by Dataset
#  - Continent pie charts
# Other 
#  - Gender disribution

outdir <- "alnstatfigs"
#for( dataset in sort(unique(stats.merged$Dataset)) ){
    #cdata <- stats.merged[stats.merged$Dataset == dataset,]
    #att <- "mapped"
for( att in c("mapped","paired","properly.paired","singletons") ) {
    print(att)
    cdata <- stats.merged[complete.cases(stats.merged[[att]]),]
    qs <- quantile( cdata[[att]] )
    irq <- qs[4] - qs[2]
    upperb <- qs[4] + irq*1.5
    lowerb <- qs[2] - irq*1.5
    print( paste("IRQ",irq,upperb,lowerb) )
    outimage <- paste(outdir,"/alignment_",att,"_violin.png",sep="")
    print(outimage)
    png(outimage)
    m <- ggplot(cdata, aes_string(x="Dataset", y=att, fill="Dataset"))
    m <- m + geom_violin(show_guide=FALSE)
    m <- m + theme(legend.position = "none")
    m <- m + theme(axis.title.y=element_blank())
    m <- m + theme(axis.title.x=element_blank())
    m <- m + theme(axis.text.x=element_text(angle = 45, hjust = 1))
    m <- m + geom_hline(yintercept=upperb, linetype="dashed")
    m <- m + geom_hline(yintercept=lowerb, linetype="dashed")
    m <- m + ng1
    fixedatt <- ifelse( att=="singletons", att, paste(att,"Reads"))
    m <- m + labs(title = simpleCap(paste( "Number of",fixedatt,"by Dataset")))
    #plotlist[[i]] <- m
#m <- m + scale_y_log10()
    m <- m + ylab("Variant rate (log)")+ xlab("All rates")
    print(m)
    dev.off()
}

# Coverage
#  - Box MeanCov
#  - Percent over 8x.... somehow
for( att in c("MeanCov","per_above_8x","per_above_15x") ) {
    print(att)
    cdata <- stats.merged[complete.cases(stats.merged[[att]]),]
    outimage <- paste(outdir,"/coverage_",att,"_attsbox.png",sep="")
    print(paste("Writing:",outimage))
    png(outimage)
    #cdata[[att]] <- as.numeric(as.character(cdata[[att]]))
    head(cdata)
    m <- ggplot(cdata, aes_string(x="Dataset", y=att, fill="Dataset"))
    m <- m + geom_boxplot()
    #m <- m + ylab("Frequency")+ xlab("Attribute")
    fixedatt <- ifelse( att=="MeanCov", "Mean Coverage", ifelse(att=="per_above_8x", "Percent Above 8x", "Percent Above 15x"))
    m <- m + labs(title = simpleCap(paste( fixedatt,"by Dataset")))
    m <- m + ng1
    print(m)
    dev.off()
}

for( dataset in unique(stats.merged$Dataset) ) {
    for( att in c("MeanCov","per_above_8x","per_above_15x") ) {
        print(dataset)
        print(att)
        cdata <- stats.merged[stats.merged$Dataset == dataset & complete.cases(stats.merged[[att]]),]
        outimage <- paste(outdir,"/batch_",dataset,"_",att,"_attsbox.png",sep="")
        print(paste("Writing:",outimage))
        png(outimage)
        #cdata[[att]] <- as.numeric(as.character(cdata[[att]]))
        head(cdata)
        m <- ggplot(cdata, aes_string(x="Batch", y=att, fill="Batch"))
        m <- m + geom_boxplot()
        #m <- m + ylab("Frequency")+ xlab("Attribute")
        fixedatt <- ifelse( att=="MeanCov", "Mean Coverage", ifelse(att=="per_above_8x", "Percent Above 8x", "Percent Above 15x"))
        m <- m + labs(title = simpleCap(paste( fixedatt,"by Batch")))
        m <- m + ng1
        m <- m + theme(axis.text.x=element_text(angle = 45, hjust = 1))
        print(m)
        dev.off()
    }
}


# Ethnicities
#  - Pie chart ethnicities by Dataset
#  - Continent pie charts
# Other 
#  - Gender disribution
for( dataset in unique(stats.merged$Dataset) ) {
    for( att in c("ethnicity","continent","country","Consang","Gender") ) {
        print(dataset)
        print(att)
        cdata <- data.frame(table(stats.merged[stats.merged$Dataset == dataset & complete.cases(stats.merged[[att]]) & stats.merged[[att]] != "",][[att]]))
        cdata <- cdata[cdata$Freq > 0,]
        cdata$Var1 <- factor(cdata$Var1, levels=unique(cdata$Var1))    
        #percents <- paste(round(cdata$Freq / sum(cdata$Freq) *100,1),"%",sep="")
        percents <- round(cdata$Freq / sum(cdata$Freq) *100,1)
        cdata$percent <- ifelse( percents > 1.0, paste(percents,"%",sep=""), "")

        percentsize <- ifelse( dim(cdata)[1] < 3, 9, 5)
        outimage <- paste(outdir,"/ident_",dataset,"_",att,"_pie.png",sep="")
        print(paste("Writing:",outimage))
        png(outimage)
        m <- ggplot(cdata, aes(x = "", y=Freq, fill = Var1)) + geom_bar(width = 1,stat="identity")
        m <- m + coord_polar(theta = "y")
        m <- m + geom_text(aes(x=1.3,y = Freq/2 + c(0, cumsum(Freq)[-length(Freq)]), label = percent), size=percentsize)
        m <- m + polar_theme
        #m <- ggplot(stats.merged, aes(x=Dataset, y=att, fill=Dataset))
        #m <- m + geom_boxplot()
        #m <- m + ylab("Frequency")+ xlab("Attribute")
        m <- m + labs(title = simpleCap(paste( dataset,"data by", att)))
        print(m)
        dev.off()
    }
}





