#!/usr/bin/env Rscript

#library(gplots)
library(ggplot2)
library(reshape)
library(grid)
library(gridExtra)

setwd("..")
source("rscripts/themes.R")

simpleCap <- function(x) {
    s <- strsplit(x, "[. ]")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
}


#allannot <- read.table("resources/annotation/patientannotation.ped",sep="\t",header=TRUE)
#festim <- read.table("results/test_inbreedings.out",sep="\t",header=TRUE)
#plink <- read.table("results/plink_ibc.het",header=TRUE)
#print( head(plink) )
#print( head(festim) )

args <- commandArgs( trailingOnly = TRUE )
prefix <- ifelse(length(args)>=1,args[1],"default")

prefix <- "merged.chimp.clean_plink"
#prefix <- "merged.chimp.regions.filt.samp.samp_plink"
ibcdatafile <- paste("results/",prefix,"_continent_ibcrates.txt",sep="")
#countryfile <- paste("results/",prefix,"_country_ibcrates.txt",sep="")

print(paste("Prefix:",prefix))
print(paste("IBC data file:",ibcdatafile))
#print(paste("Country file:",countryfile))

basename <- basename(gsub( "(.*)\\.(.*)", "\\1", ibcdatafile))
filepath <- dirname( ibcdatafile )

consang  <- read.table(ibcdatafile,sep="\t",header=TRUE)
#print(head(consang))
#print(head(melt(consang[,c("F","Origin")])))
#print("bout to cast")
consang$F <- as.numeric(as.character(consang$F))
consang$CountryConsangPercent <- as.numeric(as.character(consang$CountryConsangPercent))
consang$ContConsangPercent<- as.numeric(as.character(consang$ContConsangPercent ))

origin.mean <- aggregate(consang$F, by=list(consang$Origin), FUN=mean, na.rm=TRUE)
#print(head(origin.mean))
consang$Origin <- factor(consang$Origin, levels=as.character(origin.mean$Group.1[order(origin.mean$x)]))
#print(table(consang$Origin))
tokeep <- c("Nigeria","UK","China","Lebanon","Syria","Libya","Yemen","Egypt","Kuwait","Algeria","Iran","Palestine","Qatar","Iraq","Turkey","Saudi Arabia", "Oman","Tunisia","Iran","Pakistan")
scounts <- as.data.frame(table(consang$Origin))
scounts <- scounts[scounts$Var1 %in% tokeep & scounts$Freq > 4,]
#consang_filt <- consang[consang$Continent %in% c("Middle East","South Asia"),]
consang_filt <- consang[consang$Origin %in% scounts$Var1,]

outimage <- paste("results/figures/",prefix,"_country_ibc_violin.png",sep="")
print(outimage)
png(outimage)
m1 <- ggplot(consang_filt, aes(x=Origin, y=F))
m1 <- m1 + geom_boxplot( aes(fill=CountryConsangPercent) )
m1 <- m1 + theme(axis.text.x=element_text(angle = 45, hjust = 1))
m1 <- m1 + scale_fill_continuous(name = 'Consang\nRate')
m1 <- m1 + ng1
m1 <- m1 + labs(title="Country IBC")
m1 <- m1 + xlab("Country")+ ylab("IBC")
print(m1)
dev.off()

print(head(consang))
ethnicity.mean <- aggregate(consang$F, by=list(consang$ethnicity), FUN=mean, na.rm=TRUE)
#consang$ethnicity <- factor(consang$ethnicity, levels=as.character(ethnicity.mean[order(ethnicity.mean$F),]$ethnicity))
consang$ethnicity <- factor(consang$ethnicity, levels=as.character(ethnicity.mean$Group.1[order(ethnicity.mean$x)]))

#consang_filt <- country[country$ethnicity %in% c("Adygei","Bedouin","Druze","Mozabite","Palestinian"),]
outimage <- paste("results/figures/",prefix,"_ethnicity_ibc_violin.png",sep="")
print(outimage)
png(outimage)
m2 <- ggplot(consang_filt, aes(x=ethnicity, y=F))
m2 <- m2 + geom_violin( )
m2 <- m2 + theme(axis.text.x=element_text(angle = 45, hjust = 1))
m2 <- m2 + ng1
m2 <- m2 + labs(title="Ethnicity IBC")
m2 <- m2 + xlab("ethnicity")+ ylab("IBC")
print(m2)
dev.off()
#m2 <- m2 + geom_violin( aes(fill=ConsangPercent) )



#continent <- read.table(continentfile,sep="\t",header=TRUE)
#basename <- basename(gsub( "(.*)\\.(.*)", "\\1", continentfile))
#filepath <- dirname( continentfile )

continent.mean <- aggregate(consang$F, by=list(consang$Continent), FUN=mean, na.rm=TRUE)
print(continent.mean)
#consang$Continent <- factor(consang$Continent, levels=as.character(continent.mean[order(continent.mean$F),]$Continent))
consang$Continent <- factor(consang$Continent, levels=as.character(continent.mean$Group.1[order(continent.mean$x)]))

outimage <- paste("results/figures/",prefix,"_continent_ibc_violin.png",sep="")
print(outimage)
png(outimage)
m3 <- ggplot(consang, aes(x=Continent, y=F))
m3 <- m3 + geom_boxplot( aes(fill=ContConsangPercent) )
m3 <- m3 + theme(axis.text.x=element_text(angle = 45, hjust = 1))
m3 <- m3 + ng1
m3 <- m3 + scale_fill_continuous(name = 'Consang\nRate')
m3 <- m3 + labs(title="Continent IBC")
m3 <- m3 + xlab("Continent")+ ylab("IBC")
print(m3)
dev.off()
#m3 <- m3 + geom_boxplot( )

finalfig <- paste("results/figures/",prefix,"_ibc_final.png",sep="")
print( paste("Writing file:",finalfig) )
png(finalfig, width=7, height=4, units="in", res=400)
grid.arrange(m3,m2, widths=c(1/2,1/2), ncol=2)
dev.off()

#m <- m + geom_hline(yintercept=upperb, linetype="dashed")
#m <- m + geom_hline(yintercept=lowerb, linetype="dashed")
#m <- m + theme(legend.position = "none")
#m <- m + theme(axis.title.y=element_blank())
#m <- m + theme(axis.title.x=element_blank())

#m <- m + labs(title = simpleCap(paste( "Number of",fixedatt,"by Dataset")))
#plotlist[[i]] <- m
#m <- m + scale_y_log10()
