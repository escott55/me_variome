#!/usr/bin/env Rscript

#install.packages("maps")
#install.packages("mapproj")

library(RColorBrewer)
palette(brewer.pal(8, "Dark2"))

library(ggsubplot)
library(plyr)

#library(gplots)
library(ggplot2)
library(reshape)
library(maps)
library(mapproj)

source("themes.R")
#mapdata <- read.table("cassdata.txt", header=FALSE, sep="\t", comment.char="")
#colnames(mapdata) <- c("Country", "Freq" )
#patannot <- read.table("../processeddata/annotation/patientannotation_laser.txt", header=TRUE,sep="\t",comment.char="")
#mapdata <- patannot[patannot$Lab == "Broad" & patannot$Country != "Unknown" 
                    #& patannot$continent %in% c("Middle East","South Asia") & patannot$ethnicity != "None",
                    #colnames(patannot) %in% c("Individual.ID","Country","ethnicity","Unrelated")]
#colnames(mapdata) <- c("Individual.ID","Country","Ethnicity","keepflag")

#print( head(mapdata) )
#mapdata$Country <- as.character(mapdata$Country)
#mapdata$Country[mapdata$Country == "US"] <- "USA"
#mapdata$Country[mapdata$Country == "UAE"] <- "United Arab Emirates"
#mapdata$Country[mapdata$Country == "Korea"] <- "South Korea"
#mapdata$Country[mapdata$Country == "Palestine"] <- "West Bank"
#mapdata <- mapdata[!(mapdata$Country %in% c("India","Italy","Hungary")),]
#print("after")
#mapdata.filt <- mapdata
#mapdata.filt$Country <- factor(mapdata.filt$Country, levels=unique(mapdata.filt$Country) )
#print("after2")
#head(mapdata.filt)
#countries <- as.data.frame(table( mapdata.filt$Country  ))
#countries <- read.table("cassdata.txt", header=FALSE, sep="\t", comment.char="")
#colnames(countries) <- c("Var1", "Freq" )
#countries$Var1[countries$Var1 == "UAE"] <- "United Arab Emirates"
#countries$Var1[countries$Var1 == "Palestine Territory, Occupied"] <- "West Bank"
#countries$Var1[countries$Var1 == "Czech Republic"] <- "Czechoslovakia"
#print(countries)
#print(head(countries))

ME <- c("Egypt","Iraq","Israel","Jordan","Kuwait","Libya","Morocco","Oman","Qatar","Saudi Arabia","Syria","Turkey","United Arab Emirates","UAE","West Bank", "Yemen", "Tunisia","Algeria")
CA <- c("Iran","Pakistan","Afganistan")

mapdata <- read.table("/home/escott/workspace/variome/resources/annotation/patientannotation.ped", sep="\t", comment.char="", header=TRUE)
targetsamples <- read.table("targetsamples.txt",sep="\t")
colnames(targetsamples) <- c("Individual.ID")
mapdata.filt <- mapdata[mapdata$Individual.ID %in% targetsamples$Individual.ID & mapdata$Origin %in% ME,]
mapdata.filt$Source  <- as.character(mapdata.filt$Source)
mapdata.filt$Source[mapdata.filt$Source %in% c('Frazer','BROAD','CIDR','Murat')]  <- "Gleeson"
mapdata.filt$Source <- factor( mapdata.filt$Source, levels=sort(unique(mapdata.filt$Source)) )
mapdata.filt$Origin  <- as.character(mapdata.filt$Origin)
mapdata.filt$Origin[mapdata.filt$Origin == "UAE"] <- "United Arab Emirates"
mapdata.filt$Origin[mapdata.filt$Origin == "Palestine Territory, Occupied"] <- "West Bank"
mapdata.filt$Origin <- factor( mapdata.filt$Origin, levels=sort(unique(mapdata.filt$Origin)) )

#long <- c(-180,-50) 130
#lat <- c(10,80) 70

# Flat Middle East map
world <- map_data("world")
long <- c(-15,80)
lat <- c(12,45)

for( dsource in unique( mapdata.filt$Source)){
    mapdata.filt.filt <- mapdata.filt[mapdata.filt$Source == dsource,]
    countries <- as.data.frame(table( mapdata.filt.filt$Origin ))
    cnames <- aggregate(cbind(long, lat) ~ region, data=world, FUN=function(x)mean(range(x)))
    cnames <- cnames[cnames$region %in% countries$Var1,]
    cnames$region[cnames$region == "United Arab Emirates"] <- "UAE"

   
    missing <- countries$Var1[!(countries$Var1 %in% world$region)]
    world$variable <- countries$Freq[match(world$region,countries$Var1)]
    worldmap <- ggplot(world, aes(x=long, y=lat )) +
        geom_polygon(aes(group=group, fill=variable), colour='black') +
        geom_text(data=cnames, aes(long, lat, label = region), size=4) +
        scale_fill_gradient2( low="white", high="Red", midpoint = 1 ) +
        labs(title = paste(dsource,"Samples By Country of Origin")) +
        coord_cartesian(xlim=long, ylim=lat)

    outfile <- paste("mapfigs/",dsource,"_flatsamplecountmap.png",sep="")
    print(outfile)
    png(outfile, width=7, height=5, units="in", res=100)
    print(worldmap)
    dev.off()
}

q()
# Overlayed Pies

s <- mapdata.filt[mapdata.filt$Country=="Egypt",]
s$Ethnicity <- factor(s$Ethnicity, levels=unique(s$Ethnicity))
pie <- ggplot( s, aes( x=factor(1), fill=factor(Ethnicity) ) ) + geom_bar(width=1)
pie + coord_polar(theta = "y" )

long <- c(5,80)
lat <- c(5,55)

world <- map_data("world")
centres <- ddply(world,.(region),summarize,long=mean(long),lat=mean(lat))

ME <- c("Egypt","Iraq","Israel","Jordan","Kuwait","Libya","Morocco","Oman","Qatar","Saudi Arabia","Syria","Turkey","United Arab Emirates","West Bank", "Yemen")
CA <- c("Iran","Pakistan")
print("Overlay pies")
simdata <- merge(centres,mapdata.filt, by.x="region", by.y="Origin")
simdata.filt <- subset(simdata, region %in% c(ME,CA))
simdata.filt$Ethnicity <- factor( simdata.filt$Ethnicity, levels=unique(simdata.filt$Ethnicity) )
table(simdata.filt[,colnames(simdata.filt) %in% c("region","Ethnicity")])

world$variable <- countries$Freq[match(world$region,countries$Var1)]
worldmap <- ggplot(world, aes(x=long, y=lat )) +
    geom_polygon(aes(group=group, fill=variable), colour='black') +
    coord_cartesian(xlim=long, ylim=lat) +
    labs(title = "Ethnicities By Country") +
    #scale_color_gradient( low="white", high="blue", midpoint = 2 ) 
    #geom_text(data=cnames, aes(long, lat, label = region), size=4) +

#ggplot(simdata.filt) + geom_bar(aes(x=factor(1), fill=Ethnicity), width=1)# + coord_polar(theta="y")

#testplot <- worldmap + (geom_subplot2d(aes(long, lat, group=region,
                              #subplot = geom_bar(aes(factor(1), ..count.., fill=Ethnicity))),
                              #ref=NULL, width = rel(1.2), data = simdata.filt) + coord_polar(theta="y"))

s <- geom_subplot2d(aes(long, lat, group=region,
                              subplot = geom_bar(aes(factor(1), ..count.., fill=Ethnicity))),
                              ref=NULL, width = rel(1.2), data = simdata.filt) 

t <- geom_subplot2d(aes(long, lat, group=region,
                              subplot = geom_bar(aes(Ethnicity, ..count.., fill=Ethnicity))),
                              ref=NULL, width = rel(1.2), data = simdata.filt) 

worldmap + s
                              #subplot = geom_star(aes(0, 0,r=Ethnicity, fill=Ethnicity))),
                              

outfile <- paste("mapfigs/","cass_piecountmap.png",sep="")
print(outfile)
png(outfile)
print(worldmap + s)
dev.off()

#q()

# World 2 map
world <- map_data("world2")
unique(mapdata[!(mapdata$Country %in% unique(world$region)),]$Country)
#world_me <- world
world$variable <- countries$Freq[match(world$region,countries$Var1)]
worldmap <- ggplot(world, aes(x=long, y=lat, group=group)) +
    geom_path(colour="black") +
    scale_y_continuous(breaks=(-2:2) * 30) +
    scale_x_continuous(breaks=(-4:4) * 45) +
    geom_polygon(aes(group=group, fill=variable), colour="black", size=.5, col=NA,lwd=0) +
    scale_fill_gradient2( low="white", high="blue", midpoint = 1 ) +
    labs(title = "Number of Unrelated Samples By Country")
    #, breaks = c(1, 4, 9, 19, 29, 39, 49,100)) 

    #scale_colour_hue(h = c(120, 240)) + 
    #scale_size_manual(values = 2, guide = FALSE)


outfile <- paste("mapfigs/","cass_samplecountmap.png",sep="")
print(outfile)
png(outfile)
worldmap <- worldmap + coord_map("ortho", orientation=c(32, 25, 0))
worldmap
dev.off()

#coord_map()
#q()
# World map, using geom_path instead of geom_polygon
#worldmap <- ggplot(world, aes(x=long, y=lat, group=group)) +
    #geom_path(colour="black")
#worldmap <- worldmap + coord_map("ortho", orientation=c(32, 25, 0))
#print(worldmap)
#world_me <- subset( world, region %in% c("Egypt","Syria","Saudi Arabia","Oman","Yemen","Turkey","Iran","Iraq","Lebanon","Libya","Pakistan","Qatar") )


#world <- map("world", xlim=long, ylim=lat)
#unique(mapdata[!(mapdata$Country %in% unique(world$region)),]$Country)
#world_me <- world
#worldmap <- ggplot(world, aes(x=long, y=lat, group=group)) +
    #geom_path(colour="black") +
    #scale_y_continuous(breaks=(-2:2) * 30) +
    #scale_x_continuous(breaks=(-4:4) * 45) +
    #geom_polygon(aes(group=group, fill=variable), colour="black", size=.5, col=NA,lwd=0) +
    #scale_fill_gradient2( low="white", high="blue", midpoint = 1 ) +
    #labs(title = "Number of Unrelated Samples By Country") +
    #coord_cartesian(xlim=long, ylim=lat) +
    ##geom_text(data=cnames, aes(long, lat, label=region), size=4)
#worldmap

