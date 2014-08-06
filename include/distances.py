#!/usr/bin/python
# coding: utf-8

import os, sys
import csv
import commands
import math
import re
import pybedtools
import gzip
from pandas import *
from ggplot import *

from housekeeping import *
from localglobals import *

forceFlag=False
BUF_SIZE=30*1048576

#import rpy2.robjects as robjects
#import pandas.rpy.common as com
#from rpy2.robjects.packages import importr
#from rpy2.robjects.lib import grid
#from rpy2.robjects.lib import ggplot2

#r = robjects.r
#rprint = robjects.globalenv.get("print")
#rstats = importr('stats')
#grdevices = importr('grDevices')
#base = importr('base')

#--------------------------------------------------------------------#
#                               Annotation                           #
#--------------------------------------------------------------------#
#mytheme = {
            #'panel.background':ggplot2.element_rect(fill='white',colour='white'),
            #'panel.grid.major':ggplot2.theme_line(colour = "grey90"),
            #'axis.line':ggplot2.theme_line(size = 1.2, colour="black"),
            #'axis.ticks':ggplot2.theme_line(colour="black"),
            #'axis.text':ggplot2.element_text(colour="black",size=15),
            #'axis.title':ggplot2.element_text(colour="black",size=15),
            #'plot.title':ggplot2.element_text(face="bold", size=20,colour="black"),
            #'panel.grid.minor':ggplot2.theme_line(colour = "NA"),
            #'strip.text.y':ggplot2.element_text(colour="black",face="bold",size=15),
            #'strip.text.x':ggplot2.element_text(colour="black",face="bold",size=15)
            #}


######################################################################
# calculateDistance
######################################################################
def calculateDistance( targetfile, statsfile ):
    print "Running calculateDistance"

    #if os.path.exists( statsfile ) and forceFlag is False :
        #write2log( " - File"+statsfile+" already exists. Skipping command!", True )
        #return statsfile

    comments = []
    header = []
    
    FH = conditionalOpen( targetfile )

    distances = []
    for line in FH :
        if line[:2] == "##" or line.count("\t") < 1:
            comments.append( line )
            continue
        row = line.split("\t")
        if row[0] == "#CHROM" :
            header = row
            distances = [[0,0,0] for i in range(len(row[9:]))]
            continue
        start = int(row[1])
        ref = row[3]
        mut = row[4]
        #reftype = "human"
        #if row[7] == "Chimpref" : reftype = "chimp"
        if ref.find(",") > 0 or mut.find(",") > 0 : continue
        dbsnp = str(row[2] != ".")
        genotypes = [ x[:3] for x in row[9:] ]
        missingness = float(genotypes.count("./.")+genotypes.count(".|.")) / len(genotypes)
        totalgeno = []
        for geno in genotypes:
            if geno.find("/") > 0 : totalgeno += geno.split("/")
            elif geno.find("|") > 0 : totalgeno += geno.split("|")
        #print genotypes
        #print missingness, totalgeno.count('1'), totalgeno.count('0') 
        AF = float(totalgeno.count("1")) / (totalgeno.count("0") + totalgeno.count("1"))
        if missingness > .05 or AF < .01 : continue
        for i in range(len(genotypes)) :
            geno = genotypes[i]
            #totalgeno += geno.split("/")
            gsum = sum([int(x) if is_int(x) else 0 for x in re.split("[\/\|]", geno.strip())])
            #if gsum > 0 : print i,geno,gsum
            distances[i][gsum] += 1

    alldata = DataFrame( distances, columns=["ref","het","hom"])
    alldata["Individual.ID"] = header[9:]

    alldata.to_csv(statsfile, sep="\t",index=False)
    #write2log( "Writing file: "+statsfile, True )
    return
# End calculateDistance

######################################################################
# makeDistancePlot
######################################################################
def makeDistancePlot( alldata, figurename, feature="distance") :
    print "Running makeDistancePlot -",figurename
    alldata["distance"] = alldata.het + alldata.hom

    r_dataframe = com.convert_to_r_dataframe(alldata)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x=feature ) + \
                ggplot2.geom_density(ggplot2.aes_string(fill="factor(Continent)")) + \
                ggplot2.ggtitle("Distance from Reference by Continent") + \
                ggplot2.theme(**mytheme) #+ \
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                #ggplot2.facet_grid( robjects.Formula('RVIS_type ~ .') )

    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# END makeDistancePlot

######################################################################
# makeDistanceBox
######################################################################
def makeDistanceBox( alldata, figurename, factor="Continent", feature="distance", ptitle="Distance from Reference by Continent") :
    print "Running makeDistanceBox -",figurename
    alldata["distance"] = alldata.het + alldata.hom

    r_dataframe = com.convert_to_r_dataframe(alldata)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="factor("+factor+")", y=feature) + \
                ggplot2.geom_boxplot() + \
                ggplot2.ggtitle(ptitle) + \
                ggplot2.theme(**mytheme) + \
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})
                #ggplot2.facet_grid( robjects.Formula('RVIS_type ~ .') )

    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# END makeDistanceBox

######################################################################
# makeDistanceColoredBox
######################################################################
def makeDistanceColoredBox( alldata, figurename, feature="distance", ptitle="Distance from Reference by Continent") :
    print "Running makeDistanceColoredBox -",figurename
    alldata["distance"] = alldata.het + alldata.hom

    newdata = alldata[alldata.Continent.notnull()].sort(['Continent', 'ethnicity'], ascending=[1, 1])
    levels = newdata[["Continent","ethnicity"]].drop_duplicates().sort('Continent')
    #print levels
    levels = newdata.ethnicity.unique().tolist()
    print "Levels",levels

    r_dataframe = com.convert_to_r_dataframe(newdata[["Continent",feature,"ethnicity"]],strings_as_factors=True)
    s = robjects.FactorVector(r_dataframe.rx2("ethnicity"), levels=robjects.StrVector(levels))
    new_r_df = r_dataframe.cbind(s)
    new_r_df.colnames = robjects.StrVector(["Continent",feature,"ethnicity","Ethnicity"])
    #print robjects.r.head(new_r_df)
    p = ggplot2.ggplot(new_r_df) + \
                ggplot2.aes_string(x="Ethnicity", y=feature, fill="factor(Continent)")+ \
                ggplot2.geom_boxplot() + \
                ggplot2.ggtitle(ptitle) + \
                ggplot2.theme(**mytheme) + \
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})
                #ggplot2.facet_grid( robjects.Formula('RVIS_type ~ .') )
                #ggplot2.geom_boxplot(ggplot2.aes_string(group="Continent", order="Continent")) + \
                #ggplot2.geom_boxplot() + \

    grdevices.png(figurename, width=6.25,height=4.25,units="in",res=300)
    p.plot()
    grdevices.dev_off()
# END makeDistanceColoredBox

######################################################################
# makeDistanceColoredBoxGrid
######################################################################
def makeDistanceColoredBoxGrid( alldata, figurename, ptitle="Distance from Reference by Continent") :
    print "Running makeDistanceColoredBoxGrid -",figurename
    varsmelt = melt(alldata,id_vars=["ethnicity","Continent"], value_vars=["hom","het"])

    newdata = varsmelt[varsmelt.Continent.notnull()].sort(['Continent', 'ethnicity'], ascending=[1, 1])
    levels = newdata[["Continent","ethnicity"]].drop_duplicates().sort('Continent')
    levels = newdata.ethnicity.unique().tolist()
    print "Levels",levels

    r_dataframe = com.convert_to_r_dataframe(newdata[["Continent","variable","value","ethnicity"]],strings_as_factors=True)
    s = robjects.FactorVector(r_dataframe.rx2("ethnicity"), levels=robjects.StrVector(levels))
    new_r_df = r_dataframe.cbind(s)
    new_r_df.colnames = robjects.StrVector(["Continent","variable","value","ethnicity","Ethnicity"])
    #print robjects.r.head(new_r_df)
    p = ggplot2.ggplot(new_r_df) + \
                ggplot2.aes_string(x="Ethnicity", y="value", fill="factor(Continent)")+ \
                ggplot2.geom_boxplot() + \
                ggplot2.ggtitle(ptitle) + \
                ggplot2.theme(**mytheme) + \
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                ggplot2.facet_grid( robjects.Formula('variable ~ .'), scale="free_y" )
                #ggplot2.geom_boxplot(ggplot2.aes_string(group="Continent", order="Continent")) + \
                #ggplot2.geom_boxplot() + \

    print "Writing file:",figurename
    grdevices.png(figurename, width=6.25,height=4.25,units="in",res=300)
    p.plot()
    grdevices.dev_off()
# END makeDistanceColoredBoxGrid

######################################################################
# runDistanceAnalysis
######################################################################
def runDistanceAnalysis( targetvcf, force=False ):
    print "Running runDistanceAnalysis -",targetvcf
    sampleannot = read_csv("./resources/annotation/patientannotation.ped",sep="\t")
    hgdp = read_csv("resources/HGDP_ethnicgroups.txt",sep="\t")
    sampleannot = merge(sampleannot, hgdp[["Ethnicity","Continent","Country"]], left_on="ethnicity",right_on="Ethnicity",how="left")

    filepath, basename, suffix = getBasename( targetvcf )
    statsfile = "%s/%s_distance.txt" % (filepath, basename)
    if not os.path.exists(statsfile) : calculateDistance( targetvcf, statsfile )

    print "Plotting!!"
    alldata = read_csv( statsfile, sep="\t")
    alldata = merge(alldata,sampleannot, on="Individual.ID")
    ethcounts = alldata["ethnicity"].value_counts()
    alldata = alldata[~alldata.Continent.isin(["Oceania","NA"]) & ~(alldata.ethnicity.isin(ethcounts[ethcounts == 1].index.tolist()))]
    #print alldata.Continent.unique()

    figurename = "./results/figures/distance/%s_ethnicitybox_final.png" % basename
    makeDistanceColoredBoxGrid( alldata, figurename, ptitle="Zygosity by Ethnicity")
    sys.exit(1)

    
    figurename = "./results/figures/distance/%s_continentdistances.png" % basename
    makeDistancePlot( alldata, figurename )
    figurename = "./results/figures/distance/%s_continentbox.png" % basename
    makeDistanceBox( alldata, figurename )

    figurename = "./results/figures/distance/%s_ethnicitybox.png" % basename
    makeDistanceColoredBox( alldata, figurename, ptitle="Distance from Reference by Ethnicity")

    figurename = "./results/figures/distance/%s_hom_continentbox.png" % basename
    makeDistanceBox( alldata, figurename, feature="hom", ptitle="Homozygous Variants by Region" )

    figurename = "./results/figures/distance/%s_hom_ethnicitybox.png" % basename
    makeDistanceColoredBox( alldata, figurename, feature="hom", ptitle="Homozygosity by Ethnicity" )

    figurename = "./results/figures/distance/%s_het_continentbox.png" % basename
    makeDistanceBox( alldata, figurename, feature="het" , ptitle="Heterozygosity by Region" )

    figurename = "./results/figures/distance/%s_het_ethnicitybox.png" % basename
    makeDistanceColoredBox( alldata, figurename, feature="het", ptitle="Heterozygosity by Ethnicity" )
# END runDistanceAnalysis

######################################################################
# Main
######################################################################
if __name__ == "__main__" :
    os.chdir("..")
    filename = "./rawdata/turks/turks.names.vcf.gz"
    #outdir = "./output"
    outdir = "./liftover/tmp"

    changeLogFile( LOGDIR+"/distance_test.log" )

    targetvcf = "/home/escott/workspace/variome/rawdata/test/everything_set1.chr1.snp.vcf.gz"
    #targetvcf = "/home/escott/workspace/variome/rawdata/variome/variome_snps.vcf.gz"
    #targetvcf = "/home/escott/workspace/variome/rawdata/merged/merged.chimp.regions.filt.vcf.gz"
    #targetvcf = "/home/escott/workspace/variome/rawdata/merged/merged.regions.filt.vcf.gz"
    #targetvcf = "/home/escott/workspace/variome/rawdata/variome1/variome.regions.filt.vcf.gz"
    #targetvcf = "/home/escott/workspace/variome/rawdata/variome1/variome.regions.filt.samp.samp.vcf.gz"
    #targetvcf = "/home/escott/workspace/variome/rawdata/onekg/onekg.chimp.regions.filt.vcf.gz"
    #targetvcf = None
    if len(sys.argv) > 1 :
        targetvcf = sys.argv[1]
    else :
        print "Error: no filename submitted"

    if targetvcf is not None :
        assert os.path.exists(targetvcf)
        runDistanceAnalysis(targetvcf)
    else : 
        for tfile in filesInDir( "/home/escott/workspace/variome/rawdata/*/","*.vcf.gz" ) :
            print tfile
            runDistanceAnalysis(tfile)

    #print "Chimp file:"
    #targetvcf = "/home/escott/workspace/variome/rawdata/merged/merged.chimp.filt.vcf.gz"
    #targetvcf = "/home/escott/workspace/variome/rawdata/test/everything_set1.chr1.snp.chimp.vcf.gz"
    #filepath, basename, suffix = getBasename( targetvcf )
    #statsfile = "%s/%s_distance.txt" % (filepath, basename)
    #if not os.path.exists(statsfile) : calculateDistance( targetvcf, statsfile )

    #alldata = read_csv( statsfile, sep="\t")
    #alldata = merge(alldata,sampleannot, on="Individual.ID")
    #figurename = "./results/figures/%s_continentdistances.png" % basename
    #makeDistancePlot( alldata, figurename )
    #figurename = "./results/figures/%s_continentbox.png" % basename
    #makeDistanceBox( alldata, figurename )
    #figurename = "./results/figures/%s_hom_continentbox.png" % basename
    #makeDistanceBox( alldata, figurename, feature="hom" )
    #figurename = "./results/figures/%s_het_continentbox.png" % basename
    #makeDistanceBox( alldata, figurename, feature="het" )
   
# END MAIN
