#!/usr/bin/python
# coding: utf-8

import os, sys
import csv
import commands
import getopt
import math
import re
import pybedtools
import gzip
from pandas import *
from ggplot import *

from housekeeping import *

forceFlag=False
BUF_SIZE=30*1048576

import rpy2.robjects as robjects
import pandas.rpy.common as com
from rpy2.robjects.packages import importr
from rpy2.robjects.lib import grid
from rpy2.robjects.lib import ggplot2


r = robjects.r
rprint = robjects.globalenv.get("print")
rstats = importr('stats')
grdevices = importr('grDevices')
base = importr('base')

#--------------------------------------------------------------------#
#                               Annotation                           #
#--------------------------------------------------------------------#
mytheme = {
            'panel.background':ggplot2.element_rect(fill='white',colour='white'),
            'panel.grid.major':ggplot2.theme_line(colour = "grey90"),
            'axis.line':ggplot2.theme_line(size = 1.2, colour="black"),
            'axis.ticks':ggplot2.theme_line(colour="black"),
            'axis.text':ggplot2.element_text(colour="black",size=15),
            'axis.title':ggplot2.element_text(colour="black",size=15),
            'plot.title':ggplot2.element_text(face="bold", size=20,colour="black"),
            'panel.grid.minor':ggplot2.theme_line(colour = "NA"),
            'strip.text.y':ggplot2.element_text(colour="black",face="bold",size=15),
            'strip.text.x':ggplot2.element_text(colour="black",face="bold",size=15)
            }


######################################################################
# makeDistancePlot
######################################################################
def makeDistancePlot( alldata, figurename, feature="distance") :
    alldata["distance"] = alldata.het + alldata.hom

    r_dataframe = com.convert_to_r_dataframe(alldata)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x=feature ) + \
                ggplot2.geom_density(ggplot2.aes_string(fill="factor(continent)")) + \
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
def makeDistanceBox( alldata, figurename, feature="distance") :
    alldata["distance"] = alldata.het + alldata.hom

    r_dataframe = com.convert_to_r_dataframe(alldata)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="factor(continent)", y=feature) + \
                ggplot2.geom_boxplot() + \
                ggplot2.ggtitle("Distance from Reference by Continent") + \
                ggplot2.theme(**mytheme) #+ \
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                #ggplot2.facet_grid( robjects.Formula('RVIS_type ~ .') )

    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# END makeDistanceBox

######################################################################
# readVariantSet1
######################################################################
def readVariantSet1( statsfile, prefix, sampleannot, aflimit=1, vclasses=["HIGH"], force=False ):
    print "Running readVariantSet1"
    limit = 10
    if os.path.exists(statsfile) and not force: return read_csv(statsfile, sep="\t")
    #for var in csv.reader(open(targetfile),delimiter="\t"):
    knownsamples = sampleannot.index.tolist()
    allvars = []
    for vclass in vclasses :
        regioncounts = {}
        targetfile = "./results/classes/"+prefix+".annot_"+vclass+".tsv"
        print "Reading:",targetfile
        FM = open(targetfile,"r")
        assert os.path.exists(targetfile)

        for line in FM :
            if len(line) > 200 : continue
            var = line.rstrip().split("\t")
            if var[0] == "Sample" or var[0] not in knownsamples : continue
            if len(var) == 1 : print limit,var; sys.exit(1)
            region = sampleannot.loc[var[0]].continent
            sample = var[0]
            seen.append(key)
            af = var[7]
            if float(af) >= aflimit : continue
            #print var
            #print region, sample
            geno = 0 if var[1] == "het" else 1
            if not regioncounts.has_key(region) : regioncounts[region] = {}
            if not regioncounts[region].has_key(sample) : regioncounts[region][sample] = [0,0]

            regioncounts[region][sample][geno] += 1
            #limit += 1
            #if limit == 10 : break

        classvars = []
        for region in regioncounts :
            for samp in regioncounts[region] :
                print [samp,vclass,region,regioncounts[region][samp][0],regioncounts[region][samp][1],sum(regioncounts[region][samp])]
                classvars.append([samp,vclass,region,regioncounts[region][samp][0],regioncounts[region][samp][1],sum(regioncounts[region][samp])])
        allvars.append( DataFrame( classvars, columns=["Individual.ID","class","region","het","hom","total"]) )
    print "Length:",len(allvars)
    if len(allvars) == 1 : allvars = allvars[0]
    else : allvars = concat( allvars )
    #print allvars.head(10)
    print "Writing:",statsfile
    allvars.to_csv(statsfile, sep="\t", index=False)
    return allvars
# END readVariantSet1

######################################################################
# readVariantSet2
######################################################################
def readVariantSet2( statsfile, prefix, sampleannot, vclasses=["HIGH"], force=False ):
    print "Running readVariantSet2"
    limit = 10
    if os.path.exists(statsfile) and not force: return read_csv(statsfile, sep="\t")
    #for var in csv.reader(open(targetfile),delimiter="\t"):

    knownsamples = sampleannot.index.tolist()
    regioncounts = {}
    seen = []
    for vclass in vclasses :
        targetfile = "./results/classes/"+prefix+".annot_"+vclass+".tsv"
        print "Reading:",targetfile
        FM = open(targetfile,"r")
        assert os.path.exists(targetfile)

        for line in FM :
            if len(line) > 200 : continue
            var = line.rstrip().split("\t")
            if var[0] == "Sample" or var[0] not in knownsamples : continue
            if len(var) == 1 : print limit,var; sys.exit(1)
            region = sampleannot.loc[var[0]].continent
            af = float(var[7])
            key = var[2]+":"+var[3]
            print "key",key
            if key in seen : continue
            if not regioncounts.has_key(region) : regioncounts[region] = {}
            if not regioncounts[region].has_key(af) : regioncounts[region][af] = 0
            regioncounts[region][af] += 1

    alldata = []
    for region in regioncounts:
        for af in sorted(regioncounts[region]) :
            alldata.append([region,af,regioncounts[region][af]])
    alldata = DataFrame( alldata, columns=["Region","AF","count"] )
    return alldata
# END readVariantSet2

######################################################################
# Main
######################################################################
if __name__ == "__main__" :
    os.chdir("..")
    #filename = "./rawdata/turks/turks.names.vcf.gz"
    #outdir = "./output"
    outdir = "./liftover/tmp"

    changeLogFile( LOGDIR+"/distance_test.log" )

    #targetvcf = "/home/escott/workspace/variome/rawdata/test/everything_set1.chr1.snp.vcf.gz"
    #targetvcf = "/home/escott/workspace/variome/rawdata/variome/variome_snps.vcf.gz"
    #targetvcf = "/home/escott/workspace/variome/rawdata/merged/merged.chimp.regions.filt.vcf.gz"
    #targetvcf = "/home/escott/workspace/variome/rawdata/merged/merged.regions.filt.vcf.gz"
    #targetvcf = "/home/escott/workspace/variome/rawdata/variome1/variome.regions.filt.vcf.gz"
    #targetvcf = "/home/escott/workspace/variome/rawdata/variome1/variome.regions.filt.samp.samp.vcf.gz"
    #targetvcf = "/home/escott/workspace/variome/rawdata/onekg/onekg.chimp.regions.filt.vcf.gz"

    #prefix = "onekg.chimp.regions.recode.recode"
    prefix = "everything_set1.chr1.snp.chimp.regions.filt.samp.samp"
    targetvcf = "./results/classes/"+prefix+".annot_HIGH.tsv"

    if len(sys.argv) > 1 :
        targetvcf = sys.argv[1]
    else :
        print "Error: no filename submitted"

    sampleannot = read_csv("./resources/annotation/patientannotation.ped",sep="\t")
    sampleannot.index = sampleannot["Individual.ID"].tolist()
    #print sampleannot.head(4)

    filepath, basename, suffix = getBasename( targetvcf )
    statsfile = "%s/%s_summary.txt" % (filepath, prefix)
    #if not os.path.exists(statsfile) : calculateDistance( targetvcf, statsfile )

    #if not os.path.exists(statsfile) : 
    #else : alldata = readVariantSet1( statsfile, prefix, ["HIGH"] )
    #alldata = readVariantSet1( statsfile, prefix, ["HIGH","MODERATE","LOW"] )
    #alldata = readVariantSet1( statsfile, prefix, sampleannot, .1, ["HIGH","MODERATE","LOW","MODIFIER"], force=False )
    alldata = readVariantSet1( statsfile, prefix, sampleannot, .1, ["HIGH","MODERATE","LOW","MODIFIER"], force=False )
    
    print alldata.shape
    alldata = alldata[alldata["class"] == "LOW"]
    print alldata.shape
    print "Plotting!!"
    #alldata = read_csv( statsfile, sep="\t")
    alldata = merge(alldata,sampleannot, on="Individual.ID")
    #figurename = "./results/figures/%s_synondistances.png" % basename
    #makeDistancePlot( alldata, figurename )
    figurename = "./results/figures/%s_synonbox.png" % basename
    print "Making plot:",figurename
    makeDistanceBox( alldata, figurename, feature="total" )

    # Cumulative proportion of variance
    alldata = readVariantSet2( statsfile, prefix, sampleannot, ["HIGH","MODERATE","LOW","MODIFIER"], force=False )
    print alldata.head(10)
# END MAIN

