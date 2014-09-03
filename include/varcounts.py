#!/usr/bin/env python

#-------------------------------------------------------------------#
#                       Import Packages                             #
#-------------------------------------------------------------------#
import os, sys
import getopt
import csv
import re
import subprocess

import housekeeping as hk
import popgencommands as pop
from localglobals import *

def plotHomozygosity( samplecounts, outdir, basename ): 
    vclasses = ["Class1","Class2","Class3","Class4"]#,"Class5"
    for vclass in vclasses :
        samplecounts[vclass] = (samplecounts[vclass+"_het"] +
                                samplecounts[vclass+"_hom"])
    homcols = [vc+"_hom" for vc in vclasses]
    hetcols = [vc+"_het" for vc in vclasses]
    samplecounts["homtotal"] = samplecounts[homcols].sum(axis=1)
    samplecounts["hettotal"] = samplecounts[hetcols].sum(axis=1)
    samplecounts["allvars"] = samplecounts[homcols+hetcols].sum(axis=1)
    samplecounts["homozygosity"] = samplecounts["homtotal"] / samplecounts["allvars"].map(float)

    r_dataframe = com.convert_to_r_dataframe(samplecounts)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(Region)", y="homozygosity", fill="factor(Region)") +
                ggplot2.geom_boxplot() +
                ggplot2.ggtitle("Boxplot of Homozygosity") +
                ggplot2.theme(**mytheme) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.scale_y_continuous("Homozygosity") +
                ggplot2.scale_x_discrete("Geographic Region") )

    figurename = "%s/%s_homozygosity.png" % (outdir,basename)
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
    print samplecounts[["Sample","Region","Country","GeographicRegions","homtotal","homozygosity"]].head(10)
    scounts = melt(samplecounts,id_vars=["Sample","Region","GeographicRegions","homozygosity","Source"],
                  value_vars=vclasses)

    print "Basename:", basename
    if basename.find( "meceu" ) >= 0:
        #scounts = scounts[scounts.Region.isin(["Europe","Middle East"])] #,"South Asia"])]
        #scounts_sub = scounts[scounts.Source.isin(["BROAD","Daly"])]
        scounts_sub = scounts[((scounts.Region=="Europe") & (scounts.Source=="Daly")) |
                          ((scounts.Region=="Middle East") & (scounts.Source=="BROAD"))] #,"South Asia"
    elif basename.find( "1000G" ) >= 0 :
        print "Using 1000 Genomes pops!"
        scounts_sub = scounts[((scounts.Region=="Europe") & (scounts.Source=="1KG")) |
                              ((scounts.Region=="Africa") & (scounts.Source=="1KG")) |
                          ((scounts.Region=="Middle East") & (scounts.Source=="BROAD"))] #,"South Asia"
    elif basename.find( "variome" ) >= 0 :
        scounts_sub = scounts[scounts.Region.isin(["Europe","Africa","Middle East"])] #,"South Asia"
    elif basename.find( "test2" ) >= 0 :
        scounts_sub = scounts[scounts.Region.isin(["Europe","Africa","Middle East"])] #,"South Asia"
    else :
        print "Error: Unknown basename",basename; sys.exit(1)

    r_dataframe = com.convert_to_r_dataframe(scounts_sub)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="homozygosity", y="value" ) +
                ggplot2.geom_point(ggplot2.aes_string(colour="factor(Source)")) +
                ggplot2.ggtitle("Boxplot of Variants by Variant Class") +
                ggplot2.theme(**mytheme) +
                ggplot2.scale_y_continuous("Class Burden") +
                ggplot2.scale_x_continuous("Sample Homozygosity") +
                ggplot2.facet_grid( robjects.Formula('variable ~ Region') , scale="free") )
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                #ggplot2.scale_x_discrete("Geographic Region") +
                #ggplot2.scale_y_log10() + ggplot2.xlim(0., 1.) + \
                #ggplot2.geom_histogram(breaks=robjects.FloatVector([0,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])) + \
                #ggplot2.aes_string(fill="factor(continent)")
                #fill="factor(group)"

    figurename = "%s/%s_homozburden.png" % (outdir,basename)
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# END plotHomozygosity

def plotVariantClassCounts( scounts, outdir, basename ) :
    r_dataframe = com.convert_to_r_dataframe(scounts)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(Region)", y="value", fill="factor(Region)") +
                ggplot2.geom_boxplot() +
                ggplot2.ggtitle("Boxplot of Variants by Variant Class") +
                ggplot2.theme(**mytheme) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.scale_y_continuous("Variant Burden") +
                ggplot2.scale_x_discrete("Geographic Region") +
                ggplot2.facet_grid( robjects.Formula('variable ~ .') , scale="free") )
                #ggplot2.scale_y_log10() + ggplot2.xlim(0., 1.) + \
                #ggplot2.geom_histogram(breaks=robjects.FloatVector([0,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])) + \
                #ggplot2.aes_string(fill="factor(continent)")
                #fill="factor(group)"

    figurename = "%s/%s_sampcounts.png" % (outdir,basename)
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# ENd plotVariantClassCounts

def plotKnockOutVarCounts( samplecounts, outdir, basename ) :
    vclasses = ["ClassNMD","ClassLOF"]#,"Class5"
    for vclass in vclasses :
        samplecounts[vclass] = (samplecounts[vclass+"_het"] +
                                samplecounts[vclass+"_hom"])
    scounts = melt(samplecounts,id_vars=["Sample","Region","GeographicRegions"],
                  value_vars=vclasses)
    scount = scounts[~scounts.Region.isin(["NA","Oceania"])]
    if basename.find( "meceu" ) >= 0:
        scounts = scounts[scounts.Region.isin(["Europe","Middle East","South Asia"])]
    elif basename.find( "1000G" ) >= 0 :
        scounts = scounts[scounts.Region.isin(["Europe","Africa","Middle East"])] #,"South Asia"
    #print scounts.head(10)

    scounts = scounts[scounts.value < 50]

    r_dataframe = com.convert_to_r_dataframe(scounts)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(Region)", y="value", fill="factor(Region)") +
                ggplot2.geom_boxplot() +
                ggplot2.ggtitle("High Impact Variants by Region") +
                ggplot2.theme(**mytheme) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.scale_y_continuous("Variant Burden" ) +#, limits=robjects.IntVector((0,50))) +
                ggplot2.scale_x_discrete("Geographic Region") +
                ggplot2.facet_grid( robjects.Formula('variable ~ .'), scale="free") )
    figurename = "%s/%s_lofcounts.png" % (outdir,basename)
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# ENd plotKnockOutVarCounts

################################################################################
# plotSampleCounts
################################################################################
def plotSampleCounts( samplecounts, outdir="./results/figures/varcounts", force=False) :

    plotHomozygosity( samplecounts, outdir, basename )

    vclasses = ["Class1","Class2","Class3","Class4"]#,"Class5"
    scounts = melt(samplecounts,id_vars=["Sample","Region","GeographicRegions","Source"],
                  value_vars=vclasses)

    plotVariantClassCounts( scounts, outdir, basename ) 

    plotKnockOutVarCounts( samplecounts, outdir, basename )
# END plotSampleCounts

######################################################################
if __name__ == "__main__":

    #optlist, args = getopt.getopt( sys.argv[1:], "bk")
    #optlist = dict(optlist)
    #bedfile = optlist.get("-b",bedfile)

    # Change running directory
    os.chdir("..")

    #vcffile = "./rawdata/merge1kg/main/me1000G.clean.vcf.gz"
    #vcffile = "./rawdata/mergedaly/meceu.X.vcf.gz"
    #vcffile = "./rawdata/merge1kg/me1000G.X.vcf.gz"
    samplecountsfile = "./rawdata/mevariome/main/variome.clean_indiv.tsv"
    samplecountsfile = "./rawdata/mergedaly/main/meceu.clean_indiv.tsv"

    filepath, basename, ext = hk.getBasename( samplecountsfile )
    samplecounts = read_csv(samplecountsfile, sep="\t")

    samplecounts = addSampleAnnotation( samplecounts, "Sample" )
    print samplecounts.head()

    plotSampleCounts( samplecounts )

# END MAIN
