#!/usr/bin/python

import os, sys
import string
from pandas import *
import pandas.rpy.common as com
from collections import Counter
from scipy import stats
import math
import csv

from localglobals import *

forceFlag = False

#--------------------------------------------------------------------#
#                               Annotation                           #
#--------------------------------------------------------------------#

######################################################################
# parseRegionsFile
######################################################################
def parseRegionsFile():
    #refseqgenesbed = "./bedfiles/refseq_genes_simple.bed"
    bedfile = "./resources/regionbeds/collab_percentcoverage.txt"

    regions = {}
    for row in csv.reader(open(bedfile), delimiter="\t") :
        if not regions.has_key(row[0]) : regions[row[0]] = []
        regions[row[0]].append([int(row[1]),int(row[2])])
        #print [int(row[1]),int(row[2]),row[3]]
    return regions
# End parseRegionsFile

######################################################################
# intersectRegions
######################################################################
def intersectRegions( chrom, pos, regions ):
    pos = int(pos)
    if len(chrom) < 4 : chrom = "chr"+chrom
    if regions.has_key(chrom) :
        for region in regions[chrom] :
            if region[0] < pos and region[1] > pos :
                #genes.append(region[2])
                return True
    return False
# End intersectRegions

def regioncountPlots( regioncounts, factor="region", figureprefix="results/figures/lof", subtitle="") :
    if len(subtitle) > 0 : subtitle = "\n" + subtitle
    r_dataframe = com.convert_to_r_dataframe(regioncounts)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="factor(%s)"%factor,y="total" ) + \
                ggplot2.geom_boxplot() + \
                ggplot2.ggtitle("Loss of Function counts"+subtitle) + \
                ggplot2.theme(**mytheme) #+ \
                #ggplot2.scale_y_continuous("Internal RVIS") + \
                #ggplot2.stat_smooth(method="lm", se=False)
                #ggplot2.scale_x_continuous("External RVIS")+ \
    loffile = "%s_lof.png" %figureprefix
    print "Writing file %s"%loffile
    grdevices.png(loffile)
    p.plot()
    grdevices.dev_off()

    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="factor(%s)"%factor,y="het" ) + \
                ggplot2.geom_boxplot() + \
                ggplot2.ggtitle("Heterozygous Loss of Function counts"+subtitle) + \
                ggplot2.theme(**mytheme) #+ \
                #ggplot2.scale_y_continuous("Internal RVIS") + \
                #ggplot2.stat_smooth(method="lm", se=False)
                #ggplot2.scale_x_continuous("External RVIS")+ \
    lofhetfile = "%s_lof_het.png"%figureprefix
    print "Writing file %s" % lofhetfile
    grdevices.png(lofhetfile)
    p.plot()
    grdevices.dev_off()

    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="factor(%s)"%factor,y="hom" ) + \
                ggplot2.geom_boxplot() + \
                ggplot2.ggtitle("Homozygous Loss of Function counts"+subtitle) + \
                ggplot2.theme(**mytheme) #+ \
                #ggplot2.scale_y_continuous("Internal RVIS") + \
                #ggplot2.stat_smooth(method="lm", se=False)
                #ggplot2.scale_x_continuous("External RVIS")+ \
    lofhomfile = "%s_lof_hom.png"%figureprefix
    print "Writing file %s" %lofhomfile
    grdevices.png(lofhomfile)
    p.plot()
    grdevices.dev_off()
# end regioncountPlots

def singleDatasetAnalysis(prefix, sampleannot) :
    statsfile = "./results/%s_regioncounts.txt"%prefix
    print "Writing statsfile:",statsfile
    regioncounts = readVariantSet1( statsfile, prefix, sampleannot, vclasses=[5], maxaf=.1, force=True )

    regioncounts.index = regioncounts["Individual.ID"]
    regioncounts = merge(regioncounts,sampleannot[["Individual.ID","Source","Consang"]],on="Individual.ID")
    regioncounts.to_csv(statsfile,sep="\t",index=False)

    print "Total sum:",sum(regioncounts.total.tolist())
    tfactor = "Source"
    regioncountPlots( regioncounts, figureprefix="results/figures/lof/%s_%s"%(prefix,tfactor), factor=tfactor)
    tfactor = "region"
    regioncountPlots( regioncounts, figureprefix="results/figures/lof/%s_%s"%(prefix,tfactor), factor=tfactor)
# END singleDatasetAnalysis

def anotherSingleDatasetAnalysis( varclasses=[4], figureprefix="results/figures/variome.clean" ) :
    prefix = "variome.clean"
    vclass_str = "(%s)" % ",".join([str(x) for x in varclasses])
    statsfile = "./results/%s_regioncounts.txt"%prefix
    print "Writing statsfile:",statsfile
    regioncounts = readVariantSet1( statsfile, prefix, sampleannot, 
                                   vclasses=varclasses, maxaf=.1, force=True )
    regioncounts.index = regioncounts["Individual.ID"]
    regioncounts = merge(regioncounts,sampleannot[["Individual.ID","Source","Consang"]],
                         on="Individual.ID")
    regioncounts.to_csv(statsfile,sep="\t",index=False)
    ME_counts = regioncounts[regioncounts.region == "Middle East"]

    SA_counts = regioncounts[regioncounts.region == "South Asia"]
    CEU_counts = regioncounts[regioncounts.region == "Europe"]
    testdata = concat([CEU_counts,ME_counts, SA_counts]).reset_index()
    badsamps = ["JGG-2079-3-1","JL0557_CAGATC_L005_001","JGG-2087-5-3",
                "JGG-1476-4-1","JGG-1232-2-2","JGG-821-I-3-7","JGG-681-3-3","JGG-1698-3-4"]
    testdata = testdata[~testdata["Individual.ID"].isin(badsamps)]

    print testdata.sort('total',ascending=1).head(10) 
    print testdata.shape
    
    r_dataframe = com.convert_to_r_dataframe(testdata)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="factor(region)" ) + \
                ggplot2.geom_bar() + \
                ggplot2.ggtitle("Sample Distribution\n"+vclass_str) + \
                ggplot2.theme(**mytheme) #+ \
                #ggplot2.scale_y_continuous("Internal RVIS") + \
                #ggplot2.stat_smooth(method="lm", se=False)
                #ggplot2.scale_x_continuous("External RVIS")+ \
    #barplot ="results/figures/%s_pops.png" %(prefix)
    barplot ="%s_pops.png" %(figureprefix)
    print "Writing file %s"%barplot
    grdevices.png(barplot)
    p.plot()
    grdevices.dev_off()


    tfactor = "Source"
    regioncountPlots( testdata, figureprefix="%s_%s"%(figureprefix,tfactor), 
                     factor=tfactor, subtitle=vclass_str)
    tfactor = "region"
    regioncountPlots( testdata, figureprefix="%s_%s"%(figureprefix,tfactor), 
                     factor=tfactor, subtitle=vclass_str)

    r_dataframe = com.convert_to_r_dataframe(testdata)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="hom",y="het", colour="factor(region)" ) + \
                ggplot2.geom_point() + \
                ggplot2.ggtitle("Hom het regressions\n"+vclass_str) + \
                ggplot2.theme(**mytheme) + \
                ggplot2.theme(**{'legend.position':"none"}) + \
                ggplot2.stat_smooth(method="lm") + \
                ggplot2.scale_y_continuous("Heterozygous Burden") + \
                ggplot2.scale_x_continuous("Homozygous Burden")+ \
                ggplot2.facet_grid( robjects.Formula('region ~ .') )
    figureprefix="%s_reg"%(figureprefix)
    loffile = "%s_lof.png" %figureprefix
    print "Writing file %s"%loffile
    grdevices.png(loffile)
    p.plot()
    grdevices.dev_off()
# END anotherSingleDatasetAnalysis

def main_old(sampleannot, prefix) :
    prefix = "daily.clean"
    statsfile = "./results/%s_regioncounts.txt"%prefix
    print "Writing statsfile:",statsfile
    regioncounts = readVariantSet1( statsfile, prefix, sampleannot, vclasses=[5,4], maxaf=.1, force=True )
    regioncounts.index = regioncounts["Individual.ID"]
    regioncounts = merge(regioncounts,sampleannot[["Individual.ID","Source","Consang"]],
                         on="Individual.ID")

    regioncounts.to_csv(statsfile,sep="\t",index=False)

    euro_counts = regioncounts[regioncounts.region == "Europe"]

    for vclass in [1,2,3,4,5] :
        figureprefix = "./results/figures/variome_%d"%vclass
        anotherSingleDatasetAnalysis( [vclass], figureprefix)
# END main_old

################################################################################
# calcRegionAF
################################################################################
def calcRegionAF( vardf, region ) :
    region_totals = vardf[region].apply(lambda x: 
                                       1. if x is np.nan else
                                       sum([int(y) for y in x.split(":")])).map(float)

    region_burden = vardf[region].apply(lambda x: 
                                       0 if x is np.nan else
                                       sum([int(y) for y in x.split(":")[1:]]))
    return region_burden/region_totals
# End calcRegionAF

################################################################################
# Main
################################################################################
if __name__ == "__main__" :

    os.chdir("..")
    #sampleannot = read_csv("resources/annotation/patientannotation.ped", sep="\t")
    #sampleannot.index = sampleannot["Individual.ID"].tolist()
    sampleannot = sampleAnnotation()
    # Make X
    prefix = "everything_set1.chr1.snp.clean"
    prefix = "merged.clean"

    prefix = "test"
    varfile = "/media/data/workspace/variome/rawdata/test/everything_set1.chr1.snp_genes.tsv"
    varfile = "/media/data/workspace/variome/rawdata/test/everything_set1.chr1.snp.clean_genes.tsv"

    vardf = read_csv(varfile,sep="\t")

    print len(vardf["Middle East"].isnull())
    #print vardf[["chrom","pos","FunctionGVS","Priority"]][vardf["Middle East"].isnull()].head()
    vardf["me_af"] = calcRegionAF(vardf,"Middle East")
    vardf["ceu_af"] = calcRegionAF(vardf,"Europe")


    keepcols = ["chrom","pos","FunctionGVS","vclass","me_af","ceu_af"]
    vardf_filt = vardf[keepcols].drop_duplicates()

    r_dataframe = com.convert_to_r_dataframe(vardf_filt)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="me_af", y="ceu_af") +
                ggplot2.geom_point(ggplot2.aes_string(colour="factor(vclass)")) +
                ggplot2.ggtitle("AF comparison between CEU and ME") +
                ggplot2.theme(**mytheme) +
                ggplot2.stat_smooth(method="lm", se=False)+
                ggplot2.scale_y_continuous("CEU AF")+
                ggplot2.scale_x_continuous("ME AF"))

    figname = "results/figures/lof/%s_afcompare.png" % prefix
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()


    #print vardf["Middle East"].head(10)
    #print in_me[:10]

    is_lof = vardf.LOF.notnull()
    is_nmd = vardf.NMD.notnull()
    print vardf["LOF"][is_lof].head(10)
    print is_lof[:10]
    print sum(is_lof)
    print sum(is_nmd)

        
    #singleDatasetAnalysis( prefix, sampleannot )
    #sys.exit(1)

# END Main
    #regioncounts = merge(highimpact, sampleannot[["Individual.ID","Origin","Ethnicity","Continent"]], on="Individual.ID")
    #print "Unfiltered count:",len(highimpact)
    #highimpact = highimpact[highimpact.AF <= .1]
    #print "After AF filter count:",len(highimpact)
    
    #print sampleannot.head(3)
    #print highimpact.head(3)
    #highimpact_annot = merge(highimpact, sampleannot[["Individual.ID","Origin","Ethnicity","Continent"]], left_on="Sample", right_on="Individual.ID")

