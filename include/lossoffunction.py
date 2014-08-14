#!/usr/bin/python

import os, sys
import string
#from pandas import *
#import pandas.rpy.common as com
#from collections import Counter
import random
from scipy import stats
from scipy.stats import ttest_ind
import math
import csv
import pybedtools


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

def calcPairwisePvalues( tdf, factcol, groupcol=None, onetailed=False ) :
    pvals = []
    if groupcol is not None :
        for grp in tdf[groupcol].unique():
            lof = tdf[(tdf.Vclass=="LoF") & (tdf[groupcol] ==grp)][factcol].tolist()
            rand = tdf[(tdf.Vclass=="Random") & (tdf[groupcol] ==grp)][factcol].tolist()
            maxY = max(lof+rand)
            ts,pv = ttest_ind(lof, rand)
            #p/2 < alpha and t < 0
            #if onetailed : print "Todo"
            pvals.append( Series(data=[grp,pv,ts,"-2",maxY],
                                 index=[groupcol,"pvalue","T-statistic","x","y"]))
    else :
        lof = tdf[(tdf.Vclass=="LoF")][factcol].tolist()
        rand = tdf[(tdf.Vclass=="Random")][factcol].tolist()
        maxY = max(lof+rand)
        ts,pv = ttest_ind(lof, rand)
        #p/2 < alpha and t < 0
        #if onetailed : print "Todo"
        pvals.append( Series(data=[pv,ts,"-2",maxY],
                                 index=["pvalue","T-statistic","x","y"]))
    pvals = DataFrame(pvals)
    pvals["Pvalue"] = ["P-value: %.3g\nT-statistic: %.2f" %(x,y) 
                       for x,y in pvals[["pvalue","T-statistic"]].values]

    return pvals
# END calcPairwisePvalues

def rvisAnalysis( lofgenes, randgenes, targetdir, prefix ) :
    rvisgenes = read_csv("resources/rvis_data.txt", sep="\t")
    
    lofrvis = merge( lofgenes, rvisgenes, left_on="Gene",right_on="gene")
    lofrvis["Vclass"] = "LoF"
    randrvis = merge( randgenes, rvisgenes, left_on="Gene",right_on="gene")
    randrvis["Vclass"] = "Random"
    tcols = ["Gene","RVIS","RVIS_percent","Vclass"]
    allrvis = concat( [lofrvis[tcols], randrvis[tcols]] ).reset_index()
    
    pvals = calcPairwisePvalues( allrvis, "RVIS" )
    pvals["x"] = -5
    pvals["y"] = 0.5

    r_dataframe = com.convert_to_r_dataframe(allrvis)
    r_pvals = com.convert_to_r_dataframe(pvals)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(Vclass)", y="RVIS", group="factor(Vclass)") +
                ggplot2.geom_jitter(ggplot2.aes_string(colour="factor(Vclass)")) +
                ggplot2.ggtitle("RVIS distribution") +
                ggplot2.scale_y_continuous("RVIS") +
                #ggplot2.scale_x_discrete("ME AF") +
                ggplot2.theme(**mytheme) )
                #ggplot2.stat_smooth(method="lm", se=False)+
    figname = "%s/%s_rvisjitter.png" % (targetdir,prefix)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()

    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="RVIS") +
                ggplot2.geom_density(ggplot2.aes_string(colour="factor(Vclass)")) +
                ggplot2.ggtitle("RVIS") +
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", x="x", y="y"),
                                    hjust=0, vjust=1, size=7, data = r_pvals ) +
                ggplot2.scale_y_continuous("RVIS Density") +
                ggplot2.xlim(-5,5) +
                ggplot2.theme(**mytheme) )
    figname = "%s/%s_rvisdensity.png" % (targetdir,prefix)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()

    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string( x="RVIS" ) +
                ggplot2.stat_density(ggplot2.aes_string(ymax="..density..", ymin="-..density..",
                                                        fill="factor(Vclass)"),
                                     geom="ribbon", position="identity") + #,notch=True 
                #ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", y="0.0", x="y"),
                                 #hjust=0, data = r_pvals ) +
                ggplot2.ggtitle("RVIS distribution") +
                ggplot2.facet_grid(robjects.Formula('Vclass ~ .'), scales="free") +
                ggplot2.coord_flip() +
                ggplot2.theme(**mytheme))

    figname = "%s/%s_rvisribbon.png" % (targetdir,prefix)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END rvisAnalysis 

def psiIntersect_bedtools( vardata, bedfile ) :
    psibedfile = "resources/isoformexpress/globalpsi.bed"

    tvars = vardata[["chrom","pos"]].copy()
    tvars.rename(columns={'chrom':'#chrom'}, inplace=True)
    tvars["#chrom"] = ["chr"+str(x) for x in tvars["#chrom"]]
    tvars["end"] = tvars["pos"]+1
    
    tvars.to_csv(bedfile,sep="\t",index=False)
    psibed = pybedtools.BedTool(psibedfile)
    #resfile = "results/intersections.bed"
    #a.intersect(testexprfile).saveas(resfile)
    inter = psibed.intersect(bedfile)
    print inter
    print type(inter)
    print inter[:3]
    res = DataFrame(inter, 
                    columns=["chrom","start","end","Gene","isocov","isocount","PSI"])
    print res.head()
    return res
# END psiIntersect_bedtools

def psiIntersect( vardata, psiannotatedfile, force=False ) :

    #psiannotatedfile = os.path.join(targetdir, basename+"_psi.tsv")
    if os.path.exists(psiannotatedfile) and not force:
        print "Warning: file exists. skipping..."
        final = read_csv( psiannotatedfile, sep="\t" )

    psibedfile = "resources/isoformexpress/globalpsi.bed"
    psiranges = read_csv(psibedfile, sep="\t")
    psichrs = psiranges.groupby("#chrom")
    psivals = []
    for chrom,pos,gene in vardata[["chrom","pos","Gene"]].values :
        chrdata = psichrs.get_group("chr"+str(chrom)) 
        tregions = chrdata[( psiranges["start"] < pos ) & 
                 (psiranges["end"] > pos )]
        #print chrom,pos
        #print tregions.head()
        if len(tregions) == 1 :
            psivals.append( tregions["Global-PSI"].tolist()[0] )
        elif len(tregions) > 1 :
            print "warning: more than one match"
            if gene in tregions.Gene.tolist() :
                psivals.append( tregions[tregions.Gene == gene]["Global-PSI"].tolist()[0] )
            else :
                psivals.append( tregions["Global-PSI"].tolist()[0] )
        else :
            psivals.append( "-" )

    vardata["PSI"] = psivals
    vardata.to_csv( psiannotatedfile, sep="\t", index=False )
    return vardata
# END psiIntersect

def psiAnalysis( lofvariants, randvariants, targetdir, prefix ) :
    #lofbed = "results/tmp/%s_lofsnps.bed" %prefix
    #randbed = "results/tmp/%s_randsnps.bed" %prefix
    #lofpsi = psiIntersect_bedtools( lofvariants, lofbed )
    #randpsi = psiIntersect_bedtools( randvariants, randbed )

    print lofpsi.shape
    print randpsi.shape
# END psiAnalysis

################################################################################
# Main
################################################################################
if __name__ == "__main__" :

    os.chdir("..")
    #sampleannot = sampleAnnotation()
    # Make X
    #varfile = "/media/data/workspace/variome/rawdata/test/everything_set1.chr1.snp_genes.tsv"
    varfile = "/media/data/workspace/variome/rawdata/test2/main/test2.clean_genes.tsv"
    targetvarfile = hk.copyToSubDir( varfile, "lof" )

    targetdir, basename, suffix = hk.getBasename( targetvarfile )
    prefix = basename[:basename.find("_genes")]

    lofpsiannotatedfile = os.path.join(targetdir, basename+"_lof_psi.tsv")
    randpsiannotatedfile = os.path.join(targetdir, basename+"_rand_psi.tsv")
    print "LoF file:",lofpsiannotatedfile
    print "Rand file:",randpsiannotatedfile
    force = True

    vardf = read_csv(targetvarfile,sep="\t")
    # make Lof Variant Set
    if not os.path.exists( lofpsiannotatedfile ) :
        lofvariants = vardf[(vardf.LOF.notnull()) & (vardf.Priority != "MODIFIER")]
        print lofvariants.shape
        lofvariants["me_af"] = calcRegionAF(lofvariants,"Middle East")
        lofvariants["ceu_af"] = calcRegionAF(lofvariants,"Europe")
        psiIntersect( lofvariants, lofpsiannotatedfile )
    else :
        lofvariants = read_csv( lofpsiannotatedfile, sep="\t" )

    # make Rand Variant Set
    if not os.path.exists( randpsiannotatedfile ) or force:
        randrows = random.sample(vardf.index, len(lofvariants))
        randvariants = vardf.ix[randrows]
        randvariants["me_af"] = calcRegionAF(randvariants,"Middle East")
        randvariants["ceu_af"] = calcRegionAF(randvariants,"Europe")
        psiIntersect( randvariants, randpsiannotatedfile, force )
    else :
        randvariants = read_csv( randpsiannotatedfile, sep="\t" )

    randvariants["Vclass"] = "Random"
    lofvariants["Vclass"] = "LoF"

    lofgenes = DataFrame({'Gene':lofvariants.Gene.unique()})
    randgenes = DataFrame({'Gene':randvariants.Gene.unique()})


    # Intersect with Synthetic Lethal list
    sldata = read_csv("resources/sl_genelist.txt", sep="\t",
                      header=None, names=["Gene","interactions"])

    # Intersect with omim list
    omimgenes = read_csv("resources/omimgenes.txt", sep="\t")

    # Intersect with Kegg Pathway
    kegggenes = read_csv("resources/keggPathway.txt", sep="\t")

    # Intersect with RVIS
    #rvisAnalysis( lofgenes, randgenes, targetdir, prefix )
  
    # Intersect with PSI
    #psiAnalysis( lofvariants, randvariants, targetdir, prefix )
    tcols = ["Gene","PSI","Vclass"]
    allpsi = concat( [lofvariants[tcols], randvariants[tcols]] ).reset_index()
    print allpsi.head()
    allpsi = allpsi[allpsi.PSI != "-"]
    allpsi.PSI = allpsi.PSI.astype("float")
    r_dataframe = com.convert_to_r_dataframe(allpsi)
    #r_pvals = com.convert_to_r_dataframe(pvals)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="PSI") +
                ggplot2.geom_density(ggplot2.aes_string(colour="factor(Vclass)")) +
                ggplot2.ggtitle("PSI distribution") +
                ggplot2.scale_y_continuous("Density") +
                #ggplot2.scale_x_discrete("ME AF") +
                ggplot2.theme(**mytheme) )
                #ggplot2.stat_smooth(method="lm", se=False)+
    figname = "%s/%s_psidensity.png" % (targetdir,prefix)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()



    sys.exit(1)
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

