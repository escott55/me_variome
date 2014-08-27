#!/usr/bin/python

import os, sys
import string
#from pandas import *
#import pandas.rpy.common as com
#from collections import Counter
import random
from scipy import stats
from scipy.stats import ttest_ind
from scipy.stats import gaussian_kde
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
    statsfile = "./results/lof/%s_regioncounts.txt"%prefix
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
    statsfile = "./results/lof/%s_regioncounts.txt"%prefix
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
        finaldf = read_csv( psiannotatedfile, sep="\t" )
        return finaldf

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

def plotPSIClasses( allpsi, targetdir, prefix ):
    print allpsi.head()
    allpsi = allpsi[allpsi.PSI != "-"].drop_duplicates()
    allpsi.PSI = allpsi.PSI.astype("float")
    r_dataframe = com.convert_to_r_dataframe(allpsi)
    #r_pvals = com.convert_to_r_dataframe(pvals)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="PSI") +
                ggplot2.geom_density(ggplot2.aes_string(colour="factor(vclass)"),size=2) +
                #ggplot2.ggtitle("PSI distribution") +
                ggplot2.scale_y_continuous("Density") +
                ggplot2.scale_colour_brewer("Variant Type",palette="Set1") +
                #ggplot2.scale_x_discrete("ME AF") +
                ggplot2.theme(**mytheme) )
                #ggplot2.stat_smooth(method="lm", se=False)+
    figname = "%s/%s_psiclassdensity.pdf" % (targetdir,prefix)
    print "Writing file:",figname
    grdevices.pdf(figname)
    p.plot()
    grdevices.dev_off()
# END plotPSIClasses

def plotPSIByClass(allpsi, targetdir, prefix) :
    allpsi = allpsi[allpsi.PSI != "-"]
    allpsi.PSI = allpsi.PSI.astype("float")
    r_dataframe = com.convert_to_r_dataframe(allpsi)
    #r_pvals = com.convert_to_r_dataframe(pvals)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="PSI") +
                ggplot2.geom_density(ggplot2.aes_string(colour="factor(vclass)"),size=2) +
                ggplot2.ggtitle("PSI distribution") +
                ggplot2.scale_y_continuous("Density") +
                ggplot2.scale_colour_brewer("Variant Type",palette="Set1") +
                #ggplot2.scale_x_discrete("ME AF") +
                ggplot2.theme(**mytheme) )
                #ggplot2.stat_smooth(method="lm", se=False)+
    figname = "%s/%s_psidensity.pdf" % (targetdir,prefix)
    print "Writing file:",figname
    grdevices.pdf(figname)
    p.plot()
    grdevices.dev_off()
# END plotPSIByClass

def flyAnalysis( targetgenes ) :
    flygenes = read_csv( "resources/flybase/gene_summaries.tsv", sep="\t",
                       header=None,names=["flygene","summary"],skiprows=[0])
    flygenes["lethal"] = ["Lethal" if x.find("lethal") >= 0 else "Non"
                          for x in flygenes["summary"]]
    flygenes = (flygenes.sort("lethal")
                    .groupby("flygene").first().reset_index())
    ensembl = (read_csv( "resources/flybase/ensembl2flybase.tsv", sep="\t")
               .groupby(["Drosophila Ensembl Gene ID","Associated Gene Name"])
               .size().reset_index())

    flygenes_annot = merge( flygenes, ensembl, left_on="flygene",
                           right_on="Drosophila Ensembl Gene ID" )
    #print "Fly! all homologues"
    #print hk.dfTable(flygenes.lethal)

    loffly = merge( targetgenes, flygenes_annot, left_on="Gene",
                   right_on="Associated Gene Name")

    justfly = loffly.groupby(["flygene","lethal"]).size().reset_index()
    print "Fly!", len(justfly)
    print hk.dfTable(justfly.lethal)
# END flyAnalysis

def mouseAnalysis( targetgenes ): 
    mgigenes = read_csv("resources/mgi/mgi_allclasses.txt", sep="\t",
                        header=None, names=["Gene","GID","Phenotype","MGI"],
                        dtype={'Gene':str})

    mgigenes["lethal"] = ["Lethal" if x.find("lethal") >= 0 else "Non"
                          for x in mgigenes["Phenotype"]]
    mgigenes_sub = (mgigenes.sort("lethal")
                    .groupby("Gene").first().reset_index())

    lofmgi = merge( targetgenes, mgigenes_sub, on="Gene" )
    print lofmgi.head()

    print "Mouse", len(lofmgi)
    print hk.dfTable(lofmgi.lethal)
    return lofmgi
    #print lofgenes.shape
# END mouseAnalysis

def varsum( gcounts1, gcounts2 ) :
    ref1,het1,hom1 = [int(x) if x!="." else 0  for x in gcounts1.split(":")]
    ref2,het2,hom2 = [int(x) if x!="." else 0  for x in gcounts2.split(":")]
    return "%d:%d:%d" % (ref1+ref2, het1+het2, hom1+hom2)
# END varsum

def bootStrap( lofvars, homvars, tlength, targetdir, prefix ) :
    samplings = []
    ind = np.linspace(0,100,512)
    kde = gaussian_kde( homvars[homvars.PSI != "-"].PSI.map(float).tolist() )
    for boot in range(0,1000) :
        kdesub = gaussian_kde( kde.resample(tlength) )
        kdedf = DataFrame( {"subsetname":"Random%d" % boot, 
                          "Density":kdesub.evaluate(ind), "PSI":ind} )
        samplings.append( kdedf )

    samplings = concat( samplings ).reset_index(drop=True)

    quants = samplings.groupby(["PSI"])["Density"].quantile([.025,.5,.975]).reset_index()
    quants.rename(columns={'level_1':"Quantile",0:"Density"}, inplace=True )
    quants["linetype"] = ["Mean" if x == .5 else "95% threshold" for x in quants.Quantile]

    if "vclass" not in lofvars.columns.tolist() : 
        lofvars["vclass"] == "LoF"

    lofvars_sub = lofvars[lofvars.PSI != "-"].copy()
    lofvars_sub.PSI = lofvars_sub.PSI.astype(float)
    lofdf = []
    for vclass,lofclass in lofvars_sub.groupby("vclass") :
        kde = gaussian_kde( lofvars_sub[lofvars_sub.vclass == vclass].PSI.tolist() )
        tmpdf = DataFrame({"vclass":vclass, "Density":kde.evaluate(ind), "PSI":ind})
        lofdf.append( tmpdf )

    lofdf = concat( lofdf ).reset_index(drop=True)

    rsamplings = com.convert_to_r_dataframe(samplings)
    rlofvars = com.convert_to_r_dataframe(lofdf)
    rquants = com.convert_to_r_dataframe(quants)
    rquants = fixRLevels( rquants,"linetype", ["Mean","95% threshold"] )
    #r_pvals = com.convert_to_r_dataframe(pvals)
    p = (ggplot2.ggplot(rlofvars) +
                ggplot2.aes_string(x="PSI",y="Density") + #,group="vclass"
                ggplot2.geom_line( ggplot2.aes_string(x="PSI", y="Density", group="factor(subsetname)"),
                                  color="grey", data=rsamplings ) +
                ggplot2.geom_line( ggplot2.aes_string(x="PSI",y="Density",linetype="factor(linetype)", 
                                                      group="factor(Quantile)"),
                                  color="black", data=rquants ) +
                ggplot2.geom_line( ggplot2.aes_string(color="factor(vclass)") ) +
                #ggplot2.geom_density(ggplot2.aes_string(colour="factor(vclass)"),size=1.5,color="blue") +
                ggplot2.scale_y_continuous("Density") +
                #ggplot2.scale_x_continuous("PSI") +
                ggplot2.scale_linetype("Confidence Interval") +
                ggplot2.scale_colour_brewer("Variant Type",palette="Set1") +
                #ggplot2.theme(**{'legend.position':"none"}) +
                #ggplot2.ggtitle("PSI distribution") +
                #ggplot2.scale_colour_brewer("Variant Type",palette="Set1") +
                #ggplot2.scale_x_discrete("ME AF") +
                ggplot2.theme(**mytheme) )
                #ggplot2.stat_smooth(method="lm", se=False)+
    figname = "%s/%s_psibootstrap2.pdf" % (targetdir,prefix)
    print "Writing file:",figname
    grdevices.pdf(figname)
    p.plot()
    grdevices.dev_off()
# END bootStrap

def bootStrap2( homvars, tlength, targetdir, prefix ):
    bootvars = []
    for boot in range(0,100) :
        randrows = random.sample(homvars.index, tlength ) #len(lofvariants)
        bootset = homvars.ix[randrows]
        bootset["subsetname"] = "Random%d" % boot
        bootvars.append( bootset )

    bootvars = concat( bootvars ).reset_index(drop=True)
    allpsi = bootvars[bootvars.PSI != "-"].drop_duplicates()
    allpsi.PSI = allpsi.PSI.astype("float")
    r_dataframe = com.convert_to_r_dataframe(allpsi)
    #r_pvals = com.convert_to_r_dataframe(pvals)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="PSI") +
                ggplot2.geom_density(ggplot2.aes_string(group="factor(subsetname)"),
                                     color="grey",size=.8) +
                #ggplot2.ggtitle("PSI distribution") +
                ggplot2.scale_y_continuous("Density") +
                ggplot2.theme(**{'legend.position':"none"}) +
                #ggplot2.scale_colour_brewer("Variant Type",palette="Set1") +
                #ggplot2.scale_x_discrete("ME AF") +
                ggplot2.theme(**mytheme) )
                #ggplot2.stat_smooth(method="lm", se=False)+
    figname = "%s/%s_psibootstrap.pdf" % (targetdir,prefix)
    print "Writing file:",figname
    grdevices.pdf(figname)
    p.plot()
    grdevices.dev_off()
# END bootStrap

def randomSample( vardf, nvars, classname="Random" ) :
    randrows = random.sample(homvars.index, nvars)
    randvariants = homvars.ix[randrows]
    randvariants["subsetname"] = "Random"
    #randvariants["nhomcarriers"] = [sum([int(y) for y in x.split(":")[2:]]) 
                                        #for x in randvariants["MidEastNew"]]

    #randgenes = DataFrame({'Gene':randvariants.Gene.unique()})
    randgenes = randvariants[["Gene","nvars","nhomcarriers"]].groupby("Gene").sum().reset_index()
    return randvariants, randgenes
# END randomSample

def prepareVars( vardf, targetdir, prefix, classname="LoF", tvclass=[4], region="MidEastNew" ) :
    tcolumns = ["chrom","pos","vclass","Gene","PSI","FunctionGVS",region]
    varset = vardf[(vardf.vclass.isin(tvclass)) & (vardf.Priority != "MODIFIER")][tcolumns]
    #varset = vardf[(vardf.LOF.notnull()) & (vardf.Priority != "MODIFIER")]
    varset["ncarriers"] = [sum([int(y) for y in x.split(":")[1:]]) 
                                for x in varset[region]]
    varset["nhomcarriers"] = [sum([int(y) for y in x.split(":")[2:]]) 
                                for x in varset[region]]

    #varset = varset[(varset.ncarriers >= 1)]
    varset = varset[(varset.nhomcarriers > 1)]
    varset["subsetname"] = classname

    varset["nvars"] = 1
    vargenes = varset[["Gene","nvars","nhomcarriers"]].groupby("Gene").sum().reset_index()

    # Write out genes to file
    refseqgenes = read_csv("resources/sqlgenes.txt",sep="\t")
    refseqgenes["Gene"] = [x.strip() for x in refseqgenes.geneSymbol]
    print refseqgenes.head()
    #tmp = merge(lofgenes, refseqgenes[["geneSymbol","entrez"]], left_on="Gene",
                #right_on="entrez", how="left")

    tmp = (merge(vargenes, refseqgenes[["Gene","entrez"]], 
                 on="Gene", how="left").drop_duplicates()
           .sort(["nvars","nhomcarriers"], ascending=[False,False]))
    print tmp[tmp.entrez.isnull()].head(20)
    print "Null entrez:",sum(tmp.entrez.isnull())
    print "Not null entrez:",sum(tmp.entrez.notnull())

    genesfile = os.path.join(targetdir, prefix+"_"+classname+"_genes.tsv")
    tmp.to_csv(genesfile, sep="\t", index=False)

    return varset, vargenes
# END prepareVars

################################################################################
# Main
################################################################################
if __name__ == "__main__" :

    os.chdir("..")
    #sampleannot = sampleAnnotation()
    # Make X
    #varfile = "/media/data/workspace/variome/rawdata/test/everything_set1.chr1.snp_genes.tsv"
    #varfile = "/media/data/workspace/variome/rawdata/test2/main/test2.clean_genes.tsv"
    varfile = "/media/data/workspace/variome/rawdata/mevariome/main/variome.clean_genes.tsv"
    #varfile = "/media/data/workspace/variome/rawdata/merge1kg/main/me1000G.clean_genes.tsv"
    targetvarfile = hk.copyToSubDir( varfile, "lof" )

    targetdir, basename, suffix = hk.getBasename( targetvarfile )
    prefix = basename[:basename.find("_genes")]

    psiannotatedfile = os.path.join(targetdir, prefix+"_psi.tsv")
    print "Writing file:", psiannotatedfile
    vardf = read_csv(targetvarfile,sep="\t")
    vardf = psiIntersect( vardf, psiannotatedfile )
    vardf = vardf[ vardf.vclass > 0 ]
    #print vardf.head(10)
    plotPSIClasses( vardf, targetdir, prefix )

    vardf["MidEastNew"] = [varsum(x,y) for x,y in vardf[["Middle East","South Asia"]].values]
    vardf["meceu"] = [varsum(x,y) for x,y in vardf[["MidEastNew","Europe"]].values]
    
    lofvars, lofgenes = prepareVars( vardf, targetdir, prefix, "LoF", [4], region="MidEastNew" )

    homvars, homgenes = prepareVars( vardf, targetdir, prefix, "allhom", [1,2,3,4], region="MidEastNew" )
    plotPSIClasses( homvars, targetdir, prefix+"_hom" )

    benignvars, benigngenes = prepareVars( vardf, targetdir, prefix, "benign", [1], region="MidEastNew" )

    randvariants, randgenes = randomSample( homvars, len(lofvars), "Random" )

    bootStrap( lofvars, homvars, len(lofvars), targetdir, prefix )

    print "Homvars:",
    print hk.dfTable(homvars.FunctionGVS)
    print "LoFvars:",
    print hk.dfTable(lofvars.FunctionGVS)
    # Intersect with Synthetic Lethal list
    sldata = read_csv("resources/sl_genelist.txt", sep="\t",
                      header=None, names=["Gene","interactions"])

    # Intersect with omim list
    omimgenes = read_csv("resources/omimgenes.txt", sep="\t")

    # Intersect with Kegg Pathway
    kegggenes = read_csv("resources/keggPathway.txt", sep="\t")

    # Intersect with RVIS
    #rvisAnalysis( lofgenes, randgenes, targetdir, prefix )
  
    mouselofgenes = mouseAnalysis( lofgenes )
    #mouserandgenes = mouseAnalysis( randgenes )
    flyAnalysis( lofgenes )
    #flyAnalysis( randgenes )

    tcols = ["Gene","PSI","Lethal"]
    lofmerge = merge( lofvars[["Gene","PSI"]], mouselofgenes[["Gene","lethal"]], on="Gene" )
    lofmerge["Lethal"] = ["Mouse lethal" if x == "Lethal" else "Mouse non-lethal" 
                          for x in lofmerge.lethal]
    print lofmerge.shape
    randvariants["Lethal"] = "Random Set"
    benignvars["Lethal"] = "Benign"
    allpsi = (concat( [lofmerge[tcols], randvariants[tcols], benignvars[tcols]] )
              .drop_duplicates().reset_index())
    allpsi.rename(columns={'Lethal':'vclass'}, inplace=True)
    print allpsi.head()
    plotPSIByClass(allpsi, targetdir, prefix+"_mouse")


    # plot Fly
    tcols = ["Gene","PSI","Lethal"]
    lofmerge = merge( lofvars[["Gene","PSI"]], mouselofgenes[["Gene","lethal"]], on="Gene" )
    lofmerge["Lethal"] = ["Mouse lethal" if x == "Lethal" else "Mouse non-lethal" 
                          for x in lofmerge.lethal]
    lofmerge["vclass"] = lofmerge["Lethal"]
    bootStrap( lofmerge, homvars, len(lofmerge), targetdir, prefix+"_lethal" )
    print lofmerge.shape
    randvariants["Lethal"] = "Random Set"
    allpsi = concat( [lofmerge[tcols], randvariants[tcols]] ).drop_duplicates().reset_index()
    allpsi.rename(columns={'Lethal':'vclass'}, inplace=True)
    print allpsi.head()
    plotPSIByClass(allpsi, targetdir, prefix)

    sys.exit(1)
# END MAIN
    #lofpsiannotatedfile = os.path.join(targetdir, prefix+"_lof_psi.tsv")
    #randpsiannotatedfile = os.path.join(targetdir, prefix+"_rand_psi.tsv")
    #print "LoF file:",lofpsiannotatedfile
    #print "Rand file:",randpsiannotatedfile
    #force = True
#
    ## make Lof Variant Set
    #if not os.path.exists( lofpsiannotatedfile ) :
        #lofvariants = vardf[(vardf.LOF.notnull()) & (vardf.Priority != "MODIFIER")]
        #print lofvariants.shape
        #lofvariants["me_af"] = calcRegionAF(lofvariants,"Middle East")
        #lofvariants["ceu_af"] = calcRegionAF(lofvariants,"Europe")
        #psiIntersect( lofvariants, lofpsiannotatedfile )
    #else :
        #lofvariants = read_csv( lofpsiannotatedfile, sep="\t" )
#
    ## make Rand Variant Set
    #if not os.path.exists( randpsiannotatedfile ) or force:
        #randrows = random.sample(vardf.index, len(lofvariants))
        #randvariants = vardf.ix[randrows]
        #randvariants["me_af"] = calcRegionAF(randvariants,"Middle East")
        #randvariants["ceu_af"] = calcRegionAF(randvariants,"Europe")
        #psiIntersect( randvariants, randpsiannotatedfile, force )
    #else :
        #randvariants = read_csv( randpsiannotatedfile, sep="\t" )
#
    #randvariants["Vclass"] = "Random"
    #lofvariants["Vclass"] = "LoF"
#

    #mgilethal = read_csv("resources/mgi/mgilethal.txt", sep="\t",
                        #header=None, names=["Gene","GID","Phenotype","MGI"],
                        #dtype={'Gene':str})

    # Intersect with PSI
    #psiAnalysis( lofvariants, randvariants, targetdir, prefix )

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

