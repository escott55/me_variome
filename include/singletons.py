#!/usr/bin/env python

#-------------------------------------------------------------------#
#                       Import Packages                             #
#-------------------------------------------------------------------#
import os, sys
import getopt
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import ttest_ind
from scipy.stats import gaussian_kde

#import gzip
#from random import randint

from housekeeping import *
from localglobals import *
import popgencommands as pop
import patientInfo
import vcftoolsanalyze as vcftools

################################################################################
def getLevels( annotcol, sampleannot ):
    print "Running getLevels"
    if annotcol == "GeographicRegions2" :
        levels = target_geographic_regions2
    elif annotcol == "Continent" or annotcol == "Continent2" :
        #levels = target_continents
        #levels = target_continents + ["Africa","Anglo-American"]
        levels = ["Middle East","Anglo-American"]
    elif annotcol == "GeographicRegions" :
        levels = target_geographic_regions
        #if "Oceania" in levels : levels.remove("Oceania")
        #if "America" in levels : levels.remove("America")
        print levels
    elif annotcol == "GeographicRegions3" :
        levels = target_geographic_regions2 + ["Anglo-American"]
    else :
        print "Error: unknown levels for column -",annotcol
        sys.exit(1)

    levels = [x.strip().replace("."," ") for x in levels]
    print "Levels:",levels
    #sampleannot = sampleAnnotation(targetpats)
    print sampleannot.head()
    sampleannot = sampleannot[(sampleannot[annotcol].notnull()) &
                                (sampleannot[annotcol] != "Unknown")]
    sampleannot.loc[:,annotcol] = [x.strip().replace("."," ") for x in sampleannot[annotcol]]
    #print "Before bs filtering:", hk.dfTable( sampleannot[annotcol] )
    finallevels = [x for x in levels if x in sampleannot[annotcol].unique()]
    print "Final levels:",finallevels
    sampleannot = sampleannot[sampleannot[annotcol].isin(finallevels)]
    return finallevels, sampleannot
# END getLevels

################################################################################
def plotSingletons( singletons, sampleannot, targetdir, prefix="test", 
                   factor="GeographicRegions" ) :

    print "Running plotSingletons"
    singletons.loc[:,"count"] = singletons["count"].astype(int)

    print singletons.head()

    pvals = calcPairwisePvalues_ttest( singletons, "Vtype", 
                                      regioncol="Continent2", valuecol="count" )
    pvals.x = "Middle East"
    print "Pvalues:"
    print pvals

    tlevels, sannot = getLevels( factor, sampleannot )
    singletons_sub = singletons[singletons[factor].isin(tlevels)] 
    r_dataframe = com.convert_to_r_dataframe(singletons_sub)
    r_dataframe = fixRLevels( r_dataframe, factor, tlevels )
    r_pvals = com.convert_to_r_dataframe(pvals)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor("+factor+")",y="count" ) +
                ggplot2.geom_boxplot(notch=True) +
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", x="x", y="y"),
                                    hjust=0, vjust=1, data = r_pvals ) + #size=7,
                #ggplot2.ggtitle("Singletons by "+factor.capitalize()) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})  +
                ggplot2.scale_y_continuous("Variant count") +
                ggplot2.scale_x_discrete("") +
                ggplot2.theme(**pointtheme_nolegend) +
                ggplot2.facet_grid( robjects.Formula('Vtype ~ .'), scale="free") )
                #ggplot2.geom_boxplot(ggplot2.aes_string(fill="Continent")) \
                 #, notch=True ) + \ string.capitalize(factor)

    figname = "%s/%s_%s_singletons.pdf" % (targetdir,prefix,factor)
    print "Writing file:",figname
    #grdevices.png(figname)
    grdevices.pdf(figname,width=4,height=7)
    p.plot()
    grdevices.dev_off()
# END plotSingletons

################################################################################
def readSingletonFile( singletonfile, genesfile ):
    singles = read_csv( singletonfile, sep="\t" )
    
    vinfo = read_csv( genesfile, sep="\t" )
    vinfo = (vinfo[["chrom","pos","vclass"]]
                       .groupby(["chrom","pos"])
                       .max().reset_index())

    singletons = merge( singles, vinfo, left_on=["CHROM","POS"],
              right_on=["chrom","pos"], how="left" )
    singletons.rename(columns={'SINGLETON/DOUBLETON':"Vtype"}, inplace=True)
    #singlesannot = addSampleAnnotation( singlesannot, "INDV" )

    return singletons
    #print singlesannot.head()
    #tcols = ["CHROM","POS","SINGLETON/DOUBLETON","INDV","vclass",
             #"Continent2","GeographicRegions","GeographicRegions2"]
    #return singlesannot[tcols]
# END readSingletonFile

################################################################################
def sumRegions( grpdat ):
    if len( grpdat ) == 1: return grpdat
    summed = [0,0,0]
    for row in grpdat :
        a,b,c = [int(x) for x in row.split(":")]
        summed[0] += a
        summed[1] += b
        summed[2] += c

    return "%d:%d:%d" % tuple(summed)
# END sumRegions

def sumAF( regiondat ) :
    ref,het,hom = (0,0,0)
    for rcounts in regiondat :
        if type(rcounts) != str : continue
        a,b,c = [int(x) for x in rcounts.split(":")]
        ref += a
        het += b
        hom += c
    return (2.*hom+het) / (2.*(ref+het+hom))
# END sumAF

def mergeGeneFiles( tfiles, force=False ) :
    print "Running mergeGeneFiles"
    allregions = [y for x in tfiles.regions 
                  for y in x ]
    
    targetfile = "results/other/merged_genes.tsv"
    if os.path.exists(targetfile) and not force : 
        print "Warning: file already exists skipping"
        return targetfile

    print "All target Regions:",allregions
    tcols = ["chrom","pos","FunctionGVS","Gene","vclass","Priority","LOF","NMD"]
    genedata = []
    for grp,row in tfiles.iterrows() :
        #if grp == "onekg" : continue
        regions = row["regions"]
        genes = read_csv( row["genes"], sep="\t" )
        genes["Dataset"] = grp
        subregions = [x for x in genes.columns if x in allregions]
        print subregions
        genedata.append( genes[tcols+subregions] )

    genedata = concat( genedata )

    print genedata.head()
    #genedata.to_csv( targetfile, sep="\t", index=False )
    merged = genedata.groupby(tcols).max().reset_index()
    merged["LOF"] = [x if type(x) == str else "" for x in merged.LOF]
    merged["NMD"] = [x if type(x) == str else "" for x in merged.NMD]

    merged["AF"] = [ sumAF(x) for x in merged[allregions].values ]
    merged = merged[merged.AF > 0.0]
    merged.fillna("0:0:0", inplace=True)
    #merged = genedata.groupby(tcols).agg( sumRegions ).reset_index()

    print "Writing file:",targetfile
    merged.to_csv( targetfile, sep="\t", index=False )

    return targetfile
# END mergeGeneFiles
    #genes1 = read_csv( genesfile1, sep="\t" )
    #genes2 = read_csv( genesfile2, sep="\t" )

    #print genes1.head() 
    #print genes2.head() 
    #regions = ['Africa','America', 'East Asia', 'Europe', 'Middle East', 'South Asia']
    #tcols = ['chrom', 'pos', 'FunctionGVS', 'Gene', 'vclass', "Missing"] + regions
    #'Priority', 'penetrance', 'LOF', 'NMD', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred','Missing',
    #print genes1[tcols].head()
    #genesmerge = merge(genes1[tcols],genes2[tcols], 
                     #on=['chrom', 'pos', 'FunctionGVS', 'Gene', 'vclass'], how="outer")

    #for idx,row in genesmerge :
        #for region in regions :
            #print idx, region


def calcPairwisePvalues_ttest( tdf, groupcol, regioncol="region", valuecol="AF" ) :
    print "Running calcPairwisePvalues"

    pvals = []
    for grp, subdf in tdf.groupby([groupcol]):
        #s = pairwise_tukeyhsd(subdf[factcol], subdf.region)
        data1 = subdf[subdf[regioncol] == "Middle East"][valuecol]

        data2 = subdf[subdf[regioncol] == "Anglo-American"][valuecol]
        maxX = subdf[valuecol].max()
        maxY = subdf[valuecol].max()
        #f_value, p_value = stats.f_oneway(data1, data2)
        ts, pv = ttest_ind( data1, data2 )
        pvals.append( Series(data=[grp,pv,ts,0,maxY],
                             index=[groupcol,"pvalue","T-statistic","x","y"]))
    pvals = DataFrame(pvals)
    pvals["Pvalue"] = ["P-value: %.3g\nT: %.2f" %(x,y) 
                       for x,y in pvals[["pvalue","T-statistic"]].values]
    print pvals.head()
    return pvals
# END calcPairwisePvalues_ttest

def calcPairwisePvalues_density( tdf, groupcol, valuecol="cdf" ) :
    print "Running calcPairwisePvalues"

    pvals = []
    for grp, subdf in tdf.groupby([groupcol]):
        #s = pairwise_tukeyhsd(subdf[factcol], subdf.region)
        data1 = subdf[subdf.region == "Middle East"][valuecol]
        data2 = subdf[subdf.region == "Anglo-American"][valuecol]
        maxX = subdf["AF"].max()
        ks, p_value = stats.ks_2samp(data1, data2)
        pvals.append( Series(data=[grp,p_value,ks,maxX,200],
                             index=[groupcol,"pvalue","D-value","x","y"]))

    pvals = DataFrame(pvals)
    pvals["Pvalue"] = ["P-value: %.3g\nD: %.2f" %(x,y) 
                       for x,y in pvals[["pvalue","D-value"]].values]
    print pvals.head()
    return pvals
# END calcPairwisePvalues_density

def calcAF(genecounts):
    if type(genecounts) != str and math.isnan(genecounts) : return 0
    ref,het,hom = [int(x) for x in genecounts.split(":")]
    if ref+het+hom == 0 : return 0. #print "Err:",genecounts; 
    return (2.*hom+het) / (2.*(ref+het+hom))
# END calcAF

def mergeIndivFiles( tfiles, force=False ):
    print "Running mergeIndivFiles"
    mergedifile = "results/other/merged_indiv.tsv"

    if not force: 
        return mergedifile

    allindiv = []
    for ifile in tfiles.indiv :
        idata = read_csv( ifile, sep="\t" )
        allindiv.append( idata )

    allindiv = concat( allindiv )
    allindiv.to_csv( mergedifile, sep="\t", index=False )
    return mergedifile
# END mergeIndivFiles

def mergeAnnotFiles( tfiles, force ) :
    clustfile = "results/other/merged.annot"

    if not os.path.exists(clustfile) or force :
        allannot = []
        for afile in tfiles.annot :
            adata = read_csv( afile, sep="\t" )
            allannot.append( adata )

        allannot = concat( allannot )
        allannot.to_csv( clustfile, sep="\t", index=False )
    else : allannot = read_csv( clustfile, sep="\t" )
        
    filepats = allannot["Individual.ID"].tolist()

    if os.path.exists(clustfile) :
        sampleannot = read_csv( clustfile, sep="\t" )
        sampleannot.index = sampleannot["Individual.ID"]
        sampleannot = sampleannot.ix[filepats,:]
        mapping_regions = {"NWA":"Middle East", "NEA":"Middle East", "AP":"Middle East",
                           "SD":"Middle East", "TP":"Middle East", "CA":"Middle East",
                           "LWK":"Africa", "YRI":"Africa", "IBS":"Europe", "CEU":"Europe",
                           "TSI":"Europe", "FIN":"Europe", "GBR":"Europe",
                           "CHB":"East Asia", "CHS":"East Asia", "JPT":"East Asia",
                           "Anglo-American":"Anglo-American" }
        sampleannot["Continent2"] = [mapping_regions[x] if mapping_regions.has_key(x)
                                      else "Unknown"
                                      for x in sampleannot.GeographicRegions3]
    else :
        sampleannot = sampleAnnotation( filepats )

    print hk.dfTable( sampleannot.Continent2 )
    print hk.dfTable( sampleannot.GeographicRegions3 )

    return sampleannot
# END mergeIndivFiles

def processGeneFile( genefile, regionlist, outdir, force=False ):
    filepath, basename, ext = hk.getBasename( genefile )
    distfile = "%s/%s_dist.tsv" %(outdir,basename)
    if os.path.exists(distfile) and not force: 
        print "Warning dist file already exists. skipping...",distfile
        return read_csv(distfile,sep="\t")

    print "Reading file:",genefile
    genedata = read_csv(genefile, sep="\t")
    regionlist = [x for x in regionlist if x in genedata.columns]
    print "New regionlist:",regionlist
    afdist = genedata[["chrom","pos","FunctionGVS","vclass"]+regionlist]

    keepcols = ["chrom","pos","vclass"] + regionlist
    afdist_filt = (afdist[keepcols].groupby(["chrom","pos"]).max().reset_index())

    print afdist_filt.shape
    print afdist_filt.head()

    #figurename = "%s/%s_AFdistclasses.png" % (outdir,basename)
    #plotAFDistClass( afdist_filt, figurename )

    afdistmelt = melt(afdist_filt,id_vars=["chrom","pos","vclass"],value_vars=regionlist)
    afdistmelt.rename(columns={"variable":"region","value":"genocounts"},inplace=True)
    afdistmelt["AF"] = [calcAF(x) for x in afdistmelt["genocounts"] ]

    afdistmelt.to_csv(distfile, sep="\t", index=False)
    return afdistmelt
# END processGeneFile

def calculateDensities( distdf ) :
    # Calculate densities
    afdensity = []
    for grp, subdf in distdf.groupby(["vclass","region"]):
        vclass,region = grp
        maxX = subdf.AF.max()
        kde = gaussian_kde( subdf.AF )
        ind = np.linspace(0.,.01,512)
        maxY = max(kde.evaluate(ind))
        tmpdf = DataFrame({"vclass":vclass, "region":region,
                           "Density":kde.evaluate(ind), "AF":ind})
        afdensity.append(tmpdf)

    afdensity = concat(afdensity).reset_index(drop=True)
    return afdensity
# END calculateDensities

################################################################################
# plotClassDist
################################################################################
def plotClassDist( genefile, regionlist, outdir="./results/figures/singletons", force=False ) :
    print "Running plotClassDist",genefile

    if genefile.find("meceu") >= 0 :
        regionlist = ["Anglo-American","Europe","Middle East"]
    elif genefile.find( "1000G" ) >= 0 :
        regionlist = ["Anglo-American","Europe", "Middle East", "Africa"]

    filepath, basename, ext = hk.getBasename( genefile )
    afdist = processGeneFile( genefile, regionlist, outdir, force )

    #figurename = "%s/%s_AFdistclasses_grp.png" % (outdir,basename)
    #plotAFDistClassGroups( afdistmelt, figurename )

    print "Running plotAFDistClassDensity"
    print afdist.head(30)
    distdf = afdist[(afdist.vclass.isin([1,2,3,4])) 
                             & (afdist.AF <.01) & (afdist.AF >0.)]

    distdf = distdf[distdf.region != "Africa"]
    pvals = calcPairwisePvalues_ttest( distdf, "vclass" )
    pvals.index = pvals.vclass

    distdf = distdf[["region","vclass","AF"]].sort(["vclass","AF"])

    afcumsum = (distdf.groupby(by=['region','vclass','AF'])
                .size().groupby(level=[0,1]).cumsum().reset_index())
    afcumsum.rename(columns={0:"cumsum"}, inplace=True)

    #distdf["cumsum"] = distdf.groupby("vclass").AF.cumsum()

    figurename = "%s/%s_AFdistclasses_cumsum.png" % (outdir,basename)
    r_dataframe = com.convert_to_r_dataframe(afcumsum)
    r_pvals = com.convert_to_r_dataframe(pvals)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="AF", y="cumsum") + 
                ggplot2.geom_line( ggplot2.aes_string(colour="factor(region)"), size=1.5 ) +
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", x="x", y="y"),
                                    hjust=1, vjust=0, data = r_pvals ) + #size=7,
                ggplot2.ggtitle("AF Density by Variant Class") +
                ggplot2.facet_grid( robjects.Formula('vclass ~ .'), scale="free" ) +
                ggplot2.theme(**mytheme) )
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()


    afdensity = calculateDensities( distdf )

    afdensity["cdf"] = (afdensity.groupby(['region','vclass'])["Density"].cumsum())
    pvals = calcPairwisePvalues_density( afdensity, "vclass", "Density" )
    #cdf.rename(columns={0:"cumsum"}, inplace=True)

    r_dataframe = com.convert_to_r_dataframe(afdensity)
    r_pvals = com.convert_to_r_dataframe(pvals)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="AF", y="cdf") + 
                ggplot2.geom_line( ggplot2.aes_string(colour="factor(region)"), size=1.5 ) +
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", x="x", y="y"),
                                    hjust=1, vjust=0, data = r_pvals ) + #size=7,
                ggplot2.scale_colour_brewer("Population",palette="Set1") +
                ggplot2.scale_y_continuous("CDF") +
                ggplot2.facet_grid( robjects.Formula('vclass ~ .'), scale="free" ) +
                ggplot2.theme(**{'axis.text.y': ggplot2.element_blank(),
                          'legend.position':'top'}) +
                ggplot2.theme(**mytheme) )
                #ggplot2.ggtitle("CDF by Variant Class") +
    figurename = "%s/%s_AFdistclasses_cdf.pdf" % (outdir,basename)
    print "Making figure:",figurename
    grdevices.pdf(figurename, height=8, width=6)
    p.plot()
    grdevices.dev_off()


    maxY = afdensity.groupby("vclass")["Density"].max().reset_index()
    maxY.rename(columns={"Density":"y"},inplace=True)
    maxY.index = maxY.vclass
    pvals.update(maxY)

    figurename = "%s/%s_AFdistclasses_density2.png" % (outdir,basename)
    r_dataframe = com.convert_to_r_dataframe(afdensity)
    r_pvals = com.convert_to_r_dataframe(pvals)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="AF", y="Density") + 
                ggplot2.geom_line( ggplot2.aes_string(colour="factor(region)"), size=1.5 ) +
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", x="x", y="y"),
                                    hjust=1, vjust=1, data = r_pvals ) + #size=7,
                ggplot2.ggtitle("AF Density by Variant Class") +
                ggplot2.facet_grid( robjects.Formula('vclass ~ .'), scale="free" ) +
                ggplot2.theme(**mytheme) )
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()


    figurename = "%s/%s_AFdistclasses_density.png" % (outdir,basename)
    #hist,bin_edges=np.histogram(alldata.,bins=bins)
    r_dataframe = com.convert_to_r_dataframe(distdf)
    r_pvals = com.convert_to_r_dataframe(pvals)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="AF") + #group="factor(region)", colour="factor(region)"
                ggplot2.geom_density( ggplot2.aes_string(colour="factor(region)") ) +
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", x="x", y="y"),
                                    hjust=1, vjust=0, data = r_pvals ) + #size=7,
                ggplot2.ggtitle("AF Density by Variant Class") +
                ggplot2.facet_grid( robjects.Formula('vclass ~ .'), scale="free" ) +
                ggplot2.theme(**mytheme) )

                #ggplot2.scale_y_log10() + ggplot2.xlim(0., 1.) + \
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                #ggplot2.aes_string(fill="factor(continent)")

    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
    sys.exit(1)
# END plotClassDist

def parseRegionCounts( rlist ) :
    #allref,allhet,allhom = (0,0,0)
    allcounts = []
    for region, regc in rlist.iteritems() :
        ref,het,hom = regc.split(":")
        s = Series({"Region":region,"ref":ref,"het":het,"hom":hom})
        allcounts.append(s)

    allcounts = DataFrame(allcounts)
    return allcounts
# END parseRegionCounts

################################################################################
def singletonWorkflow2( genefile, sampleannot ) :
    print "Running singletonWorkflow"

    tregions = ["Africa","Anglo-American","Middle East"]
    singletons = []
    nline = 0
    header = []
    for row in csv.reader( open(genefile), delimiter="\t" ) :
        if nline == 0 : header = row
        nline += 1
        s = Series( data=rwo, index=header )
        if s.AF > .1 : continue
        allcounts = parseRegionCounts( s[tregions] )
        if allcounts.het.sum() + allcounts.hom.sum() > 1 : continue
        if allcounts.het.sum() == 1 :
            print "found singletons!"


    #for grp, row in tfiles.iterrows() :
        #print grp
        ##if grp == "onekg" : continue
        #singles = readSingletonFile( row["singletons"], row["genes"] )
        #singletons.append( singles )

    #singletons = concat( singletons )#, singles3

    singleton_counts = DataFrame({'count': singletons.groupby(["INDV","Vtype"])
                                  .size()}).reset_index()

    #singleton_counts = addSampleAnnotation( singleton_counts , "INDV" )
    singleton_counts = merge( singleton_counts, sampleannot, left_on="INDV",
                             right_on="Individual.ID" )

    singleton_counts["Vtype"] = ["Singleton" if x == "S" else "Doubleton" 
                                 for x in singleton_counts.Vtype]
    
    tregions = ["Anglo-American","Europe","Northeast Africa","Northwest Africa",
                "Arabian Peninsula","Turkish Peninsula","Syrian Desert",
                "Central Asia","Africa"]

    print hk.dfTable(singleton_counts.Continent2)
    print singleton_counts.columns
    print singleton_counts.head()

    scounts = singleton_counts[(singleton_counts.Continent2.notnull()) &
                                    (singleton_counts.Continent2 != "Unknown")]
                                    #(singleton_counts.GeographicRegions.isin(tregions))]

    outliers = scounts.sort("count",ascending=False).head(15)["INDV"]
    scounts = scounts[~(scounts.INDV.isin(outliers))]
    targetdir = "results/figures/singletons"
    prefix = "merged"

    plotSingletons( scounts, sampleannot, targetdir, prefix, factor="Continent2" )

    #targetfactors = ["Continent2","GeographicRegions3", "GeographicRegions2"]
    targetfactors = ["Continent2"]
    #for annotcol in targetfactors :
        #plotSingletons( scounts, sampleannot, targetdir, prefix, factor=annotcol )
# END singletonWorkflow2

################################################################################
def singletonWorkflow( tfiles, sampleannot ) :
    print "Running singletonWorkflow"

    singletons = []
    for grp, row in tfiles.iterrows() :
        print grp
        #if grp == "onekg" : continue
        singles = readSingletonFile( row["singletons"], row["genes"] )
        singletons.append( singles )

    singletons = concat( singletons )#, singles3
    singleton_counts = DataFrame({'count': singletons.groupby(["INDV","Vtype"])
                                  .size()}).reset_index()

    #singleton_counts = addSampleAnnotation( singleton_counts , "INDV" )
    singleton_counts = merge( singleton_counts, sampleannot, left_on="INDV",
                             right_on="Individual.ID" )

    singleton_counts["Vtype"] = ["Singleton" if x == "S" else "Doubleton" 
                                 for x in singleton_counts.Vtype]
    
    tregions = ["Anglo-American","Europe","Northeast Africa","Northwest Africa",
                "Arabian Peninsula","Turkish Peninsula","Syrian Desert",
                "Central Asia","Africa"]

    print hk.dfTable(singleton_counts.Continent2)
    print singleton_counts.columns
    print singleton_counts.head()

    scounts = singleton_counts[(singleton_counts.Continent2.notnull()) &
                                    (singleton_counts.Continent2 != "Unknown")]
                                    #(singleton_counts.GeographicRegions.isin(tregions))]

    outliers = scounts.sort("count",ascending=False).head(15)["INDV"]
    scounts = scounts[~(scounts.INDV.isin(outliers))]
    targetdir = "results/figures/singletons"
    prefix = "merged"

    plotSingletons( scounts, sampleannot, targetdir, prefix, factor="Continent2" )

    #targetfactors = ["Continent2","GeographicRegions3", "GeographicRegions2"]
    targetfactors = ["Continent2"]
    #for annotcol in targetfactors :
        #plotSingletons( scounts, sampleannot, targetdir, prefix, factor=annotcol )
# END singletonWorkflow
    #singles1 = readSingletonFile( dalysingletonfile, dalygenesfile )
    #singles2 = readSingletonFile( variomesingletonfile, variomegenesfile )
    #singles3 = readSingletonFile( onekgsingletonfile, onekggenesfile )

################################################################################
# Main
################################################################################
if __name__ == "__main__" :
    os.chdir("..")

    tfiles = DataFrame(
        {"singletons":
              ["./rawdata/daily/main/vstats/daily.clean.singletons", 
               "./rawdata/mevariome/main/vstats/variome.clean.singletons",
               "./rawdata/onekg/main/vstats/onekg.clean.singletons"
              ],
         "genes":
              ["./rawdata/daily/main/classify/daily.clean.sfilt2_genes.tsv",
               "./rawdata/mevariome/main/classify/variome.clean.sfilt2_genes.tsv",
               "./rawdata/onekg/main/classify/onekg.clean.sfilt2_genes.tsv"
              ],
         "indiv":
              ["./rawdata/daily/main/classify/daily.clean.sfilt2_indiv.tsv",
               "./rawdata/mevariome/main/classify/variome.clean.sfilt2_indiv.tsv",
               "./rawdata/onekg/main/classify/onekg.clean.sfilt2_indiv.tsv"
              ],
         "annot":
              ["./rawdata/daily/main/clust/daily.clean.annot",
               "./rawdata/mevariome/main/clust/variome.clean.annot",
               "./rawdata/onekg/main/clust/onekg.clean.annot"
              ],
         "regions":
              [["Anglo-American"],
               ["Middle East"],
               ["Africa"]]},
           index=["daily","variome","onekg" ] )

    force=False
    sampleannot = mergeAnnotFiles( tfiles, False )  
    ifile = mergeIndivFiles( tfiles, True )

    mfile = mergeGeneFiles( tfiles, True )
    #print sampleannot.head()
    #print hk.dfTable(sampleannot.Continent2)
    #print hk.dfTable(sampleannot.GeographicRegions3)
    singletonWorkflow( tfiles, sampleannot )
    #singletonWorkflow2( mfile, sampleannot )

    plotClassDist( mfile, ["Middle East","Anglo-American","Africa"], force=force )
    sys.exit(1)
    for grp, row in tfiles :
        if grp == "onekg" : continue
        plotClassDist( row["genes"], ["Middle East","Europe"] )
       
# END MAIN
    #dalysingletonfile = "./rawdata/daily/main/vstats/daily.clean.singletons"
    #dalygenesfile = "./rawdata/daily/main/daily.clean_genes.tsv"
    #variomesingletonfile = "./rawdata/mevariome/main/vstats/variome.clean.singletons"
    #variomegenesfile = "./rawdata/mevariome/main/variome.clean_genes.tsv"
#
    #onekgsingletonfile = "./rawdata/onekg/main/vstats/onekg.clean.singletons"
    #onekggenesfile = "./rawdata/onekg/main/onekg.clean_genes.tsv"


