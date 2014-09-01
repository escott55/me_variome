#!/usr/bin/env python

#-------------------------------------------------------------------#
#                       Import Packages                             #
#-------------------------------------------------------------------#
import os, sys
import getopt
#import gzip
#from random import randint

from housekeeping import *
from localglobals import *
import popgencommands as pop
import patientInfo
import vcftoolsanalyze as vcftools

################################################################################
def getLevels( annotcol, targetpats ):
    print "Running getLevels"
    if annotcol == "GeographicRegions2" :
        levels = target_geographic_regions2
    elif annotcol == "Continent" or annotcol == "Continent2" :
        #levels = target_continents
        levels = target_continents
    elif annotcol == "GeographicRegions" :
        levels = target_geographic_regions
        #if "Oceania" in levels : levels.remove("Oceania")
        #if "America" in levels : levels.remove("America")
        print levels
    else :
        print "Error: unknown levels for column -",annotcol
        sys.exit(1)

    levels = [x.strip().replace("."," ") for x in levels]
    sampleannot = sampleAnnotation(targetpats)
    sampleannot = sampleannot[(sampleannot[annotcol].notnull()) &
                                (sampleannot[annotcol] != "Unknown")]
    sampleannot[annotcol] = [x.strip().replace("."," ") for x in sampleannot[annotcol]]
    finallevels = [x for x in levels if x in sampleannot[annotcol].unique()]
    sampleannot = sampleannot[sampleannot[annotcol].isin(finallevels)]
    return finallevels, sampleannot
# END getLevels

################################################################################
def plotSingletons( singletons, targetdir, prefix="test", factor="GeographicRegions" ) :
    print "Running plotSingletons"
    singletons.loc[:,"count"] = singletons["count"].astype(int)
    singletons_sub = singletons[(singletons[factor].notnull()) 
                            & (singletons[factor] != "Unknown")]
    tlevels, sannot = getLevels( factor, singletons_sub.INDV.tolist() )
    print tlevels
    r_dataframe = com.convert_to_r_dataframe(singletons_sub)
    r_dataframe = fixRLevels( r_dataframe, factor, tlevels )
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor("+factor+")",y="count",fill="factor(Continent2)" ) +
                ggplot2.geom_boxplot() +
                ggplot2.ggtitle("Singletons by "+factor.capitalize()) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})  +
                ggplot2.scale_y_continuous("Variant count") +
                ggplot2.scale_x_discrete(string.capitalize(factor)) +
                ggplot2.theme(**pointtheme_nolegend) +
                ggplot2.facet_grid( robjects.Formula('Vtype ~ .'), scale="free") )
                #ggplot2.geom_boxplot(ggplot2.aes_string(fill="Continent")) \
                 #, notch=True ) + \
    figname = "%s/%s_%s_singletons.png" % (targetdir,prefix,factor)
    print "Writing file:",figname
    grdevices.png(figname)
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

def mergeGeneFiles( tfiles ) :
    print "Running mergeGeneFiles"
    allregions = [y for x in tfiles.regions 
                  for y in x ]
    
    targetfile = "results/other/merged.tsv"
    if os.path.exists(targetfile) : 
        return targetfile

    print allregions
    tcols = ["chrom","pos","FunctionGVS","Gene","vclass"]
    genedata = []
    for grp,row in tfiles.iterrows() :
        if grp == "onekg" : continue
        regions = row["regions"]
        genes = read_csv( row["genes"], sep="\t" )
        genes["Dataset"] = grp
        genedata.append( genes[tcols+allregions] )

    genedata = concat( genedata )

    limit = 100
    merged = genedata.groupby(tcols).agg( sumRegions ).reset_index()

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

################################################################################
# plotClassDist
################################################################################
def plotClassDist( genefile, regionlist, outdir="./results/figures/singletons" ) :

    if genefile.find("meceu") >= 0 :
        regionlist = ["Europe","Middle East"]
    elif genefile.find( "1000G" ) >= 0 :
        regionlist = ["Europe", "Middle East", "Africa"]

    filepath, basename, ext = hk.getBasename( genefile )
    print "Reading file:",genefile
    genedata = read_csv(genefile, sep="\t")
    regionlist = [x for x in regionlist if x in genedata.columns]
    afdist = genedata[["chrom","pos","FunctionGVS","vclass","AF","Polyphen2_HVAR_pred"]+regionlist]
    figurename = "%s/%s_AFdistclasses.png" % (outdir,basename)
    #afdist = DataFrame( afdist, columns=["chrom","pos","FunctionGVS","vclass","AF"]+regionlist )
    keepcols = ["chrom","pos","AF","vclass"] + regionlist
    print afdist.shape
    afdist_filt = (afdist[keepcols]
                   .groupby(["chrom","pos","AF"])
                   .max().reset_index())

    print afdist_filt.shape
    print afdist_filt.head()

    plotAFDistClass( afdist_filt, figurename )

    afdistmelt = melt(afdist_filt,id_vars=["chrom","pos","vclass"],value_vars=regionlist)
    afdistmelt.to_csv("%s/%s_dist.tsv" %(outdir,basename), sep="\t", index=False)

    figurename = "%s/%s_AFdistclasses_grp.png" % (outdir,basename)
    plotAFDistClassGroups( afdistmelt, figurename )

    figurename = "%s/%s_AFdistclasses_density.png" % (outdir,basename)


    print "Running plotAFDistClassDensity"
    print distdf.head(30)
    print distdf[distdf.variable == "Africa"].head()
    distdf_filt = distdf[(distdf.vclass.isin([1,2,3,4])) & (distdf.AF <.05) & (distdf.AF >.0)]

    print distdf_filt.groupby("variable").size().head()

    #hist,bin_edges=np.histogram(alldata.,bins=bins)
    r_dataframe = com.convert_to_r_dataframe(distdf_filt)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="AF", group="factor(variable)", colour="factor(variable)") +
                ggplot2.geom_density() +
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

def singletonWorkflow( singletonfiles ) :
    singletons = []
    for grp, row in singletonfiles.iterrows() :
        print grp
        if grp == "onekg" : continue
        singles = readSingletonFile( row["singletons"], row["genes"] )
        singletons.append( singles )

    #singles1 = readSingletonFile( dalysingletonfile, dalygenesfile )
    #singles2 = readSingletonFile( variomesingletonfile, variomegenesfile )
    #singles3 = readSingletonFile( onekgsingletonfile, onekggenesfile )

    singletons = concat( singletons )#, singles3
    singleton_counts = DataFrame({'count': singletons.groupby(["INDV","Vtype"]).size()}).reset_index()
    singleton_counts = addSampleAnnotation( singleton_counts , "INDV" )
    singleton_counts["Vtype"] = ["Singleton" if x == "S" else "Doubleton" 
                                 for x in singleton_counts.Vtype]
    
    tregions = ["Europe","Northeast Africa","Northwest Africa","Arabian Peninsula","Turkish Peninsula","Syrian Desert","Central Asia"]#,"Africa"
    print hk.dfTable(singleton_counts.GeographicRegions)

    scounts = singleton_counts[(singleton_counts.GeographicRegions.notnull()) &
                                        (singleton_counts.GeographicRegions.isin(tregions))]

    outliers = scounts.sort("count",ascending=False).head(15)["INDV"]
    scounts = scounts[~(scounts.INDV.isin(outliers))]
    targetdir = "results/figures/singletons"
    prefix = "meceu"


    targetfactors = ["Continent2","GeographicRegions", "GeographicRegions2"]

    for annotcol in targetfactors :
        plotSingletons( scounts, targetdir, prefix, factor=annotcol )
# END singletonWorkflow

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
              ["./rawdata/daily/main/daily.clean_genes.tsv",
               "./rawdata/mevariome/main/variome.clean_genes.tsv",
               "./rawdata/onekg/main/onekg.clean_genes.tsv"
              ],
         "regions":
              [["Europe"],
               ["Middle East"],
               ["Africa"]]},
           index=["daily","variome","onekg" ] )

    #singletonWorkflow( tfiles )
  
    mfile = mergeGeneFiles( tfiles )

    plotClassDist( mfile, ["Middle East","Europe","Africa"] )

    #for grp, row in tfiles :
        #plotClassDist( row["genes"], ["Middle East","Europe","Africa"] )
       
# END MAIN
    #dalysingletonfile = "./rawdata/daily/main/vstats/daily.clean.singletons"
    #dalygenesfile = "./rawdata/daily/main/daily.clean_genes.tsv"
    #variomesingletonfile = "./rawdata/mevariome/main/vstats/variome.clean.singletons"
    #variomegenesfile = "./rawdata/mevariome/main/variome.clean_genes.tsv"
#
    #onekgsingletonfile = "./rawdata/onekg/main/vstats/onekg.clean.singletons"
    #onekggenesfile = "./rawdata/onekg/main/onekg.clean_genes.tsv"


