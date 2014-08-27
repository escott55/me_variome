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
    print singletons.head()
    return singletons
# END readSingletonFile

def mergeGenefiles( genesfile1, genefiles2 ) :
    genes1 = read_csv( genesfile1, sep="\t" )
    genes2 = read_csv( genesfile2, sep="\t" )

    print genes1.head() 
    print genes2.head() 
    regions = 'Africa','America', 'East Asia', 'Europe', 'Middle East', 'South Asia']
    tcols = ['chrom', 'pos', 'FunctionGVS', 'Gene', 'vclass', "Missing"] + regions
    #'Priority', 'penetrance', 'LOF', 'NMD', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred','Missing',
    print genes1[tcols].head()
    genesmerge = merge(genes1[tcols],genes2[tcols], 
                     on=['chrom', 'pos', 'FunctionGVS', 'Gene', 'vclass'], how="outer")

    for idx,row in genesmerge :
        for region in regions :
        print idx
# END mergeGenefiles

################################################################################
# Main
################################################################################
if __name__ == "__main__" :
    os.chdir("..")

    dalysingletonfile = "./rawdata/daily/main/vstats/daily.clean.singletons"
    dalygenesfile = "./rawdata/daily/main/daily.clean_genes.tsv"
    variomesingletonfile = "./rawdata/mevariome/main/vstats/variome.clean.singletons"
    variomegenesfile = "./rawdata/mevariome/main/variome.clean_genes.tsv"

    singles1 = readSingletonFile( dalysingletonfile, dalygenesfile )
    singles2 = readSingletonFile( variomesingletonfile, variomegenesfile )

    singletons = concat( [singles1, singles2] )
    singletons.rename(columns={'SINGLETON/DOUBLETON':"Vtype"}, inplace=True)
    singleton_counts = DataFrame({'count': singletons.groupby(["INDV","Vtype"]).size()}).reset_index()
    singleton_counts = addSampleAnnotation( singleton_counts , "INDV" )
    singleton_counts["Vtype"] = ["Singleton" if x == "S" else "Doubleton" 
                                 for x in singleton_counts.Vtype]
    
    tregions = ["Europe","Northeast Africa","Northwest Africa","Arabian Peninsula","Turkish Peninsula","Syrian Desert","Central Asia","Africa"]
    scounts = singleton_counts[(singleton_counts.GeographicRegions.notnull()) &
                                        (singleton_counts.GeographicRegions.isin(tregions))]

    outliers = scounts.sort("count",ascending=False).head(15)["INDV"]
    scounts = scounts[~(scounts.INDV.isin(outliers))]
    targetdir = "results/figures/singletons"
    prefix = "meceu"


    targetfactors = ["Continent2","country", "GeographicRegions", "GeographicRegions2"]

    for annotcol in targetfactors :
        plotSingletons( scounts, targetdir, prefix, factor=annotcol )

    #optlist, args = getopt.getopt( sys.argv[1:], "r:ot")
    #optlist = dict(optlist)
#
    #dataset = optlist.get("-r",None)
    #if dataset is None and not optlist.has_key("-t"):
        #print "Error: no dataset provided"
        #print " -r <string> with the name of the dataset"
        #print " -o flag to overwrite everything"
        #print " -t run a test"
        #sys.exit(1)

    #path = os.path.abspath("./rawdata/")
    #if dataset == "test2" :
        #singletonfile = "./rawdata/test2/main/vstats/test2.clean.singletons"
        #genesfile = "./rawdata/test2/main/test2_genes.tsv"
        #vcffile = "./rawdata/test2/main/test2.clean.vcf.gz"
#
    #elif dataset == "mevariome" :
        #singletonfile = "./rawdata/mevariome/main/vstats/variome.clean.singletons"
        #genesfile = "./rawdata/mevariome/main/variome.clean_genes.tsv"
        #vcffile = "./rawdata/mevariome/main/variome.clean.vcf.gz"
#
    #elif dataset == "mergedaly" :
        #singletonfile = "./rawdata/mergedaly/main/vstats/mecue.clean.singletons"
        #genesfile = "./rawdata/mergedaly/main/meceu.clean_genes.tsv"
        #vcffile = "./rawdata/mergedaly/main/meceu.clean.vcf.gz"
#
    #elif dataset == "daily" :
        #singletonfile = "./rawdata/daily/main/vstats/daily.clean.singletons"
        #genesfile = "./rawdata/daily/main/daily.clean_genes.tsv"
        #vcffile = "./rawdata/daily/main/daily.clean.vcf.gz"
#
    #elif dataset == "merge1kg" :
        #singletonfile = "./rawdata/merge1kg/main/vstats/me1000G.clean.singletons"
        #genesfile = "./rawdata/merge1kg/main/me1000G.clean_genes.tsv"
        #vcffile = "./rawdata/merge1kg/main/me1000G.clean.vcf.gz"
#
    #else :
        #print "Error: unknown dataset",dataset
        #sys.exit(1)
   
       
# END MAIN
