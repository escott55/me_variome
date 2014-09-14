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


def plotHWEafBinned(hwedf, outdir, basename, tregions=["Europe","Middle_East"] ):
    bins=[0.,.005,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]

    print hk.dfTable(hwedf.Region)
    hwebinned = hwedf[hwedf.Region.isin(tregions)][["Region","AF","P"]]
    hwebinned["bins"] = cut(hwebinned.AF, bins)

    print hwebinned.head(10)
    r_dataframe = com.convert_to_r_dataframe(hwebinned)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(bins)", y="P", fill="factor(Region)") +
                ggplot2.geom_boxplot() + \
                ggplot2.ggtitle("HWE Probability by AF") +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45),
                                 'legend.position':"none"}) +
                ggplot2.scale_x_discrete("HWE Binned") +
                ggplot2.scale_y_continuous("HWE") +
                ggplot2.facet_grid( robjects.Formula('Region ~ .') ) +
                ggplot2.theme(**mytheme))
                #ggplot2.scale_y_log10() + ggplot2.xlim(0., 1.) + \
                #ggplot2.geom_histogram(breaks=robjects.FloatVector([0,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])) + \
                #ggplot2.aes_string(fill="factor(continent)")

    figurename = "%s/%s_hweaf.png" % (outdir, basename)
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# END plotHWEafBinned

def plotHWEafRibbon(hwedf, outdir, basename, tregions=["Europe","Middle_East"] ):
    bins=[0.,.005,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]

    print hk.dfTable(hwedf.Region)
    hwebinned = hwedf[hwedf.Region.isin(tregions)][["Region","AF","P"]]
    hwebinned["bins"] = cut(hwebinned.AF, bins)

    print hwebinned.head(10)
    r_dataframe = com.convert_to_r_dataframe(hwebinned)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string( x="P" ) +
                ggplot2.stat_density(ggplot2.aes_string(ymax="..density..", ymin="-..density..",
                                                        fill="factor(Region)"),
                                     geom="ribbon", position="identity") + #,notch=True 
                ggplot2.ggtitle("HWE Probability by AF") +
                ggplot2.facet_grid(robjects.Formula('Region ~ bins'), scales="free") +
                ggplot2.scale_x_discrete("HWE Binned") +
                ggplot2.scale_y_continuous("HWE") +
                ggplot2.coord_flip() +
                ggplot2.theme(**mytheme))


    figurename = "%s/%s_hweribbon.png" % (outdir, basename)
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# END plotHWEafRibbon

def calcAF( vcounts ) :
    if vcounts.find(".") >= 0 :
        return 0
    elif vcounts.find("/") >= 0 : 
        ref,het,hom = [int(x) for x in vcounts.split("/")]
    elif vcounts.find(":") >= 0 :
        ref,het,hom = [int(x) for x in vcounts.split(":")]
    else : 
        print "Error: unknown delimiter",vcounts; sys.exit(1)
    return (2*hom+het) / (2.*(ref+het+hom))
# END calcAF

def hwTernaryPlot( targetdf, targetdir, prefix, vmax=10000 ):
    print "Running hwTernaryPlot"
    tmpfile = "%s/%s.geno" % (targetdir, prefix)
    print "Tempfile:",tmpfile
    ternaryplot = "%s/%s_ternary.png" % (targetdir, prefix)
    qqplot = "%s/%s_qq.png" % (targetdir, prefix)
    genoplot = "%s/%s_geno.png" % (targetdir, prefix)
    clrplot = "%s/%s_geno.png" % (targetdir, prefix)
    #if len(targetdf) > vmax :
        #targetdf[["AA","AB","BB"]].head(vmax).to_csv(tmpfile, index=False)
    #else :
    targetdf[["AA","AB","BB"]].drop_duplicates().to_csv(tmpfile, index=False)

    rcmd = ("library(HardyWeinberg)\n"+
            'genodata<-read.csv("'+tmpfile+'")\n'+
            'print(head(genodata))\n'
            'png("'+ternaryplot+'")\n'+
            'HWTernaryPlot(genodata)\n'+
            'dev.off()\n' +
            'png("'+genoplot+'")\n'+
            'HWTernaryPlot(as.matrix(genodata), nsim=30)\n'+
            'dev.off()\n' +
            'png("'+clrplot+'")\n'+
            'HWClrPlot(as.matrix(genodata))\n'+
            'dev.off()\n' +
            'png("'+qqplot+'")\n'+
            'HWGenotypePlot(genodata)\n'+
            'dev.off()\n' )
            #'genodata <- UniqueGenotypeCounts(genodata)\n'+
    out = robjects.r( rcmd )
# END hwTernaryPlot

def plotHardy( filedf, genesfile, prefix, targetdir ) :
    variantinfo = read_csv( genesfile, sep="\t" )
    variantinfo_sub = (variantinfo[["chrom","pos","vclass"]]
                       .groupby(["chrom","pos"])
                       .max().reset_index())
    allpvals = []
    for index, row in filedf.iterrows() :
        print row
        print row["pop"], row["hardyfile"]
        hardydata = read_csv(row["hardyfile"], sep="\t" )

        hardyfilt = hardydata[(hardydata.P > 0) & (hardydata.ChiSq.notnull())]
        hardyfilt["AF"]  = [calcAF(x) for x in hardyfilt["OBS(HOM1/HET/HOM2)"]]
        hardyannot = merge(hardyfilt, variantinfo_sub, left_on=["CHR","POS"], 
                           right_on=["chrom","pos"], how="left")
        hardyannot_sub = hardyannot[hardyannot.vclass.isin( [1,2,3,4] )]
        print "Before:",len(hardydata), "After:",len(hardyfilt)
        print "annot",len(hardyannot), "Annot filt:",len(hardyannot_sub)
        print hardyannot_sub.head(10)
        hardyannot_sub["Region"] = row['pop']
        allpvals.append(hardyannot_sub[["chrom","pos","AF","vclass","P","Region",
                                        "OBS(HOM1/HET/HOM2)","E(HOM1/HET/HOM2)"]])
        #allpvals.append(DataFrame({'pval':hardyannot_sub.P, 'grp':row['pop']}))
        #allpvals.append(DataFrame({'pval':hardyfilt.P, 'grp':row['pop']}))
        #qqplot( hardyfilt.P, figname=prefix+"_"+row["pop"]+".png",
               #ptitle="QQ of HWE in "+row["pop"]+" variants" )

    allpvals = concat( allpvals ).reset_index(drop=True)

    allpvals["AA"] = [int(x.split("/")[0]) for x in allpvals["OBS(HOM1/HET/HOM2)"]]
    allpvals["AB"] = [int(x.split("/")[1]) for x in allpvals["OBS(HOM1/HET/HOM2)"]]
    allpvals["BB"] = [int(x.split("/")[2]) for x in allpvals["OBS(HOM1/HET/HOM2)"]]
       
    for region, rdata in allpvals.groupby("Region") :
        print rdata.head()
        hwTernaryPlot( rdata, targetdir, prefix+"_"+region )

    sys.exit(1)
    singletons = allpvals[ (allpvals.nhet == 1) & (allpvals.nhom == 0)]
    doubletons = allpvals[ (allpvals.nhet == 0) & (allpvals.nhom == 1)]
    #singletons.groupby( )

    allpvals.to_csv(os.path.join(targetdir,prefix+"_alldata.tsv"),sep="\t")
    print allpvals.head(10)

    plotHWEafBinned(allpvals, targetdir, prefix )
    plotHWEafRibbon(allpvals, targetdir, prefix )
    #qqplot( allpvals.pval, pfactor=allpvals.grp, 
            #figname=os.path.join(targetdir,prefix+"_qq.png"),
            #ptitle="QQ of HWE for all pops" )
# END plotHardy


################################################################################
# Main
################################################################################
if __name__ == "__main__" :
    os.chdir("..")

    optlist, args = getopt.getopt( sys.argv[1:], "r:ot")
    optlist = dict(optlist)

    dataset = optlist.get("-r",None)
    if dataset is None and not optlist.has_key("-t"):
        print "Error: no dataset provided"
        print " -r <string> with the name of the dataset"
        print " -o flag to overwrite everything"
        print " -t run a test"
        sys.exit(1)

    #hk.changeLogFile( hk.LOGDIR+"/classifyVars_test.log", True )

    print "Using dataset:", dataset
    path = os.path.abspath("./rawdata/")
    if dataset == "test" :
        #vcffile = path+"/test/everything_set1.chr1.snp.clean.vcf.gz"
        rohfile = "./rawdata/test/main/everything_set1.chr1.snp.clean.plink.filt.hom"
        genesfile = "./rawdata/test/main/everything_set1.chr1.snp_genes.tsv"
        vcffile = "./rawdata/test/main/everything_set1.chr1.snp.clean.vcf.gz"

        #vcffile = path+"/test/everything_set1.chr1.snp.vcf.gz"
    elif dataset == "mevariome" :
        rohfile = "./rawdata/mevariome/main/variome.clean.plink.filt.hom"
        genesfile = "./rawdata/mevariome/main/variome.clean_genes.tsv"
        vcffile = "./rawdata/mevariome/main/variome.clean.vcf.gz"

    elif dataset == "mergedaly" :
        rohfile = "./rawdata/mergedaly/main/meceu.clean.plink.filt.hom"
        genesfile = "./rawdata/mergedaly/main/meceu.clean_genes.tsv"
        vcffile = "./rawdata/mergedaly/main/meceu.clean.vcf.gz"

    elif dataset == "daily" :
        rohfile = "./rawdata/daily/main/daily.clean.plink.filt.hom"
        genesfile = "./rawdata/daily/main/daily.clean_genes.tsv"
        vcffile = "./rawdata/daily/main/daily.clean.vcf.gz"

    elif dataset == "test2" :
        rohfile = "./rawdata/test2/main/test2.clean.plink.filt.hom"
        genesfile = "./rawdata/test2/main/test2.clean_genes.tsv"
        vcffile = "./rawdata/test2/main/test2.clean.vcf.gz"

    elif dataset == "merge1kg" :
        rohfile = "./rawdata/merge1kg/main/me1000G.clean.plink.filt.hom"
        genesfile = "./rawdata/merge1kg/main/me1000G.clean_genes.tsv"
        vcffile = "./rawdata/merge1kg/main/me1000G.clean.vcf.gz"

    else :
        print "Error: unknown dataset",dataset
        sys.exit(1)

    filepath, filename = os.path.split( vcffile )                    
    targetdir = os.path.join(filepath,"hwe")
    hk.makeDir( targetdir )
    vtoolsoutfiles = vcftools.run_vcftools_analyze( vcffile )

    print vtoolsoutfiles.head()
    print vtoolsoutfiles.hardyfile.tolist()

    plotHardy( vtoolsoutfiles, genesfile, dataset, targetdir )
# END MAIN    
#elif dataset == "daily" :
        #vcffile = path+"/daily/daily.clean.vcf.gz"
    #elif dataset == "merge1kg" :
        #vcffile = path+"/merge1kg/main/me1000G.clean.vcf.gz"
    #elif dataset == "mergedaly" :
        #vcffile = path+"/mergedaly/main/meceu.clean.vcf.gz"
    #elif dataset == "casanova" :
        #vcffile = path+"/casanova/casanova.snp.recal.clean.vcf.gz"
    #elif dataset == "CEU" :
        #vcffile = path+"/fixed/daily.clean.Europe.vcf.gz"
    #elif dataset == "ME" :
        #vcffile = path+"/fixed/variome.clean.Middle_East.vcf.gz"
    #elif not optlist.has_key("-t") :

