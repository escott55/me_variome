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


def plotHardy( filedf, prefix ) :
    allpvals = []
    for index, row in filedf.iterrows() :
        print row
        print row["pop"], row["hardyfile"]
        hardydata = read_csv(row["hardyfile"], sep="\t" )
        print hardydata.head(10)
        hardyfilt = hardydata[hardydata.P > 0]
        print "Before:",len(hardydata), "After:",len(hardyfilt)
        allpvals.append(DataFrame({'pval':hardyfilt.P, 'grp':row['pop']}))
        qqplot( hardyfilt.P, figname=prefix+"_"+row["pop"]+".png",
               ptitle="QQ of HWE in "+row["pop"]+" variants" )
    allpvals = concat( allpvals ).reset_index(drop=True)
    print allpvals.head(10)
    qqplot( allpvals.pval, pfactor=allpvals.grp, figname=prefix+"_qq.png",
               ptitle="QQ of HWE for all pops" )
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

    elif dataset == "test2" :
        rohfile = "./rawdata/test2/main/test2.clean.plink.filt.hom"
        genesfile = "./rawdata/test2/main/test2.clean_genes.tsv"
        vcffile = "./rawdata/merged/main/test2.clean.vcf.gz"

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
    else :
        print "Error: unknown dataset",dataset
        sys.exit(1)

                        
    vtoolsoutfiles = vcftools.run_vcftools_analyze( vcffile )

    print vtoolsoutfiles.head()
    print vtoolsoutfiles.hardyfile.tolist()

# END MAIN
