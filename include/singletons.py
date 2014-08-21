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

    path = os.path.abspath("./rawdata/")
    if dataset == "test2" :
        #vcffile = path+"/test/everything_set1.chr1.snp.clean.vcf.gz"
        singletonfile = "./rawdata/test2/main/vstats/test2.clean.singletons"
        genesfile = "./rawdata/test2/main/test2_genes.tsv"
        vcffile = "./rawdata/test2/main/test2.clean.vcf.gz"

        #vcffile = path+"/test/everything_set1.chr1.snp.vcf.gz"
    elif dataset == "mevariome" :
        singletonfile = "./rawdata/mevariome/main/vstats/variome.clean.singletons"
        genesfile = "./rawdata/mevariome/main/variome.clean_genes.tsv"
        vcffile = "./rawdata/mevariome/main/variome.clean.vcf.gz"

    elif dataset == "mergedaly" :
        singletonfile = "./rawdata/mergedaly/main/vstats/mecue.clean.singletons"
        genesfile = "./rawdata/mergedaly/main/meceu.clean_genes.tsv"
        vcffile = "./rawdata/mergedaly/main/meceu.clean.vcf.gz"

    elif dataset == "daily" :
        singletonfile = "./rawdata/daily/main/vstats/daily.clean.singletons"
        genesfile = "./rawdata/daily/main/daily.clean_genes.tsv"
        vcffile = "./rawdata/daily/main/daily.clean.vcf.gz"

    elif dataset == "merge1kg" :
        singletonfile = "./rawdata/merge1kg/main/vstats/me1000G.clean.singletons"
        genesfile = "./rawdata/merge1kg/main/me1000G.clean_genes.tsv"
        vcffile = "./rawdata/merge1kg/main/me1000G.clean.vcf.gz"

    else :
        print "Error: unknown dataset",dataset
        sys.exit(1)

    print "something"
    singles = read_csv( singletonfile, sep="\t" )
    
    vinfo = read_csv( genesfile, sep="\t" )
    vinfo = (vinfo[["chrom","pos","vclass"]]
                       .groupby(["chrom","pos"])
                       .max().reset_index())

    singletons = merge( singles, vinfo, left_on=["CHROM","POS"],
              right_on=["chrom","pos"] )
    
    print singletons.head()
# END MAIN
