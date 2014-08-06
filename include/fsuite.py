i#!/usr/bin/env python

#-------------------------------------------------------------------#
#                       Import Packages                             #
#-------------------------------------------------------------------#
import os, sys
import getopt
import csv
import re
import subprocess

import housekeeping as hk
from localglobals import *

def runFSuite( bedfile, targetdir=None ) :
    print "Running runFSuite"
    filepath,basename,ext = hk.getBasename(bedfile)
    targetdir = targetdir if targetdir is not None else filepath
    print "Using targetdir:",targetdir
    hk.makeDir(targetdir)

    bimfile = "%s/%s.bim" % (filepath,basename)
    famfile = "%s/%s.fam" % (filepath,basename)

    #1st step : not executed as reference allele frequencies are given.

    #2nd step : create submaps
    #perl ../FSuite/fsuite.pl --create-submaps --map example.bim --n-submaps 5 --coverage > log_example.txt
    step2 = ["fsuite.pl","--create-submaps","--map",bimfile,"--n-submaps","5","--coverage"]

    #from now, we work on submaps in example/submaps
    cleanup = ["rm","-rf","submaps"]

    #3rd step : run FEstim and perform mating-type
    #perl ../FSuite/fsuite.pl --FEstim --bfile example --freq reference.frq --sub-folder 
    #submaps_example --mating-type  >> log_example.txt

    step3 = ["fsuite.pl","--FEstim","--bfile","example","--freq","reference.frq",
             "--sub-folder", "submaps_example","--mating-type"]

    #4th step : calculate FLOD on inbred cases
    #perl ../FSuite/fsuite.pl --FLOD --bfile example --freq reference.frq --sub-folder 
    #submaps_example --inbred #example >> log_example.txt

    step4 = ["fsuite.pl","--FLOD","--bfile","example","--freq","reference.frq","--sub-folder",
            "submaps_example","--inbred","example"]

    #5th step : calculate HFLOD on inbred cases
    #perl ../FSuite/fsuite.pl --HFLOD >> log_example.txt
    step5 = ["fsuite.pl","--HFLOD"]
 
    #6th step: plot HBD segments 
    step6_1 = ["fsuite.pl --plot-HBD --chr 20  >> log_example.txt"
    #perl ../FSuite/fsuite.pl --plot-HBD --fid fam1 --iid ind1  >> log_example.txt
    #perl ../FSuite/fsuite.pl --plot-HBD --circos  >> log_example.txt
# END workflow
#example of others options
#perl ../FSuite/fsuite.pl --estimate-frequencies --bfile example --freq-all
#perl ../FSuite/fsuite.pl --create-submaps --map example.bim --hotspots hg18 --n-submaps 5 --coverage


######################################################################
if __name__ == "__main__":

    optlist, args = getopt.getopt( sys.argv[1:], "bk")
    optlist = dict(optlist)

    # Change running directory
    os.chdir("..")

    #bedfile = "./rawdata/variome1/variome.chimp.regions.filt.samp.plink.mod.filt.uniq.bed"
    #bedfile = "./rawdata/onekg/onekg.chimp.regions.filt.samp.plink.filt.bed"
    bedfile = "./rawdata/test/everything_set1.chr1.snp.clean.plink.filt.bed"

    bedfile = optlist.get("-b",bedfile)
    assert os.path.exists(bedfile)

    filepath, filename = os.path.split(bedfile)
    targetdir = os.path.join(filepath,"pca")
    outliers = runFSuite( bedfile, targetdir )
# END Main


