#!/usr/bin/python

import os, sys
import subprocess
import glob
import pybedtools
#import resource
#resource.setrlimit(resource.RLIMIT_NOFILE, (1000,-1))

from scipy.stats import binom_test

from localglobals import *


def getAnnotCallFiles( dataset, targetdir ) :
    print "Running getAnnotCallFiles"
    #path = os.path.join("rawdata",dataset,"main/calls/annotcalls")
    path = os.path.join( targetdir, "annotcalls" )
    callfiles = hk.filesInDir( path, "tsv" )

    cfiledict = {}
    for cfile in callfiles :
        fpath,basename = os.path.split(cfile)
        samp = basename[:basename.find("_annotcalls")]
        cfiledict[samp] = cfile

    return cfiledict
# END getAnnotCallFiles

def identifySingletons( callfile ):
    #print "Running identifySingletons"
    scalls = read_csv( callfile, sep="\t" ) 
    genotypes = DataFrame( [Series(data=x.split(":"),index=["ref","het","hom"])
                                   for x in scalls["Vcounts"] ] )
    scalls = scalls.join( genotypes )

    #print scalls.head()
    scalls["singles"] = "common"
    scalls.loc[(scalls["het"] == "1") & (scalls["hom"] == "0"),"singles"] = "Singleton"
    scalls.loc[(scalls["het"] == "0") & (scalls["hom"] == "1"),"singles"] = "Doubleton"

    print "Singletons!", len(scalls[scalls.singles == "Singleton"])
    print "Doubletons!"
    print scalls[scalls.singles == "Doubleton"].head()
    singlecounts  = scalls.groupby(["Region","sample","vclass","singles"]).size().reset_index()
    return singlecounts
# END identifySingletons

################################################################################
def countAllSingletons( dataset, force=False ) :
    targetdir = os.path.join("rawdata",dataset,"main/calls")

    singlecountfile = os.path.join(targetdir,dataset+"_commonness.txt")
    #print singlecountfile
    if os.path.exists(singlecountfile) and not force :
        print "Warning: file found"
        singlecounts = read_csv( singlecountfile, sep="\t" )
        return singlecounts

    callfiles = getAnnotCallFiles( dataset, targetdir )

    # Identify singletons and doubletons
    singlecounts = []
    for samp in callfiles :
        print samp, "file:",callfiles[samp]
        sscounts = identifySingletons( callfiles[samp] )
        singlecounts.append( sscounts )
        
    singlecounts = concat( singlecounts )

    singlecounts.to_csv( singlecountfile, sep="\t", index=False )
    return singlecounts
# END countAllSingletons

################################################################################
if __name__ == "__main__" :
    os.chdir("..")

    dataset = "mevariome"

    singles = countAllSingletons( dataset )

################################################################################
