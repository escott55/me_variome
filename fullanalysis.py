#!/bin/env python

# Starting with a VCF file 
# VCFtools - Convert to plink file with some filtering
# Plink - Filter variants and convert to bed file
# King - Calculate Kinship and unrelated individual list
# Admixture - Determine optimal number of clusters
# python - figure out if clusters makes sense

#-------------------------------------------------------------------#
#                       Import Packages                             #
#-------------------------------------------------------------------#
import os, sys
import commands
from time import clock,time
import glob
import getopt
import csv
import re

from housekeeping import *

#-------------------------------------------------------------------#
#                            GLOBALS                                #
#-------------------------------------------------------------------#
DEBUGFLAG = True
forceFlag = False
#DBSNPLIST = "/local/gleeson/analysis/inbreed/dbsnp_137_idlist.hg19.txt"
DBSNPLIST = "/home/escott/workspace/inbreed/dbsnp_137_idlist.hg19.txt"

#-------------------------------------------------------------------#
#                            Commands                               #
#-------------------------------------------------------------------#

######################################################################
# runAdmixtureClustering
######################################################################
def runAdmixtureClustering(bedfile_uniq,annotfile,currK,forceFlag=forceFlag) :
    write2log( " - Running:"+whoami(), True )
    print "Target File:",bedfile_uniq, "Annot:",annotfile,"CurrK:",currK
    targetdir, filename, suffix = getBasename(bedfile_uniq)

    Qfile = targetdir + filename+"."+str(currK)+".Q"
    groupsfile = "%s/%s.%d.groups.txt" % (targetdir,filename,currK)
    logfile = "%s/%s_admixclustering.%d.log" % (LOGDIR,filename,currK)
    print groupsfile
    if os.path.exists( groupsfile ) and forceFlag is False:
        write2log(" - File"+groupsfile+" already exists. Skipping command!",False)
        groups = parseGroupsFile( groupsfile )
        return groups

    command = "./admix.R %s %s %d" % (bedfile_uniq, annotfile, currK)
    runout = runCommand( command )
    OUT = open(logfile,"wb")
    OUT.write(runout)
    OUT.close()
    groups = parseGroupsFile( groupsfile )
    return groups
# runAdmixtureClustering

######################################################################
# plotInbreedingCoefficient
######################################################################
def plotInbreedingCoefficient( bedfile_uniq, currK ) :
    write2log( " - Running:"+whoami(), True )
    print "Target File:",bedfile_uniq
    targetdir, filename, suffix = getBasename(bedfile_uniq)

    filepattern = "%s.%s" % (filename, currK)
    command = "./inbreedCoeff.R %s %s" % ( targetdir, filepattern )
    print "Command",command
    output = runCommand( command )

    return groups
# plotInbreedingCoefficient

#-------------------------------------------------------------------#
#                            Workflows                              #
#-------------------------------------------------------------------#

######################################################################
# preprocess_vcf
######################################################################
def preprocess_vcf( targetdir, annotfile, fprefix=None, outdir=None ) :
    write2log( " - Running:"+whoami(), True )

    if outdir is None :
        outdir = targetdir
    print "Targetdir:", targetdir
    print "Annotfile:", annotfile
    print "Prefix:", fprefix
    print "Outdir:", outdir

    patientdata = parsePatientPed( annotfile )

    filetype = ".vcf.gz"
    allpedfiles = []
    for vcffile in filesInDir(targetdir, filetype, prefix=fprefix ) :
        print "VCF:",vcffile
        pedfile = run_vcftools( vcffile, outdir )
        print "PED:",pedfile
        modped = modify_ped( pedfile, patientdata )
        allpedfiles.append( modped )

    print allpedfiles
    #sys.exit(1)

    filetype = "plink.mod.ped"
    #filetype = "plink.ped"
    bestK = -1
    #for pedfile in filesInDir(outdir, filetype, prefix=fprefix ) :
    allpedfiles = filesInDir(outdir, filetype, prefix=fprefix )
    maxK = 10
    if len(allpedfiles) > 1 :
        print "Bizzare Number of pedfiles :",len(allpedfiles)
        sys.exit(1)
    elif len(allpedfiles) == 0  :
        print "Error: No pedfiles found!"
        sys.exit(1)

    pedfile = allpedfiles[0]
    bedfile = run_plink_convert(pedfile)
    bedfile_filt = run_plink_filter(bedfile)
    patientlist, uniqfile, kinfile = run_king(bedfile_filt)
    print "Kin file:",
    make_kinshipplots( kinfile )
    print "Remaining Patients:",len(patientlist)
    bedfile_uniq = run_plink_remove(bedfile_filt, uniqfile)
    fixedpatientlist, fixeduniqfile, fixedkinfile = run_king(bedfile_uniq)
    make_kinshipplots( fixedkinfile )
    groups = determinegroups( uniqfile, patientdata )
    print "bed",bedfile,"uniq",uniqfile,"filt",bedfile_filt
    if len(groups) < maxK :
        maxK = len(groups)

    bestK = findBestK( bedfile_uniq, maxK )
    run_admixture( bedfile_uniq, bestK, True, forceFlag )
    run_admixture( bedfile_uniq, bestK+1, True, forceFlag )
    run_admixture( bedfile_uniq, bestK+2, True, forceFlag )
    run_plink_het( bedfile_uniq )
    print "BestK:",bestK
    if bestK < 0 :
        print "Error: No best K found"
        sys.exit(1)

    for K in range(2,maxK) :
        patientgroups = annotatePatients( patientlist, bedfile_uniq, patientdata, K )

    return bedfile_uniq, bestK, groups
# End preprocess_vcf
    #bedfile = run_plink_filter(pedfile)
    #patientlist, uniqfile = run_king(bedfile)
    #bedfile_filt = run_plink_remove(bedfile, uniqfile)
    #groups = determinegroups( uniqfile, patientdata )
    #bestK = findBestK( bedfile_filt, len(groups)+1 )

######################################################################
# inbreedingCoefficient
######################################################################
def inbreedingCoefficient( bedfile, annotfile, bestK, groups ) :
    targetdir, filename, suffix = getBasename(bedfile)
    print "bedfile:",targetdir,filename,suffix
    patientdata = parsePatientPed( annotfile )

    # GRAPH ADMIXTURE
    admixgroups = {}
    #maxK = 24
    #if len(groups) < maxK :
        #maxK = len(groups)
    #for K in range(2,maxK) :
        #admixgroups = runAdmixtureClustering(bedfile,annotfile,K)

    #forceReRun = True
    admixgroups = runAdmixtureClustering(bedfile,annotfile,bestK, forceFlag)


    #sys.exit(1)
    #for p in patientdata :
        #print p
    for gnum in sorted(admixgroups) :
        print "Group #",gnum,":",len(admixgroups[gnum])
        keepfile = "%s/%s.%d.grp%s.txt" % (LOGDIR,filename,bestK,gnum)
        print "Keepfile:",keepfile
        OUT = open( keepfile, "wb" )
        for patient in admixgroups[gnum] :
            if patientdata.has_key(patient) :
                OUT.write( "\t".join(patientdata[patient][:2]) +"\n")
            else :
                OUT.write( patient+"\t"+patient +"\n" )
        OUT.close()
        hetfile = run_plink_het( bedfile, bestK, gnum, keepfile, forceFlag )
    plotInbreedingCoefficient( bedfile, bestK )
# inbreedingCoefficient

######################################################################
# Main
######################################################################
if __name__ == "__main__":

    optlist, args = getopt.getopt( sys.argv[1:], "iegfacmb")

    if len(optlist) == 0:
        print "Error: no options given: -t -w"
        print " -i run inhouse"
        print " -e run 1000 exomes"
        print " -g run 1000 genomes"
        print " -c run 1000 genomes chrom1"
        print " -a run all internal and external"
        print " -m run bigmerged"
        print " -b run bigmerged2"
        print " -f boolean Force overwrite"
        sys.exit(1)

    # Set boolean options
    for opttype in optlist :
        if opttype[0] == "-f" :
            toggleForceFlag()

    #print runCommand( "ulimit -n unlimited" )
    for opttype in optlist :
        if opttype[0] == "-i" :
            basedir = "inhouse/"
            changeLogDir( basedir )
            targetdir = basedir+"cohort/"
            #prefix = "plateII_"
            #prefix = None
            #prefix = "merged_snps.trim"
            targetdir = basedir+"allplates/"
            outdir = basedir+"cohort/"
            prefix = "merged.percentcoverage"
            patientannotation = basedir+"inhouse_samples_annotation.ped"
            #preprocess_vcf(targetdir, patientannotation, prefix )
            bedfile, bestK, groups = preprocess_vcf(targetdir, patientannotation, prefix, outdir )
            for currK in range(3,7) :
                inbreedingCoefficient( bedfile,patientannotation,currK,groups )
        elif opttype[0] == "-e" :
            basedir = "onekg/"
            changeLogDir( basedir )
            targetdir = basedir+"data"
            outdir = basedir+"cohort/"
            patientannotation = basedir+"phase1_samples_integrated_20101123.ped"
            #preprocess_vcf(targetdir, patientannotation, "ALL" )
            prefix = "ALL.genome"
            preprocess_vcf(targetdir, patientannotation, prefix, outdir )
        elif opttype[0] == "-c" :
            basedir = "onekg/"
            changeLogDir( basedir )
            targetdir = basedir+"data"
            outdir = basedir+"cohort/"
            patientannotation = basedir+"phase1_samples_integrated_20101123.ped"
            prefix = "chr1.phase1"
            preprocess_vcf(targetdir, patientannotation, prefix, outdir )
        elif opttype[0] == "-a" :
            basedir = "merged/"
            changeLogDir( basedir )
            #targetdir = basedir+"data/"
            #targetdir = "/local/gleeson/analysis/inbreed/rawdata/"
            targetdir = "/home/escott/workspace/inbreed/rawdata/"
            outdir = basedir+"unaffonly/"
            patientannotation = basedir+"combined.ped"
            prefix = "merged.1kg.region"
            bedfile, bestK, groups = preprocess_vcf(targetdir, patientannotation, prefix, outdir )
            print "Finished VCF:",bedfile, bestK
            for currK in range(bestK,bestK+3) :
                inbreedingCoefficient( bedfile,patientannotation,currK,groups )

        elif opttype[0] == "-m" :
            basedir = "bigmerged/"
            changeLogDir( basedir )
            targetdir = basedir+"data/"
            outdir = basedir+"cohort/"
            patientannotation = basedir+"combined.ped"
            prefix = "merged.eichler"
            bedfile, bestK, groups = preprocess_vcf(targetdir, patientannotation, prefix, outdir )
            print "Finished VCF:",bedfile, bestK
            for currK in range(bestK,bestK+3) :
                inbreedingCoefficient( bedfile,patientannotation,currK,groups )

        elif opttype[0] == "-b" :
            basedir = "bigmerge2/"
            changeLogDir( basedir )
            targetdir = basedir+"data/"
            outdir = basedir+"cohort/"
            patientannotation = basedir+"combined.ped"
            prefix = "merged.daly"
            bedfile, bestK, groups = preprocess_vcf(targetdir, patientannotation, prefix, outdir )
            print "Finished VCF:",bedfile, bestK
            for currK in range(bestK,bestK+3) :
                inbreedingCoefficient( bedfile,patientannotation,currK,groups )


        elif opttype[0] == "-g" :
            basedir = "onekgenomes/"
            changeLogDir( basedir )
            targetdir = basedir+"cohort/"
            patientannotation = basedir+"phase1_samples_integrated_20101123.ped"
            preprocess_vcf(targetdir, patientannotation, "plateIV_2" )
        else :
            print "Error: Unknown operation"

# End Main
######################################################################

