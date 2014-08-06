#!/usr/bin/env python

import os,sys
import csv
import re
import getopt
from pandas import *

sys.path.append('include')
from housekeeping import *
import popgencommands as pop
import runfestim
import changeref
import patientInfo
import admixture
import vcftoolsanalyze

######################################################################
# preprocessFiles
######################################################################
def preprocessFiles( vcffile, keepfile=None ) :
    write2log( "Running "+whoami() )
    printCurrentTime()
    bedfile = "./resources/regionbeds/percentcoverage.txt"
    recodefiles = []

    originalvcf = vcffile
    if originalvcf[-3:] != ".gz" :
        vcffile = compress( originalvcf )

    assert(os.path.exists(vcffile))
    keeppats = patientInfo.getPats(vcffile, keepfile)

    print "Keeping:",len(keeppats),"Pats"
    pseqfile, outliers = pop.pseq_istats( vcffile )

    # remove indels, non-"PASS", Coverage, and multi-allelic sites
    # Depth and Quality
    patfiltvcf = pop.run_vcftools(vcffile, keeppats )
    pseqfile, outliers = pop.pseq_istats( patfiltvcf )

    chimpvcf = changeref.forceRef( patfiltvcf, "chimp" ) 
    print chimpvcf
    assert os.path.exists( chimpvcf ), "Uhh.... shit" 
    pseqfile, outliers = pop.pseq_istats( chimpvcf )

    # Convert vcf to plink
    pedfile = pop.run_vcftools_plink( chimpvcf )
    print "Pedfile :",pedfile
    bedfile_filt = pop.run_plink_filter( pedfile )

    # Assess relatedness
    patientlist, uniqfile, kinfile = pop.run_king(bedfile_filt)
    make_kinshipplots( kinfile )
    print "Remaining Patients:",len(patientlist)

    # keep only samples in uniqfile
    bedfile_uniq = pop.run_plink_remove(bedfile_filt, uniqfile)
    
    # merge vcf and bed file?
    trimmedfile = pop.intersectBED( chimpvcf, bedfile )
    pseqfile, outliers = pop.pseq_istats( trimmedfile )

    # get pats again
    keeppats = patientInfo.getPats(uniqfile)
    uniqvcf = pop.run_vcftools(trimmedfile, keeppats)

    # Make plink, then filter by MAF and missingness, recode12
    modped = runfestim.calculateIBC( uniqvcf )

    recodefiles.append(uniqvcf)
    return recodefiles
# END preprocessFiles

######################################################################
# chimpifyFiles
######################################################################
def chimpifyFiles( vcffile, dataset="test", keepfile=None, rerun=False, allpseq=False, stayhuman=False ) :
    write2log( "Running "+whoami() )
    printCurrentTime()
    bedfile = os.path.abspath("./resources/regionbeds/collab_percentcoverage.txt")

    #annotfile = "./patientannot/annotation/patientannotation_laser.txt"
    #patientdata = pop.parsePatientPed( annotfile )
    #print len(patientdata)

    sampleannot = pop.parseAnnotation()
    sys.exit(1)
    print sampleannot[["Family.ID","Individual.ID","Paternal.ID","Maternal.ID","Gender","Phenotype"]]
    sampleannot.index = sampleannot["Individual.ID"]
    assert(os.path.exists(vcffile))

    originalvcf = vcffile
    if originalvcf[-3:] != ".gz" : vcffile = compress( originalvcf )
    if allpseq : pseqfile, outliers = pop.pseq_istats( vcffile )
    #print "Pseq:",pseqfile, "Outliers:",len(outliers)

    #vcfnames = patientInfo.fixPatNames( vcffile )
    #print "Finished fixPatNames:", vcfnames

    if stayhuman : chimpvcf = vcffile
    else : 
        print "Running Force chimp reference"
        chimpvcf = changeref.forceRef( vcffile, "chimp" ) 
        statsfile = changeref.referenceStats( chimpvcf )
        changeref.make_refplots( statsfile )

    #assert os.path.exists( chimpvcf )
    if allpseq : pseqfile, outliers = pop.pseq_istats( chimpvcf )
    #print "Pseq:",pseqfile, "Outliers:",len(outliers)

    print "Trim to a set of regions"
    #trimmedvcf = pop.intersectBED( chimpvcf, bedfile )

    #if dataset == "hgdp" : trimmedvcf = chimpvcf
    #else :
        #trimmedvcf = pop.intersectBED( chimpvcf, bedfile )
        #if allpseq : pseqfile, outliers = pop.pseq_istats( trimmedvcf )
        #print "Pseq:",pseqfile, "Outliers:",len(outliers)

    #vcftoolsanalyze.run_vcftools_analyze( trimmedvcf )

    # remove indels, non-"PASS", Coverage, and multi-allelic sites
    # Filter Depth and Quality
    filtvcf = pop.run_vcftools( trimmedvcf, force=rerun )
    pseqfile, outliers = pop.pseq_istats( filtvcf )

    # identify samples to keep
    #print sampleannot.Source.unique().tolist()
    #print sum(sampleannot.Source == "1000g")
    outliers = [x for x in outliers if x not in sampleannot["Individual.ID"][sampleannot.Source == "1000g"].tolist()]
    mixedsamples = sampleannot[sampleannot.Ethnicity.isin(["ASW","Mixed"])]["Individual.ID"].tolist()
    #print len(outliers)
    #print len(mixedsamples)

    #if vcffile.find("onekg") >= 0 : print "skipping outliers"; outliers=[]
    allrm = sampleannot.loc[outliers+mixedsamples][["Individual.ID","Source","continent","Ethnicity","ethnicity"]]
    allrm.to_csv( "results/"+dataset+"removedoutliers.txt", sep="\t" )
    
    keeppats = patientInfo.getPats(filtvcf, keeplist=keepfile, removelist=outliers+mixedsamples)
    print "keeppats:",keeppats[:10]
    print "Keeping:",len(keeppats),"Pats"

    # Filter vcf samples
    filtvcfsamp = pop.run_vcftools_filtsamples(filtvcf, keeppats, force=True )
    pseqfile, outliers = pop.pseq_istats( filtvcfsamp )
    print "Pseq:",pseqfile, "Outliers:",len(outliers)

    vcftoolsanalyze.run_vcftools_analyze( filtvcfsamp )

    # Convert vcf to plink
    tpedfile = pop.run_vcftools_plink( filtvcfsamp, force=True )
    print "Pedfile :",tpedfile

    # Modify Ped File annotation
    pop.modify_tped( tpedfile, sampleannot, force=True )

    # LD vars, MAF >=.05, Missingness >=.05
    bedfile_filt = pop.run_plink_filter( tpedfile, force=True )

    roh_calls = pop.run_plink_homozyg( bedfile_filt, force=True )

    # Assess relatedness
    patientlist, uniqfile, kinfile = pop.run_king(bedfile_filt)

    pop.make_kinshipplots( kinfile )
    print "Remaining Patients:",len(patientlist)

    # keep only samples in uniqfile
    bedfile_uniq = pop.run_plink_remove(bedfile_filt, uniqfile, force=True)
    
    # get pats again
    print "#"*30,"- run_vcftools_filtsamples, File:",uniqfile
    keeppats = patientInfo.getPats(uniqfile)
    print "keeppats:",keeppats[:10]
    uniqvcf = pop.run_vcftools_filtsamples( filtvcfsamp, keeppats, force=True )
    print "#"*30
    print keeppats
    print uniqvcf

    # Make plink, then filter by MAF and missingness, recode12
    #print "#"*30,"- calculateIBS, File:",uniqvcf
    #modped = runfestim.calculateIBC( uniqvcf )

    # Admixture
    bestK = admixture.admixtureAnalysis( bedfile_filt, force=True )

    return uniqvcf
# END chimpifyFiles
    #filepath, basename, suffix = getBasename( bedfile_uniq )
    #admixdir = filepath+"/admixture/"
    #print "Admixdir:",admixdir
    #makeDir(admixdir)
    #print "Bedfile uniq:",bedfile_uniq
    #runCommand("mv %s/%s* %s" % (filepath,basename,admixdir))
    #bedfile_uniq = admixdir+basename+suffix
    #print keeppats
    #for pat in keeppats1 :
        #if pat not in keeppats :
            #print pat

######################################################################
# processVCF
######################################################################
def processVCF( ):
   print "Nothing" 
# END processVCF

######################################################################
# Main
######################################################################
if __name__ == "__main__" :

    optlist, args = getopt.getopt( sys.argv[1:], "r:oh")
    optlist = dict(optlist)

    dataset = optlist.get("-r",None)
    if dataset is None : 
        print "Error: no dataset provided"
        print " -r <string> with the name of the dataset"
        print " -o flag to overwrite everything"
        sys.exit(1)

    changeUlimit(4046)

    print "Using dataset:", dataset
    path = os.path.abspath("./rawdata/")
    #vcffile = path+"ciliopathies/ciliopathies.unfilt.vcf.gz"
    if dataset == "ciliopathies" :
        vcffile = path+"/ciliopathies/ciliopathies.vcf.gz"
    elif dataset == "onekg" :
        vcffile = path+"/onekg/onekg.vcf.gz"
    elif dataset == "daily" :
        vcffile = path+"/daily/daily.vcf.gz"
    elif dataset == "eichler" :
        vcffile = path+"/eichler/eichler.vcf.gz"
    elif dataset == "turks" :
        vcffile = path+"/turks/turks.vcf.gz"
    elif dataset == "test" :
        vcffile = path+"/test/everything_set1.chr1.snp.vcf.gz"
    elif dataset == "fowzan" :
        vcffile = path+"/fowzan/fowzan.snp.recal.vcf.gz"
    elif dataset == "casanova" :
        vcffile = path+"/casanova/casanova.snp.recal.vcf.gz"
    elif dataset == "variome" :
        vcffile = path+"/variome/variome_snps.vcf.gz"
    elif dataset == "variome1" :
        vcffile = path+"/variome1/variome.vcf.gz"
    elif dataset == "hgdp" :
        vcffile = path+"/hgdp/HGDP_938.vcf.gz"
    elif dataset == "merged" :
        vcffile = path+"/merged/merged.vcf.gz"

    bedfile = "./resources/regionbeds/collab_percentcoverage.txt"

    filepath, basename, suffix = getBasename( vcffile )
    changeLogFile( LOGDIR+"/process"+basename+".log" )

    #if vcffile[-3:] != ".gz" :
        #cfile = compress( vcffile )

    #preprocessFiles( vcffile )

    stayhuman = optlist.has_key("-h")
    print "Should we stay human??", stayhuman

    chimpifyFiles( vcffile, dataset, rerun=optlist.has_key("-o"), stayhuman=stayhuman )

    targetdir = path+"merged"

# END Main

    #path = "./rawdata/"
    #bedfile = "./regionbeds/percentcoverage.bed"
    #targetfile = path+"merged.vcf"
    #targetdir = path+"turks"

    #allfiles = filesInDir( path, "vcf.gz" )
    #cfile = compress( path+"daily.vcf" )

    #allfiles = [path+"ciliopathies.unfilt.vcf.gz", \
         #path+"daily.vcf.gz", \
         #path+"eichler.vcf.gz" ]
    #allfiles = [path+"ciliopathies.unfilt.vcf.gz", \
         #path+"ALL.genome.phase1.vcf.gz" ]

    #allfiles = [ path+"eichler.vcf.gz",
                 #path+"daily.vcf.gz" ]


