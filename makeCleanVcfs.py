#!/usr/bin/env python

#-------------------------------------------------------------------#
#                       Import Packages                             #
#-------------------------------------------------------------------#
import os, sys
import getopt
import csv
import re
import gzip
import subprocess
from pandas import *
#from random import randint

sys.path.append('include')
from localglobals import *
import popgencommands as pop
import king
import patientInfo
import classifyVars
import roh_analysis
import vcftoolsanalyze as vcftools
import roh_more
import admixture
import runfestim
import changeref
import math
import eigensoft

################################################################################
# imputeVars
################################################################################
def imputeVars( vcffile, pedfile, prefix="test" ) :
    removefile = "toremove.txt"
    cpats = patientInfo.currentFilePats( vcffile )
    removeset = [x.rstrip() for x in open(removefile,"r").readlines()]
    print removeset
    toremove = [x for x in removeset if x in cpats]
    newremovefile = "tmpremove.txt"
    if len(toremove) > 0 :
        OUT = open(newremovefile,"wb")
        OUT.write("\n".join(toremove))
        OUT.close()

               #" ref=/home/escott/Packages/Beagle/onekg.vcf.gz" +
    command = ["java","-Xmx4000m","-jar","/home/escott/Packages/Beagle/b4.r1274.jar",
               "gtgl=%s" % vcffile,
               "ped=%s" % pedfile,
               "out=%s" % prefix,
               "nthreads=4"]
    if len(toremove) > 0 : command.append(" excludesamples=%s" % newremovefile)
    print " ".join(command)
    stdout = subprocess.check_output(command, stderr=subprocess.STDOUT)
# END imputeVars

def vcftools_to_plink( vcffile, keeplist, force=False):
    targetdir, filename, suffix = hk.getBasename(vcffile)
    plinkfile = os.path.join(targetdir,filename+".plink")
        
    finalped = plinkfile+".tped"
    if os.path.exists( finalped )  and not force :
        print "Warning file already exists:",finalped,"skipping!"
        return finalped

    command = ["vcftools","--gzvcf",vcffile,
                "--keep",keeplist,
                "--plink-tped","--out",plinkfile]
        #" --maf .002" \
        #" --remove-filtered-geno-all " \
        #" --FILTER-summary"\
        #" --plink-tped --out %s" \
        #" --snps %s" \
        #" --min-meanDP 4 --minQ 20" \
        #" --thin 10" \

    print " ".join(command)
    runout = subprocess.check_output( command )
    
    return finalped
# END vcftools_to_plink

################################################################################
# mad
################################################################################
def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))

def isOutlier( x, mean, std ):
    return abs((x-mean)/std) > 5

################################################################################
# popOutliers
################################################################################
def popOutliers( vcffile, outprefix=None ) : #, sampleannot
    print "Running popOutliers"
    originaldir, vfilename, vsuffix = hk.getBasename(vcffile)
    vtargetdir = outprefix if outprefix is not None else originaldir
    #print vtargetdir, outprefix, vfilename
    pseqfile = os.path.join(vtargetdir, "%s_istats.txt"%vfilename)
    pseqerr = os.path.join(vtargetdir, "%s_istats.err"%vfilename)
    outlierfile = os.path.join(vtargetdir, "%s.qcoutliers"%vfilename)
    if not os.path.exists( pseqfile ) or hk.fileIsEmpty( pseqfile ):
        with open(pseqfile,"wb") as out, open(pseqerr,"wb") as err:
            print " ".join(["pseq",vcffile,"i-stats"])
            subprocess.call(["pseq",vcffile,"i-stats"],stdout=out,stderr=err)
            out.close(); err.close()

    #print pseqerr, pseqfile
    istats = read_csv( pseqfile, sep="\t" )
    alloutliers = []
    for col in istats.columns : 
        if col == "ID" : continue
        if col in ["DP"] : continue
        mean = istats[col].mean()
        std = istats[col].std()
        outliers = istats[["ID",col]][[isOutlier(x,mean,std) for x in istats[col].tolist()]]
        for out,val in outliers[["ID",col]].values :
            alloutliers.append([out,val,col,mean,std])

    if len(alloutliers) > 0 :
        alloutliers = DataFrame( alloutliers, columns=["ID","Value","Col","Mean","STD"])
    else :
        alloutliers = DataFrame( columns=["ID","Value","Col","Mean","STD"] )

    #print "Shape:",alloutliers.shape
    print "Going to remove",len(alloutliers.ID.unique()),"samples"
    print "Writing outlier file:",outlierfile
    alloutliers.to_csv(outlierfile,sep="\t")
    #print alloutliers[alloutliers.ID == "JL0688_ACAGTG_BD2GJBACXX_L006_001"]
    return alloutliers.ID.unique().tolist()
# END popOutliers


################################################################################
################################################################################
def seriousCleanVcf( vcffile, keepfile, outprefix, rerun=False ):
    print "Running seriousCleanVcf"

    #outprefix = "%s/%s"%(targetdir,basename)
    scleanfile = outprefix+'.sclean.vcf'
    if os.path.exists( scleanfile+".gz" ) and not rerun :
        return scleanfile+".gz"

    command = ["vcftools","--gzvcf", vcffile,
               "--remove-filtered-all",
               "--remove-indels",
               "--keep",keepfile,
               "--maf",".01",
               "--min-alleles","2",
               "--max-alleles","2",
               "--recode","--out",outprefix]

    print " ".join(command)
    cleanvcf = outprefix+".recode.vcf"
    print "Cleanvcf should exist:",cleanvcf
    if not os.path.exists( cleanvcf ) or rerun:
        out = subprocess.check_output( command )

    assert os.path.exists( cleanvcf ), out

    subprocess.call( ['mv',cleanvcf,scleanfile] )
    subprocess.call( ["bgzip","-f",scleanfile] )
    subprocess.call( ["tabix","-p","vcf",scleanfile+".gz"] )

    return scleanfile+".gz"
    #print vcffile
# END seriousCleanVcf

################################################################################
# Calculate Relationships
# Select unique
# Make pedfile 
# Identify Outliers
################################################################################
def identifySamplesToKeep( vcffile,toremove="./toremove.txt", force=False) :
    print "Running identifySamplesToKeep"
    filepath, filename, suffix = hk.getBasename(vcffile)
    targetdir = os.path.join(filepath,"qc")
    hk.makeDir( targetdir )
    outprefix = "%s/%s" % (targetdir, filename)

    finalkeepfile = outprefix+".finalpats"
    if os.path.exists(finalkeepfile) and not force:
        print "Finalpats file exists! skipping....", finalkeepfile
        return finalkeepfile

    autoremove = []
    if os.path.exists(toremove):
        autoremove = ([x.strip() for x in open(toremove).readlines()
                       if len(x.strip()) > 0])

    qcoutliers = popOutliers( vcffile, targetdir )
    print "Finished popOutliers"

    filepats = patientInfo.currentFilePats( vcffile )
    beforepatsfile = outprefix+".beforepats"
    fp_annot = addSampleAnnotation( DataFrame({'filepats':filepats}), mergecol="filepats" )

    print fp_annot.head()
    print "Writing filter file", beforepatsfile
    fp_annot[["filepats","Continent","GeographicRegions2"]].to_csv( 
        beforepatsfile, sep="\t", index=False )

    firstkeepfile = outprefix+".keep"
    OUT = open(firstkeepfile,"wb")
    for samp in filepats :
        if samp in autoremove : print "autoremove",samp; continue
        if samp in qcoutliers : print "qcoutlier",samp; continue # skip outliers
        OUT.write(samp+"\n")

    OUT.close()
    
    #cleanvcffile = pop.run_vcftools_clean( vcffile, firstkeepfile, outprefix, force=force )
    cleanvcffile = seriousCleanVcf( vcffile, firstkeepfile, outprefix, force )

    # unrelated samples
    #unrelatedfile = vcftools.identifyUnrelatedSamples( cleanvcffile, targetdir, force=force )
    unrelatedfile, famrelfile, interrelfile = king.king_workflow( cleanvcffile, force )

    unrelatedkeepfile = outprefix+".unrelpats"
    OUT = open( unrelatedkeepfile, "w")
    for row in csv.reader(open(unrelatedfile), delimiter="\t") :
        if len(row) <= 1 : continue
        OUT.write(row[1]+"\n") 
    OUT.close()

    # to Tped, filter out indels, non-bi-allelic, maf < .002
    #tpedfile = pop.run_vcftools_plink( cleanvcffile, force=force, 
                                      #keeplist=unrelatedfile )
    tpedfile = vcftools_to_plink( cleanvcffile, unrelatedkeepfile, force=force)
    # convert to bed file
    print "$"*50
    print "got to here"
    bedfile = pop.run_plink_convert(tpedfile, force=force)

    # PCA outlier analysis
    outliers = eigensoft.pcaOutlierAnalysis( bedfile, force=force )

    finalpats = []
    print "Writing final keep file:",finalkeepfile
    OUT = open(finalkeepfile,"wb")
    for line in [x.strip() for x in open(unrelatedfile).readlines()] :
        fam, samp = line.split("\t")
        if samp not in outliers : 
            finalpats.append(samp)
            OUT.write(samp+"\n")
    OUT.close()

    afterfilterfile = outprefix+".afterpats"
    print "finalpats length:",len(finalpats)
    fp_annot = addSampleAnnotation( DataFrame({'filepats':finalpats}), mergecol="filepats" )
    print "After annot:", len(fp_annot)
    #print fp_annot.head()
    print "Writing filter file", afterfilterfile
    fp_annot[["filepats","Continent","Continent2","GeographicRegions2","Origin",
              "ethnicity","Country","country"]].to_csv( 
        afterfilterfile, sep="\t", index=False )

    sys.exit(1)
    return finalkeepfile
# END identifySamplesToKeep

################################################################################
# makeCleanVCF
################################################################################
def makeCleanVCF( vcffile, keepfile=None, toforce=False ) :
    print "Running makeCleanVCF",vcffile
    targetdir, filename, suffix = hk.getBasename(vcffile)
    outdir = os.path.join(targetdir,"main")
    hk.makeDir(outdir)
    outprefix = "%s/%s" % (outdir, filename)
    #popsprefix = "%s/vstats/pops/%s" % (outdir, filename)

    filepats = patientInfo.currentFilePats( vcffile )
    sampleannot = sampleAnnotation( filepats )
    regionlist = sorted(sampleannot.Continent.unique().tolist())

    # identify samples to keep 
    # remove related
    # remove outliers
    keepfile = identifySamplesToKeep( vcffile, force=True ) 
    print "Keep file is:",keepfile

    #if keepfile is not None :
    print "Fixing vcffile:",vcffile
    cleanvcffile = pop.run_vcftools_clean( vcffile, keepfile, outprefix, force=toforce )
    print "New vcffile:", cleanvcffile

    pseqfile, outliers = pop.pseq_istats( cleanvcffile )
    print "Pseq:",pseqfile, "Outliers:",len(outliers)

    # to Tped, filter out indels, non-bi-allelic, maf < .002
    tpedfile = pop.run_vcftools_plink( cleanvcffile, force=toforce )
    # LD filtering
    bedfile_filt = pop.run_plink_filter( tpedfile, force=toforce )
 
    # ROH Analysis
    figprefix = "./results/figures/roh/%s" % (filename)
    roh_calls = pop.run_plink_homozyg( bedfile_filt, force=toforce )
    roh_analysis.makeRohBoxPlots( roh_calls, sampleannot, figprefix )
    roh_analysis.makeRohCumPlot( roh_calls, sampleannot, figprefix )

    # Admixture
    print "Run admixture on your own!"
    #admixture.admixtureAnalysis( bedfile_filt, sampleannot )

    # ROH exome coverage
    roh_more.exomeCoverage( roh_calls, sampleannot, outprefix=filename )

    # IBC 
    runfestim.calculateIBC( cleanvcffile, rerun=toforce )

    # VCF tools analysis
    vcftools.run_vcftools_analyze( cleanvcffile, keepfile )

    # classify vars
    classifyVars.runClassifyWorkflow( cleanvcffile, sampleannot, regionlist, "./results/classes")
  
    # Chimpify
    chimpvcf = changeref.forceRef( cleanvcffile, "chimp" )
    statsfile = changeref.referenceStats( chimpvcf )
    changeref.make_refplots( statsfile )

    # classify again
    classifyVars.runClassifyWorkflow( chimpvcf, sampleannot, regionlist, "./results/classes")

    return 
# END makeCleanVCF
    #smallvcf = classifyVars.makeSmallVCF( cleanvcffile, force=False )
    #allvcffiles = classifyVars.annotateVCF( smallvcf, force=False )

    #print allvcffiles
    #geneannotfile = classifyVars.classifyVars( allvcffiles["scores"], 
                                              #cleanvcffile, sampleannot,
                                              #regionlist, 
                                              #filepats,
                                              #"./results/classes", 
                                              #force=True )
    #classifyVars.plotClassDist( geneannotfile, regionlist )
 
################################################################################
# Main
################################################################################
if __name__ == "__main__" :

    optlist, args = getopt.getopt( sys.argv[1:], "r:ot")
    optlist = dict(optlist)

    dataset = optlist.get("-r",None)
    if dataset is None and not optlist.has_key("-t"):
        print "Error: no dataset provided"
        print " -r <string> with the name of the dataset"
        print " -o flag to overwrite everything"
        sys.exit(1)

    #filename = "/home/escott/workspace/variome/rawdata/onekg/onekg.chimp.regions.filt.samp.samp.vcf.gz"
    #filename = "/home/escott/workspace/variome/rawdata/test/everything_set1.chr1.snp.chimp.regions.filt.samp.samp.vcf.gz"
    #filename = "/home/escott/workspace/variome/rawdata/merged/merged.chimp.regions.filt.samp.samp.vcf.gz"
    #filename = "/home/escott/workspace/variome/rawdata/daily/daily.chimp.regions.filt.samp.samp.vcf.gz"
    #filename = "/home/escott/workspace/variome/rawdata/variome1/variome.chimp.regions.filt.samp.samp.vcf.gz"
    path = os.path.abspath("./rawdata/")
    if dataset == "test" :
        #vcffile = path+"/test/everything_set1.chr1.snp.chimp.vcf.gz"
        vcffile = path+"/test/everything_set1.chr1.snp.vcf.gz"
    elif dataset == "test2" :
        vcffile = path+"/test2/test2.vcf.gz"
    elif dataset == "onekg" :
        #vcffile = path+"/onekg/onekg.chimp.vcf.gz"
        vcffile = path+"/onekg/onekg.vcf.gz"
    elif dataset == "daily" :
        #vcffile = path+"/daily/daily.chimp.vcf.gz"
        vcffile = path+"/daily/daily.vcf.gz"
    elif dataset == "mevariome" :
        vcffile = path+"/mevariome/variome.vcf.gz"
    elif dataset == "mergedaly" :
        vcffile = path+"/mergedaly/meceu.vcf.gz"
    elif dataset == "merge1kg" :
        vcffile = path+"/merge1kg/me1000G.vcf.gz"
    elif dataset == "merge1kgimp" :
        vcffile = path+"/merge1kg/me1000G_imp.vcf.gz"
    elif dataset == "mergedalyimp" :
        vcffile = path+"/mergedalyimp/meceu_imp.vcf.gz"
    elif dataset == "casanova" :
        vcffile = path+"/casanova/casanova.snp.recal.chimp.vcf.gz"
    else :
        print "Error: dataset not found", dataset
        sys.exit(1)

    print "Using dataset:", dataset
    keepfile = "%s/%s/tokeep.txt" % (path, dataset)
    patientannotation = "./resources/annotation/patientannotation.ped"
    filepath, basename, ext = pop.getBasename( vcffile )
    makeCleanVCF( vcffile, toforce=False )

# END MAIN
    #elif dataset == "variome" :
        #vcffile = path+"/variome/variome.vcf.gz"
    #elif dataset == "variome1" :
        ##vcffile = path+"/variome1/variome.chimp.vcf.gz"
        #vcffile = path+"/variome1/variome.vcf.gz"
    #elif dataset == "variome2" :
        #vcffile = path+"/variome2/variome.recal.vcf.gz"
    #elif dataset == "hgdp" :
        #vcffile = path+"/hgdp/HGDP_938.chimp.vcf.gz"
    #elif dataset == "merged" :
        #vcffile = path+"/merged/merged.vcf.gz"
    #elif dataset == "merged1" :
        #vcffile = path+"/merged1/merged.vcf.gz"

