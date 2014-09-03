#!/usr/bin/env python

#-------------------------------------------------------------------#
#                       Import Packages                             #
#-------------------------------------------------------------------#
import os, sys
import getopt
import csv
import re
import subprocess

import housekeeping as hk
import popgencommands as pop
from localglobals import *

################################################################################
# modifyFamFile
################################################################################
def modifyFamFile( famfile, targetfile, annotcol="Continent" ):
    print "Using annotation column:",annotcol
    famdata = read_csv(famfile, delim_whitespace=True, header=None,
                               names=["Fam","IID","MID","PID","Gender","AFF"])
    famdata = famdata.fillna(0)
    #famdata[famdata.MID == "nan"] = 0
    #famdata[famdata.PID == "nan"] = 0
    famdata.index = famdata.IID
    #sampleannot = sampleAnnotation()
    famdata[annotcol] = "Unk"
    #famdata.update(sampleannot)
    famdata = addSampleAnnotation( famdata, update=True )
    famdata[annotcol] = famdata[annotcol].apply(lambda x:
                                                x.strip().replace(' ','.'))
    famdata["Gender"] = famdata.Gender.astype(int)
    famdata.to_csv(targetfile, sep=" ", header=False, index=False,
                   columns=["Fam","IID","MID","PID","Gender",annotcol])
    return famdata[annotcol].unique().tolist()
# END modifyMapFile

################################################################################
def runClustering( bedfile, force=False, mdsplot=20 ):
    print "Running runClustering"
    filepath,basename,suffix = hk.getBasename(bedfile)
    mydata = os.path.join(filepath,basename)
    genomefile = mydata+".genome"

    annotfamfile = mydata+"_annot.fam"
    mdistfile = mydata+".mdist"
    mdsfile = mydata+".mds"

    if os.path.exists(mdistfile) and os.path.exists(mdsfile) and not force : 
        return annotfamfile, mdistfile, mdsfile

    famdata = read_csv( mydata+".fam", delim_whitespace=True, header=None,
                       names=["FID","IID","MID","PID","SEX","AFF"])

    famdata = addSampleAnnotation( famdata, "IID" )
    print famdata.head()
    print "writing file :",annotfamfile
    famdata.to_csv( annotfamfile, sep="\t",index=False )

    if not os.path.exists(genomefile) or force :
        command = ["plink","--bfile",mydata,"--genome","--out",mydata]
        print " ".join(command)
        out = subprocess.check_output(command)

    # Produce NxN matrix
    print "NxN cluster"
    command = ["plink",
               "--bfile",os.path.join(filepath,basename),
               "--read-genome",genomefile,
               "--cluster",
               "--out",os.path.join(filepath,basename)]
    print " ".join(command)
    out = subprocess.check_output(command)

    # Produce IBS distance Matrix
    print "IBS clustering"
    command = ["plink","--bfile",os.path.join(filepath,basename),
               "--read-genome",genomefile,
               "--cluster","--distance-matrix",
               "--out",os.path.join(filepath,basename)]
    print " ".join(command)
    out = subprocess.check_output(command)

    # Run MDS clustering on IBS distance matrix
    print "MDS cluster"
    command = ["plink","--bfile",os.path.join(filepath,basename),
               "--read-genome",genomefile,
               "--cluster","--mds-plot",str(mdsplot),
               "--out",os.path.join(filepath,basename)]
    print " ".join(command)
    out = subprocess.check_output(command)

    return annotfamfile, mdistfile, mdsfile
# END runClustering

def seriousClean( vcffile, rerun=False ):
    print "Running seriousClean"
    patientdata = sampleAnnotation()
    filepath, basename, suffix = hk.getBasename(vcffile)
    targetdir = "%s/clust/" % filepath
    hk.makeDir(targetdir)
    cptargetfile = targetdir+basename+suffix+".gz"
    cleanped = os.path.join(targetdir,basename)

    if not os.path.exists(cptargetfile) :
        out = hk.runCommand( "cp %s/%s%s* %s" % (filepath,basename,suffix, targetdir) )
        print "Copy targetfile:",cptargetfile
    
    print "Making sample missingness estimates", cleanped
    #command = ["vcftools","--gzvcf", cptargetfile,
               #"--missing-indv","--out",cleanped]
    #out = subprocess.check_output( command )

    command = ["vcftools","--gzvcf", cptargetfile,
               "--remove-filtered-all",
               "--remove-indels",
               "--maf","0.001",
               "--hwe","0.001",
               "--min-alleles","2",
               "--max-alleles","2",
               "--plink-tped","--out",cleanped]
    
    print " ".join(command)
    tped = cleanped+".tped"
    if not os.path.exists( tped ) or rerun:
        out = subprocess.check_output( command )
    
    #recodefile = plink_recode12( tped, force=rerun )
    # Modify ped
    #modped = pop.modify_ped( recodefile, patientdata, calcdistances=True, shorten=True, force=rerun )
    #pop.modify_tped( tped, patientdata, force=rerun )
    #print modped
    #frqfile = plink_makefrqfile( recodefile, force=rerun )
    #newfrq, excludebase = pop.excludeSnps( recodefile, force=rerun )
    #bedfile = pop.run_plink_convert(excludebase+".ped", force=rerun)
    bedfile = pop.run_plink_convert(tped, force=rerun)
    return bedfile
    #print vcffile
# END seriousClean
 
######################################################################
if __name__ == "__main__":

    #optlist, args = getopt.getopt( sys.argv[1:], "bk")
    #optlist = dict(optlist)
    #bedfile = optlist.get("-b",bedfile)
    
    # Change running directory
    os.chdir("..")

    #bedfile = "./rawdata/variome1/variome.chimp.regions.filt.samp.plink.mod.filt.uniq.bed"
    #bedfile = "./rawdata/onekg/onekg.chimp.regions.filt.samp.plink.filt.bed"
    #bedfile = "./rawdata/test/everything_set1.chr1.snp.clean.plink.filt.bed"
    #bedfile = "./rawdata/mevariome/main/ibc/variome.clean.Middle_East.recode12.mod.bed"
    #bedfile = "./rawdata/mevariome/main/variome.clean.plink.filt.bed"
    #bedfile = "./rawdata/mevariome/main/variome.clean.plink.filt.bed"
    #vcffile = "./rawdata/merge1kg/main/me1000G.clean.vcf.gz"
    #vcffile = "./rawdata/merge1kg/main/me1000G.clean.vcf.gz"
    vcffile = "./rawdata/mergedaly/meceu.X.vcf.gz" 
    vcffile = "./rawdata/merge1kg/me1000G.X.vcf.gz"
    vcffile = "./rawdata/mevariome/variome.X.vcf.gz"
    vcffile = "./rawdata/mevariome/main/variome.clean.vcf.gz"
    #vcffile = "./rawdata/merge1kg/main/me1000G.clean.vcf.gz"
    #vcffile = "./rawdata/onekg/main/onekg.clean.vcf.gz"
    #vcffile = "./rawdata/test2/main/test2.clean.vcf.gz"
    
    #bedfile = seriousClean2( vcffile, True )
    bedfile = seriousClean( vcffile, False )

    famfile, distfile, mdsfile = runClustering( bedfile, force=True )

# END MAIN
