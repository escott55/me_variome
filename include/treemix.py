#!/usr/bin/env python

#-------------------------------------------------------------------#
#                       Import Packages                             #
#-------------------------------------------------------------------#

import os, sys
import getopt
import csv
import re
import subprocess
import gzip

import housekeeping as hk
import popgencommands as pop
from localglobals import *

#1. outstem.cov.gz. The covariance matrix between populations estimated from the data
#2. outstem.covse.gz. The standard errors for each entry in the covariance matrix
#3. outstem.modelcov.gz. The fitted covariance (W in Pickrell and Pritchard [2012]) according
#   to the model
#4. outstem.treeout.gz. The fitted tree model and migration events
#5. outstem.vertices.gz. This and the following file (outstem.edges.gz) contain the internal
#   structure of the inferred graph. Modifying these files will cause issues if you try 
#   to read the graph back in, so we recommend against this.
#6. outstem.edges.gz.

################################################################################
# plink_makefrqfile
################################################################################
def plink_makefrqfile( bedfile, clusterfile, force=False ) :
    print "Running plink_makefrqfile"
    filepath,basename,suffix = hk.getBasename(bedfile)
    # Make frqfile
    frqfile = filepath+"/"+basename+".frq.strat"
    gzfrqfile = frqfile+".gz"
    #tmpfile = filepath+"/"+basename+".frq"
    #print "Checking frq file:",frqfile, os.path.exists( frqfile )
    # --1
    if not os.path.exists( gzfrqfile ) or force :
        command = ["plink","--noweb",
                   "--bfile", bedfile[:-4],
                   "--freq","--allow-no-sex",
                   "--within", clusterfile,
                   "--out",filepath+"/"+basename]
        print "Command:"," ".join(command)
        out = subprocess.check_output( command )
        ret = subprocess.call( ["gzip",frqfile] )
        assert ret == 0 
        assert os.path.exists( gzfrqfile )
    return gzfrqfile
# END plink_makefrqfile
    #print out
    #command = "awk '$1 !~ /CHR/ {print $1,$2,$3,$4,$5,$6}' %s > %s" % (tmpfile, frqfile )
    #print command
    #out = hk.runCommand( command )

# functionsfolder = "/home/escott/Packages/treemix-1.12/src"
# source("/home/escott/Packages/treemix-1.12/src/plotting_funcs.R")
#>source("src/plotting funcs.R")
#>plot tree("outstem")
#>plot resid("outstem", "poporder")
def plotTree( prefix, pop_order ) :
    print "Running plotTree -",prefix
    #print robjects.r("3*4")
    treefigfile = prefix+"_treetest.png"
    covfigfile = prefix+"_treecov.png"
    residfigfile = prefix+"_treeresid.png"

    poporderfile = prefix+".poporder"
    OUT = open(poporderfile, "w")
    for pid in pop_order :
        print >> OUT, pid
    OUT.close()

    print "Making figure:",treefigfile
    rcmd = ("source('/home/escott/Packages/treemix-1.12/src/plotting_funcs.R')\n"+
        'png("'+treefigfile+'")\n'+
        'plot_tree("'+prefix+'",mbar=F)\n'+
        'dev.off()\n' +
        'png("'+covfigfile+'")\n'+
        'plot_cov("'+prefix+'", "'+poporderfile+'")\n'+
        'dev.off()\n' +
        'png("'+residfigfile+'")\n'+
        'plot_resid("'+prefix+'", "'+poporderfile+'")\n'+
        'dev.off()\n')
        #'pop_order <- c("'+'","'.join(pop_order)+'")\n'

    #print "R command", rcmd
    out = robjects.r( rcmd )
    "Make figures out:"
    #print out
    assert os.path.exists(treefigfile)
# END plotTree

def makeClusterFile( famfile, pfactor="Continent" ) :
    filebase, ext = os.path.splitext( famfile )
    patientdata = sampleAnnotation()

    famdata = read_csv( famfile, header=None, delim_whitespace=True,
                       names=["FID","IID","PID","MID","gender","aff"] )
    #print famdata.head(10)
    famdata["IID"] = [x[:x.find(".")] if x.find(".") > 0 else x
                      for x in famdata["IID"].tolist()]
    tids = famdata["IID"].tolist()
    patientdata["IID"] = [x[:x.find("_")] if x[:x.find("_")] in tids 
                          and x.find("JL") > 0
                          else x
                          for x in patientdata["Individual.ID"]]
    alldata = merge(famdata, patientdata, on="IID", how="left")
                    #left_on="IID", right_on="Individual.ID",how="left")
    #print alldata.head()
    assert len(alldata) == len(famdata)

    alldata["cluster"] = factorize( alldata[pfactor] )[0]
    #print alldata[["FID","IID","cluster",pfactor]].head()
    clusterfile = filebase+".clust"

    print "Writing file:", clusterfile
    alldata[["FID","IID","cluster",pfactor]].to_csv(clusterfile, 
                                                    index=False, header=False,sep="\t")
    return clusterfile
# END makeClusterFile

def treemixRegions( treemixfile ) :
    FH = gzip.open( treemixfile )
    head = FH.readline()
    FH.close()
    regions = head.rstrip().split(" ")
    return regions
# END treemixRegions

# plink2treemix.py 
# first compress, then convert
# plink --bfile data --freq --missing --within data.clust
# gzip plink.frq
# plink2treemix.py plink.frq.gz treemix.frq.gz
def convertPlink2Treemix( bedfile, pfactor="Continent", force=False ):
    print "Running convertPlink2Treemix -",bedfile
    filepath, basename, suffix = hk.getBasename(bedfile)
    famfile = "%s/%s.fam" %(filepath,basename)
    mapfile = "%s/%s.map" %(filepath,basename)

    treemixfile = "%s/%s.tree.frq.gz" % (filepath,basename)
    if os.path.exists(treemixfile) and not force: 
        print "Warning: treemix file already exists -",treemixfile
        regions = treemixRegions( treemixfile )
        return treemixfile, regions

    clusterfile = makeClusterFile( famfile, pfactor )
    frqfile = plink_makefrqfile( bedfile, clusterfile )
    subprocess.call(["plink2treemix.py",frqfile,treemixfile,clusterfile])

    regions = treemixRegions( treemixfile )
    return treemixfile, regions
# END convertPlink2Treemix

# treemix -i input_file.gz -o out_stem
# -root outgroup
# -k blocks of snps 1000
# -m migration events (2)
# -cor_mig -climb Known migration events
# -noss sample size correction
def runTreeMix( treemixfile, migrationevents=0, root=None, annotcol=None, force=False ) :
    filepath, basename, suffix = hk.getBasename(bedfile)
    targetbase = os.path.join(filepath,basename)+"_m"+str(migrationevents)
    if annotcol is not None  : targetbase = targetbase+"_"+annotcol
    verticesfile = targetbase+".verticies.gz"
    #if os.path.exists(verticesfile) : return targetbase
    treemixout = "%s.treemixout" % (targetbase)
    if os.path.exists(treemixout) and not force : 
        return targetbase

    print "Running runTreeMix -",treemixfile
    treemixcommand = ["treemix",
                      "-i",treemixfile,
                      "-bootstrap",
                      "-m",str(migrationevents),
                      "-o", targetbase] #"-se",

    if root is not None: treemixcommand = treemixcommand+["-root",root]
    print "Treemix command"
    print " ".join(treemixcommand)
    out = subprocess.check_output( treemixcommand, stderr=subprocess.STDOUT )
    open(treemixout,"w").writelines(out)
    return targetbase
# END runTreeMix

################################################################################
def seriousClean( vcffile, keepfile=None, rerun=False ):
    print "Running seriousClean"
    filepath, basename, suffix = hk.getBasename(vcffile)

    cleanped = "%s/%s"%(filepath,basename)
    #bedfile = cleanped+".bed"
    bedfile = cleanped+".filt.bed"
    if os.path.exists(bedfile) and not rerun : 
        return bedfile

    command = ["vcftools","--gzvcf", vcffile,
               "--remove-filtered-all",
               "--remove-indels",
               "--maf",".005",
               "--min-alleles","2",
               "--max-alleles","2",
               "--plink-tped","--out",cleanped]
    
    if keepfile is not None : 
        assert os.path.exists(keepfile)
        command = command+["--keep",keepfile]

    print " ".join(command)
    tped = cleanped+".tped"
    if not os.path.exists( tped ) or rerun:
        out = subprocess.check_output( command )

    #bedfile = pop.run_plink_convert(tped, force=rerun)
    bedfile = pop.run_plink_filter(tped, force=rerun)
    return bedfile
# END seriousClean
    #cptargetfile = os.path.join(targetdir,basename+suffix+".gz")
    #if not os.path.exists(cptargetfile) :
        #out = hk.runCommand( "cp %s/%s%s* %s" % (filepath,basename,suffix, targetdir) )
        #print "Copy targetfile:",cptargetfile
    #recodefile = plink_recode12( tped, force=rerun )
    #print "#"*50
    #print "Made recodefile:",recodefile
    ## Modify ped
    #modped = pop.modify_ped( recodefile, patientdata, calcdistances=True, shorten=True, force=rerun )
    #frqfile = plink_makefrqfile( modped, force=rerun )
    #newfrq, excludebase = pop.excludeSnps( modped, force=rerun )

################################################################################
def getLevels( annotcol, targetvcf ):
    print "Running getLevels"
    if annotcol == "GeographicRegions2" :
        levels = target_geographic_regions2
    elif annotcol == "Continent" :
        levels = target_continents
    elif annotcol == "GeographicRegions" :
        levels = target_geographic_regions
        #if "Oceania" in levels : levels.remove("Oceania")
        #if "America" in levels : levels.remove("America")
        print levels
    else :
        print "Error: unknown levels for column -",annotcol
        sys.exit(1)

    levels = [x.strip().replace(" ",".") for x in levels]
    targetpats = patientInfo.getPats( targetvcf )
    sampleannot = sampleAnnotation(targetpats)

    sampleannot = sampleannot[(sampleannot[annotcol].notnull()) &
                                (sampleannot[annotcol] != "Unknown")]
    sampleannot[annotcol] = [x.strip().replace(" ",".") for x in sampleannot[annotcol]]
    
    finallevels = [x for x in levels if x in sampleannot[annotcol].unique()]
    sampleannot = sampleannot[sampleannot[annotcol].isin(finallevels)]
    #keepfiles = sampleannot["Individual.ID"].tolist()
    return finallevels, sampleannot#, keepfiles
# END getLevels

################################################################################
if __name__ == "__main__":

    #optlist, args = getopt.getopt( sys.argv[1:], "bk")
    #optlist = dict(optlist)
    # Change running directory
    os.chdir("..")

    #bedfile = "./rawdata/mevariome/main/pca/variome.clean.recode12.mod.exclude.bed"
    #vcffile = "./rawdata/mevariome/main/variome.clean.vcf.gz"
    vcffile = "./rawdata/merge1kg/main/me1000G.clean.vcf.gz"
    #vcffile = "./rawdata/onekg/main/onekg.clean.vcf.gz"
    #vcffile = "./rawdata/test2/main/test2.clean.vcf.gz"

    #filepath, filename = os.path.split(vcffile)
    #targetdir = os.path.join(filepath,"tree")
    #hk.makeDir(targetdir)

    annotcol = "GeographicRegions2"
    targetvcf = hk.copyToSubDir( vcffile, "tree/"+annotcol )
    filepath, filename, suffix = hk.getBasename(targetvcf)

    levels, sampleannot = getLevels( annotcol, targetvcf )
    print levels
    print hk.dfTable(sampleannot[annotcol])
    
    # Make Keep file
    keepfile = os.path.join(filepath, filename+".keeppats")
    sampleannot[["Individual.ID"]].to_csv(keepfile, header=None,index=None)

    bedfile = seriousClean( targetvcf, keepfile, rerun=False )

    treemixfile, regions = convertPlink2Treemix( bedfile, annotcol, force=False )
    #print treemixfile
    print "Regions:",regions
    
    root = None
    if "Africa" in regions : root = "Africa"
    if "YRI" in regions : root = "YRI,LWK"
    else : None

    #for numevents in [0] :
    for numevents in [0,1,2,3,4] :
        treemixoutbase = runTreeMix( treemixfile, numevents, root, annotcol, force=True )
        treemixsuffixes = ["cov","covse","edges","modelcov","treeout","verticies"]
        for suffix in treemixsuffixes :
            cfile = "%s.%s.gz" %(treemixoutbase,suffix)
            print "for file:",suffix,"file exists?",os.path.exists(cfile)

        plotTree( treemixoutbase, regions )

    #bedfile = optlist.get("-b",bedfile)
    #assert os.path.exists(bedfile)

    #filepath, filename = os.path.split(bedfile)
    #targetdir = os.path.join(filepath,"pca")

    #annotname = ["GeographicRegions","Origin","Continent","ethnicity","Source"]
    #for annot in annotname :
        #outliers = pcaOutlierAnalysis( bedfile, targetdir, force=True, annotname=annot )
        #print outliers

# END Main

