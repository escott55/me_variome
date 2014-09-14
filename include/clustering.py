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

def plotRDendrograms( bedfile, famfile, distfile, mdsfile, nclust=10, force=False ) :
    print "Running plotRDendrograms"
    filepath,basename,suffix = hk.getBasename(bedfile)
    fbase = basename[:basename.find(".")]
    #mydata = os.path.join(filepath,basename)
    #genomefile = mydata+".genome"

    clustfile = os.path.join(filepath,"ibs_%s_ward.txt"%fbase)
    if os.path.exists(clustfile) and not force :
        return clustfile

    rcmd = ('source("rscripts/dendrogram.R")\n'+
            'famfile<-"%s"\n'%famfile+
            'distfile<-"%s"\n'%distfile+
            'mdsfile<-"%s"\n'%mdsfile+
            'nclust<-%d\n'%nclust+
            'targetdir<-"%s"\n'%filepath+
            'plotIBSdendrogram(famfile,distfile,"ward",nclust,targetdir)\n'+
            'plotMDSdendrogram(famfile,mdsfile,"ward",targetdir)\n'+
            '\n'
           )
    print rcmd
    out = robjects.r( rcmd )
    assert os.path.exists(clustfile), "No clust file produced!"
    return clustfile
# END plotRDendrograms

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
        try :
            out = subprocess.check_output(command)
        except subprocess.CalledProcessError as e :
            print "Plink error:",e; sys.exit(1)

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

################################################################################
def parseClustFile( bedfile, clustfile ) :
    print "Running parseClustFile"
    filepath,basename,suffix = hk.getBasename(bedfile)
    #fbase = basename[:basename.find(".")]
    newannotfile = os.path.join(filepath,basename+".annot")
    famfile = os.path.join(filepath,basename+".fam")
    clustdata = read_csv( clustfile, sep="\t" )

    #print hk.dfTable(clustdata[clustdata.clust == 8].group)
    
    grdict = {"Arabian Peninsula":"AP","Central Asia":"CA",
              "Northeast Africa":"NEA","Northwest Africa":"NWA",
              "Syrian Desert":"SD","Turkish Peninsula":"TP"}
    grpcounts = []
    for grp, cdata in clustdata.groupby(["clust"]) :
        gtbl = hk.dfTable( cdata["group"] ).reset_index(drop=True)
        grpname = "None"
        if ((gtbl.ix[0,"Variable"] == "Unknown") and
            (gtbl.ix[1,"Variable"] == "Africa")) :
            grpname = "Africa"
        elif ((gtbl.ix[0,"Variable"] == "Northeast Africa") and
           (len(cdata) == 74)) : 
            grpname = "Syrian Desert"
        else : 
            grpname = gtbl.ix[0,"Variable"]

        if grdict.has_key( grpname ) : grpname = grdict[grpname]
        grpcounts.append(Series({"grpname":grpname,"clust":grp,
                                 "nsamp":gtbl["Count"].sum()}))

    grpcounts = DataFrame(grpcounts)
    #print grpcounts.head()

    #newgrps = []
    grpcounts["GeographicRegions3"] = grpcounts.grpname
    for grp, data in grpcounts.groupby("grpname") :
        if len(data) > 1 and grp != "Unknown":
            cnt = 1
            for idx,row in data.iterrows() :
                grpcounts.ix[idx,"GeographicRegions3"] = grp+str(cnt)
                cnt += 1 

    print grpcounts.head()

    keepcols = ["label","clust"]
    clustannot = merge( clustdata[keepcols], grpcounts, on="clust" )
    clustannot.rename(columns={"label":"Individual.ID"},inplace=True)

    famdata = read_csv(famfile, delim_whitespace=True, header=None, 
                       names=["FID","Individual.ID","MID","PID","Gender","AFF"] )
    clustannot = merge( famdata[["FID","Individual.ID"]],clustannot, 
                       on="Individual.ID",how="left") 

    print clustannot.head()
    clustannot = addSampleAnnotation( clustannot, "Individual.ID" )
    clustannot.loc[clustannot.GeographicRegions2.isnull(),"GeographicRegions2"] = "Unknown"
    tofix = clustannot.GeographicRegions3.isnull()
    clustannot.loc[tofix,"GeographicRegions3"] = clustannot.loc[tofix,"GeographicRegions2"].tolist()
    tofix = clustannot.GeographicRegions3.isin(["America","East Asia","Oceania"])
    clustannot.loc[tofix,"GeographicRegions3"] = "Unknown"
    print "Writing file:",newannotfile
    clustannot.to_csv( newannotfile, sep="\t", index=False )

    return newannotfile
# END parseClustFile

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
 
def clusterWorkflow( vcffile, nclust, force=False ) :
    bedfile = seriousClean( vcffile, False )

    famfile, distfile, mdsfile = runClustering( bedfile, force=False )

    clustfile = plotRDendrograms( bedfile, famfile, distfile, mdsfile, nclust, True )

    newannotfile = parseClustFile( bedfile, clustfile )
# END clusterWorkflow

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
    #vcffile = "./rawdata/mevariome/variome.vcf.gz"
    vcffile = "./rawdata/merge1kg/main/me1000G.clean.vcf.gz"
    #vcffile = "./rawdata/merge1kg/me1000G.vcf.gz"
    vcffile = "./rawdata/onekg/main/onekg.clean.vcf.gz"
    #vcffile = "./rawdata/mergedaly/main/meceu.clean.vcf.gz"
    vcffile = "./rawdata/daily/main/daily.clean.vcf.gz"
    #vcffile = "./rawdata/test2/main/test2.clean.vcf.gz"
    
    #bedfile = seriousClean2( vcffile, True )

    if vcffile.find("daily") : nclust = 2
    else : nclust = 10
    clusterWorkflow( vcffile, nclust )
# END MAIN
