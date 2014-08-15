#!/usr/bin/env python

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
import getopt
import csv
import re
import gzip
import subprocess
from pandas import *
#from random import randint

from housekeeping import *


#-------------------------------------------------------------------#
#                            GLOBALS                                #
#-------------------------------------------------------------------#
DBSNPLIST = "rawdata/dbsnp_137_idlist.hg19.txt"
DEBUG = False
#-------------------------------------------------------------------#
#                            Commands                               #
#-------------------------------------------------------------------#

######################################################################
# plink_fileinput
######################################################################
def plink_fileinput( pedfile ) :
    filepath, basename, suffix = getBasename(pedfile)
    if suffix == ".ped" :
        pfile = "--file %s/%s" % (filepath,basename)
    elif suffix == ".bed" :
        pfile = "--bfile %%s/s" % (filepath,basename)
    elif suffix == ".tped" :
        pfile = "--tfile %%s/s" % (filepath,basename)
    else : print "Error: unknown suffix", suffix; sys.exit(1)
    return pfile
# End plink_fileinput

######################################################################
# parsePatientPed
# Input - ped
# Output - dictionary patid -> [patdata]
######################################################################
def parsePatientPed( pedfile ) :
    write2log( " - Running:"+whoami()+" file:"+pedfile, True )
    assert( os.path.exists(pedfile) )
    FILE = open(pedfile, "r")
    reader = csv.reader(FILE,delimiter="\t")
    data = {}
    header = []
    for row in reader :
        #print "Row:",row[:6]
        if len(row) < 1 :
            continue
        if len(header) == 0 :
            header = row
            continue
        #if row[0] == "1884" :
            #print row
        for i in range(len(header)) :
            if header[i] in ["Individual ID","original","Individual.ID"] :
                data[row[i].strip()] = row
                #print "Patient:",row[i]
                patid = parse_patid(row[i].strip())
                data[patid] = row
        #print "\t".join(row)
    FILE.close()
    print "Num Patient entries:", len(data)
    data["header"] = header
    return data
# End parsePatientPed

######################################################################
# parseAnnotation
######################################################################
def parseAnnotation( ) :
    print "Running parseAnnotation"
    alnstats = read_csv("resources/annotation/alignmentstats.txt", sep="\t")
    alnstats = (alnstats[alnstats["bamfile"].notnull()]
                [["Patient","ethnicity","ethnicity_prob","continent","continent_prob",
                  "country","country_prob"]])

    daly = read_csv("resources/annotation/daly_annot.ped",sep="\t")
    eichler = read_csv("resources/annotation/eichler_annot.ped",sep="\t")
    onekg = read_csv("resources/annotation/onekg_annot.ped",sep="\t")
    hgdp = read_csv("resources/annotation/hgdp_annot.ped",sep="\t")
    casanova = read_csv("resources/annotation/casanova_annot.ped",sep="\t")
    fowzan = read_csv("resources/annotation/fowzan_annot.ped",sep="\t")
    gleeson = read_csv("resources/annotation/gleeson_annot.ped",sep="\t")

    basiccols = ["Family.ID","Individual.ID","Gender","Phenotype","Source"]
    targetcolumns = [ "Family.ID", "Individual.ID", "Paternal.ID", "Maternal.ID", "Gender", "Phenotype",
               "Other.info", "Origin", "Consang", "Ethnicity", "Plate", "Source", "ethnicity",
               "ethnicity_prob", "continent", "continent_prob", "country", "country_prob"]
    if DEBUG : print casanova[basiccols].head(10)
    #print "Daly",[x for x in targetcolumns if x not in daly.columns.tolist()]
    #print "Onekg",[x for x in targetcolumns if x not in onekg.columns.tolist()]
    #print "gleeson",[x for x in targetcolumns if x not in gleeson.columns.tolist()]
    #print "fowzan",[x for x in targetcolumns if x not in fowzan.columns.tolist()]
    #print "casanova",[x for x in targetcolumns if x not in casanova.columns.tolist()]
    newdata = concat([casanova,fowzan,gleeson])
    allannot = merge( newdata, alnstats, how="left",left_on="Individual.ID",right_on="Patient" )
    allannot = concat([allannot,onekg,daly,hgdp])
    allannot.index = allannot["Individual.ID"]
    if DEBUG : print allannot.shape
    if DEBUG : print allannot[allannot.Patient.isnull()].head(10)
    allannot["Paternal.ID"] = allannot["Paternal.ID"].fillna("0")
    allannot["Maternal.ID"] = allannot["Maternal.ID"].fillna("0")
    allannot["Gender"] = allannot["Gender"].fillna("-9").map(str)
    allannot["Phenotype"] = allannot["Phenotype"].fillna("-9").map(str)
    allannot.to_csv( "resources/annotation/patientannotation.ped", sep="\t", 
                    index=False, cols=targetcolumns )
    if DEBUG : print allannot[basiccols].head(20)
    return allannot
# End parseAnnotation

######################################################################
# run_vcftools_plink
# Input - vcffile
# Process - if toremove.txt, remove samples, if keeplist keepsamples
#    filter indels, non-bi-allelic
# Output - ped
######################################################################
def run_vcftools_plink( vcffile, outdir=None, keeplist=None, force=False ) :
    write2log( " - Running:"+whoami(), True )

    targetdir, filename, suffix = getBasename(vcffile)
    if outdir is None : outdir = targetdir
    if outdir[-1] == "/" :
        plinkfile = "%s%s.plink" % (outdir, filename)
    else :
        plinkfile = "%s/%s.plink" % (outdir, filename)
    #print plinkfile

    logfile = "%s/%s_vcftools.log" % (LOGDIR,filename)

    finalped = plinkfile+".tped"
    if os.path.exists( finalped )  and not force :
        write2log(" - File"+plinkfile+".tped already exists. Skipping command!",True)
        return plinkfile+".tped"

    otheroptions = ""
    removefile = targetdir+"/toremove.txt"
    if keeplist is not None :
        for patient in open(keeplist).readlines() :
            if len(patient) == 0 or patient[0] == "#" :
                continue
            otheroptions += ' --indv "%s"' % patient.strip()
        #otheroptions += ' --keep %s' % keeplist
    elif os.path.exists( removefile ):
        print "Remove File:",removefile
        f = open( removefile, "r" )
        for patient in f.readlines():
            if len(patient) == 0 or patient[0] == "#" :
                continue
            otheroptions += ' --remove-indv "%s"' % patient.strip()
        #otheroptions = "--remove %s/toremove.txt" % targetdir
        f.close()

    # thin <int> filter variants within a given proximity
    # remove-indels
    # mix/max-alleles - Preven n allelic sites
    # snps <filename> include snps based on list of ids
    command = "vcftools --gzvcf %s" \
        " --remove-indels --min-alleles 2 --max-alleles 2" \
        " %s" \
        " --plink-tped --out %s" \
        % ( vcffile, otheroptions, plinkfile )
        #" --maf .002" \
        #" --remove-filtered-geno-all " \
        #" --FILTER-summary"\
        #" --plink-tped --out %s" \
        #" --snps %s" \
        #% ( vcffile, otheroptions, DBSNPLIST, plinkfile )
        #" --min-meanDP 4 --minQ 20" \
        #" --thin 10" \

    print command
    runout = runCommand( command )
    print "Writing logfile:",logfile
    OUT = open(logfile,"wb")
    OUT.write(runout)
    OUT.close()
    print runout

    #cfile = compress( finalped )

    return finalped
    #command = "%s/bwa aln -q %d -f %s -t %d %s %s" % (BINDIR, qualcutoff, saifile, NUMTHREADS, BWAIDX, fastqfile )
    #runCommand( command )
# End run_vcftools_plink

######################################################################
# run_vcftools
# Input - vcffile, keeplist
# Output - vcffile recoded
######################################################################
def run_vcftools_filtsamples( vcffile, keeplist, force=False ) :
    write2log( " - Running:"+whoami(), True )

    targetdir, filename, suffix = getBasename(vcffile)
    print "tdir:",targetdir,"fname",filename,"suffix:",suffix

    if targetdir[-1] == "/" :
        newvcffile = "%s%s" % (targetdir, filename)
    else :
        newvcffile = "%s/%s" % (targetdir, filename)

    print "looking for ",newvcffile+".samp.vcf.gz"
    if os.path.exists( newvcffile+".samp.vcf.gz" ) and not force :
        print " - File"+newvcffile+".samp.vcf.gz : already exists. Skipping command!"
        return newvcffile+".samp.vcf.gz"

    if keeplist is None or len(keeplist) == 0 :
        print "Error: keeplist invalid:", keeplist; sys.exit(1)

    otheroptions = " ".join(['--indv "%s"' % x for x in keeplist])

    command = "vcftools --gzvcf %s" \
        " %s --recode" \
        " --out %s" \
        % ( vcffile, otheroptions, newvcffile )

    print command
    out = runCommand( command )
    print out

    assert os.path.exists( newvcffile+".recode.vcf" ), out
    runCommand( "mv %s.recode.vcf %s.samp.vcf" % (newvcffile,newvcffile) )
    cfile = compress( newvcffile+".samp.vcf" )
    return cfile
# END run_vcftools_filtsamples
 
################################################################################
# run_vcftools_clean
# Input - vcffile, keepfile
################################################################################
def run_vcftools_clean( vcffile, keepfile, outprefix, force=False ) :
    print "Running run_vcftools_clean",vcffile
    assert os.path.exists(vcffile)
    
    if os.path.exists( outprefix+".clean.vcf.gz" ) and not force :
        print " - File"+outprefix+".clean.vcf.gz : already exists. Skipping command!"
        return outprefix+".clean.vcf.gz"

    otheroptions = ''
    if os.path.exists(keepfile) :
        for patient in [x.strip() for x in open(keepfile).readlines()]:
            if len(patient) == 0 or patient[0] == "" : continue
            otheroptions += ' --indv "%s"' % patient

    # remove-indels 
    # mix/max-alleles - Preven n allelic sites
    # snps <filename> include snps based on list of ids
    # Geno - removes genotypes without sufficient coverage
    command = "vcftools --gzvcf %s" \
        " --min-alleles 2 --max-alleles 2" \
        " --max-missing 0.95" \
        " %s" \
        " --recode" \
        " --out %s" \
        % ( vcffile, otheroptions, outprefix )
        #" --remove-filtered-all" \
        #" --min-meanDP 4 --minQ 20" \
    
    print command
    out = runCommand( command )
        
    assert os.path.exists( outprefix+".recode.vcf" ), out
    runCommand( 'mv %s.recode.vcf %s.clean.vcf' % (outprefix,outprefix) )
    cfile = compress( outprefix+".clean.vcf", force )
    print " - End run_vcftools"
    return cfile
# END run_vcftools_clean

######################################################################
# run_vcftools
# Input - vcffile, keeplist
# Process -  Filter by allelism, missingness, hwe, dp, Qual, 
# Notes - regions crashes due to a lack of memory
# Output - vcffile recoded
######################################################################
def run_vcftools( vcffile, keeplist=None, maffilter=False, force=False ) :
    write2log( " - Running:"+whoami(), True )

    targetdir, filename, suffix = getBasename(vcffile)
    print "tdir:",targetdir,"fname",filename,"suffix:",suffix

    assert os.path.exists(vcffile)

    if vcffile[-3:] != ".gz" :
        vcffile = compress( vcffile )

    if targetdir[-1] == "/" :
        newvcffile = "%s%s" % (targetdir, filename)
    else :
        newvcffile = "%s/%s" % (targetdir, filename)
    #print newvcffile

    if os.path.exists( newvcffile+".filt.vcf.gz" ) and not force :
        print " - File"+newvcffile+".filt.vcf.gz : already exists. Skipping command!"
        return newvcffile+".filt.vcf.gz"

    otheroptions = " "
    #print "$"*70
    #print filename[:5], filename[:4]
    if filename[:5] != "onekg" and filename[:4] != "HGDP" :
        otheroptions = "--min-meanDP 8 --minQ 30 "

    if keeplist is not None :
        otheroptions += " ".join(['--indv "%s"' % x for x in keeplist])

    if maffilter :
        otheroptions += " --maf 0.1" \
                " --max-maf 1.0" \
                " --remove-indels" \
                " --remove-filtered-all" 

    if os.path.exists("./toremove.txt") :
        removefile = read_csv("./toremove.txt",header=False,comment='#').dropna(how='all').reset_index(drop=True)
        removefile.columns = ["Badids"]
        otheroptions += " "+" ".join(['--remove-indv "%s"' % x for x in removefile["Badids"].tolist()])

    print "VCF file:",filename
    print "targetdir:",targetdir
    print "targetfile:",newvcffile+".filt.vcf.gz"
    #print "Other options:",otheroptions

    # thin <int> filter variants within a given proximity
    # remove-indels
    # mix/max-alleles - Preven n allelic sites
    # snps <filename> include snps based on list of ids
    # Geno - removes genotypes without sufficient coverage
    command = "vcftools --gzvcf %s" \
        " --min-alleles 2 --max-alleles 2" \
        " --max-missing 0.95" \
        " --hwe 0.001" \
        " %s --recode" \
        " --out %s" \
        % ( vcffile, otheroptions, newvcffile )
        #"--snps %s" \
        #" --min-meanDP 4 --minQ 20" \
        #" --thin 10" \
        #" --chr 11"\

    print command
    out = runCommand( command )

    assert os.path.exists( newvcffile+".recode.vcf" ), out
    runCommand( 'mv %s.recode.vcf %s.filt.vcf' % (newvcffile,newvcffile) )
    cfile = compress( newvcffile+".filt.vcf" )
    print " - End run_vcftools"
    return cfile
    #command = "%s/bwa aln -q %d -f %s -t %d %s %s" % (BINDIR, qualcutoff, saifile, NUMTHREADS, BWAIDX, fastqfile )
    #runCommand( command )
# End run_vcftools
    #elif os.path.exists( removefile ):
        #print "Remove File:",removefile
        #f = open( removefile, "r" )
        #for patient in f.readlines():
            #patient = patient.strip()
            #if len(patient) == 0 or patient[0] == "#" :
                #continue
            #otheroptions += ' --remove-indv "%s"' % patient
        #f.close()

######################################################################
# run_plink_convert
# Input - ped
# Process - Makes a bed file, no filtering
# Output - bed
######################################################################
def run_plink_convert( pedfile, force=forceFlag ) :
    write2log( " - Running:"+whoami(), True )
    targetdir, filename, suffix = getBasename(pedfile)
    basename = "%s/%s" % (targetdir, filename)
    if os.path.exists( basename+".bed" ) and not force:
        write2log(" - File"+basename+" already exists. Skipping command!",True)
        return basename+".bed"
    tfile = "--file "+basename
    if suffix == ".tped" : tfile = "--tfile "+basename
    # --1
    command = "plink --noweb %s --make-bed --out %s" % ( tfile, basename )
    out = runCommand( command )
    print "run_plink_convert:",out
    return basename+".bed"
# End run_plink_convert

######################################################################
# run_plink_homozyg
# Input - ped
# Process - Makes a bed file, no filtering
# Output - bed
######################################################################
def run_plink_homozyg( plinkfile, force=forceFlag ) :
    write2log( " - Running:"+whoami(), True )
    targetdir, filename, suffix = getBasename(plinkfile)

    basename = "%s/%s" % (targetdir, filename)

    if suffix == ".ped" : pfile = "--file %s" % basename
    elif suffix == ".bed" : pfile = "--bfile %s" % basename
    elif suffix == ".tped" : pfile = "--tfile %s" % basename
    else : print "Error: Suffix unrecognized:",filename,suffix

    if os.path.exists( basename+".hom" ) and not force:
        write2log(" - File"+basename+".hom already exists. Skipping command!",True)
        return basename+".hom"

    command = "plink --noweb %s --homozyg --out %s" % ( pfile, basename )
    out = runCommand( command )
    print "run_plink_convert:",out
    return basename+".hom"
# End run_plink_homozyg

######################################################################
# run_plink_filter
# Input - ped/bed
# Process - Identify variants that fall within LD blocks
#    Filter LD vars, MAF >= .5, Missingness >= .05, add hwe...
# Output - bed
######################################################################
def run_plink_filter( targetfile, force=False ) :
    write2log( " - Running:"+whoami(), True )
    targetdir, filename, suffix = getBasename(targetfile)
    print targetdir, filename, suffix

    basename = "%s/%s" % (targetdir, filename )
    bedfile = basename+".filt.bed"
    if os.path.exists( bedfile ) and not force :
        write2log(" - File"+bedfile+" already exists. Skipping command!",True)
        return bedfile

    if suffix == ".ped" : pfile = "--file %s" % basename
    elif suffix == ".bed" : pfile = "--bfile %s" % basename
    elif suffix == ".tped" : pfile = "--tfile %s" % basename
    else : print "Error: Suffix unrecognized:",filename,suffix

    # a) consider a window of 50 SNPs, 
    # b) calculate LD between each pair of SNPs in the window, 
    # c) remove one of a pair of SNPs if the LD is greater than 0.5,
    # d) shift the window 5 SNPs forward and repeat the procedure
    # --1 
    command = "plink --noweb %s" \
            " --indep-pairwise 50 5 0.5" % ( pfile )
            #" --geno 0.05" \
            #" --maf 0.05" \
            #" %s" \
    print "Command:",command

    out = runCommand( command )
    print "run_plink_filter:",out

    # Mind - exclude individuals with too much missing genotype data
    # Geno - exclude variants not covered across patients
    # MAF - filter variants based on minor allele frequency
    command = "plink --noweb %s" \
        " --extract plink.prune.in --make-bed" \
        " --out %s" \
        % ( pfile, basename+".filt" )
        #" --maf 0.05" \
        #" --filter-cases "\
    print "Command:",command
    out = runCommand( command )
        #" --mind 0.1 "\
    print "run_plink_filter:",out

    return bedfile
# End run_plink_filter

######################################################################
# excludeSnps
# Goal - remove variants from dataset that fall 
# within the same centimorgan
################################################################################
################################################################################
def excludeSnps( pedfile, force=False ) :
    write2log(" - Running "+whoami(), True)
    filepath, basename, suffix = getBasename(pedfile)
    frqfile = filepath+"/"+basename+".frq"
    mapfile = filepath+"/"+basename+".map"
    excludefile = filepath+"/"+basename+".exclude.txt"
    pfile = plink_fileinput( pedfile )
    print "Fixing mapfile:",mapfile

    newfrq = filepath+"/"+basename+".exclude.frq"
    #newmap = filepath+"/"+basename+".exclude.map"
    thinped = filepath+"/"+basename+".exclude"

    if os.path.exists(thinped+".ped") and os.path.exists(newfrq) and not force :
        return newfrq, thinped

    print "Reading:",frqfile
    # VID - variant id
    freqdata = read_csv(frqfile,sep=" ",header=None,
                        names=["chrom","vid","ref","alt","freq","count"],
                        dtype={"vid":str})
    #freqdata["pos"] = [int(x.split(".")[1]) for x in freqdata.vid]
    zerofreq = freqdata.vid[freqdata.freq.isin([0,1])].tolist()
    print "Zerofreq:",zerofreq
    
    #evars = {}
    #for row in csv.reader(open(mapfile),delimiter=" ") :
    # http://www.nature.com/nature/journal/v409/n6822/fig_tab/409951a0_T1.html
    # calculate Centimorgans!!!!! 1.30 cM/Mb
    print "Reading mapfile:",mapfile
    mapdata = read_csv( mapfile, header=None, delim_whitespace=True,
                        names=["chrom","vid","centimorgans","pos"],
                        dtype={"vid":str})

    mapdata["centimorgans"] = ["%.1f" % (float(pos) / 1000000 * 1.3)
                               for pos in mapdata.pos]
    alldata = merge(mapdata, freqdata, on="vid" )
    print "mapdata:",mapdata.shape
    print "Alldata:",alldata.shape
    for i in range(1,1000) :
        if freqdata.ix[i,"vid"]  != mapdata.ix[i,"vid"] :
            print "Error",i,freqdata.ix[i,"vid"], mapdata.ix[i,"vid"]

    freqdata["pos"] = mapdata["pos"]
    maxfreqbycm = alldata[alldata.freq > 0.0].groupby("centimorgans")["freq"].idxmax()
    newfrqdata = freqdata.ix[maxfreqbycm,:].sort(["chrom","pos"])
    newfrqdata[["chrom","vid","ref","alt","freq","count"]].to_csv(
        newfrq,sep=" ",index=False,header=False)

    # Write out exclude variants
    newfrqdata["vid"].to_csv(excludefile, index=False,header=False)
    filetype = "--file" if pedfile[-3:] == "ped" else "--bfile"
    command = ["plink","--noweb",filetype,pedfile[:-4],"--allow-no-sex",
               "--extract",excludefile,"--recode","--out",thinped]
    print " ".join(command)
    try :
        out = subprocess.check_output( command )
    except subprocess.CalledProcessError as e :
        print "Error:",e
        sys.exit(1)

    print "Thin ped!:",thinped
    assert os.path.exists(thinped+".ped")
    return newfrq, thinped
# END excludeSnps

######################################################################
# make_pcaplots
######################################################################
def make_pcaplots( bedfile ) :
    write2log( " - Running:"+whoami(), True )

    targetdir, filename, suffix = getBasename(bedfile)
    prefix = filename
    if filename.find(".") > 0 : prefix = filename[:filename.find(".")]
    # MAKE RELATEDNESS PLOTS
    command = "./rscripts/runPCA.R %s %s" % bedfile, prefix
    runout = runCommand(command)
    print runout
    return
# End make_pcaplots

######################################################################
# run_plink_remove
# Input - ped/bed, samplelist
# Process - filter samples by keeplist
# Output - bed
######################################################################
def run_plink_remove( targetfile, keeplist, force=True ) :
    write2log( " - Running:"+whoami(), True )
    targetdir, filename, suffix = getBasename(targetfile)

    basename = "%s/%s" % (targetdir, filename)
    bedfile = basename+".uniq.bed"
    if os.path.exists( bedfile ) and not force:
        write2log(" - File"+bedfile+" already exists. Skipping command!",True)
        return bedfile

    if suffix == ".ped" :
        pfile = "--file %s" % basename
    elif suffix == ".bed" :
        pfile = "--bfile %s" % basename

    # --1
    command = "plink --noweb %s" \
            " --keep %s --make-bed --out %s.uniq" \
            % ( pfile, keeplist, basename )

    runCommand( command )
    return bedfile
# End run_plink_remove

######################################################################
# run_plink_het
# Input - ped/bed
# Process - LD filtering, make het file
# Output - het (inbreeding coefficients)
######################################################################
def run_plink_het( targetfile, currK=None, groupnum=None, keeplistfile=None,force=forceFlag ) :
    write2log( " - Running:"+whoami(), True )
    print "Targetfile:",targetfile,"Group:",groupnum,"keeplist:",keeplistfile
    targetdir, filename, suffix = getBasename(targetfile)

    basename = "%s/%s" % (targetdir, filename)
    outfile = basename
    if currK is not None :
        outfile = outfile+"."+str(currK)
    if groupnum is not None :
        outfile = outfile+".grp"+str(groupnum)

    print outfile
    if os.path.exists( outfile+".het" ) and not force:
        write2log(" - File"+outfile+" already exists. Skipping command!",True)
        return outfile+".het"

    if suffix == ".ped" :
        pfile = "--file %s" % basename
    elif suffix == ".bed" :
        pfile = "--bfile %s" % basename

    otheroptions = ""
    if keeplistfile is not None :
        otheroptions += ' --keep %s ' % keeplistfile
    #if keeplist is not None :
        #for patient in keeplist:
            #if len(patient) == 0 or patient[0] == "#" :
                #continue
            #otheroptions += ' --indv "%s"' % patient.strip()

    # --1
    command = "plink --noweb %s" \
            " --indep-pairwise 50 5 0.5" \
            " --out %s --het %s" \
            % ( pfile, outfile, otheroptions )

    print "Command:",command
    runCommand( command )
    return outfile+".het"
# End run_plink_het

######################################################################
# modify_ped_old
######################################################################
def modify_ped_old( pedfile, patientdata ) :
    write2log( " - Running:"+whoami(), True )
    assert( pedfile[-3:] == "ped" )
    targetdir, filename, suffix = getBasename(pedfile)
    pedmodoutfile = "%s/%s.mod.ped" % (targetdir,filename)
    mapmodoutfile = "%s/%s.mod.map" % (targetdir,filename)

    if os.path.exists( pedmodoutfile ) and forceFlag is False:
        write2log(" - File"+pedmodoutfile+" already exists. Skipping command!",True)
        return pedmodoutfile

    FILE = open( pedfile, "r" )
    reader = csv.reader( FILE, delimiter="\t" )
    OUT = open( pedmodoutfile, "wb" )
    #for row in reader :
    patlist = []
    for line in FILE :
        row = line.strip().split("\t",8)
        patid = re.sub( r"[\/\|]","-",row[0] )
        cleanpatid = parse_patid( patid )
        #patid = string.replace( r"_wex1","", patid )
        if patientdata.has_key(patid) :
            currpatdata = patientdata[patid][:6]
            if currpatdata[1] in patlist :
                print "Duplicate:",currpatdata[1]
                continue
            patlist.append(currpatdata[1])
            OUT.write( "%s\t%s\n" % ("\t".join(currpatdata),"\t".join(row[6:])) )
            #print "Prev",row[:6]
            #print "New",currpatdata
        elif patientdata.has_key(cleanpatid) :
            #if cleanpatid in patlist :
                #print "Duplicate:", cleanpatid
                #continue
            currpatdata = patientdata[cleanpatid][:6]
            if currpatdata[1] in patlist :
                print "Duplicate:",currpatdata[1]
                continue
            patlist.append(currpatdata[1])
            OUT.write( "%s\t%s\n" % ("\t".join(currpatdata),"\t".join(row[6:])) )
            #print "Prev",row[:6]
            #print "New",currpatdata
        #else :
            #row[0] = re.sub( "[\/\|]|_wex1","-",row[0] )
            #row[1] = re.sub( "[\/\|]|_wex1","-",row[1] )
            #OUT.write( "%s\n" % "\t".join(row) )
    #print patlist
    FILE.close()
    OUT.close()

    FILE = open("%s/%s.map" %(targetdir,filename),"r" ) # open map file
    OUT = open( mapmodoutfile, "wb" )
    reader = csv.reader( FILE, delimiter="\t" )
    for row in reader :
        if len(row) < 0 :
            continue
        if row[1][:2] != "rs" :
            row[1] = row[0]+"."+row[3]
        OUT.write( "\t".join(row) +"\n" )
    OUT.close()
    FILE.close()
    return pedmodoutfile
# End modify_ped_old

######################################################################
# modify_ped
# Input - ped, patientdata
# Process - Fix sample names, modify var ids, calculate cMs 
# Output - ped,map
######################################################################
def modify_ped( pedfile, patientdata, shorten=False, calcdistances=False, force=forceFlag ) :
    write2log( " - Running:"+whoami(), True )
    assert( pedfile[-4:] == ".ped" )
    targetdir, filename, suffix = getBasename(pedfile)
    pedmodoutfile = "%s/%s.mod.ped" % (targetdir,filename)
    mapmodoutfile = "%s/%s.mod.map" % (targetdir,filename)

    if os.path.exists( pedmodoutfile ) and os.path.exists(mapmodoutfile) and not force :
        write2log(" - File"+pedmodoutfile+" already exists. Skipping command!",True)
        return pedmodoutfile

    FILE = open( pedfile, "r" )
    reader = csv.reader( FILE, delimiter="\t" )
    OUT = open( pedmodoutfile, "wb" )
    #for row in reader :

    print patientdata.head(4)
    patientdata.index = patientdata["Individual.ID"].tolist()
    patlist = []
    for line in FILE :
        row = re.split("[ \t]", line.strip(), 8 )
        #gcalls = [x.replace('0','1') for x in row[6:]]
        gcalls = row[6:]
        patid = re.sub( r"[\/\|]","-",row[0] )
        cleanpatid = parse_patid( patid )
        #patid = string.replace( r"_wex1","", patid )
        #print "Patid:",cleanpatid
        if patid in  patientdata["Individual.ID"].tolist():
            currpatdata = patientdata.loc[patid][["Family.ID","Individual.ID","Paternal.ID","Maternal.ID","Gender","Phenotype"]].map(str).tolist()
            if shorten and len(currpatdata[1]) > 21 and currpatdata[1].find("_") >= 0: 
                print "Shorten!!"
                currpatdata[1] = currpatdata[1][:currpatdata[1].find('_')]
            #print "Currpat data:",currpatdata
            if currpatdata[1] in patlist :
                print "Duplicate:",currpatdata[1]
                continue
            patlist.append(currpatdata[1])
            #OUT.write( "%s\t%s\n" % ("\t".join(currpatdata),"\t".join(row[6:])) )
            #print "Gcalls:",gcalls
            OUT.write( "%s\t%s\n" % ("\t".join(currpatdata),"\t".join(gcalls)) )
        elif cleanpatid in  patientdata["Individual.ID"].tolist():
            #if cleanpatid in patlist :
                #print "Duplicate:", cleanpatid
                #continue
            #currpatdata = patientdata[cleanpatid][:6]
            currpatdata = patientdata.loc[cleanpatid][["Family.ID","Individual.ID","Paternal.ID","Maternal.ID","Gender","Phenotype"]].tolist()
            if len(currpatdata[1]) > 21 and currpatdata[1].find("_") >= 0: 
                print "Shorten!!"
                currpatdata[1] = currpatdata[1][:currpatdata[1].find('_')]
            print "Currpat data:",currpatdata
            if currpatdata[1] in patlist :
                print "Duplicate:",currpatdata[1]
                continue
            patlist.append(currpatdata[1])
            #print "New",currpatdata #,[x[:8] for x in row]
            #print "Gcalls:",gcalls
            #OUT.write( "%s\t%s\n" % ("\t".join(currpatdata)," ".join(row[6:])) )
            OUT.write( "%s\t%s\n" % ("\t".join(currpatdata),"\t".join(gcalls)) )
            #print "Prev",row[:6],"New",currpatdata
        else :
            OUT.write( "%s\n" % ("\t".join(row)) )
    FILE.close()
    OUT.close()

    #print "#"*70

    FILE = open("%s/%s.map" %(targetdir,filename),"r" ) # open map file
    print "Making file:",mapmodoutfile
    OUT = open( mapmodoutfile, "wb" )
    #reader = csv.reader( FILE, delimiter=" " )
    for line in FILE :
        if line.find('\t') >= 0 : row = line.rstrip().split('\t')
        else: row = line.rstrip().split(' ')
        if len(row) < 0 :
            continue
        if calcdistances and int(row[2]) == 0 :
            row[2] = str(round(float(row[3]) / 1000000 * 1.3,1))
            #print "New Dist", row[2], row[3]
        if row[1][:2] != "rs" :
            row[1] = row[0]+"."+row[3]
        #OUT.write( "\t".join(row) +"\n" )
        OUT.write( " ".join(row) +"\n" )
    OUT.close()
    FILE.close()
    return pedmodoutfile
# End modify_ped

######################################################################
# modify_tped
# Input - tped, tfam, patientdata
# Process - Fix sample names, modify var ids, calculate cMs 
# Output - tped,tfam
######################################################################
def modify_tped( tpedfile, patientdata, force=forceFlag ) :
    write2log( " - Running:"+whoami(), True )
    print tpedfile
    assert( tpedfile[-4:] == "tped" )
    targetdir, filename, suffix = getBasename(tpedfile)

    targetcolumns = ["Family.ID","Individual.ID","Paternal.ID","Maternal.ID","Gender","Phenotype"]
    infocp = patientdata[targetcolumns].copy()
    infocp.index = infocp["Individual.ID"].tolist()
    infocp[["Gender","Phenotype"]] = infocp[["Gender","Phenotype"]].astype(int).astype(str)

    tfamfile = "%s/%s.tfam" %(targetdir,filename)
    famdata = read_csv(tfamfile,sep="\t",header=None, names=targetcolumns)
    famdata.index = famdata["Individual.ID"].tolist()

    famdata.update(infocp)
    famdata.to_csv(tfamfile,sep="\t",index=False,header=False)
# End modify_tped
    #tfamcopy = "%s/%s.1.tfam" %(targetdir,filename)
    #runCommand( "cp %s %s" % (tfamfile,tfamcopy) )

######################################################################
# determinegroups
######################################################################
def determinegroups( patientlistfile, patientdata ):
    write2log( " - Running:"+whoami(), True )
    FILE = open(patientlistfile,"r")
    reader = csv.reader(FILE,delimiter="\t")
    groups = {}
    for row in reader :
        if patientdata.has_key(row[1]) :
            #print patientdata[row[1]]
            if not groups.has_key(patientdata[row[1]][6]):
                groups[patientdata[row[1]][6]] = 0
            groups[patientdata[row[1]][6]] += 1
    FILE.close()
    #print groups
    return groups
# End determinegroups

######################################################################
# run_admixture
# cv - cross validation
# B - bootstrapping
######################################################################
def run_admixture( targetfile, K, bootstrap=False, forceFlag=False ) :
    write2log( " - Running:"+whoami(), True )
    targetdir, filename, suffix = getBasename(targetfile)

    #if bootstrap :
        #Qfile = targetdir+"/"+filename+"."+str(K)+".best.Q"
        #logfile = "%s/%s_admix.%d.best.log" % (LOGDIR,filename,K)
    #else :
    Qfile = targetdir+"/"+filename+"."+str(K)+".Q"
    logfile = "%s/%s_admix.%d.log" % (LOGDIR,filename,K)
    print "Qfile:",Qfile
    if os.path.exists( Qfile ) and forceFlag is False:
        write2log(" - File"+Qfile+" already exists. Skipping command!",False)
        with open(logfile, 'r') as f:
            runout = f.read()
        return runout

    if bootstrap :
        command = "admixture -B %s %s" % ( targetfile, K )
    else :
        command = "admixture --cv %s %s" % ( targetfile, K )

    runout = runCommand( command )
    OUT = open(logfile,"wb")
    OUT.write(runout)
    OUT.close()
    command = "mv %s*.Q %s ; mv %s*.P %s" % (filename, targetdir, filename, targetdir)
    runCommand( command )
    return runout
# End run_admixture

######################################################################
# findBestK
######################################################################
def findBestK( targetfile, maxK ) :
    write2log( " - Running:"+whoami(), True )
    targetdir, filename, suffix = getBasename(targetfile)
    logfile = "%s/%s_findBestK.%d.log" % (LOGDIR,filename,maxK)
    bestk = 0
    lowesterr = 1.0
    errorrates = {}
    for K in range(2,maxK+1) :
        runout = run_admixture( targetfile, K )
        cverrorsearch = re.search( "CV error \WK=\d+\W+([\.\d]+)", runout )
        cverror = 0
        if cverrorsearch is not None :
            cverror = cverrorsearch.group(1)
        print "K:",K,"Error:",cverror
        errorrates[K] = cverror
        if float(cverror) < lowesterr :
            lowesterr = float(cverror)
            bestk = K
    errlogfile = "%s/%s_error.log" % (LOGDIR, filename)
    OUT = open(errlogfile, "wb")
    OUT.write("K\tError\n")
    for K in sorted(errorrates) :
        OUT.write( "%s\t%s\n" % (K, errorrates[K]) )
    OUT.close()
    # MAKE ERROR RATE GRAPH HERE
    command = "./errorrate.R %s" % errlogfile
    runout = runCommand(command)
    OUT = open(logfile,"wb")
    OUT.write(runout)
    OUT.close()
    return bestk
# End findBestK

######################################################################
# annotatePatients
######################################################################
def annotatePatients( groups, targetfile, patientdata, bestK ):
    write2log( " - Running:"+whoami(), True )
    targetdir, filename, suffix = getBasename(targetfile)
    Qfile = targetdir +"/"+ filename+"."+str(bestK)+".Q"

    groupfile = LOGDIR+"/"+filename+"_patientgroups."+str(bestK)+".txt"
    OUT = open(groupfile, "wb")

    #FILE = open(uniquelist,"r")
    #reader = csv.reader(FILE,delimiter="\t")
    #groups = []
    #for row in reader :
        #groups.append( row[1] )
    #FILE.close()

    FILE = open(Qfile, "r")
    reader = csv.reader(FILE, delimiter=" ")
    patientgroups = {}
    count = 0
    separation = {}
    #print "Groups:",groups
    OUT.write("Index\tProb\tGroup\tPatient\tEthnicity\n")
    for row in reader :
        if len(row) == 0 : continue
        maxval = 0

        for i in range(len(row)) :
            if float(row[i]) > maxval :
                patgroup = i
                maxval = float(row[i])
        ethnicity = "unknown"
        #print "Count",count, "Max",len(groups),"Row:",row
        #print "curr group",groups[count]
        if patientdata.has_key(groups[count]) :
            ethnicity = patientdata[groups[count]][6]
        OUT.write("%d\t%f\t%d\t%s\t%s\n" % (count, max([float(x) for x in row]),patgroup+1, groups[count],ethnicity) )
        patientgroups[groups[count]] = patgroup+1
        if not separation.has_key(ethnicity) :
            separation[ethnicity] = [0 for x in range(len(row))]
        separation[ethnicity][patgroup] += 1
        count += 1
    OUT.write("#Group"+"\t"+"\t".join("groups"+str(x) for x in range(1,len(row)+1))+"\n")
    for eth in separation :
        OUT.write( eth+"\t"+"\t".join([str(x) for x in separation[eth]]) +"\n")
    return patientgroups
# annotatePatients

######################################################################
# runAdmixtureClustering
######################################################################
def runAdmixtureClustering( bedfile_uniq,annotfile,currK,forceFlag=forceFlag ):
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
#def plotInbreedingCoefficient( bedfile_uniq) :
    #write2log( " - Running:"+whoami(), True )
    #print "Target File:",bedfile_uniq
    #targetdir, filename, suffix = getBasename(bedfile_uniq)
    #command = "./inbreedCoeff.R %s %s" % ( targetdir, filename )
    #print "Command",command
    #output = runCommand( command )
    #return
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
        pedfile = run_vcftools_plink( vcffile, outdir )
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
    write2log( " - Running:"+whoami(), True )
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
    plotInbreedingCoefficient( bedfile )
# inbreedingCoefficient

######################################################################

#-------------------------------------------------------------------#
#                            Commands                               #
#-------------------------------------------------------------------#

######################################################################
# printCurrentTime
######################################################################
#def printCurrentTime() :
    #now = datetime.datetime.now()
    #print "Current date and time using strftime:"
    #print now.strftime("%Y-%m-%d %H:%M")
# END printCurrentTime()

######################################################################
# getPats1
######################################################################
def getPats1( vcffile, keeplist=None ):
    #pfile = "./merged/combined.ped"
    pfile = "./patientannot/annotation/patientannotation_laser.txt"

    targetdir, filename, suffix = getBasename(vcffile)
    allpats = currentFilePats( vcffile )

    print allpats
    removelist = []
    removefile = targetdir+"/toremove.txt"
    if os.path.exists( removefile ):
        print "Removelist found:", removefile
        print "Remove File:", removefile
        f = open( removefile, "r" )
        for patient in f.readlines():
            if len(patient) == 0 or patient[0] == "#" :
                continue
            removelist.append(patient.strip())
        f.close()

    fulllist = []
    for pat in allpats :
        if pat not in removelist :
            fulllist.append(pat)

    subfulllist = []
    for pat in fulllist :
        newpat = "-".join(pat.split("-")[1:])
        subfulllist.append(newpat)

    # If we have a list...
    keeppats = []
    if keeplist is not None :
        print "Using keeplist:", keeplist
        keeppats = [x.strip() for x in open(keeplist).readlines()]
        #print len(keeppats)
        keeppats = [x for x in keeppats if x in allpats]
        #print len(keeppats)
        return keeppats

    hd = {}
    for row in csv.reader(open(pfile),delimiter="\t") :
        if len(hd) == 0 :
            for i in range(len(row)) :
                hd[row[i]] = i
            print hd
            continue
        #print row[hd["Lab"]]
        if (row[hd["Lab"]] == "Broad" and \
            row[hd["Source"]] not in ["Pilot"]) :
            print "Broad"
            #["PlateI", "PlateII","Pilot"]) :
            #and row[hd["continent"]] in ["Middle East", "South Asia"]) \
            #or (row[hd["Lab"]] == "Daly" and row[hd["continent"]] == "Europe") \
            #or (row[hd["Lab"]] == "Eichler" and row[hd["continent"]] == "Europe"):
            #patname = row[hd["Individual.ID"]]
            patname = row[hd["original"]]
            if fulllist is not None :
                fixpatname1 = "-".join(patname.split("-")[1:])
                if patname in fulllist :
                    keeppats.append( patname )
                elif fixpatname1 in subfulllist :
                    for i in range(len(subfulllist)) :
                        if fixpatname1 == subfulllist[i] :
                            break
                    keeppats.append( fulllist[i] )
                #else :
                    #print "Error: patname not in fulllist:", patname
            else :
                print patname
                keeppats.append( patname )
    print len(keeppats), len(fulllist)
    #for pat in fulllist :
        #if pat not in keeppats :
            #print pat
    return keeppats
# END getPats1

######################################################################
# currentFilePats1
######################################################################
def currentFilePats1( cfile ):
    targetdir, filename, suffix = getBasename(cfile)

    patlist = []
    if suffix == ".ped" :
        command = "cut -f 2 %s" % cfile
        #print command
        out = runCommand( command ).strip()
        patlist = out.split("\n")

    elif suffix == ".txt" :
        command = "cut -f 2 %s" % cfile
        #print command
        out = runCommand( command ).strip()
        patlist = out.split("\n")
    elif suffix == ".fam" :
        command = '''awk '{print $2}' %s''' % cfile
        #print command
        out = runCommand( command ).strip()
        patlist = out.split("\n")
    elif cfile[-3:] == ".gz" and suffix == ".vcf" :
        headline = None
        #print "Opening:",cfile
        FH = gzip.open( cfile )
        while( headline is None ):
            line = FH.readline()
            #print line[:10]
            if line[:6] == "#CHROM" :
                headline = line.strip()
        FH.close()
        patlist = headline.split("\t")[9:]
    elif suffix == ".vcf" :
        command = ''' head -200 %s | awk '$1 == "#CHROM"' ''' % cfile
        #print command
        out = runCommand( command ).strip()
        #print out[:40]
        patlist = out.split("\t")[9:]
    else :
        print "Error: Unknown filetype:",suffix
        print "filebase:",filename
        sys.exit(1)
    return patlist
# END currentFilePats1


######################################################################
# mergeVcfs
######################################################################
def mergeVcfs( path, filelist, targetfile ) :
    write2log( " - Running:"+whoami(), True )
    print "filelist:",filelist

    if os.path.exists( targetfile ) and forceFlag is False:
        print " - File"+targetfile+" already exists. Skipping command!"
        return

    for eachfile in filelist :
        assert eachfile[-3:] == ".gz", "Error: file isnt bgzipped"
    #for vfile in allfiles :
        #newfile = compress(vfile)

    assert( existsInPath("vcf-merge") )
    command = "vcf-merge --ref-for-missing '0/0' -d %s > %s" \
             % (" ".join( filelist ), targetfile)
    print command
    runCommand( command, True )
    return
# END mergeVcfs

######################################################################
# mergeAllVcfs
######################################################################
def mergeAllVcfs( path, filename, additionalfile=None ) :
    print "Running mergeAllVcfs. Path:",path
    allfiles = filesInDir( path, "recode.vcf" )
    if additionalfile is not None  and additionalfile[-3:] == "vcf" :
        allfiles = allfiles + [additionalfile]

    if os.path.exists( filename ) and forceFlag is False:
        print " - File"+filename+" already exists. Skipping command!"
        return

    #for vfile in allfiles :
        #newfile = compress(vfile)

    allfiles = filesInDir( path, "recode.vcf.gz" )
    if additionalfile is not None  and additionalfile[-3:] == ".gz" :
        allfiles = allfiles + [additionalfile]

    print "Allfiles:",allfiles
    assert( existsInPath("vcf-merge") )
    command = "vcf-merge %s > %s" \
             % (" ".join( allfiles ), filename)
    print command
    runCommand( command, True )
    return
# END mergeAllVcfs
    #command = "bgzip %s; tabix -p vcf %s.gz" % (outfile, outfile )
    #for i in range(len(allfiles)) :
        #if allfiles[i].find("merged") > -1 :
            #del allfiles[i]
        #command = "bgzip %s; tabix -p vcf %s.gz" % (vfile,vfile)
        #runCommand( command, True )

######################################################################
# compress
######################################################################
def compress( file, force=False ) :
    print " - Running compress"
    if file[-3:] == ".gz" :
        return file

    outfile = file+".gz"
    if os.path.exists( outfile ) and not force :
        print " - File"+outfile+" already exists. Skipping command!"
        return outfile
    elif os.path.exists( outfile ) and force :
        runCommand( "rm %s" % outfile )
    command = "bgzip %s; tabix -p vcf %s" % (file,outfile)
    runCommand( command, True )
    print " - End compress"
    return outfile
# END compress

######################################################################
# parseOutlierFile
######################################################################
def parseOutlierFile( outlierfile ) :
    alldata = {}
    head = []
    for row in csv.reader(open(outlierfile), delimiter="\t") :
        if len(head) == 0 :
            head = row
            continue
        pat,val,att = row
        if not alldata.has_key(pat) : alldata[pat] = []
        alldata[pat].append(att)

    return [x for x in alldata if len(alldata[x]) > 1]
#END parseOutlierFile

######################################################################
# pseq_istats
######################################################################
def pseq_istats( vcffile, force=forceFlag ) :
    print "VCF:",vcffile

    vtargetdir, vfilename, vsuffix = getBasename(vcffile)
    if len(vtargetdir) > 0 :
        pseqfile = "%s/%s_istats.txt" % (vtargetdir,vfilename)
        outlierfile = "%s/%s_istats_outliers.txt" % (vtargetdir,vfilename)
    else :
        pseqfile = "%s_istats.txt" % (vfilename)
        outlierfile = "%s_istats_outliers.txt" % (vfilename)

    if os.path.exists( pseqfile ) and not force :
        print " - File"+pseqfile+" already exists. Skipping command!"
    else :
        print "Writing to file:",pseqfile
        command = "pseq %s i-stats > %s" % (vcffile,pseqfile)
        out = runCommand( command, True )

    fixedbase = vfilename.replace(".","_")
    outdir = os.path.abspath(vtargetdir+"/"+fixedbase+"_pseq")
    print "Outdir:",outdir
    makeDir(outdir)
    command = "rscripts/pseqgraphs.R %s %s" % (pseqfile, outdir)
    write2log(command)
    out = runCommand(command)
    write2log(out)
    command = "rscripts/pseqgraphs2.R %s %s" % (pseqfile, outdir)
    write2log(command)
    out = runCommand(command)
    write2log(out)

    outliers = None
    if fileIsOk(outlierfile) :
        outliers = parseOutlierFile( outlierfile )
    return pseqfile, outliers
# End pseq_istats

######################################################################
# intersectBED
######################################################################
def intersectBED( vcffile, bedfile ) :
    print "VCF:",vcffile, "BED:",bedfile
    vtargetdir, vfilename, vsuffix = getBasename(vcffile)
    btargetdir, bfilename, bsuffix = getBasename(bedfile)
    #headfile = "%s/%s.head" % (btargetdir,vfilename)
    #outfile = "%s/%s.%s%s" % (vtargetdir,vfilename,bfilename,vsuffix)
    outfile = "%s/%s.regions%s" % (vtargetdir,vfilename,vsuffix)

    if os.path.exists( outfile ) and forceFlag is False:
        print " - File"+outfile+" already exists. Skipping command!"
        return outfile

    print "Writing to file:",outfile
    command = "gunzip -c %s | head -300 | awk '$1 ~ /#/' > %s" % (vcffile,outfile)
    runCommand( command, True )
    command = "intersectBed -wa -a %s -b %s >> %s" % (vcffile,bedfile,outfile)
    runCommand( command, True )

    cfile = compress( outfile )

    return cfile
# End intersectBED

######################################################################
# Main
######################################################################
if __name__ == "__main__":

    optlist, args = getopt.getopt( sys.argv[1:], "c:tp:")
    optlist = dict(optlist)
#
    if len(optlist) == 0:
        print "Error: no options given: -c"
        print " -c functionname: run check"
        print " -t run a test"
        sys.exit(1)

    os.chdir("..")

    vcffile = "testdata/testdata.vcf"
    gzvcffile = "testdata/testdata.vcf.gz"
    pedfile = "testdata/testdata.plink.ped"
    famfile = "testdata/testdata.fam"
    txtfile = "testdata/testdata.unrelated.txt"

    outdir = "testdata"
    functionname = optlist.get("-c",None)

    DEBUG = True

    if optlist.has_key("-t") : 
        annot = parseAnnotation()

    elif optlist.has_key("-p") : 
        vcffile = optlist.get("-p",None)
        if vcffile is None or not os.path.exists(vcffile) : print "Error:",vcffile; sys.exit(1)
        pseq_istats( vcffile )
    elif functionname == "currentFilePats" :
        print "Vcf file"
        plist = currentFilePats(vcffile)
        print plist[:20]

        print "GZIP VCF file"
        plist = currentFilePats(gzvcffile)
        print plist[:20]

        print "PED file"
        plist = currentFilePats(pedfile)
        print plist[:20]

        print "Fam file"
        plist = currentFilePats(famfile)
        print plist[:20]

        print "Txt file"
        plist = currentFilePats(txtfile)
        print plist[:20]
    elif functionname == "run_vcftools" :
        print run_vcftools( gzvcffile )
        print run_vcftools_plink( gzvcffile, outdir )
    elif functionname == "modify_tped" :
        print "Got here!"
        patientannotation = "./resources/annotation/patientannotation.ped"
        patientdata = read_csv( patientannotation, sep="\t" )
        tpedfile = "/scratch/escott/workspace/variome/rawdata/variome1/variome.chimp.regions.filt.samp.plink.tped"
        modify_tped( tpedfile, patientdata,  force=True)
        #shorten=False, calcdistances=False,
# END Main

