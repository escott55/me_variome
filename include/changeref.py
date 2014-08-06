#!/usr/bin/python
# coding: utf-8

import os, sys
import csv
import commands
import math
import re
import pybedtools
import gzip

from housekeeping import *

forceFlag=False
BUF_SIZE=30*1048576

######################################################################
# liftOver
######################################################################
def liftOver( targetfile, outfile=None ) :
    assert( existsInPath("liftOver") )

    #chainfile = "liftover/hg19ToHg18.over.chain"
    chainfile = "/home/escott/resources/liftover/hg19ToPanTro4.over.chain"

    filepath, basename, suffix = getBasename( targetfile )
    if outfile is None : outfile = "%s.out" % targetfile
    #unmapped = "unmapped.out"
    unmapped = "%s/%s_unmapped.bed" % (filepath, basename)
    command = "liftOver %s %s %s %s" % (targetfile, chainfile, outfile, unmapped)
    runCommand( command, True )
    return outfile
# END liftOver

######################################################################
# getSequence
######################################################################
def getSequence( chrom, start, end, tref="human" ) :
    assert( existsInPath("samtools") )
    assert( existsInPath("grep") )
    if chrom[:3] != "chr" :
        chrom = "chr"+chrom
    if chrom == "chr23" :
        chrom = "chrX"
    elif chrom == "chr24" :
        chrom = "chrY"

    if tref == "human" :
        ref = "/home/escott/resources/hg19/hg19.fa"
    elif tref == "chimp" :
        ref = "/home/escott/resources/chimp/panTro4.fa"
    command = 'samtools faidx \
        %s %s:%s-%s |\
        grep -v ">" | tr -d "\n"' % (ref, chrom, start, end)
    runout = runCommand( command, False )
    return runout
    #command = 'samtools faidx \
        #/home/Gleeson/escott/resources/hg19all/hg19.fa %s:%s-%s |\
        #grep -v ">" | tr -d "\n" | tr "atgcn" "ATGCN"' % (chrom, start, end)
# END getSequence

######################################################################
# getSeqHead
######################################################################
def getSeqHead( row ) :
    #seqhead = "fam:%s|chr%s:%s|%s/%s|%s|%s|%s" % (row[0],row[1],row[2],row[4],
            #row[5], row[7], row[10],row[12])
    tmp = re.match( r"\w+", row[7] )
    genename = "Unk"
    if tmp is not None :
        genename = tmp.group(0)
    seqhead = "%s_c%sp%s_%s" % (row[0],row[1],row[2],genename)
    return seqhead
#END getSeqHead

#G/C→S, A/T→W, G/A→R, T/C→Y, G/T→K, A/C→M
def iupac( ref, mut ) :
    assert( ref != mut )
    if len(ref)  > 1 or len(mut) > 1 :
        #print "Indel!"
        return "X"
    if ref in ["G","C"] and mut in ["G","C"] :
        return "S"
    elif ref in ["A","T"] and mut in ["A","T"] :
        return "W"
    elif ref in ["G","A"] and mut in ["G","A"] :
        return "R"
    elif ref in ["T","C"] and mut in ["T","C"] :
        return "Y"
    elif ref in ["G","T"] and mut in ["G","T"] :
        return "K"
    elif ref in ["A","C"] and mut in ["A","C"] :
        return "M"
    else :
        print "Error: ref/mut",ref,"/",mut
        sys.exit()
    return
# END iupac

######################################################################
# parseRegionsFile
######################################################################
def parseRegionsFile():
    refseqgenesbed = "./bedfiles/refseq_genes_simple.bed"
    regions = {}
    for row in csv.reader(open(refseqgenesbed), delimiter="\t") :
        if not regions.has_key(row[0]) : regions[row[0]] = []
        regions[row[0]].append([int(row[1]),int(row[2]),row[3]])
        #print [int(row[1]),int(row[2]),row[3]]
    return regions
# End parseRegionsFile

######################################################################
# intersectWithGenes
######################################################################
def intersectWithGenes( chrom, pos, regionsfile ):
    #refseqgenesbed = "./bedfiles/refseq_fullgenes.bed"
    #a = pybedtools.BedTool( refseqgenesbed )
    #a.intersect(bedfile).saveas("test.bed")
    pos = int(pos)
    genes = []
    if len(chrom) < 4 : chrom = "chr"+chrom
    if regionsfile.has_key(chrom) : 
        for region in regionsfile[chrom] :
            if region[0] < pos and region[1] > pos :
                genes.append(region[2])
    return ";".join(genes)
# End intersectGenes

######################################################################
# findPops
######################################################################
def findPops( samples ) :
    #annotfile = "/home/escott/consang/processeddata/annotation/patientannotation_laser.txt"
    annotfile = "./patientannot/annotation/patientannotation_laser.txt"
    head = {}
    popassoc = {}
    for row in csv.reader(open(annotfile), delimiter="\t") :
        if len(head) == 0 :
            for i in range(len(row)) :
                head[row[i]] = i
            continue
        popassoc[ row[head["Individual.ID"]] ] = row[head["ethnicity"]]
        popassoc[ row[head["original"]] ] = row[head["ethnicity"]]
    pops = []
    for sam in samples :
        if popassoc.has_key( sam ) :
            #print sam, popassoc[sam]
            pops.append( popassoc[sam] )
        else :
            #print "Error! no ethnicity found"
            #print sam
            pops.append( "Unk" )
    return pops
# END findPops

######################################################################
# calculatePopFreqs
######################################################################
def calculatePopFreqs( calls, pops ) :
    assert len(pops) == len(calls), "Error: different lengths"
    reffreqs = {}
    altfreqs = {}

    for p in pops :
        reffreqs[p] = 0
        altfreqs[p] = 0

    for i in range(len(calls)) :
        mycall = calls[i][:3]
        #print mycall
        if mycall == "1/1" :
            altfreqs[ pops[i] ] += 2
        elif mycall == "0/1" or mycall == "1/0" :
            altfreqs[ pops[i] ] += 1
            reffreqs[ pops[i] ] += 1
        elif mycall == "0/0" :
            reffreqs[ pops[i] ] += 2
        #else : 
            #print "Error: missing"
    return [reffreqs[x] for x in sorted(reffreqs)], [altfreqs[x] for x in sorted(altfreqs)]
# End calculatePopFreqs

######################################################################
# flipatts
######################################################################
def flipatts( atts ) :
    alldata = atts.split(":")
    for i in range(len(alldata)) :
        catt = alldata[i]
        if catt.find(",") >= 0 :
            subatts = catt.split(",")
            subatts.reverse()
            alldata[i] = ",".join(subatts)
    #print "first:", atts, "new",":".join(alldata)
    return ":".join(alldata)
# END flipatts

################################################################################
# vcf2Bed
################################################################################
def vcf2Bed( datalines, targetbed ):
    OUT = open( targetbed, "wb" )
    alldata = {}
    for line in datalines :
        if line[0] == "#" : continue
        row = line.rstrip().split("\t")
        start = int(row[1])
        ref = row[3]
        mut = row[4]
        end = start + max(len(ref), len(mut))
        key = "chr%s:%s:%s:%s" % (row[0], row[1], ref, mut )
        OUT.write( "chr%s\t%d\t%d\t%s\n" % ( row[0], start, end, key ))
        alldata[key] = row
    OUT.close()
    return alldata
# END vcf2Bed

################################################################################
# processVCFLines
################################################################################
def processVCFLines( genodata, humanbed, chimpbed, fixedvcf, refname="chimp" ) :

    # Get Chimp positions for all variants
    liftOver( humanbed, chimpbed )

    write2log( "Writing to:"+fixedvcf, True )
    #if group == 0 :
        #OUT = open(fixedvcf, "wb")
    OUT = open(fixedvcf, "ab")

    limit = -1
    for line in open(chimpbed):
        chrom, start, end, key = line.rstrip().split("\t")
        humchrom, humpos, humref, mut =  key.split(":")
        seq = getSequence( chrom, start, int(end)-1, refname ) # chimp
        chimpref = seq.upper()
        if len(chimpref) > 10 :
            #print "chrom",chrom,"start", start, "end", end
            write2log( "Err: hum"+humref[:10]+"chimp"+chimpref[:10].strip()+" len:"+str(len(chimpref)), True )
            #print chimpref.strip()
        else :
            #print ":",humref, chimpref, mut,":"
            newdata = genodata[key]
            if humref == chimpref :
                OUT.write( "\t".join(newdata) +"\n" )
            elif chimpref == mut :
                newdata[3] = chimpref
                newdata[4] = humref
                newdata[7] = "Chimpref"
                lineformat = newdata[8].split(":")
                for i in range(9,len(newdata)) :
                    if newdata[i].find(":") > -1 : 
                        geno, atts = newdata[i].split(":",1)
                    else :
                        geno  = newdata[i]
                        atts = ":".join(["." for x in lineformat])
                        #print geno+":"+atts
                    newatts = flipatts( atts )
                    #print geno
                    if geno == "0/0" :
                        newdata[i] = "1/1:"+newatts
                    elif geno == "1/0" :
                        newdata[i] = "0/1:"+newatts
                    elif geno == "0/1" :
                        newdata[i] = "1/0:"+newatts
                    elif geno == "1/1" :
                        newdata[i] = "0/0:"+newatts
                OUT.write( "\t".join(newdata) +"\n" )
            else : 
                write2log( "Error: multi refs %s %s %s" % (humref, chimpref, mut ))

        limit -= 1
        if limit == 0 : break
    OUT.close()
    return
# END processVCFLines

######################################################################
# forceRef
######################################################################
def forceRef( targetfile, refname="chimp" ):
    write2log( "#"*60, True )
    write2log( "Running forceRef:",targetfile, True )
    filepath, basename, suffix = getBasename( targetfile )
    #regionsdata = parseRegionsFile()
    fixedvcf = "%s/%s.%s%s" % ( filepath, basename, refname, suffix )
    zippedvcf = "%s/%s.%s%s.gz" % ( filepath, basename, refname, suffix )
    humanbed = "%s/human_%s.bed" % ( filepath, basename )
    chimpbed = "%s/%s_%s.bed" % ( filepath, refname, basename )
    #humanbed = "rawdata/turks/human_"+basename+".bed"
    #chimpbed = "rawdata/turks/"+refname+"_"+basename+".bed"

    if os.path.exists( zippedvcf ) and not forceFlag :
        write2log( " - File"+zippedvcf+" already exists. Skipping command!", True )
        return zippedvcf

    if os.path.splitext( targetfile )[1] == ".gz" : 
        FH = gzip.open( targetfile )
    else :
        FH = open(targetfile)

    header = []
    comments = []
    while peekChar(FH) == "#" :
        line = FH.readline().rstrip()
        if line[:2] == "##" : comments.append(line)
        elif line[0] == "#" : header = line.split("\t")
        else : print "Error?:",line; sys.exit(1)
        
    OUT = open(fixedvcf, "wb")
    # write out comment lines
    if len(comments) > 0 : OUT.write( "\n".join(comments) +"\n" )
    if len(header) > 0 : OUT.write( "\t".join(header)+"\n" )
    OUT.close()

    tmp_lines = FH.readlines(BUF_SIZE)

    group = 0
    while tmp_lines:
        print "Group:",group, "len:",len(tmp_lines)
        genodata = vcf2Bed( tmp_lines, humanbed )
        processVCFLines( genodata, humanbed, chimpbed, fixedvcf, refname )
        for line in tmp_lines :
            row = line.rstrip().split("\t")
        tmp_lines = FH.readlines(BUF_SIZE)
        group +=1

    zippedvcf = tabixbgzip( fixedvcf )
    return zippedvcf
# END forceRef
 
######################################################################
# referenceStats
######################################################################
def referenceStats( targetfile ):
    write2log("#"*70+"\n# Running"+whoami()+"\n"+"#"*70, True )

    filepath, basename, suffix = getBasename( targetfile )
    statsfile = "%s/%s_refstats.txt" % (filepath, basename)

    if os.path.exists( statsfile ) and forceFlag is False :
        write2log( " - File"+statsfile+" already exists. Skipping command!", True )
        return statsfile

    comments = []
    header = []
    #if os.path.splitext( targetfile )[1] == ".gz" :
        #FH = gzip.open( targetfile )
    #else :
        #FH = open(targetfile)

    FH = conditionalOpen( targetfile )

    vardata = {}
    limit = 0
    passstats = {"chimp":{}, "human":{}}
    dbsnpstats = {"chimp":{}, "human":{}}
    bins = [0.0,.05] + [float(x)/10 for x in range(1,11)]
    missingbins = {"chimp":{}, "human":{}}
    missingbins["chimp"] = {x:0 for x in bins}
    missingbins["human"] = {x:0 for x in bins}

    AFbins = {"chimp":{}, "human":{}}
    AFbins["chimp"] = {x:0 for x in bins}
    AFbins["human"] = {x:0 for x in bins}
    #for row in csv.reader( FH, delimiter="\t" ) :
    for line in FH :
        if line[:2] == "##" or line.count("\t") < 1:
            comments.append( line )
            continue
        row = line.split("\t")
        if row[0] == "#CHROM" :
            header = row
            continue
        start = int(row[1])
        passed = row[6]
        reftype = "human"
        if row[7] == "Chimpref" : reftype = "chimp"
        dbsnp = str(row[2] != ".")
        genotypes = [ x[:3] for x in row[9:] ]
        missingness = float(genotypes.count("./.")+genotypes.count(".|.")) / len(genotypes)
        totalgeno = []
        for geno in genotypes : 
            #totalgeno += geno.split("/")
            totalgeno += re.split("[\/\|]", geno)
        AF = (float(totalgeno.count("1")) / 
              (totalgeno.count("0") + totalgeno.count("1")))
        #print passed, reftype, dbsnp, missingness, AF 
        if not passstats[reftype].has_key( passed ) : 
            passstats[reftype][passed] = 0
        passstats[reftype][passed] += 1
        if not dbsnpstats[reftype].has_key( dbsnp ) : 
            dbsnpstats[reftype][dbsnp] = 0
        dbsnpstats[reftype][dbsnp] += 1

        for i in range(len(bins)) :
            #print missingness, i, bins[i], bins[i+1], bins[i+1] > missingness
            if bins[i+1] >= missingness :
                missingbins[reftype][bins[i]] += 1
                break
        #print "AF -",AF
        for i in range(len(bins)) :
            if bins[i+1] >= AF :
                AFbins[reftype][bins[i]] += 1
                break

    write2log( "Writing file: "+statsfile, True )
    OUT = open(statsfile, "wb")
    OUT.write("/\n")
    OUT.write( "Missingness\n" )
    OUT.write( "Rtype\t"+"\t".join([ str(x) for x in bins ]) +"\n" )
    for rtype in missingbins :
        OUT.write( rtype+"\t"+"\t".join([ str(missingbins[rtype][x]) for x in bins ]) +"\n")

    OUT.write("/\n")
    OUT.write( "Allele Frequency\n" )
    OUT.write( "Rtype\t"+"\t".join([ str(x) for x in bins ]) +"\n" )
    for rtype in AFbins :
        OUT.write(rtype+"\t"+"\t".join([ str(AFbins[rtype][x]) for x in bins ])+"\n")

    OUT.write("/\n")
    OUT.write( "DBsnp\n" )
    OUT.write("Rtype\t"+"\t".join([ str(x) for x in dbsnpstats[rtype]])+"\tRatio\n")
    for rtype in dbsnpstats :
        ratio = float(dbsnpstats[rtype]["True"])/sum(dbsnpstats[rtype].values())
        OUT.write( rtype+"\t"+"\t".join([ str(dbsnpstats[rtype][x]) for x in dbsnpstats[rtype] ])+"\t"+"%.3f"%ratio +"\n" )

    OUT.write("/\n")
    OUT.write( "Passed\n" )
    OUT.write( "Rtype\t"+"\t".join([ str(x) for x in passstats[rtype] ]) +"\n" )
    for rtype in passstats :
        OUT.write( rtype+"\t"+"\t".join([ str(passstats[rtype][x]) for x in passstats[rtype] ]) +"\n" )
    OUT.close()
    return statsfile
# End referenceStats

######################################################################
# make_refplots
######################################################################
def make_refplots( statsfile ) :
    write2log( " - Running:"+whoami(), True )
    targetdir, filename, suffix = getBasename(statsfile)
    # MAKE REF PLOTS
    command = "./rscripts/refstats.R %s" % os.path.abspath(statsfile)
    runout = runCommand(command)
    print runout
    return
# End make_refplots

######################################################################
# Main
######################################################################
if __name__ == "__main__" :
    os.chdir("..")
    filename = "./rawdata/turks/turks.names.vcf.gz"
    #outdir = "./output"
    outdir = "./liftover/tmp"

    changeLogFile( LOGDIR+"/changeref_test.log" )

    if len(sys.argv) > 1 :
        filename = sys.argv[1]
    else :
        write2log( "Error: no filename submitted" )

    chimpvcf = forceRef( filename, "chimp" )

    statsfile = referenceStats( chimpvcf )

    make_refplots( statsfile )

    #print "Direct",direct,"File",file
    #outputfile = "%s/%s.%s" % (outdir, basename,"txt")
    #print "Outfile:",outputfile

    #humanbed = "tmp/humanpos.bed"
    #chimpbed = "tmp/chimppos.bed"
    #alldata = {}
    #samples = []
    #pops = []
    #OUTBED = open(humanbed, "wb")
    #for row in csv.reader(open(filename),delimiter="\t") :
        #if row[0][:2] == "##" :
            #continue
        #if row[0] == "#CHROM" :
            #samples = row[9:]
            #pops = findPops( samples )
            #continue
        #elif row[0][0] == "c" :
            #print "Error: different reference being used"
            #sys.exit(1)
        #if len(row) < 4 : print row
        #chrom = row[0]
        #start = int(row[1])-1
        #end = int(row[1])+1
        #pos = row[1]
        #ref = row[3]
        #alt = row[4]
        #key = "chr%s:%s:%s-%s" % (chrom,pos,ref,alt)
        #genes = intersectWithGenes( chrom, pos, regionsdata )
        #refcounts, altcounts = calculatePopFreqs( row[9:], pops )
        #alldata[key] = ["","",ref]+refcounts+[alt]+[altcounts,genes,pos]
        #OUTBED.write("chr%s\t%d\t%d\t%s\n" % (chrom,start,end, key))
    #OUTBED.close()

    #liftOver( humanbed, chimpbed )

    #for row in csv.reader(open(humanbed),delimiter="\t") :
        #chrom, start, end, key = row
        #seq = getSequence( chrom, start, end, "human" )
        #alldata[key][0] = seq
        ##print "%s\t%s\t%s\t%s" % (chrom, start, end, seq)
#
    #for row in csv.reader(open(chimpbed),delimiter="\t") :
        #chrom, start, end, key = row
        #seq = getSequence( chrom, start, end, "chimp" )
        #alldata[key][1] = seq
        ##alldata[key].append(seq)
        ##print "%s\t%s\t%s\t%s" % (chrom, start, end, seq)
#
    #OUT = open( "dadifile.txt", "wb")
    #for entry in sorted(alldata) :
        #OUT.write( alldata[entry] +"\n" )
    #OUT.close()
    ##regionlength = 600
    ##snpBasedFasta( filename, regionlength )
    ##genericBasedFasta( filename, regionlength, outputfile )
#
#
## END
