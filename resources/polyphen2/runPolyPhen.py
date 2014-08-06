#!/usr/bin/python

import os.path
import csv
from time import clock, time
import sys
import commands
import getopt

sys.path.append('/home/database/src')
import db
import queries

######################################################################
# any
######################################################################
def any(iterable):
    for element in iterable:
        if element:
            return True
    return False
# END any

######################################################################
# existsInPath 
######################################################################
def existsInPath(executable) :
    allpaths = [os.path.join(p, executable) for p in os.environ["PATH"].split(os.pathsep)]
    #print [os.path.exists(p) for p in allpaths ]
    if not any([os.path.exists(p) for p in allpaths ]) :
        #print "can't find %s" % executable
        return False
    return True
# END existsInPath

######################################################################
# runCommand
######################################################################
def runCommand( targetcommand, debugflag=False ):
    if debugflag :
        print targetcommand
    runout = commands.getoutput( targetcommand )
    if debugflag :
        print runout
    return runout
# END runCommand


######################################################################
# makeVariantFile
######################################################################
def makeVariantFile( conn, limit=30, gene=None ) :
    if gene is not None :
        allvars = retrieveVariants(conn, limit)
    else :
        allvars = retrieveVariantsByGene(conn, gene, limit)
    tmpfile = "./tmp/variants.input"
    OUT = open(tmpfile, "wb")
    seen = []
    for var in allvars :
        if len(str(var[0])) < 4 :
            varid = "chr%d:%d\t%s/%s" % (var)
            #print "Var:chr%d:%d\t%s/%s" %(var)
            #print "Var:%s" %(varid)
            if varid not in seen :
                OUT.write( "chr%d:%d\t%s/%s\n" %(var) )
                seen.append(varid)
    OUT.close()
    return tmpfile
# End makeVariantFile

######################################################################
# makeTmpVariantFile
######################################################################
def makeTmpVariantFile( allvars ) :
    tmpfile = "./tmp/variants.input"
    OUT = open(tmpfile, "wb")
    seen = []
    for var in allvars :
        if len(str(var[0])) < 4 :
            varid = "chr%d:%d\t%s/%s" % (var)
            #print "Var:chr%d:%d\t%s/%s" %(var)
            #print "Var:%s" %(varid)
            if varid not in seen :
                OUT.write( "chr%d:%d\t%s/%s\n" %(var) )
                seen.append(varid)
    OUT.close()
    return tmpfile
# End makeTmpVariantFile

######################################################################
# mapSnps
######################################################################
def mapSnps( variantfile ):
    print "Varfile:",variantfile
    pphinputfile = "./tmp/subs.pph.input"
    featuresfile = "./tmp/snps.features"
    logfile = "./tmp/mapsnps.log"
    #command = '~/bin/mapsnps.pl -g hg19 -m -U -y subs.pph.input snps.list 1>snps.features 2>mapsnps.log &'
    command = '~/bin/mapsnps.pl -g hg19 -m -U -y %s %s 1>%s 2>%s' \
        %(pphinputfile, variantfile, featuresfile, logfile )
    print command
    runout = runCommand( command )
    print runout
    return pphinputfile
# End mapSnps

######################################################################
# runPph
######################################################################
def runPph( pphinputfile ):
    print "PPH file:",pphinputfile
    featuresfile = "./tmp/pph.features"
    logfile = "./tmp/run_pph.log"
    command = 'run_pph.pl %s 1>%s 2>%s' \
        %(pphinputfile, featuresfile, logfile )
    print command
    runout = runCommand( command )
    print runout
    return featuresfile
# End runPph

######################################################################
# runWeka
######################################################################
def runWeka( featuresfile ):
    print "Features file:",featuresfile
    predfile = "./tmp/pph.predictions"
    command = 'run_weka.pl %s 1>%s' \
        %( featuresfile, predfile )
    print command
    runout = runCommand( command )
    print runout
    return predfile
# End runWeka

######################################################################
# parsePredictions
######################################################################
def parsePredictions( predfile ):
    FILE = open(predfile, "r")
    reader = csv.reader( FILE, delimiter="\t" )
    header = []
    allpreds = []
    for row in reader :
        keys = [x.strip() for x in row]
        if len(header) == 0 :
            keys.append("GenomicPos")
            values = range(0,len(keys))
            header = dict(zip(keys,values))
            #print header
        else :
            notes = keys[header["GenomicPos"]]
            varid, mut, ucnum, gen, acc = notes[2:].split("|")
            chrom, pos = varid.split(":")
            chrnum = chrom[3:]
            #print "chr:",chrom,"num:",chrnum,"pos:",pos,
            #print "Pred:",keys[header["prediction"]],
            #print "Prob:",keys[header["pph2_prob"]]
            allpreds.append( [chrnum,pos,keys[header["o_pos"]],keys[header["acc"]],acc,keys[header["prediction"]],keys[header["pph2_prob"]]] )
    FILE.close()
    return allpreds
# End parsePredictions

######################################################################
# runPolyPhen
######################################################################
def runPolyPhen( chrom, pos, ref, mut ) :
    variantfile = "./tmp/variants.input"
    if len(str(var[0])) < 4 :
    varid = "chr%d:%d\t%s/%s" % (var)
    OUT = open(tmpfile, "wb")
    OUT.write( "chr%d:%d\t%s/%s\n" %(chrom,pos,ref,mut) )
    OUT.close()

    pphinputfile = mapSnps( variantfile )
    featuresfile = runPph( pphinputfile )
    predfile = runWeka( featuresfile )
    print "predfile:",predfile
    predictions = parsePredictions(predfile)
    print predictions
# End runPolyPhen

######################################################################
# Main
######################################################################
if __name__ == "__main__" :
    #conn = db.Conn("localhost")

   # get options
    if len(sys.argv) < 4:
        print "Error: no options given: chrom pos ref mut"
        sys.exit(1)


    chrom, pos, ref, mut = sys.argv
    runPolyPhen( chrom, pos, ref, mut )

# End Main
