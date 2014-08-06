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
import commands
from time import clock,time
import datetime
import subprocess
import glob
import getopt
import string
import inspect
import csv
import re
import gzip
from sets import Set
import pandas

#-------------------------------------------------------------------#
#                            GLOBALS                                #
#-------------------------------------------------------------------#
DEBUGFLAG = True
LOGDIR = "rawdata/logs"
LOGFILE="%s/default.txt" % LOGDIR
forceFlag = False

#-------------------------------------------------------------------#
#                       LOG FUNCTIONS                               #
#-------------------------------------------------------------------#

######################################################################
# makeDir
######################################################################
def makeDir( targetdir ):
    if not os.path.exists(targetdir):
        print "Making dir:",targetdir
        os.makedirs(targetdir)
# END makeDir

######################################################################
# copyToSubDir
######################################################################
def copyToSubDir( tfile, foldername, force=False ) :
    filepath, filename = os.path.split( tfile )
    tdir = os.path.join(filepath, foldername)
    makeDir( tdir )

    targetfile = os.path.join(tdir,filename)
    if not os.path.exists(targetfile) or force:
        command = ["cp",os.path.join(filepath,filename), tdir]
        retval = subprocess.call( command )
    assert os.path.exists(targetfile)
    return targetfile
# End copyToSubDir

######################################################################
# changeLogDir
######################################################################
def changeLogDir( intermediatedir ):
    global LOGDIR
    newlogdir = "rawdata/%s/logs" % intermediatedir
    LOGDIR = newlogdir
    makeDir( newlogdir )
# End changeLogDir

######################################################################
# changeLogFile
######################################################################
def changeLogFile( filename, overwrite=False ):
    global LOGFILE
    global LOGDIR
    path, basename, suffix = getBasename(filename)
    LOGFILE = filename
    if os.path.exists(path):
        LOGDIR = path
    else :
        print "Error: path does not exist!",path, basename, suffix
        sys.exit(1)
    if overwrite :
        open( filename, "w" ).write("")
# End changeLogFile

######################################################################
# whoami
######################################################################
def whoami( ):
    # return function name
    return inspect.stack()[1][3]
#End whoami

######################################################################
# write2log
######################################################################
def write2log( text, echo=False, overwrite=False, filename=None ):
    #print "Echo", echo, "Overwrite",overwrite,"filename", filename
    if echo :
        print text
    if filename is None :
        filename = LOGFILE
    else :
        changeLogFile(filename)
    if os.path.exists( filename ) and overwrite == False :
        open( filename, "a" ).write("\n"+text+"\n")
    else :
        open( filename, "w" ).write("\n"+text+"\n")
    return
# End write2log

#-------------------------------------------------------------------#
#                       Helper Functions                            #
#-------------------------------------------------------------------#

######################################################################
# conditionalOpen
######################################################################
def conditionalOpen( targetfile ) :
    if os.path.splitext( targetfile )[1] == ".gz" :
        FH = gzip.open( targetfile )
    else :
        FH = open(targetfile)
    return FH
# End conditionalOpen

######################################################################
# toggleForceFlag
######################################################################
def toggleForceFlag( ):
    global forceFlag
    forceFlag = True
# End toggleForceFlag

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
# uniqify
######################################################################
def uniqify(seq, idfun=None): 
   # order preserving
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result
# End uniqify

######################################################################
# runCommand
######################################################################
def runCommand( targetcommand, debugflag=DEBUGFLAG ):
    if debugflag :
        write2log( targetcommand )
    status, runout = commands.getstatusoutput( targetcommand )
    if status != 0 :
        print "Error: Status -", status
        print "Command:",targetcommand
        print "Runout:", runout
        sys.exit(1)
    if debugflag :
        write2log( runout )
    return runout
# END runCommand

######################################################################
# is_int
#"""Returns true if a can be an interger"""
######################################################################
def is_int(a) :
    try:
        int (a)
        return True
    except:
        return False
# END is_int

def is_float(a) :
    try:
        float (a)
        return True
    except:
        return False
# END is_float

######################################################################
# changeUlimit
# 8192, 4096
######################################################################
def changeUlimit( newlimit=4096 ): 
    #write2log( " - Running "+whoami(), True )
    if not is_int( newlimit) :
        print "Error: limit is not an int", newlimit
        sys.exit(1)
    command = " ulimit -n %d" % newlimit
    try :
        out = runCommand( command )
    except:
        print "Warning: could not change ulimit"
    return
# END changeUlimit

######################################################################
# existsInPath
######################################################################
def existsInPath( executable ):
    allpaths = [os.path.join(p, executable) for p in os.environ["PATH"].split(os.pathsep)]
    if not any([os.path.exists(p) for p in allpaths ]) :
        #print "can't find %s" % executable
        return False
    return True
# END existsInPath

######################################################################
# filesInDir
######################################################################
def filesInDir( path, filetype, prefix=None ):
    filelist = []
    if prefix is not None :
        filelist = glob.glob( os.path.join( path, "%s*%s" %(prefix, filetype)))
    else :
        filelist = glob.glob( os.path.join( path, "*."+filetype ) )
    if filetype[-5:] == "fastq" :
        prefixdict = {}
        for file in sorted(filelist) :
            prefix = getBasename( file )
            if not prefixdict.has_key(prefix[:-3]) :
                prefixdict[prefix[:-3]] = []
            prefixdict[prefix[:-3]].append( file )
        return prefixdict.values()
    return filelist
# End filesInDir

######################################################################
# getBasename
######################################################################
def getBasename( originalfilename ) :
    filepath, filename = os.path.split( originalfilename )
    #print "path",filepath,"name", filename
    if originalfilename[-3:] == ".gz" :
        basename, suffix = os.path.splitext( filename[:-3] )
        #print "Base:", basename, "Suff:",suffix
    else :
        basename, suffix = os.path.splitext( filename )
    return filepath, basename, suffix
# End getBasename

######################################################################
# parse_patid
######################################################################
def parse_patid( originalpatientid ) :
    res = re.search( "\d+\-\d+\-\d+", originalpatientid )
    if res is not None :
        #print res.group(0)
        return res.group(0)
    res = re.search( "\d+\-[IV]+\-\d+\-\d+", originalpatientid )
    if res is not None :
        #print res.group(0)
        return res.group(0)
    return originalpatientid
# parse_patid

######################################################################
# iupac
#G/C->S, A/T->W, G/A->R, T/C->Y, G/T->K, A/C->M
######################################################################
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
# printCurrentTime
######################################################################
def printCurrentTime() :
    now = datetime.datetime.now()
    print "Current date and time using strftime:"
    print now.strftime("%Y-%m-%d %H:%M")
# END printCurrentTime()

######################################################################
# tabixbgzip
######################################################################
def tabixbgzip( tfile, force=forceFlag ) :
    write2log( "Running "+whoami() )
    assert os.path.exists( tfile )
    print " - Running tabixbgzip"
    if tfile[-3:] == ".gz" :
        return tfile

    outfile = tfile+".gz"
    if os.path.exists( outfile ) and not force :
        write2log( " - File"+outfile+" already exists. Skipping command!" )
        return outfile
    command = "bgzip -f %s; tabix -p vcf %s" % (tfile,outfile)
    runCommand( command, True )
    return outfile
# END tabixbgzip

######################################################################
# fileIsEmpty
######################################################################
def fileIsEmpty( tfile ):
    print "fileIsEmpty:",tfile
    assert os.path.exists(tfile) 
    statinfo = os.stat(tfile)
    return statinfo.st_size == 0
# END fileIsEmpty

######################################################################
# fileIsOk
######################################################################
def fileIsOk( tfile ):
    print "Tfile:",tfile
    return os.path.exists(tfile) and not fileIsEmpty(tfile) and forceFlag is False
# END fileIsOk

######################################################################
# fixfilepath
######################################################################
def fixfilepath( tfile, shouldexist=True ) :
    newfile = None
    if os.path.isabs( tfile ) :
        newfile = tfile
    elif tfile.find( "~" ) > -1 :
        #print "Expand user"
        newfile = os.path.expanduser( os.path.normpath(tfile) )
    else :
        print "Correct relative path"
        newfile = os.path.realpath( os.path.normpath(tfile) )

    if shouldexist and not os.path.exists(newfile):
        print "Error: File doesnt exist and should? -", tfile
        sys.exit(1)
    return newfile
# End fixfilepath

def readColumn( inputfile, columnnum=1 ) :
    column = subprocess.check_output(["cut", "-f%d" %columnnum, inputfile])
    #print column[:10]
    return [x.rstrip() for x in column.split("\n") if len(x) > 0]
# END readColumn

######################################################################
# CommentStripper
######################################################################
def CommentStripper (iterator):
    for line in iterator:
        if line [:1] == '#':
            continue
        if not line.strip ():
            continue
        yield line #.replace('\x00', '')

######################################################################
# peekChar
######################################################################
def peekChar(fileobj):
    ch = fileobj.read(1)
    fileobj.seek(-1,1)
    return ch

######################################################################
# dfTable
######################################################################
def dfTable( targetlist ):
    if type(targetlist) is pandas.core.series.Series :
        targetlist = targetlist.tolist()
    return pandas.DataFrame( [[x,targetlist.count(x)] for x in Set(targetlist)],
                     columns=["Variable","Count"]).sort("Count",ascending=0)
# END dfTable
    
######################################################################
# dfProportion
######################################################################
def dfProportion( targetlist ):
    from sets import Set
    if type(targetlist) is pandas.core.series.Series :
        targetlist = targetlist.tolist()
    datapts = float(len(targetlist))
    return pandas.DataFrame( [[x,targetlist.count(x)/datapts] for x in Set(targetlist)],
                     columns=["Variable","Count"]).sort("Count",ascending=0)
# END dfProportion

######################################################################
# MAIN
######################################################################
if __name__ == "__main__" :
    print "I should really write a main section"

    #if len(sys.argv) > 1 :
        #tfile = sys.argv[1]
    #else :
        #print "Error: not arguments provided"
        #sys.exit(1)

    printCurrentTime()
    sys.exit(1)
    print "Orig:",tfile
    print "New:",fixfilepath(tfile)

# MAIN
