#!/usr/bin/env python

#-------------------------------------------------------------------#
#                       Import Packages                             #
#-------------------------------------------------------------------#
import os, sys
import getopt
import csv
import re
import gzip
#from random import randint

from housekeeping import *

######################################################################
# readAnnotTable
######################################################################
def readAnnotTable( ):
    pfile = "./patientannot/annotation/patientannotation_laser.txt"
    hd = {}
    annottable = {}
    for row in csv.reader(open(pfile),delimiter="\t") :
        if len(hd) == 0 :
            for i in range(len(row)) :
                hd[row[i]] = i
        annottable[row[hd["Normalized"]]] = row
    return annottable
# End readAnnotTable

######################################################################
# batcheffectremove
######################################################################
def batcheffectremove(plates=["Pilot"]) :
    #pfile = "./patientannot/annotation/patientannotation_laser.txt"
    pfile = "./patientannot/annotation/patientannotation.ped"
    hd = {}
    toremove = {}
    for row in csv.reader(open(pfile),delimiter="\t") :
        if len(hd) == 0 :
            for i in range(len(row)) :
                hd[row[i]] = i
            continue
        #print row[hd["Lab"]]
        if (row[hd["Lab"]] == "Broad" and \
            row[hd["Source"]] in plates) :
            patname = row[hd["original"]]
            normpatname = row[hd["Individual.ID"]]
            toremove[patname] = normpatname
    return toremove
# END batcheffectremove
            #["PlateI", "PlateII","Pilot"]) :
            #and row[hd["continent"]] in ["Middle East", "South Asia"]) \
            #or (row[hd["Lab"]] == "Daly" and row[hd["continent"]] == "Europe") \
            #or (row[hd["Lab"]] == "Eichler" and row[hd["continent"]] == "Europe"):
            #if fulllist is not None :
                #fixpatname1 = "-".join(patname.split("-")[1:])
                #if patname in fulllist :
                    #keeppats.append( patname )
                #elif fixpatname1 in subfulllist :
                    #for i in range(len(subfulllist)) :
                        #if fixpatname1 == subfulllist[i] :
                            #break
                    #keeppats.append( fulllist[i] )
               #else :
                    #print "Error: patname not in fulllist:", patname
 
######################################################################
# removePatsList
######################################################################
def removePatsList( targetdir ) :
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
    return removelist
# END removePatsList

######################################################################
# getPats
######################################################################
def getPats( vcffile, keeplist=None, removelist=None ):
    write2log(" - Running "+whoami()+" - "+vcffile, True)
    #pfile = "./merged/combined.ped"

    targetdir, filename, suffix = getBasename(vcffile)
    allpats = currentFilePats( vcffile )
    #print allpats

    keeppats = []
    # If we have a keep list...
    if keeplist is not None :
        if not os.path.exists(keeplist) : 
            print "Error: no keeplist found:",keeplist
            sys.exit(1)
        print "Using keeplist:", keeplist
        for pat in open(keeplist).readlines() :
            pat = pat.strip()
            if pat[0] == "#" : continue
            if pat in allpats : 
                if pat not in removelist :
                    keeppats.append(pat)
    else : 
        # Filter if in removelist
        if removelist is None : removelist = removePatsList(targetdir)
        #removelist2 = batcheffectremove()
        for pat in allpats :
            if pat not in removelist : #and pat not in removelist2.keys() :
                keeppats.append(pat)

    print "Original:",len(allpats), "Kept:",len(keeppats)
    #for pat in fulllist :
        #if pat not in keeppats :
            #print pat
    return keeppats
# END getPats
        #keeppats = [x.strip() for x in open(keeplist).readlines()]
        #print len(keeppats)
        #keeppats = [x for x in keeppats if x in allpats]
        #print len(keeppats)
    #subfulllist = []
    #for pat in fulllist :
        #newpat = "-".join(pat.split("-")[1:])
        #subfulllist.append(newpat)

######################################################################
# parsePatDictionary
######################################################################
def parsePatDictionary( allpats ):
    # Parse Dictionary
    dictfile = "./resources/dict_alldata.txt"
    pdict = {}
    for row in csv.reader(open(dictfile), delimiter="\t") :
        if len(row) == 0  or row[0][0] == "#" : continue
        if len(row) == 8 and row[-1] != "Exome" : continue
        if row[-1] == "Discard" : continue
        pdict[row[1]] = row[0]
    return pdict
# END parsePatDictionary

######################################################################
# getFixedPatList
######################################################################
def getFixedPatList( vcffile ) :
    allpats = currentFilePats( vcffile )
    annot = readAnnotTable()

    pdict = parsePatDictionary(allpats)

    fixed = [ None for x in allpats ]
    for i in range(len(allpats)) :
        pat = allpats[i]
        if pdict.has_key(pat) and annot.has_key(pdict[pat]) :
            #print "Found:",pat,":",pdict[pat]
            fixed[i] = pdict[pat].strip().replace(" ","_")

    #print "Fixed:","F:"+"\nF:".join(fixed)
    #print "Original:","O:"+"\nO:".join(allpats)
    #print "Fixed:","F:"+"\nF:".join([str(x) for x in fixed])
    #print set([x for x in allpats if allpats.count(x) > 1])
    tset = set([x for x in fixed if fixed.count(x) > 1])
    for x in tset:
        if x is None : continue
        print "Replicate:",x
        print [ i for i in range(len(fixed)) if fixed[i] == x]
        print [ fixed[i] for i in range(len(fixed)) if fixed[i] == x]
        print [ allpats[i] for i in range(len(fixed)) if fixed[i] == x]

    return fixed
# END getFixedPatList
    #tset = [x for x in allpats if x.find("brain") > -1 ]
    #for x in tset:
        #print [ i for i in range(len(fixed)) if allpats[i] == x]
        #print [ fixed[i] for i in range(len(fixed)) if allpats[i] == x]
        #print [ allpats[i] for i in range(len(fixed)) if allpats[i] == x]

######################################################################
# fixPatNames
######################################################################
def fixPatNames( vcffile, force=forceFlag ) :
    write2log(" - Running "+whoami(), True)
    targetdir, filename, suffix = getBasename(vcffile)

    newvcf = "%s/%s.names%s" % (targetdir, filename, suffix)
    compnewvcf = "%s/%s.names%s.gz" % (targetdir, filename, suffix)
    if os.path.exists( compnewvcf )  and not force :
        write2log(" - File"+compnewvcf+" already exists. Skipping command!",True)
        return compnewvcf

    fixed = getFixedPatList( vcffile )

    print "Fixed:",fixed

    FH = conditionalOpen(vcffile)
    header = []
    original = []
    line = FH.readline()
    limit = 200
    linenum = 1
    while line :
        if line[:2] == "##" : header.append(line)
        elif line[:6] == "#CHROM" :
            row = line.strip().split("\t")
            original = line
            subrow = [fixed[i] for i in range(len(fixed)) if fixed[i] != None]
            header.append( "\t".join(row[:9] + subrow) + "\n" )
            break
        #else : 
            #row = line.strip().split("\t")
            #subrow = [row[i+9] for i in range(len(row[9:])) if fixed[i] != None]
            #break
        line = FH.readline()
        linenum += 1
        limit -= 1
        if limit == 0 : break
    FH.close()

    print original[:60]
    print header[-1][:60]

    start = time()
    tmp = [i+9 for i in range(len(fixed)) if fixed[i] == None]
    colstokeep = ",".join([str(x) for x in tmp])
    OUT = open( newvcf, "wb" )
    OUT.writelines( header )
    OUT.close()

    print "Cols to keep :",colstokeep
    if len(colstokeep) > 0 and vcffile[-3:] == ".gz" :
        command = "gunzip -c %s | tail -n +%s | cut --complement -f%s >> %s" % (vcffile, linenum+1, colstokeep, newvcf)
    elif len(colstokeep) > 0 :
        command = "tail -n +%s %s | cut --complement -f%s >> %s" % (linenum+1, vcffile, colstokeep, newvcf)
    elif vcffile[-3:] == ".gz" :
        command = "gunzip -c %s | tail -n +%s >> %s" % (vcffile, linenum+1, newvcf)
    else :
        command = "tail -n +%s %s >> %s" % (linenum+1, vcffile, newvcf)

    print command
    out = runCommand( command )
    print out
 
    compnewvcf = tabixbgzip( newvcf )

    elapsed = str(datetime.timedelta(seconds=(time()-start))) 
    write2log( "Original:"+str(len(fixed)), True )
    write2log( "Failedcount:"+str(len([x for x in fixed if x is None])), True )
    return compnewvcf
# END fixPatNames
        #command = "gunzip -c %s | awk '$1 !~ /^#/' | cut --complement -f%s >> %s" % ( vcffile, colstokeep, newvcf)
        #command = "sed -e '/^#/d' %s | cut --complement -f%s >> %s" % ( vcffile, colstokeep, newvcf)
    #write2log( "Fixed:"+str(len(subrow)) )
    #print "\n".join(fixed)
    #FH = conditionalOpen(vcffile)
    #OUT = gzip.open(newvcf,"wb")
    #line = FH.readline()
    #limit = 200
    #while line :
        #if line[:2] == "##" : OUT.write(line)
        #elif line[:6] == "#CHROM" :
            #row = line.strip().split("\t")
            #subrow = [fixed[i] for i in range(len(fixed)) if fixed[i] != None]
            #OUT.write( "\t".join(row[:9]+subrow)+"\n" )
            #break
        #else : 
            #row = line.strip().split("\t")
            #subrow = [row[i+9] for i in range(len(row[9:])) if fixed[i] != None]
            #OUT.write( "\t".join(row[:9]+subrow)+"\n" )
        #line = FH.readline()
        #limit -= 1
        #if limit == 0 : break
    #FH.close()
    #OUT.close()
    #tmp = range(1,9) + [i+9 for i in range(len(fixed)) if fixed[i] != None]



######################################################################
# currentFilePats
######################################################################
def currentFilePats( cfile ):
    write2log(" - Running "+whoami(), True)
    targetdir, filename, suffix = getBasename(cfile)

    patlist = []
    if suffix == ".ped" :
        command = "cut -f 2 %s" % cfile
        out = runCommand( command ).strip()
        patlist = out.split("\n")

    elif suffix == ".txt" :
        command = "cut -f 2 %s" % cfile
        out = runCommand( command ).strip()
        patlist = out.split("\n")
    elif suffix == ".fam" :
        command = '''awk '{print $2}' %s''' % cfile
        out = runCommand( command ).strip()
        patlist = out.split("\n")
    elif cfile[-3:] == ".gz" and suffix == ".vcf" :
        headline = None
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
        out = runCommand( command ).strip()
        #print out[:40]
        patlist = out.split("\t")[9:]
    else :
        print "Error: Unknown filetype:",suffix
        print "filebase:",filename
        sys.exit(1)
    return patlist
# END currentFilePats

######################################################################
# Main
######################################################################
if __name__ == "__main__" :
    os.chdir("..")

    #vcffile = "/home/escott/inbreed/rawdata/turks/turks.vcf.gz"
    vcffile = "/home/escott/workspace/inbreed/rawdata/ciliopathies/ciliopathies.vcf.gz"
    keepfile = "./rawdata/turkstokeep.txt"
    #vcffile = "/home/escott/inbreed/rawdata/daily/daily.vcf.gz"

    changeLogFile( LOGDIR+"/patientInfo_test.log" )

    if len(sys.argv) > 1 :
        vcffile = sys.argv[1]
    else :
        write2log( "Error: no vcffile submitted" )

    pats1 = getPats( vcffile )
    #fixed = getFixedPatList( vcffile )
    #print fixed
    newvcf = fixPatNames( vcffile, True )
    pats2 = getPats( newvcf )
    print "Finished processing:",newvcf
    print len(pats1)
    print len(pats2)
    #keeppats = getPats( vcffile, keepfile )
    #print keeppats

# End Main
