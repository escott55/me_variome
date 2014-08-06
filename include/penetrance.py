#!/usr/bin/python

from localglobals import *

    if os.path.splitext( callvcf )[1] == ".gz" :
        FH = gzip.open( callvcf )
    else :
        FH = open(callvcf)

    #for row in csv.reader( FH, delimiter="\t" ) :
    linecount = 0
    for line in FH :
        linecount += 1
        if line[:2] == "##" or len(line) < 1:
            comments.append( line )
            continue
        #print "#",linecount
        row = line.rstrip().split("\t")
        #print linecount,"Row:",row[:3]
        if row[0] == "#CHROM" :
            header = row
            samples = [x.strip() for x in row[9:]]
            #print "Info:",row[:9]
            #print "Samples",samples[:10],"..."
            continue

        chrom,pos,dbsnp,ref,mut,qual,passed = row[:7]
        pos = int(pos)
        dbsnp = str(dbsnp != ".")
        end = pos + max(len(ref), len(mut))
        if len(mut)+len(ref) > 100 : continue
        key = "chr%s:%s:%s:%s" % (chrom, pos, ref, mut )
