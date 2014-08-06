#!/usr/bin/python

import os, sys
from pandas import *

def intersectRanges( grange, iranges ):
    foundvec = [0 for x in range(len(grange))]
    for start,end in iranges.values :
        if start < end :
            for i in range(grange.index(start),grange.index(end)) : foundvec[ i ] +=1 
        else :
            for i in range(grange.index(end),grange.index(start)) : foundvec[ i ] +=1 
    return foundvec
# END intersectRanges

#gtfcols = ["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"]
#gtfdata = read_csv("refseq.gtf",sep="\t",header=None,names=gtfcols)
gtfdata = read_csv("refseqgenes.txt", sep="\t" )


OUT = open("penetrance.txt","wb")
limit = -1
for index, data in gtfdata.groupby(["name2","chrom"]) :
    gene,chrom = index
    allexonranges = ",".join(data.exonStarts.tolist() + data.exonEnds.tolist())
    generanges = sorted(set([int(x) for x in allexonranges.split(",") if len(x) > 0]))
    allstarts = [int(x) for x in (",".join(data.exonStarts.tolist())).split(",") if len(x) > 0]
    allends = [int(x) for x in (",".join(data.exonEnds.tolist())).split(",") if len(x) > 0]
    irange = DataFrame({'start':allstarts,'end':allends})
    covvec = intersectRanges( generanges, irange)
    print "Gene:",gene
    print "# isoforms:",len(data)
    #print covvec
    for i in range(len(generanges)) :
        if covvec[i] == 0 : continue
        OUT.write( "%s\t%s\t%d\t%d\t%d\t%d\n"%(gene,chrom,generanges[i],generanges[i+1],covvec[i], len(data)) )
    if limit == 0 : break
    limit -= 1
OUT.close()

    #for index, isodata in data.iterrows() :
        #print isodata.name
        #iranges = DataFrame({'start':isodata.exonStarts.split(","), 'end':isodata.exonEnds.split(',')})
        #print iranges.head(10)
