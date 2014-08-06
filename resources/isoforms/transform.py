#!/usr/bin/python

import os,sys
from pandas import *


################################################################################
################################################################################
if __name__ == "__main__" :
    trgtfile = "refseqgenes.txt"

    trgtdata = read_csv(trgtfile,sep="\t")
    print trgtdata.head()

    #subset = trgtdata.head(20)
    #test = [row['name'],row['name2'],row['chrom'], row['exonStarts'].split(','), row['exonEnds'].split(',') for _, row in subset.iterrows()]

    expanded = concat([DataFrame({'acc':row['name'],'gene':row['name2'],
                               'chrom':row['chrom'], 
                               'exonStarts':row['exonStarts'].strip(',').split(','),
                               'exonEnds':row['exonEnds'].strip(',').split(',')})
                              for _, row in trgtdata.iterrows()])

    
    print row["exonStarts"]
    print expanded.head(50)
    expanded.to_csv("refseqstarts.txt",sep="\t",index=False)

