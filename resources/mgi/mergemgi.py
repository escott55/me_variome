#!/usr/bin/python

import os, sys
import csv

def removeNonAscii(s): 
    s = s.replace("'","")
    return "".join(i for i in s if ord(i) <128)

if __name__ == "__main__" :

    genefile = "HMD_HumanPhenotype.rpt"
    termfile = "MGI_GenePheno.rpt"
    annotfile = "VOC_MammalianPhenotype.csv"

    FILE = open( annotfile, "r" )
    reader = csv.reader(FILE,delimiter="\t")
    annot = {}
    for row in reader :
        if len(row[0]) == 0 : continue
        annot[row[0].strip()] = removeNonAscii(row[1].strip())

    FILE = open( genefile, "r" )
    reader = csv.reader( FILE, delimiter="\t" )
    genetable = {}
    for row in reader :
        if len(row) == 0 : continue
        #print row[3]
        gene = row[0].strip()
        entrez = row[1].strip()
        genetable[row[3].strip()] = gene+"\t"+entrez
    FILE.close()

    #sys.exit(1)
    FILE = open( termfile, "r" )
    reader = csv.reader( FILE, delimiter="\t" )
    count = 0
    for row in reader :
        pid = row[4].strip()
        mid = row[6].strip()
        if mid not in genetable : 
            count += 1
            continue
        print "%s\t%s\t%s" %( genetable[mid],annot[pid],mid)

    #print "Error Count:",count
# End Main
