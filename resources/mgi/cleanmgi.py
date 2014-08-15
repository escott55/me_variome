#!/usr/bin/python

import os, sys
import csv
import string

if __name__ == "__main__" :
    datafile = "HMD_HumanPhenotype.rpt"
    annotfile = "VOC_MammalianPhenotype.rpt"

    FILE = open( annotfile, "r" )
    reader = csv.reader(FILE,delimiter="\t")
    annot = {}
    for row in reader :
        annot[row[0]] = row[1]

    header = ["Gene","Entrez","Homologene","Mouse","MGI","Pheno","Nothing"]
    FILE = open(datafile, "r" )
    reader = csv.reader(FILE,delimiter="\t")
    for row in reader :
        #print row
        gene = row[0].strip()
        entrez = row[0].strip()
        phenocount = len(string.split(row[5].strip()," "))
        for pid in string.split(row[5]," ") :
            if len(pid) == 0 : continue
            #print pid
            #+[row[4],pid],
            print "%s\t%s\t%s\t%d" % (gene,entrez,annot[pid],phenocount)

# End MAin
