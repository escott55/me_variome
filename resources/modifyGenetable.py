#!/usr/bin/python

import os, sys
from pandas import *

if __name__ == "__main__" :

    tfile = "human_gene_info.txt"
    genes = read_csv(tfile, sep="\t" )
    print genes.head()

    genes = genes[["GeneID","Symbol","Synonyms","chromosome"]]

    s = genes['Synonyms'].str.split('|').apply(Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'Synonyms'
    del genes['Synonyms']
    synonyms = genes.join(s)

    #print genes.head()
    
    synonyms.rename(columns={"Synonyms":"Gene"},inplace=True)
    genes["Gene"] = genes["Symbol"]
    del genes["Symbol"]

    synonyms["class"] = "synonyms"
    genes["class"] = "symbol"
    newgenes = concat([synonyms,genes]).sort("class")

    newgenes = newgenes[(newgenes.Gene != "-") & ~(newgenes.duplicated("Gene"))]

    print newgenes[newgenes.Gene == "AK"]
    #print newgenes.head()
    newgenes.to_csv("human_genes.txt",sep="\t",index=False)
# END MAIN
