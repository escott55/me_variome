#!/usr/bin/python

import os,sys
import csv
import getopt
import pybedtools
from pandas import *

################################################################################
# isoformExpression
#"Brain-Nucleus_accumbens_basal_ganglia",
#"Brain-Putamen_basal_ganglia",
#"Brain-Spinal_cord_cervical_c-1",
#"Brain-Substantia_nigra",
#"Breast-Mammary_Tissue",
#"Cells-EBV-transformed_lymphocytes",
#"Cells-Leukemia_cell_line_CML",
#"Cells-Transformed_fibroblasts",
#"Colon-Transverse",
#"Esophagus-Mucosa",
#"Esophagus-Muscularis",
#"Fallopian_Tube",
#"Heart-Atrial_Appendage",
#"Heart-Left_Ventricle",
#"Kidney-Cortex",
#"Liver","Lung",
#"Muscle-Skeletal",
#"Nerve-Tibial","Ovary","Pancreas",
#"Pituitary","Stomach","Testis","Thyroid","Uterus",
#"Vagina", "Whole_Blood"
#"Adipose-Subcutaneous","Adipose-Visceral_Omentum",
#"Adrenal_Gland","Artery-Aorta","Artery-Coronary",
#"Artery-Tibial","Brain-Amygdala",
#"Brain-Anterior_cingulate_cortex_BA24", 
#"Brain-Caudate_basal_ganglia",
#"Brain-Cerebellar_Hemisphere",

#
################################################################################
def isoformExpression( isoformfile, generef, force=False ) :
    isodict = read_csv( "rawdata/CG_named_loci.tsv", sep="\t" )
    isodict.index= isodict["IsoformID"]

    if os.path.exists(isoformfile) and not force : 
        genedata = read_csv(isoformfile, sep="\t",index_col="GeneSymbol")
        return genedata

    tissuelist = [
                  "Brain-Cerebellum",
                  "Brain-Cortex"
                  #"Brain-Frontal_Cortex_BA9",
                  #"Brain-Hippocampus"
                  #"Brain-Hypothalamus"
                ]

    #tissue = "Brain-Cerebellar_Hemisphere"
    #targetfile = "rawdata/Brain-Cerebellar_Hemisphere_fpkm.tsv"
    isodata = {}
    for tissue in tissuelist :
        targetfile = "rawdata/%s_fpkm.tsv" % tissue
        print targetfile, os.path.exists(targetfile)
        limit = -1
        cdata = {}
        FH = open(targetfile)
        for row in csv.reader(FH,delimiter="\t"):
            if row[0] == "IsoformID": continue
            isoid = row[0]
            meanscore = sum([float(x) for x in row[1:]])/len(row)
            #print isoid, meanscore
            cdata[isoid] = meanscore
            limit -= 1
            if limit == 0 : break
        isodata[tissue] = Series( cdata )
        FH.close()

    isodata = DataFrame( isodata ).reset_index()

    genedata = merge(generef,isodata,left_on="IsoformID",right_on="index").drop("index",axis=1)
    #print genedata.head(10)
    genedata = genedata.groupby("GeneSymbol").agg(lambda col: 
                                              ";".join(col.map(str).tolist()))
    #print genedata.head(10)
    genedata.to_csv(isoformfile,sep="\t")
    return genedata
# END isoformExpression

def intersectRanges( grange, iranges ):
    foundvec = [0 for x in range(len(grange))]
    isovec = [[] for x in range(len(grange))]
    for isoid,start,end in iranges[["IsoformID","start","end"]].values :
        if start < end :
            for i in range(grange.index(start),grange.index(end)) : 
                foundvec[ i ] +=1
                isovec[i].append(isoid)
        else :
            for i in range(grange.index(end),grange.index(start)) :
                foundvec[ i ] +=1
                isovec[i].append(isoid)
    founddf = DataFrame({"RangeCount":foundvec, "Isoset":isovec})
    return founddf
# END intersectRanges

def parseIsoformRegions( targetfile, generef, force=False, limit=-1):
    print "Running parseIsoformRegions"
    regionfile = "CG.gtf"

    if os.path.exists(targetfile) and not force : 
        gtfdata = read_csv(targetfile, sep="\t")
        return gtfdata

    generegions = []
    for row in csv.reader(open(regionfile),delimiter="\t") :
        info = Series( dict([grp.strip().replace('"','').split(" ") 
                       for grp in row[8].rsplit(";")
                       if len(grp) > 0]) )
        generegions.append([info.transcript_id, row[0],row[3],row[4]])
        limit -= 1
        if limit == 0 : break

    gtfdata = DataFrame(generegions, columns=["IsoformID","chrom","start","end"])
    gtfdata = merge(gtfdata,generef, on="IsoformID").reset_index(drop=True)
    gtfdata.to_csv(targetfile, index=False, sep="\t")
    return gtfdata
# END parseIsoformRegions

################################################################################
# calculateRelativeExpression
################################################################################
def calculateRelativeExpression( exprfile ):
    print "Running calculateRelativeExpression"
    generef = read_csv("CG_named_loci.tsv",sep="\t")
    generef.columns = [x.strip() for x in generef.columns]

    targetfile = "results/isoexpr.tsv"
    exprdata = isoformExpression(targetfile, generef, force=True)

    gtfsummaryfile = "results/gtfsummary.tsv"
    gtfdata = parseIsoformRegions(gtfsummaryfile, generef, force=True)

    OUT = open( exprfile , "w" )
    OUT.write("chrom\tstart\tend\tGene\tisocov\tisocount\trelexpr\ttotalexpr\tIsoforms\n")
    for index, irange in gtfdata.groupby(["GeneSymbol","chrom"]) :
        gene,chrom = index
        geneexpr = DataFrame({"IsoformID":exprdata.ix[gene]["IsoformID"].split(';'),
                             'Expr':exprdata.ix[gene]["Brain-Cerebellum"].split(';')})
        isoforms = irange.IsoformID.unique()
        geneexpr = geneexpr[geneexpr.IsoformID.isin(isoforms)]
        totalexpr = sum(geneexpr.Expr.map(float))
        generanges = sorted(set(irange.start.tolist()+irange.end.tolist()))

        covdf = intersectRanges( generanges, irange)
        for i in range(len(generanges)) :
            if covdf.ix[i].RangeCount == 0 : continue
            isoset = covdf.ix[i].Isoset
            subexpr = (geneexpr[geneexpr.IsoformID.isin(isoset)].Expr.map(float))
            relexpr = sum(subexpr)/totalexpr if totalexpr > 0 else 1
            s = ("%s\t%d\t%d\t%s\t%d\t%d\t%.3f\t%.3f\t%s\n" %
                             (chrom, int(generanges[i]),
                              int(generanges[i+1]),gene, 
                              covdf.ix[i].RangeCount, 
                              len(isoforms), 
                              relexpr, totalexpr,
                              ",".join(covdf.ix[i].Isoset)) )
            OUT.write( s )

    OUT.close()
# END calculateRelativeExpression

################################################################################
if __name__ == "__main__" :
    optlist, args = getopt.getopt( sys.argv[1:], "s:c")
    optlist = dict(optlist)

    dataset = optlist.get("-r",None)
    if len(optlist) == 0 :
        print "Error: no dataset provided"
        print " -c calculate all expression values"
        print " -s <chr:pos> retrieve relative expression for query"
        sys.exit(1)

    testexprfile = "results/relativeexpr1.bed"
    exprfile = "results/relativeexpr.bed"
    if optlist.has_key("-c") : 
        calculateRelativeExpression(exprfile)
    else: 
        query = optlist.get("-s",None)
        assert query.find(":") > 0, "Error: incorrect formatted query"+ query
        chrom,pos = query.split(":")
        tmpfile = "results/snps.bed"
        OUT = open(tmpfile,"w")
        OUT.write("chr%s\t%s\t%d\n" %(chrom,pos,int(pos)+1))
        OUT.close()
        a = pybedtools.BedTool(testexprfile)
        #resfile = "results/intersections.bed"
        #a.intersect(testexprfile).saveas(resfile)
        inter = a.intersect(tmpfile)
        print inter
        #, trackline="track name='SNPs in exons' color=128,0,0"

# END MAIN
################################################################################
        #allpoints = [x for plist in irange.start+irange.end for x in plist]
        #irange = []
        #for index, isorow in data.iterrows() :
            #print isorow
            #irange.append(DataFrame({'IsoformID':isorow.IsoformID,
                                     #'start':isorow.start,'end':isorow.end}))
        #irange = concat( irange )
    #generegions = {}
    #for row in csv.reader(open(regionfile),delimiter="\t") :
        #info = Series( dict([grp.strip().replace('"','').split(" ") 
                       #for grp in row[8].rsplit(";")
                       #if len(grp) > 0]) )
        #if not generegions.has_key(info.transcript_id) :
            #generegions[info.transcript_id] = {"chr":"","start":[],"end":[]}
        #generegions[info.transcript_id]["chr"] = row[0]
        #generegions[info.transcript_id]["start"].append(row[3])
        #generegions[info.transcript_id]["end"].append(row[4])
        #limit -= 1
        #if limit == 0 : break
    #gtfdata = DataFrame(generegions).transpose()


