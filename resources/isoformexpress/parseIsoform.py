#!/usr/bin/python

import os,sys
import csv
import getopt
import pybedtools
import gzip
from pandas import *

################################################################################
# isoformExpression
AllTissues = [
    "Brain-Nucleus_accumbens_basal_ganglia",
    "Brain-Putamen_basal_ganglia",
    "Brain-Spinal_cord_cervical_c-1",
    "Brain-Substantia_nigra",
    "Breast-Mammary_Tissue",
    "Brain-Cerebellum",
    "Brain-Cortex",
    "Brain-Frontal_Cortex_BA9",
    "Brain-Hippocampus",
    "Brain-Hypothalamus",
    "Colon-Transverse",
    "Esophagus-Mucosa",
    "Esophagus-Muscularis",
    "Fallopian_Tube",
    "Heart-Atrial_Appendage",
    "Heart-Left_Ventricle",
    "Kidney-Cortex",
    "Liver","Lung",
    "Muscle-Skeletal",
    "Nerve-Tibial","Ovary","Pancreas",
    "Pituitary","Stomach","Testis","Thyroid","Uterus",
    "Vagina", "Whole_Blood",
    "Adipose-Subcutaneous","Adipose-Visceral_Omentum",
    "Adrenal_Gland","Artery-Aorta","Artery-Coronary",
    "Artery-Tibial","Brain-Amygdala",
    "Brain-Anterior_cingulate_cortex_BA24", 
    "Brain-Caudate_basal_ganglia",
    "Brain-Cerebellar_Hemisphere",
]
    #"Cells-EBV-transformed_lymphocytes",
    #"Cells-Leukemia_cell_line_CML",
    #"Cells-Transformed_fibroblasts",
#
################################################################################
def isoformExpression( isoformfile, generef, tissuelist, force=False ) :
    #isodict = read_csv( "rawdata/CG_named_loci.tsv", sep="\t" )
    #isodict.index= isodict["IsoformID"]

    if os.path.exists(isoformfile) and not force : 
        genedata = read_csv(isoformfile, sep="\t",index_col="GeneSymbol")
        return genedata

    #tissue = "Brain-Cerebellar_Hemisphere"
    #targetfile = "rawdata/Brain-Cerebellar_Hemisphere_fpkm.tsv"
    isodata = {}
    for tissue in tissuelist :
        targetfile = "rawdata/%s_fpkm.tsv.gz" % tissue
        print targetfile, os.path.exists(targetfile)
        limit = -1
        cdata = {}
        FH = gzip.open(targetfile)
        for row in csv.reader(FH,delimiter="\t"):
            if row[0] == "IsoformID": continue
            isoid = row[0]
            meanscore = sum([float(x.replace(",",".")) for x in row[1:]])/len(row)
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
    print "Running parseIsoformRegions", targetfile
    regionfile = "CG.gtf"

    if os.path.exists(targetfile) and not force : 
        print "Warning: isoformfile already exists. skipping..."
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

def makeExpressionFile( exprdata, gtfdata, targetdir, tissue, force=False ):
    print "Running makeExpressionFile:",tissue
    exprfile = "%s/%s_expr.txt" % (targetdir, tissue)
    if os.path.exists( exprfile ) and not force :
        print "Warning: Tissue found already",tissue,"skipping..."
        return exprfile

    OUT = open( exprfile , "w" )
    OUT.write("chrom\tstart\tend\tGene\tisocov\tisocount\tpsi\tsubexpr\ttotalexpr\tIsoforms\n")
    for index, irange in gtfdata.groupby(["GeneSymbol","chrom"]) :
        gene,chrom = index
        geneexpr = DataFrame({"IsoformID":exprdata.ix[gene]["IsoformID"].split(';'),
                             'Expr':exprdata.ix[gene][tissue].split(';')})
        isoforms = irange.IsoformID.unique()
        geneexpr = geneexpr[geneexpr.IsoformID.isin(isoforms)]
        totalexpr = sum(geneexpr.Expr.map(float))
        generanges = sorted(set(irange.start.tolist()+irange.end.tolist()))

        covdf = intersectRanges( generanges, irange)
        for i in range(len(generanges)) :
            if covdf.ix[i].RangeCount == 0 : continue
            isoset = covdf.ix[i].Isoset
            subexpr = sum(geneexpr[geneexpr.IsoformID.isin(isoset)].Expr.map(float))
            psi = "%.1f" % (subexpr/totalexpr*100) if totalexpr > 0 else "-"
            s = ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%.3f\t%.3f\t%s\n" %
                             (chrom, int(generanges[i]),
                              int(generanges[i+1]),gene, 
                              covdf.ix[i].RangeCount, 
                              len(isoforms), 
                              psi, subexpr, totalexpr,
                              ",".join(covdf.ix[i].Isoset)) )
            OUT.write( s )

    OUT.close()
    return exprfile
# END makeExpressionFile

################################################################################
# calculateRelativeExpression
################################################################################
def calculateRelativeExpression( targetdir="results" ):
    print "Running calculateRelativeExpression"
    generef = read_csv("CG_named_loci.tsv",sep="\t")
    generef.columns = [x.strip() for x in generef.columns]

    tissuelist = [
                  "Brain-Cerebellum",
                  "Brain-Cortex"
                  #"Brain-Frontal_Cortex_BA9",
                  #"Brain-Hippocampus"
                  #"Brain-Hypothalamus"
                ]
    #tissuelist = AllTissues 

    targetfile = "results/isoexpr.tsv"
    exprdata = isoformExpression(targetfile, generef, tissuelist, force=False)

    gtfsummaryfile = os.path.join(targetdir,"gtfsummary.tsv")
    gtfdata = parseIsoformRegions(gtfsummaryfile, generef, force=False)

    allfiles = {}
    tissueexprdir = os.path.join(targetdir,"tissueexpr")
    for tissue in tissuelist : 
        cfile = makeExpressionFile( exprdata, gtfdata, tissueexprdir, tissue )
        allfiles[tissue] = cfile
                                   
    print allfiles
    alldata = None
    psifile = os.path.join(targetdir,"psi_results.txt")
    psibedfile = "globalpsi.bed"
    for tissue in allfiles : 
        print tissue, allfiles[tissue]
        tdata = read_csv( allfiles[tissue], sep="\t" )
        print tdata.head()
        if alldata is None :
            alldata = tdata[["chrom","start","end","Gene","isocov","isocount"]]

        alldata[tissue] = tdata["psi"]

    print alldata[allfiles.keys()].head()
    #alldata["Global-PSI"] = "-"
    psi = []
    for idx, row in alldata[allfiles.keys()].iterrows() :
        subrow = [float(x) for x in row if x != "-"]
        if len(subrow) == 0 : psi.append("-")
        else : psi.append("%.1f"%(sum(subrow)/len(subrow)))

    alldata["Global-PSI"] = psi
    alldata.rename(columns={'chrom':'#chrom'}, inplace=True)
    print alldata["#chrom"].unique()
    targetchromosomes = ["chr"+str(x) for x in range(1,25)+["X","Y"]]
    cleandf = alldata[(alldata["#chrom"].isin(targetchromosomes)) &
                     (alldata["start"] < alldata["end"])]
    cleandf["chrom"] = [x[3:] if x[3:] in ["X","Y","M"] else int(x[3:])
                        for x in cleandf["#chrom"]]

    cleandf.sort(["chrom","start","end"], inplace=True)
    #alldata["Global-PSI"] = alldata[allfiles.keys()].mean(0,numeric_only=True)
    print cleandf.head()
    cleandf.to_csv(psifile, sep="\t", index=False )
    tcols = ["#chrom","start","end","Gene","isocov","isocount","Global-PSI"]
    cleandf[tcols].to_csv(psibedfile, sep="\t", index=False )
    return psibedfile
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

    targetdir = "results"
    testexprfile = "results/relativeexpr1.bed"
    exprfile = "results/relativeexpr.bed"
    psibedfile = "globalpsi.bed"
    if optlist.has_key("-c") : 
        calculateRelativeExpression()
    else: 
        query = optlist.get("-s",None)
        assert query.find(":") > 0, "Error: incorrect formatted query"+ query
        chrom,pos = query.split(":")
        tmpfile = "results/snps.bed"
        OUT = open(tmpfile,"w")
        OUT.write("chr%s\t%s\t%d\n" %(chrom,pos,int(pos)+1))
        OUT.close()
        a = pybedtools.BedTool(psibedfile)
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


