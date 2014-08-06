#!/usr/bin/python

import os, sys
import subprocess
import glob
import pybedtools
#import resource
#resource.setrlimit(resource.RLIMIT_NOFILE, (1000,-1))

from localglobals import *
import roh_figures as rohfig

################################################################################
# intersectRegions
################################################################################
def intersectRegions( chrom, pos, regions ):
    pos = int(pos)
    if len(chrom) < 4 : chrom = "chr"+chrom
    if regions.has_key(chrom) :
        for region in regions[chrom] :
            if region[0] < pos and region[1] > pos :
                #genes.append(region[2])
                return True 
    return False
# End intersectRegions

################################################################################
# calculateRegionsCovered1
################################################################################
def calculateRegionsCovered1( regionsfile ) :
    breakpoints = {}
    for row in csv.reader( open(regionsfile), delimiter="\t" ):
        if hk.is_int(row[0]) :
            chrom, start, end = [int(x) for x in row[:3]]
        elif row[0] == 'X' :
            chrom = 23
            start, end = [int(x) for x in row[1:3]]
        elif row[0] == 'Y' :
            chrom = 24
            start, end = [int(x) for x in row[1:3]]
        else : continue

        if not breakpoints.has_key( chrom ) : breakpoints[chrom] =[]
        breakpoints[chrom].append(start)
        breakpoints[chrom].append(end)
    dist = {}
    for chrom in breakpoints :
        dist[chrom] = 0
        chrombps = sorted(breakpoints[chrom])
        for i in range(len(chrombps)-1) :
            dist[chrom] += chrombps[i+1] - chrombps[i]
    return Series(dist)
# END calculateRegionsCovered1

################################################################################
# calculateRegionsCovered
################################################################################
def calculateRegionsCovered( regionsfile ) :
    breakpoints = {}
    for row in csv.reader( open(regionsfile), delimiter="\t" ):
        if hk.is_int(row[0]) :
            chrom, start, end = [int(x) for x in row[:3]]
        elif row[0] == 'X' :
            chrom = 23
            start, end = [int(x) for x in row[1:3]]
        elif row[0] == 'Y' :
            chrom = 24
            start, end = [int(x) for x in row[1:3]]
        else : continue

        if not breakpoints.has_key( chrom ) : breakpoints[chrom] = 0
        breakpoints[chrom] += abs(end-start)
    return Series(breakpoints)
# END calculateRegionsCovered

################################################################################
# sortBed
################################################################################
def sortBed( bedfile ) :
    filepath,basename,ext = hk.getBasename( bedfile )
    print basename, ext
    sortedfile = "%s/%s_sorted%s" % (filepath, basename, ext)
    print sortedfile
    with open(sortedfile, "w") as outfile:
        subprocess.call( ["sort","-V","-k1,1","-k2,2n",bedfile], stdout=outfile )
    return sortedfile
# END sortBed

################################################################################
# makeIntersectBed
################################################################################
def makeIntersectBed( data, exomeregions, intersectbed, tmprohbedfile ):
    tmproh = data.sort(["CHR","POS1","POS2"])
    tmproh.to_csv(tmprohbedfile, columns=["CHR","POS1","POS2","IID"], sep="\t",index=False,header=False)
    command = ["bedtools","intersect","-a",tmprohbedfile,"-b",exomeregions]#,"-sorted"
    print " ".join(command)
    stdout = subprocess.check_output( command )
    if len(stdout.strip()) == 0 :
        print "Error: stdout - ",stdout
        print "Command -"," ".join(command)
        return "Failed"

    intersectdata = (DataFrame( [x.split("\t") for x in stdout.split("\n")],
                              columns=["chrom","start","end","Individual.ID"] )
                     .drop_duplicates())
    print "Writing file:",intersectbed
    intersectdata.to_csv( intersectbed, header=False, index=False, sep="\t" )
    return "Success"
# END makeIntersectBed

################################################################################
# exomeCoverage
################################################################################
def exomeCoverage( rohfile, sampleannot, outprefix="test", force=False ):
    print "Running exomeCoverage", rohfile
    targetdir = "results/rohcov/%s" % outprefix
    hk.makeDir(targetdir)

    overlapfile = "results/rohcov/%s_chrom_overlap.txt" % (outprefix )
    if os.path.exists(overlapfile) and not force :
        print "Final file already exists:",overlapfile
        bedfiles = retrieveBedFiles(outprefix)
        allbases = read_csv( overlapfile, sep="\t" )
        return allbases,bedfiles

    #[["IID","CHR","POS1","POS2","KB"]]
    rohdf = read_csv(rohfile,delim_whitespace=True, skipinitialspace=True,
                     usecols=["IID","CHR","POS1","POS2","KB"])
                     ##dtype={'POS1':np.int64,'POS2':np.int64})#'CHR':np.int32,

    rohdf[['POS1', 'POS2']] = rohdf[['POS1', 'POS2']].astype(np.int64)
    print "rohdf:",rohdf.shape
    rohannot = merge( sampleannot[["Individual.ID","Continent","ethnicity"]],
                   rohdf, right_on="IID", left_on="Individual.ID", how="right") 
    rohannot = rohannot[rohannot.Continent.notnull()]
    
    rohbedfile = "results/rohcov/%s.bed" % outprefix
    print "Writing file:",rohbedfile
    rohannot.to_csv(rohbedfile, columns=["CHR","POS1","POS2","IID"], sep="\t",index=False,header=False)
    rohbedsort = sortBed( rohbedfile ) 
                
    rohbases = calculateRegionsCovered( rohbedsort )
                
    #exomeregions = "resources/regionbeds/refseqgenes_b37_sorted.bed"
    # Maybe replace this with well covered regions
    exomeregions = "resources/regionbeds/refseq_condensed.bed"
    exometotals = calculateRegionsCovered( exomeregions )
    
    #sort -V -k1,1 -k2,2n refseqgenes_b37.bed
    #sort -V -k1,1 -k2,2n in.bed
    allbases = {}
    bedfiles = {}
    limit = -1
    for samp, data in rohannot.groupby( "IID" ) :
        #print "Processing:",samp
        intersectbed = "%s/%s_exomeroh.bed" % (targetdir, samp)
        if not os.path.exists(intersectbed) :
            tmprohbedfile = "%s/%s_sort.bed" % (targetdir,samp)
            retval = makeIntersectBed( data, exomeregions, intersectbed, tmprohbedfile )
            if retval == "Failed" : continue
        totalbases = calculateRegionsCovered( intersectbed )
        #print "Total bases:",totalbases
        bedfiles[samp] = intersectbed
        allbases[samp] = totalbases
        limit -= 1
        if limit == 0 : break
                
    allbases = DataFrame( allbases ).transpose().fillna(0).astype(int).reset_index()
    allbases.rename(columns={'index':'Individual.ID'}, inplace=True)
    print "Allbases:",allbases.head(10)
    
    print "Writing file:",overlapfile
    foundcols = [x for x in exometotals.index if x in allbases.columns]
    print foundcols
    samplesums = allbases[foundcols].fillna(0).apply(sum,axis=1)

    print "Samplesums:",samplesums.head(10)
    exomesum = float(exometotals[foundcols].sum())
    print "Exomesum:",exomesum
    allbases["percent"] = (samplesums / exomesum *100)
    print "Allbases:",allbases.head(10)

    allbases.to_csv( overlapfile, sep="\t", index=False)
    rohfig.percentAnalysis( allbases, exometotals, outprefix )
    rohfig.percentAnalysisRegions( allbases, exometotals, sampleannot, outprefix )
    return allbases, bedfiles
# END exomeCoverage

################################################################################
# retrieveBedFiles
################################################################################
def retrieveBedFiles( outprefix, targetdir="./results/rohcov/",filetype="exomeroh.bed" ):
    filelist = glob.glob( os.path.join( targetdir,outprefix,"*"+filetype))
    bedfiles = {}
    for bfile in filelist :
        filepath, filename = os.path.split(bfile)
        sampname = filename[:-len(filetype)-1]
        bedfiles[sampname] = bfile

    return bedfiles
# END retrieveBedFiles

def summarizeOverlap( prefix, sampleannot, bedfiles ):
    popfiles = []
    for pop, data in sampleannot.groupby("Continent") :
        print "Pop:", pop, len(data)
        #print data["Individual.ID"].tolist()
        popregionsfile = "./results/rohcov/%s_%s_rohcounts.bed" % (prefix.replace(" ","_"),pop.replace(" ","_"))
        popregions = {}
        foundcount = 0
        for samp in data["Individual.ID"].tolist() :
            if bedfiles.has_key(samp) :
                #print "Found file", samp
                foundcount +=1
                for row in csv.reader(open(bedfiles[samp]),delimiter="\t") :
                    exonid = "\t".join(row[0:3])
                    if not popregions.has_key(exonid) : popregions[exonid] = 0
                    popregions[exonid] += 1
            else : print "Error sample not found:", samp
        flattened = []
        popfiles.append([pop, len(data), foundcount, popregionsfile])
        for exon in sorted(popregions) :
            if len(exon.strip()) == 0 : continue
            flattened.append(exon.split("\t")+[popregions[exon]])
        popregions = (DataFrame(flattened,columns=["chrom","start","end","count"])
                      .convert_objects(convert_numeric=True)
                      .sort(["chrom","start","end"]))
        popregions.to_csv(popregionsfile, header=False, index=False, sep="\t")
    popfiles =DataFrame(popfiles,columns=["Pop","len","found","file"])
    return popfiles
# END summarizeOverlap


################################################################################
# readAllRohRegions
################################################################################
def readAllRohRegions( bedfiles ):
    allregions = []
    for pat in bedfiles : 
        alldata = read_csv( bedfiles[pat], sep="\t", header=None, 
                           names=["chrom","start","end","sample"] )
        allregions.append(alldata)
    allregions = concat(allregions)
    return allregions
# END readAllRohRegions

################################################################################
# splitVcf
################################################################################
def splitVcf( vcffile, targetdir, force=False ) :
    filepath, basename, ext = hk.getBasename( vcffile )
    #print filepath
    #targetdir = "%s/calls" % filepath

    sampfiledesc =  "%s/calls/%s_sampfiles.tsv" % (filepath, basename)
    if os.path.exists(sampfiledesc) and not force :
        sampfiles = read_csv(sampfiledesc,sep="\t")
        sampfiles.index = sampfiles["Sample"]
        return sampfiles

    #vcfread = csv.reader( open( vcffile ),delimiter="\t" )
    vcffh = gzip.open( vcffile )#, "rU"
    header, comments = [], []
    while len(header) == 0 : 
        line = vcffh.next().strip()
        if line[:2] == "##" : comments.append(line)
        else : header = line[1:].split("\t")

    allcallsfile = os.path.join( targetdir, "%s_calls.tsv"%basename)
    outfh = open( allcallsfile, "wb" ) 
    outfh.write("#chrom\tstart\tend\tsample\tcall\n")
    for row in csv.reader( vcffh, delimiter="\t" ) :
        #print len(row), len(header)
        s = Series( data=row, index=header )
        end = int(s["POS"]) + len(s["REF"])
        for samp in header[9:] :
            if s[samp][:3] in ["1/0","0/1","1/1"] : 
                fline = ("%s\t%s\t%d\t%s\t%s\n" 
                         % (s["CHROM"],s["POS"],end,samp,s[samp]))
                outfh.write(fline)

    outfh.close()

    sampfiles = []
    alldata = read_csv( allcallsfile, delimiter="\t" )
    for samp, subdata in alldata.groupby("sample") :
        sfile = os.path.join( targetdir, "%s_calls.tsv"%samp)
        print "Samp:",samp,sfile
        sampfiles.append(Series(data=[samp,sfile],index=["Sample","Callfile"]))
        subdata.to_csv( sfile, sep="\t", index=False )
    sampfiles = DataFrame(sampfiles)

    sampfiles.to_csv(sampfiledesc, sep="\t", index=False)
    sampfiles.index = sampfiles["Sample"]
    return sampfiles
# END splitVcf
    #sampfh = {}
    #for samp in header[9:] : 
        #sfile = os.path.join( targetdir, "%s_calls.tsv"%samp)
        #sampfiles.append(Series(data=[samp,sfile],index=["Sample","Callfile"]))
        #sampfh[samp] = open( sfile, "wb" )
    #sampfiles = DataFrame(sampfiles)

################################################################################
# rohIntersectCalls
################################################################################
def rohIntersectCalls( tcalls, tregions ) :
    print "Running rohIntersectCalls"
    print tcalls.shape
    print tregions.shape
    tcallsgrps = tcalls.groupby("chrom")
    #print tcalls.pos >= tregions["start"]
    #print tcalls.pos >= tregions["end"]
    burden = [0 for i in range(0,len(tregions))]
    for chrom, chrregions in tregions.groupby("chrom")  :
        if chrom not in tcallsgrps.groups.keys() : continue
        chromcalls = tcallsgrps.get_group(chrom)
        for idx, region in chrregions.iterrows() :
            #print "Region:",region["chrom"],region["start"],region["end"]
            #tregions.ix[index]["burden"]
            hits = chromcalls[ (chromcalls.pos >= region["start"]) &
                            (chromcalls.pos <= region["end"] )]
            burden[idx] = len(hits)
            #if len(hits) > 0 : print hits.head()

    tregions["burden"] = burden
    calldf = tregions[tregions.burden > 0]
    print "rohblocks",calldf.shape
    return calldf
    #tmpcallsfile = "./results/sample_calls.bed"
    #tmpregionsfile = "./results/roh_regions.bed"
    #tcalls.to_csv( tmpcallsfile ,sep="\t",header=False,index=False, 
                #cols=["chrom","start","end","call"])
    #tregions["rohid"] = ["roh%d"%x for x in range(1,len(tregions)+1)]
    #tregions.to_csv( tmpregionsfile ,sep="\t",header=False,index=False)
    #a = pybedtools.BedTool( tmpregionsfile )
    #b = pybedtools.BedTool( tmpcallsfile )
    #a_and_b = a.intersect(b)
    #print "A and B head:", a_and_b.head()
    #if len(a_and_b) == 0 : 
        #print "Warning: empty dataframe!"
        #return DataFrame(tregions.columns)
    #calldf = DataFrame(data=[list(x) for x in a_and_b],columns=tregions.columns)
    #print "Tmp:",calldf.head()
    #rohburden = DataFrame({'burden':calldf.groupby("rohid").size()}).reset_index()
    #print rohburden.burden.max()
    #return calldf
# END rohIntersectCalls
        
################################################################################
# rohVarOverlap
################################################################################
def rohVarOverlap( rohfile, genesfile, vcffile, allregions, force=False ) :
    print "Running rohVarOverlap"

    filepath, basename, ext = hk.getBasename( vcffile )
    targetdir = "%s/calls" % filepath
    hk.makeDir( targetdir )

    varcountsfile = "%s/%s_varcounts.tsv" % (targetdir, basename)
    if os.path.exists(varcountsfile) and not force :
        print "Warning: final file already exists", varcountsfile
        varcounts = read_csv(varcountsfile,sep="\t")
        return varcounts

    vardata = read_csv( genesfile, sep="\t")
                           #dtype={'chrom':np.int32,'pos':np.int32})
    sampfiles = splitVcf( vcffile, targetdir )

    print sampfiles
    allregions = allregions[allregions.chrom.notnull()]

    allregions[["chrom","start","end"]] = allregions[["chrom","start","end"]].astype(int)
    allregions["RohLen"] = allregions["end"] - allregions["start"]

    delvars = vardata[vardata.vclass.isin([5,4,3])] # filter for deleterious vars
    delvars.index = delvars["chrom"].map(str)+":"+delvars["pos"].map(str)
    print "Delvars:",delvars.shape

    benignvars = vardata[vardata.vclass.isin([1,0,-1])] # filter for deleterious vars
    benignvars.index = benignvars["chrom"].map(str)+":"+benignvars["pos"].map(str)
    print "benignvars:",benignvars.shape
    
    varcounts = {'Individual.ID':[], 'delvars':[], 'delroh':[], 
                 'benignvars':[],'benignroh':[]} #,'longroh':[],'shortroh':[]
    #rohburdens = {'sample':[], 'Rohlen':[], 'delburden':[]}
    for samp in sampfiles.index :
        print "Processing sample:",samp#,sampfiles.ix[samp]
        s_calls = read_csv( sampfiles.ix[samp]["Callfile"], sep="\t", 
                           dtype={'start':int,'end':int})
        s_calls[s_calls["#chrom"] == "X"] = "23"
        s_calls["#chrom"] = s_calls["#chrom"].astype(int)
        s_calls.index = s_calls["#chrom"].map(str)+":"+s_calls["start"].map(str)
        print "s_calls Before shape:",s_calls.shape

        print "Vardata shape:",delvars.shape
        #print "Delvars","$"*60
        #print delvars.head()
        #print "S_calls","$"*60
        #print s_calls.head()
        delcalls = (merge(delvars,s_calls,left_index=True, right_index=True)
                    [["chrom","pos"]].drop_duplicates())

        benigncalls = (merge(benignvars,s_calls,left_index=True, right_index=True)
                       [["chrom","pos"]].drop_duplicates())
        s_regions = allregions[allregions.sample == samp]
        if len(s_regions) == 0 : 
            print "Warning: sample has no regions!:",samp
            continue
        print "Dels:", delcalls.shape
        print "Benign:", benigncalls.shape
        #print "S regions head:",s_regions.head(10)
        #print "S regions shape:", s_regions.shape

        print "Making benign overlaps"
        benignrohcounts = rohIntersectCalls( benigncalls, s_regions ) 
        print "Making del overlaps"
        delrohcounts = rohIntersectCalls( delcalls, s_regions ) 

        varcounts['Individual.ID'].append(samp)
        varcounts['delvars'].append(len(delcalls))
        varcounts['delroh'].append(delrohcounts.burden.sum())
        varcounts['benignvars'].append(len(benigncalls))
        varcounts['benignroh'].append(benignrohcounts.burden.sum())
        print "Intersect variants with ROH regions"

    varcounts = DataFrame(varcounts)

    print varcounts.head(10)
    varcounts.to_csv( varcountsfile, sep="\t", index=False )
    return varcounts
# END rohVarOverlap

def plotVariantSlopes( varcounts, dataset ):
    vcounts = melt( varcounts_merge, 
                   id_vars=["Individual.ID","percent"], 
                   value_vars=["benignroh","benignnonroh","delroh","delnonroh"])
    vcounts["rohclass"] = ["Outside" if x.find("non") >= 0 else "Within" for x in vcounts.variable]
    vcounts["vtype"] = ["Benign" if x.find("benign") == 0 else "Deleterious" for x in vcounts.variable]

    vcounts_annot = addSampleAnnotation( vcounts, mergecol="Individual.ID" )
    vcounts_annot = vcounts_annot[ vcounts_annot.Continent.isin(["Europe","Middle East"]) ]
    print vcounts_annot.head()
    lmdata = []
    for idx, grp in vcounts_annot.groupby(["Continent","rohclass","vtype"]) :  
        lm = stats.linregress( grp["percent"].tolist(), grp["value"].tolist() )
        lmdata.append( Series(data=idx+lm, index=["Continent","rohclass","vtype","slope", 
                                                 "intercept", "rval", "pval", "stderr"]) )

    lmdata = DataFrame(lmdata)

    lmdata["xpos"] = 15

    benignMaxY = vcounts[vcounts.vtype == "Benign"]["value"].max()
    delMaxY = vcounts[vcounts.vtype == "Deleterious"]["value"].max()
    #lmdata["ypos"] = maxY
    for idx, subdata in lmdata.groupby("vtype"): 
        maxY = vcounts[vcounts.vtype==idx]["value"].max()
        lmdata.loc[lmdata.vtype==idx,"ypos"] = [maxY*.95 if x == "Outside" 
                                                 else maxY *.87 
                                                 for x in subdata.rohclass]

    lmdata["Slope"] = ["%s Slope: %.3f" % (x,y) 
                       for x,y in lmdata[["rohclass","slope"]].values]


    r_dataframe = com.convert_to_r_dataframe(vcounts_annot)
    r_lm = com.convert_to_r_dataframe(lmdata)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="percent",
                                   y="value", 
                                   group="factor(rohclass)",
                                   colour="factor(rohclass)") +
                ggplot2.geom_point() +
                ggplot2.geom_abline(ggplot2.aes_string(intercept="intercept",slope="slope"), 
                                    color="black", data=r_lm) +
                ggplot2.geom_text(ggplot2.aes_string(label="Slope",x="xpos",y="ypos"), 
                                  size=4, hjust=0, data=r_lm) +
                #ggplot2.ggtitle("Benign Variant counts within\nRuns of Homozygosity") +
                ggplot2.scale_y_continuous("Count of Benign variants") +
                ggplot2.scale_x_continuous("Percent of Genome in RoH") +
                ggplot2.scale_colour_manual(name="Variant Location",
                                            values=robjects.StrVector(["red","blue"]), 
                                            breaks=robjects.StrVector(["benignroh","benignnonroh"]), 
                                            labels=robjects.StrVector(["Within RoH","Outside RoH"]) ) +
                ggplot2.facet_grid(robjects.Formula('vtype ~ Continent'), scale="free_y") +
                #ggplot2.stat_smooth(method="lm", se=False) +
                ggplot2.theme(**pointtheme) )
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 90)}) +

    figname = "results/figures/roh/%s_roh_varcounts.png" % (dataset)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotVariantSlopes

################################################################################
# Main
################################################################################
if __name__ == "__main__" :

    os.chdir("..")
    optlist, args = getopt.getopt( sys.argv[1:], "r:ot")
    optlist = dict(optlist)

    #dataset = optlist.get("-r","test")
    dataset = optlist.get("-r",None)
    if dataset == "test" :
        rohfile = "./rawdata/test/main/everything_set1.chr1.snp.clean.plink.filt.hom"
        #genesfile = "./rawdata/test/main/everything_set1.chr1.snp.clean_genes.tsv"
        genesfile = "./rawdata/test/main/everything_set1.chr1.snp_genes.tsv"
        vcffile = "./rawdata/test/main/everything_set1.chr1.snp.clean.vcf.gz"
    elif dataset == "test2" :
        rohfile = "./rawdata/test2/main/test2.clean.plink.filt.hom"
        genesfile = "./rawdata/test2/main/test2.clean_genes.tsv"
        vcffile = "./rawdata/merged/main/test2.clean.vcf.gz"
    elif dataset == "merged" :
        rohfile = "./rawdata/merged/main/merged.clean.plink.filt.hom"
        genesfile = "./rawdata/merged/main/merged.clean_genes.tsv"
        vcffile = "./rawdata/merged/main/merged.clean.vcf.gz"
    elif dataset == "mevariome" :
        rohfile = "./rawdata/mevariome/main/variome.clean.plink.filt.hom"
        genesfile = "./rawdata/mevariome/main/variome.clean_genes.tsv"
        vcffile = "./rawdata/mevariome/main/variome.clean.vcf.gz"

    #elif dataset == "variome" :
        #rohfile = "./rawdata/variome/main/variome.clean.plink.filt.hom"
        #genesfile = "./rawdata/variome/variome.clean_genes.tsv"
        #vcffile = "./rawdata/variome/main/variome.clean.vcf.gz"
    #elif dataset == "variome1" :
        #rohfile = "./rawdata/variome1/admixture/variome.clean.plink.filt.hom"
        #genesfile = "./rawdata/variome/variome.clean_genes.tsv"
        #vcffile = "./rawdata/variome/main/variome.clean.vcf.gz"
    else :
        print "Error: dataset not recognized -",dataset
        sys.exit(1)

    #rohfile = "./rawdata/merged/merged.clean.plink.filt.hom"

    filepath, basename, ext = hk.getBasename( vcffile )
    targetdir = "%s/calls" % filepath

    rohdata = read_csv(rohfile, delim_whitespace=True )
    targetsamples = rohdata.IID.unique().tolist()
    sampleannot = sampleAnnotation(targetsamples)

    allbases, bedfiles = exomeCoverage( rohfile, sampleannot, 
                                       outprefix=dataset, force=False )
    print allbases.head(10)

    allregions = readAllRohRegions( bedfiles ) 

    # sample, delvars, delroh, benignvars, benignroh, delnonroh
    varcounts = rohVarOverlap( rohfile, genesfile, vcffile, allregions ) 

    varcounts["delnonroh"] = varcounts["delvars"] - varcounts["delroh"]
    varcounts["benignnonroh"] = varcounts["benignvars"] - varcounts["benignroh"]

    varcounts_merge = merge( varcounts, 
                            allbases[["Individual.ID","percent"]], 
                            on="Individual.ID" )
    plotVariantSlopes( varcounts_merge, dataset )

    # plot it
    #vcounts = melt( varcounts_merge, 
                   #id_vars=["Individual.ID","percent"], 
                   #value_vars=["delroh","delnonroh"])

    #vcounts_annot = addSampleAnnotation( vcounts, mergecol="Individual.ID" )
    #vcounts_annot = vcounts_annot[ vcounts_annot.Continent.isin(["Europe","Middle East"]) ]
    #r_dataframe = com.convert_to_r_dataframe(vcounts_annot)
    #p = (ggplot2.ggplot(r_dataframe) +
                #ggplot2.aes_string(x="percent",
                                   #y="value", 
                                   #group="factor(variable)",
                                   #colour="factor(variable)") +
                #ggplot2.geom_point() +
                #ggplot2.ggtitle("Variant counts within\nRuns of Homozygosity") +
                #ggplot2.scale_y_continuous("Count of Deleterious variants") +
                #ggplot2.scale_x_continuous("Percent of Genome in RoH") +
                #ggplot2.scale_colour_manual(name="Location",
                                            #values=robjects.StrVector(["red","blue"]), 
                                            #breaks=robjects.StrVector(["delroh","delnonroh"]), 
                                            #labels=robjects.StrVector(["Within RoH","Outside RoH"]) ) +
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                #ggplot2.facet_grid(robjects.Formula('Continent ~ .')) +
                #ggplot2.stat_smooth(method="lm", se=False) +
                #ggplot2.theme(**mytheme) )
    #figname = "results/figures/roh/%s_roh_varcounts.png" % (dataset)
    #print "Writing file:",figname
    #grdevices.png(figname)
    #p.plot()
    #grdevices.dev_off()


    #slope : float #slope of the regression line
    #intercept : float #intercept of the regression line
    #r-value : float #correlation coefficient
    #p-value : float #two-sided p-value for a hypothesis test whose null hypothesis is
        #that the slope is zero.
    #stderr : float #Standard error of the estimate
    # BENIGN Variants

    #rohfig.makeIdeograms( dataset, bedfiles, sampleannot )

    #intersectbed = "results/rohcov/%s_%s_exomeroh.bed" % (outprefix, samp)

    # Make summarized rohfile file for plotting
    # run Rscript plotting software

# END MAIN
    #lmdata[lmdata["variable"] == "benignnonroh"]["ypos"] = .9 * maxY
    #lmdata[lmdata["variable"] == "benignroh"]["ypos"] = .85 * maxY
    #for idx, row in lmdata.iterrows() :
        #if row["variable"].find("non") >= 0: 
            #rohclass = "Outside RoH"
            #lmdata.loc[idx,"ypos"] = maxY * .95
        #else : 
            #rohclass = "Within RoH"
            #lmdata.loc[idx,"ypos"] = maxY * .85
        #lmdata.loc[idx,"Slope"]  = "%s Slope: %.3f" % (rohclass,row["slope"])


    #rohfile = "rawdata/variome1/roh/variome.clean.plink.filt.hom"
    #rohfile = "./rawdata/onekg/old2/roh/onekg.clean.plink.filt.hom"
    #rohfile = "./rawdata/variome1/roh/variome.clean.plink.filt.hom"
    #rohfile = "./rawdata/variome/main/variome.clean.plink.filt.hom"
    #genesfile = "./rawdata/variome/variome.clean_genes.tsv"
    #vcffile = "./rawdata/variome/variome.clean.vcf.gz"
