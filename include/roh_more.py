#!/usr/bin/python

import os, sys
import subprocess
import glob
import pybedtools
#import resource
#resource.setrlimit(resource.RLIMIT_NOFILE, (1000,-1))

from scipy.stats import binom_test 

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
def exomeCoverage( rohfile, sampleannot, targetdir, basename="test", force=False ):
    print "Running exomeCoverage", rohfile

    newtargetdir = os.path.join(targetdir, "exomeint" )
    hk.makeDir(newtargetdir)

    overlapfile = os.path.join(targetdir, basename+"_chrom_overlap.txt" )
    if os.path.exists(overlapfile) and not force :
        print "Final file already exists:",overlapfile
        bedfiles = retrieveBedFiles(newtargetdir)
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
    
    rohbedfile = os.path.join(targetdir,prefix+".bed")
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
        intersectbed = "%s/%s_exomeroh.bed" % (newtargetdir, samp)
        if not os.path.exists(intersectbed) :
            tmprohbedfile = "%s/%s_sort.bed" % (newtargetdir,samp)
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
    rohfig.percentAnalysis( allbases, exometotals, basename )
    rohfig.percentAnalysisRegions( allbases, exometotals, sampleannot, basename )
    return allbases, bedfiles
# END exomeCoverage

################################################################################
def retrieveBedFiles( targetdir="./results/rohcov/",filetype="exomeroh.bed" ):
    #print "Running retrieveBedFiles", targetdir
    filelist = glob.glob( os.path.join( targetdir,"*"+filetype))
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
# readExomeRohRegions
################################################################################
def readExomeRohRegions( bedfiles ):
    allregions = []
    for pat in bedfiles : 
        alldata = read_csv( bedfiles[pat], sep="\t", header=None, 
                           names=["chrom","start","end","sample"] )
        allregions.append(alldata)
    allregions = concat(allregions)
    allregions = allregions[allregions.chrom.notnull()]

    allregions[["chrom","start","end"]] = allregions[["chrom","start","end"]].astype(int)
    allregions["RohLen"] = allregions["end"] - allregions["start"]
   
    return allregions
# END readExomeRohRegions

################################################################################
# splitVcf
################################################################################
def splitVcf( vcffile, targetdir, force=False ) :
    print "Running splitVcf - file:", vcffile
    filepath, basename, ext = hk.getBasename( vcffile )
    #print filepath
    #targetdir = "%s/calls" % filepath

    sampfiledesc =  "%s/calls/%s_sampfiles.tsv" % (filepath, basename)
    if os.path.exists(sampfiledesc) and not force :
        print "Already exists:",sampfiledesc
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

    print "Header:",header[:10]

    newtargetdir = os.path.join( targetdir,"scalls" )
    hk.makeDir( newtargetdir )
    allcallsfile = os.path.join( newtargetdir, "%s_calls.tsv"%basename)
    outfh = open( allcallsfile, "wb" ) 
    outfh.write("#chrom\tstart\tend\tsample\tcall\n")
    indelcnt = 0
    ##for row in csv.reader( vcffh, delimiter="\t" ) :
    for line in vcffh : 
        if len(line) > 40000 : 
            indelcnt += 1
            print "Removing indel #",indelcnt
            continue
        row = line.split("\t")
        #print line, len(row), row[:10]
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
        sfile = os.path.join( newtargetdir, "%s_calls.tsv"%samp)
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

def splitRange( x, firstname="start", secondname="end" ) :
    return Series({firstname:int(x[1:x.find(",")]),
                   secondname:int(x[x.find(",")+2:len(x)-1])})
# ENd splitRange

################################################################################
def rohIntersectCalls2( tcalls, tregions ) :
    print "Running rohIntersectCalls2"
    results = []
    tcallsgrps = tcalls.groupby("chrom")
    for chrom, chrregions in tregions.groupby("chrom")  :
        chromcalls = tcallsgrps.get_group(chrom)
        chrroh = np.unique(np.sort(np.concatenate(([1],chrregions["start"],
                                                   chrregions["end"]+1,
                                                   [chromcalls["pos"].max()]))))
        print "positions:",len(chrroh)
        chrcounts = chromcalls.groupby( cut( chromcalls["pos"], chrroh ) ).size().reset_index()
        chrcounts["chrom"] = chrom
        chrcounts = chrcounts.join( chrcounts['pos'].apply(splitRange) )
        chrcounts["Autozygosity"] = chrcounts["start"].isin(chrregions["start"])
        results.append(chrcounts)

    results = concat(results)
    results.fillna(0, inplace=True)
    del results["pos"]
    results.rename(columns={0:'burden'}, inplace=True )
    #print results.head()
    #print len(results)
    return results
# END rohIntersectCalls2

################################################################################
def rohIntersectCallsClass( tcalls, tregions  ) :
    print "Running rohIntersectCallsClass"
    results = []
    tcallsgrps = tcalls.groupby("chrom")
    #print tcallsgrps.size()
    for chrom, chrregions in tregions.groupby("chrom")  :
        chromcalls = tcallsgrps.get_group(chrom)
        chrroh = np.unique(np.sort(np.concatenate(([1],chrregions["start"],
                                                   chrregions["end"]+1,
                                                   [chromcalls["pos"].max()]))))
        #print "positions:",len(chrroh)
        chrcounts = (chromcalls.groupby( [cut( chromcalls["pos"], chrroh ), "vclass"] )
                     .size().reset_index())
        chrcounts["vclass"] = "class"+chrcounts["vclass"].map(str)
        #print chrcounts.head()
        chrcounts = chrcounts.pivot( 'pos', 'vclass', 0 ).reset_index()
        chrcounts["chrom"] = chrom
        chrcounts = chrcounts.join( chrcounts['pos'].apply(splitRange) )
        chrcounts["Autozygosity"] = chrcounts["start"].isin(chrregions["start"])
        #print chrcounts.head()
        results.append(chrcounts)

    results = concat(results)
    results.fillna(0, inplace=True)
    del results["pos"]
    return results
# END rohIntersectCallsClass

################################################################################
def rohIntersectCallsAF( tcalls, tregions ) :
    print "Running rohIntersectCallsAF"
    results = []
    tcallsgrps = tcalls.groupby("chrom")
    afbins=[0.,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]

    for chrom, chrregions in tregions.groupby("chrom")  :
        chromcalls = tcallsgrps.get_group(chrom)
        chrroh = np.unique(np.sort(np.concatenate(([1],chrregions["start"],
                                                   chrregions["end"]+1,
                                                   [chromcalls["pos"].max()]))))
        #print "positions:",len(chrroh)
        #print chromcalls.groupby([cut( chromcalls["AF"], afbins)]).size().reset_index()
        #for idx, cdata in chromcalls.groupby([cut( chromcalls["AF"], afbins)]) :
            #print idx
            #print cdata.head()
        chrcounts = (chromcalls.groupby( [cut( chromcalls["pos"], chrroh ), 
                                          cut( chromcalls["AF"], afbins), "vclass"] )
                     .size().reset_index())
        chrcounts["vclass"] = "class"+chrcounts["vclass"].map(str)
        chrcounts = pivot_table(chrcounts, index=['pos','AF'], 
                                columns='vclass', values=0 ).reset_index()
        chrcounts["chrom"] = chrom
        chrcounts = chrcounts.join( chrcounts['pos'].apply(splitRange) )
        chrcounts["Autozygosity"] = chrcounts["start"].isin(chrregions["start"])
        results.append(chrcounts)

    results = concat(results)
    results.fillna(0, inplace=True)
    del results["pos"]
    #print results.head()
    #print len(results)
    return results
# END rohIntersectCallsAF

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

    tregions.loc[:,"burden"] = burden
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
        
def calcAF(genecounts):
    #print [x for x in genecounts if x.find(":") == -1]
    ref,het,hom = [int(x) for x in genecounts.split(":")]
    if ref+het+hom == 0 : print "Err:",genecounts; return 0
    return (2.*hom+het) / (2.*(ref+het+hom))
# END calcAF

def prepareGenesFile( genesfile ):

    vardata = read_csv( genesfile, sep="\t" )

    tregions = ["Africa","America","East Asia","Europe","Middle East","Anglo-American","South Asia"]
    tregions = [x for x in tregions if x in vardata.columns]
    keepcols = ["chrom","pos"] + tregions
    vardata_filt = vardata.groupby(keepcols)["vclass"].max().reset_index()
    vardata = melt( vardata_filt, id_vars=["chrom","pos","vclass"], value_vars=tregions )
    vardata.rename(columns={"variable":"Region", "value":"Vcounts"},inplace=True)
    vardata["AF"] = [calcAF(x) for x in vardata["Vcounts"].tolist()]

    vardata = vardata[(vardata.AF > 0.0) & (vardata.AF < 1.0)]

    delvars = vardata[vardata.vclass.isin([5,4,3])] # filter for deleterious vars
    delvars.index = delvars["chrom"].map(str)+":"+delvars["pos"].map(str)
    #print "Delvars:",delvars.shape

    benignvars = vardata[vardata.vclass == 1] 
    benignvars.index = benignvars["chrom"].map(str)+":"+benignvars["pos"].map(str)
    #print "benignvars:",benignvars.shape
 
    return vardata.groupby("Region"), delvars.groupby("Region"), benignvars.groupby("Region")
# END prepareGenesFile

def makeSampleCallsFiles( sampfiles, groupedvardata, sampleannot, targetdir, force=False ) :
    print "Running makeSampleCallsFiles"
    sampfiles = merge(sampfiles, sampleannot[["Continent2","GeographicRegions3"]], 
                      left_index=True, right_index=True )
    sampfiles["Annotcall"] = "none"

    newtargetdir = os.path.join(targetdir, "annotcalls" )
    hk.makeDir(newtargetdir)

    sampnum = 1
    for samp,srow in sampfiles.iterrows() :
        sampnum += 1
        targetfile = os.path.join(newtargetdir,samp+"_annotcalls.tsv")
        if os.path.exists(targetfile) and not force : 
            #print "Warning: file exists -",targetfile; 
            sampfiles.loc[samp,"Annotcall"] = targetfile
            continue

        sregion = srow["Continent2"]
        if sregion not in groupedvardata.groups.keys() : 
            #print "Error: Unknown region:",sregion; 
            continue

        print "Processing sample #",sampnum,":",samp 
        s_calls = read_csv( srow["Callfile"], sep="\t", 
                           dtype={'start':int,'end':int})

        if "X" in s_calls["#chrom"].unique().tolist() : 
            s_calls.loc[s_calls["#chrom"] == "X","#chrom"] = "23"

        #print "Unique chromosomes"
        #print ":".join([str(x) for x in s_calls["#chrom"].unique().tolist()])

        s_calls["#chrom"] = s_calls["#chrom"].astype(int)
        s_calls.index = s_calls["#chrom"].map(str)+":"+s_calls["start"].map(str)
        #print "s_calls Before shape:",s_calls.shape
    
        vardata = groupedvardata.get_group( sregion )
        if "X" in vardata["chrom"].unique().tolist() : 
            vardata.loc[vardata["chrom"] == "X","chrom"] = "23"
        vardata["chrom"] = vardata["chrom"].astype(int)

        #sampcalls = merge(vardata,s_calls,left_index=True, right_index=True)
        sampcalls = merge(vardata,s_calls,left_on=["chrom","pos"], right_on=["#chrom","start"])
        #print "Unique sampcalls chromosomes"
        #print sampcalls.groupby(["chrom"]).size()
        #print "Writing file:",targetfile
        sampcalls.to_csv( targetfile, sep="\t", index=False )
        sampfiles.loc[samp, "Annotcall"] = targetfile

    print "No Annotation files:",sampfiles[sampfiles["Annotcall"] == "none"].head()
    print "Len:",len(sampfiles[sampfiles["Annotcall"] == "none"])
    sampfiles = sampfiles[sampfiles["Annotcall"] != "none"]
    return sampfiles
# END makeSampleCallsFiles
        #print "Unique vardata chromosomes"
        #print ":".join([str(x) for x in s_calls["#chrom"].unique().tolist()])
        #print vardata.groupby(["chrom"]).size()
        #print "Vardata"
        #print vardata[vardata["chrom"] == 22].head()
        #print "S_calls"
        #print s_calls[s_calls["#chrom"] == 22].head()

def plotIndvAFDistribution( s_calls, samp, targetdir="results/figures/roh/indiv" ) :
    afbins=[0.,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]
    afsum = s_calls.groupby([cut(s_calls["AF"], afbins),"vclass"]).size().reset_index()
    afsum.rename(columns={0:"burden"}, inplace=True )

    afsum["AF_bin"] = [float(x[x.find(",")+2:len(x)-1]) for x in afsum["AF"]]

    totals = afsum.groupby(["vclass"])["burden"].sum().reset_index()
    print totals.head()

    totals.rename(columns={"burden":"total"}, inplace=True )
    afsum = merge(afsum,totals, on=["vclass"] )

    print afsum.head()
    afsum["proportions"] = afsum["burden"] / afsum["total"]

    r_dataframe = com.convert_to_r_dataframe(afsum)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(AF_bin)", y="proportions") +
                ggplot2.geom_bar(stat="identity",position="dodge") +
                ggplot2.geom_hline(yintercept=0) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45),
                                 'legend.position':'top'}) +
                ggplot2.scale_y_continuous("Proportion",
                                           expand=robjects.IntVector((0,0))) +
                ggplot2.scale_x_discrete("Allele Frequency") +
                                           #expand=robjects.IntVector((0,0))) +
                #ggplot2.scale_fill_brewer( "Autozygosity", palette="Set1" ) +
                ggplot2.facet_grid( robjects.Formula('vclass ~ .'), scale="free_y" ) +
                ggplot2.theme(**mytheme) )

    figurename = os.path.join(targetdir,samp+"_sf.png" )
    print "Making figure:",figurename
    grdevices.png(figurename) 
    p.plot()
    grdevices.dev_off()
# END plotIndvAFDistribution

def fractionHomozygotes( genecounts ) :
    ref,het,hom = [float(x) for x in genecounts.split(":")]
    if ref+het+hom == 0 : print "Err:",genecounts; return 0
    return float(hom) / (ref+het+hom)
# END fractionHomozygotes

def rohIntersectVars( tcalls, tregions ) :
    print "Running rohIntersectVars"
    #print tcalls.head()

    results = []
    tcalls["hfraction"] = [fractionHomozygotes( x ) for x in tcalls.Vcounts]
    tcalls["geno"] = [x[:3] for x in tcalls.call]

    print hk.dfTable( tcalls.geno )
    tcalls["gt"] = ["het" if x == "0/1" else "hom" for x in tcalls.geno]
    tcallsgrps = tcalls[["chrom","pos","vclass","Region","AF","gt","hfraction"]].groupby("chrom")
    #print tcallsgrps.size()

    for chrom, chrregions in tregions.groupby("chrom")  :
        chromcalls = tcallsgrps.get_group(chrom)
        chrroh = np.unique(np.sort(np.concatenate(([1],chrregions["start"],
                                                   chrregions["end"]+1,
                                                   [chromcalls["pos"].max()]))))
        #print "positions:",len(chrroh)
        chromcalls["rohbin"] = cut( chromcalls["pos"], chrroh )
        binpositions = DataFrame([ concat([Series({'rohbin':x}),splitRange(x)]) 
                                  for x in chromcalls['rohbin'].unique() ])

        chromcalls = merge( chromcalls, binpositions, on="rohbin" )
        chromcalls["Autozygosity"] = ["In" if x else "Out" 
                                      for x in chromcalls["start"].isin(chrregions["start"])]

        chromcalls["Autozygosity"] = chromcalls["Autozygosity"] +"_"+ chromcalls["gt"]

        results.append(chromcalls[["chrom","pos","vclass","Region","AF","hfraction","Autozygosity"]])

    results = concat(results)
    results.fillna(0, inplace=True)
    return results
# END rohIntersectVars

def makeCountsFiles( sampfiles, targetdir, rohregions, classcountsfile, afcountsfile, 
                    varrohfile, force=False ) :
    print "Running makeCountsFiles"
    if (os.path.exists(classcountsfile) and os.path.exists(afcountsfile) 
        and os.path.exists(varrohfile) and not force ):
        classcounts = read_csv( classcountsfile, sep="\t" )
        afcounts = read_csv( afcountsfile, sep="\t" )
        allvrohcounts = read_csv( varrohfile, sep="\t" )
        return classcounts, afcounts, allvrohcounts
 
    classcounts = []
    afcounts = []
    allvrohcounts = []

    rohcallsdir = os.path.join( targetdir, "rohcalls" )
    hk.makeDir( rohcallsdir )

    rohafdir = os.path.join( targetdir, "rohaf" )
    hk.makeDir( rohafdir )

    varrohdir = os.path.join( targetdir, "varroh" )
    hk.makeDir( varrohdir )

    sampnum = 1
    for samp,srow in sampfiles.iterrows() :
        print "Processing sample #",sampnum,":",samp 
        sampnum += 1
        targetfile = os.path.join(rohcallsdir,samp+"_rohcallcounts.tsv")
        s_calls = read_csv( srow["Annotcall"], sep="\t", 
                           dtype={'chrom':int,'start':int,'end':int})
        s_calls = s_calls[s_calls.vclass >= 0] # filter classes
        #plotIndvAFDistribution( s_calls, samp )

        s_roh = rohregions[rohregions.sample == samp]
        if len(s_roh) == 0 : 
            print "Warning: sample has no roh!:",samp;continue

        if os.path.exists(targetfile) and not force :
            classcounts.append( read_csv( targetfile, sep="\t" ) )
        else :
            rohcallcounts = rohIntersectCallsClass( s_calls, s_roh )
            rohcallcounts['Individual.ID'] = samp
            rohcallcounts.to_csv( targetfile, sep="\t", index=False )
            classcounts.append( rohcallcounts )

        targetfile = os.path.join( rohafdir,samp+"_rohafcounts.tsv" )
        if os.path.exists(targetfile) and not force :
            afcounts.append( read_csv( targetfile, sep="\t" ) )
        else :
            rohafcounts = rohIntersectCallsAF( s_calls, s_roh )
            rohafcounts['Individual.ID'] = samp
            rohafcounts.to_csv( targetfile, sep="\t", index=False )
            afcounts.append( rohafcounts )

        targetfile = os.path.join( varrohdir,samp+"_varrohcounts.tsv" )
        if os.path.exists(targetfile) and not force :
            allvrohcounts.append( read_csv( targetfile, sep="\t" ) )
        else :
            varrohcounts = rohIntersectVars( s_calls, s_roh )
            #print varrohcounts.head()
            varrohcounts['Individual.ID'] = samp
            varrohcounts.to_csv( targetfile, sep="\t", index=False )
            allvrohcounts.append( varrohcounts )


    classcounts = concat(classcounts)
    afcounts = concat(afcounts)
    allvrohcounts = (concat(allvrohcounts)
                     .groupby(["chrom","pos","vclass","Region","AF","hfraction","Autozygosity"])
                     .size().reset_index())

    allvrohcounts = pivot_table(allvrohcounts, index=['chrom','pos','vclass',"hfraction",'Region','AF'], 
                            columns='Autozygosity', values=0 ).reset_index()
    allvrohcounts.fillna( 0, inplace=True )
    allvrohcounts["In_hom"] = allvrohcounts["In_hom"].astype(int)
    allvrohcounts["In_het"] = allvrohcounts["In_het"].astype(int)
    allvrohcounts["Out_hom"] = allvrohcounts["Out_hom"].astype(int)
    allvrohcounts["Out_het"] = allvrohcounts["Out_het"].astype(int)

    classcounts.to_csv( classcountsfile , sep="\t", index=False )
    afcounts.to_csv( afcountsfile , sep="\t", index=False )
    allvrohcounts.to_csv( varrohfile, sep="\t", index=False )

    return classcounts, afcounts, allvrohcounts
# END makeCountsFiles

# sample, delvars, delroh, benignvars, benignroh, delnonroh
def makeVarcountsFile( classcounts, varcountsfile, force=False ) :
    varcounts = melt( classcounts, id_vars=["Individual.ID","Autozygosity"], 
                                 value_vars=["class1","class2","class3","class4"] ).reset_index()
    #print "Varcounts melted"
    print varcounts.head()
    varcounts["impact"] = ["del" if x in ["class3","class4"] else "benign" 
                           for x in varcounts["variable"] ]
    #print varcounts["Autozygosity"].unique()
    varcounts["loc"] = ["roh" if x else "nonroh" for x in varcounts["Autozygosity"] ]

    varcounts["label"] = varcounts["impact"]+varcounts["loc"]
    varcounts = (varcounts.groupby(["Individual.ID","Autozygosity","label"])["value"].sum().reset_index())
    varcounts = varcounts.pivot('Individual.ID','label','value').reset_index()
    varcounts["delvars"] = varcounts["delroh"] + varcounts["delnonroh"]
    varcounts["benignvars"] = varcounts["benignroh"] + varcounts["benignnonroh"]
    #print varcounts.head(10)
    varcounts.to_csv( varcountsfile, sep="\t", index=False )
    return varcounts
# END makeVarcountsFile

################################################################################
def rohVarOverlap( rohfile, genesfile, vcffile, rohregions, sampleannot, force=False ) :
    print "Running rohVarOverlap"

    filepath, basename, ext = hk.getBasename( vcffile )
    targetdir = "%s/calls" % filepath
    hk.makeDir( targetdir )

    classcountsfile = "%s/%s_classcounts.tsv" % (targetdir, basename)
    afcountsfile = "%s/%s_afcounts.tsv" % (targetdir, basename)
    varcountsfile = "%s/%s_varcounts.tsv" % (targetdir, basename)
    varrohfile = "%s/%s_varroh.tsv" % (targetdir, basename)

    if (os.path.exists(varcountsfile) and os.path.exists(afcountsfile) and 
        os.path.exists(classcountsfile) and os.path.exists(varrohfile) and not force) :
        print "Warning: final file already exists", varcountsfile
        varcounts = read_csv(varcountsfile,sep="\t")
        classcounts = read_csv( classcountsfile, sep="\t" )
        afcounts = read_csv( afcountsfile, sep="\t" )
        vrohcounts = read_csv( varrohfile, sep="\t" )
        return varcounts, classcounts, afcounts, vrohcounts

    groupedvardata, alldelvars, allbenignvars = prepareGenesFile( genesfile )
    print "Grouped Var size:\n", groupedvardata.size()

    sampfiles = splitVcf( vcffile, targetdir )

    sampfiles = makeSampleCallsFiles( sampfiles, groupedvardata, sampleannot, targetdir )
    assert len(sampfiles) > 0

    classcounts, afcounts, vrohcounts = makeCountsFiles( sampfiles, targetdir, rohregions,
                                            classcountsfile, afcountsfile, varrohfile, force=force )

    varcounts = makeVarcountsFile( classcounts, varcountsfile, force=force )
    return varcounts, classcounts, afcounts, vrohcounts
# END rohVarOverlap

def plotVariantSlopes( varcounts, percentcov, dataset ):
    print "Running plotVariantSlopes"
    print varcounts.head()
    print percentcov.columns()
    varcounts_merge = merge( varcounts, 
                            percentcov[["Individual.ID","percent"]], 
                            on="Individual.ID" )

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

def localGetSampleAnnot( clustfile, filepats ):
    #filepats = patientInfo.currentFilePats( vcffile )
    if os.path.exists(clustfile) :
        sampleannot = read_csv( clustfile, sep="\t" )
        sampleannot.index = sampleannot["Individual.ID"]
        sampleannot = sampleannot.ix[filepats,:]
        mapping_regions = {"NWA":"Middle East", "NEA":"Middle East", "AP":"Middle East",
                           "SD":"Middle East", "TP":"Middle East", "CA":"Middle East",
                           "LWK":"Africa", "YRI":"Africa", "IBS":"Europe", "CEU":"Europe",
                           "TSI":"Europe", "FIN":"Europe", "GBR":"Europe",
                           "CHB":"East Asia", "CHS":"East Asia", "JPT":"East Asia",
                           "Anglo-American":"Anglo-American"
                          }
        sampleannot["Continent2"] = [mapping_regions[x] if mapping_regions.has_key(x)
                                      else "Unknown"
                                      for x in sampleannot.GeographicRegions3]
    else :
        sampleannot = sampleAnnotation( filepats )

    return sampleannot
# END localGetSampleAnnot

def plotRohLengthBurden( classcounts, dataset ) :

    autoz = classcounts[classcounts.Autozygosity]

    autoz["rohlen"] = (autoz["end"] - autoz["start"]) / 1000
    autoz = melt( autoz, id_vars=["rohlen","Individual.ID"], 
                 value_vars=["class1","class2","class3","class4"] )

    autoz["rohlen_bin"] = (np.digitize( autoz.rohlen, np.arange( 1000, autoz.rohlen.max(), 1000 ))+1) *1000
    print autoz.head()
    print "number of points:",autoz.shape
    
    r_dataframe = com.convert_to_r_dataframe(autoz)

    p = (ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x = "factor(rohlen_bin)",y="value" ) +
                ggplot2.geom_boxplot() + 
                ggplot2.geom_hline(yintercept=0) +
                ggplot2.ggtitle("Variant Burden by ROH by Length") +
                ggplot2.scale_x_discrete("ROH Length (KB)", 
                                           expand=robjects.IntVector((0,0))) +
                ggplot2.scale_y_continuous("Exome Variant Burden",
                                           expand=robjects.IntVector((0,0))) +
                ggplot2.stat_smooth(method="loess", se=False) +
                #ggplot2.scale_colour_brewer(palette="Set1") +
                ggplot2.facet_grid( robjects.Formula('variable ~ .'), scale="free_y" ) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                ggplot2.theme(**mytheme) )

    figname = "results/figures/roh/%s_roh_burden.png" % (dataset)
    print "Writing file", figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotRohLengthBurden

def plotRohSiteFreq( afcounts, dataset ):

    print afcounts.head()

    afsum = afcounts.groupby(["AF","Autozygosity"])["class1","class2",
                                                    "class3","class4"].sum().reset_index()
    afsum["AF"] = [float(x[x.find(",")+2:len(x)-1]) for x in afsum["AF"]]

    afsum = melt( afsum, id_vars=["AF","Autozygosity"], 
                 value_vars=["class1","class2","class3","class4"])
    afsum.rename( columns={"variable":"Class","value":"burden"}, inplace=True)

    print afsum.head()
    sys.exit(1)

    totals = afsum.groupby(["Autozygosity","Class"])["burden"].sum().reset_index()
    print totals.head()
    totals.rename(columns={"burden":"total"}, inplace=True )
    afsum = merge(afsum,totals, on=["Autozygosity","Class"] )
    print afsum.head()
    afsum["proportions"] = afsum["burden"] / afsum["total"]
    afsum["Autozygosity"] = ["RoH" if x else "Non-RoH" for x in afsum["Autozygosity"]]

    r_dataframe = com.convert_to_r_dataframe(afsum)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(AF)", y="proportions", fill="factor(Autozygosity)") +
                ggplot2.geom_bar(stat="identity",position="dodge") +
                ggplot2.geom_hline(yintercept=0) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45),
                                 'legend.position':'top'}) +
                ggplot2.scale_y_continuous("Proportion",
                                           expand=robjects.IntVector((0,0))) +
                ggplot2.scale_x_discrete("Allele Frequency") +
                                           #expand=robjects.IntVector((0,0))) +
                ggplot2.scale_fill_brewer( "Autozygosity", palette="Set1" ) +
                ggplot2.facet_grid( robjects.Formula('Class ~ .'), scale="free_y" ) +
                ggplot2.theme(**mytheme) )
                #ggplot2.ggtitle("Proportions of Auto by AF") + \

    figurename = "results/figures/roh/%s_sf_autoz.png" % dataset
    print "Making figure:",figurename
    grdevices.png(figurename) 
    p.plot()
    grdevices.dev_off()
# END plotRohSiteFreq

def plotVarRohSiteFreq( varafcounts, dataset ):
    print "Running plotVarRohSiteFreq"
    print varafcounts.head()

    varafcounts["Autozygosity"] = [x == 0 for x in varafcounts["Out"]]
    varafcounts.loc[:,"vclass"] = ["class%s"%x for x in varafcounts.vclass]

    afbins=[0.,.005,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]
    afsum = (varafcounts.groupby([cut(varafcounts["AF"], afbins),
                                  "Autozygosity","vclass"]).size().reset_index())
    afsum.rename(columns={0:"burden"}, inplace=True )

    afsum["AF_bin"] = [float(x[x.find(",")+2:len(x)-1]) for x in afsum["AF"]]


    #afsum = varafcounts.groupby(["AF","Autozygosity","vclass"]).size().reset_index()
    #afsum["AF"] = [float(x[x.find(",")+2:len(x)-1]) for x in afsum["AF"]]

    #afsum = melt( afsum, id_vars=["AF","Autozygosity"], 
                 #value_vars=["class1","class2","class3","class4"])
    #afsum.rename( columns={"variable":"Class","value":"burden"}, inplace=True)

    totals = afsum.groupby(["Autozygosity","vclass"])["burden"].sum().reset_index()
    print totals.head()
    totals.rename(columns={"burden":"total"}, inplace=True )
    afsum = merge(afsum,totals, on=["Autozygosity","vclass"] )
    print afsum.head()
    afsum["proportions"] = afsum["burden"] / afsum["total"]
    afsum["Autozygosity"] = ["Within RoH" if x else "Outside RoH" for x in afsum["Autozygosity"]]

    r_dataframe = com.convert_to_r_dataframe(afsum)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(AF_bin)", y="proportions", fill="factor(Autozygosity)") +
                ggplot2.geom_bar(stat="identity",position="dodge") +
                ggplot2.geom_hline(yintercept=0) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45),
                                 'legend.position':'top'}) +
                ggplot2.scale_y_continuous("Proportion",
                                           expand=robjects.IntVector((0,0))) +
                ggplot2.scale_x_discrete("Allele Frequency") +
                                           #expand=robjects.IntVector((0,0))) +
                ggplot2.scale_fill_brewer( "Autozygosity", palette="Set1" ) +
                ggplot2.facet_grid( robjects.Formula('vclass ~ .'), scale="free_y" ) +
                ggplot2.theme(**mytheme) )
                #ggplot2.ggtitle("Proportions of Auto by AF") + \

    figurename = "results/figures/roh/%s_sf_autoz3.png" % dataset
    print "Making figure:",figurename
    grdevices.png(figurename) 
    p.plot()
    grdevices.dev_off()
# END plotRVarohSiteFreq

################################################################################
# Main
################################################################################
if __name__ == "__main__" :

    os.chdir("..")
    optlist, args = getopt.getopt( sys.argv[1:], "r:ot")
    optlist = dict(optlist)

    #dataset = optlist.get("-r","test")
    dataset = optlist.get("-r","test2")
    if dataset == "test2" :
        rohfile = "./rawdata/test2/main/test2.clean.plink.filt.hom"
        genesfile = "./rawdata/test2/main/test2.clean_genes.tsv"
        vcffile = "./rawdata/test2/main/test2.clean.vcf.gz"
        clustfile = "./rawdata/mevariome/main/clust/variome.clean.annot"
    elif dataset == "daily" :
        rohfile = "./rawdata/daily/main/daily.clean.plink.filt.hom"
        genesfile = "./rawdata/daily/main/classify/daily.clean.sfilt_genes.tsv"
        vcffile = "./rawdata/daily/main/daily.clean.vcf.gz"
        clustfile = "./rawdata/daily/main/clust/daily.clean.annot"
    elif dataset == "onekg" :
        rohfile = "./rawdata/onekg/main/onekg.clean.plink.filt.hom"
        genesfile = "./rawdata/onekg/main/classify/onekg.clean.sfilt_genes.tsv"
        vcffile = "./rawdata/onekg/main/onekg.clean.vcf.gz"
        clustfile = "./rawdata/onekg/main/clust/onekg.clean.annot"
    elif dataset == "mevariome" :
        rohfile = "./rawdata/mevariome/main/variome.clean.plink.filt.hom"
        genesfile = "./rawdata/mevariome/main/classify/variome.clean.sfilt_genes.tsv"
        vcffile = "./rawdata/mevariome/main/variome.clean.vcf.gz"
        clustfile = "./rawdata/mevariome/main/clust/variome.clean.annot"
    else :
        print "Error: dataset not recognized -",dataset
        sys.exit(1)

    assert os.path.exists(vcffile)
    assert os.path.exists(genesfile)
    assert os.path.exists(rohfile)
    #rohfile = "./rawdata/merged/merged.clean.plink.filt.hom"

    filepath, basename, ext = hk.getBasename( vcffile )
    targetdir = "%s/calls" % filepath

    rohdata = read_csv(rohfile, delim_whitespace=True )
    #targetsamples = rohdata.IID.unique().tolist()
    #sampleannot = sampleAnnotation(targetsamples)

    sampleannot = localGetSampleAnnot( clustfile, rohdata.IID.unique().tolist() )

    percentcov, bedfiles = exomeCoverage( rohfile, sampleannot, targetdir, 
                                          basename=dataset, force=False )

    exomeregions = readExomeRohRegions( bedfiles ) 
    #tcolumns = ["chrom","start","end","sample","RohLen"]
    allrohregions = rohdata[["IID","CHR","POS1","POS2","KB","NSNP"]]
    allrohregions.rename(columns={"IID":"sample","CHR":"chrom","POS1":"start","POS2":"end"}, inplace=True)
    print allrohregions.head()
    allrohregions["RohLen"] = allrohregions["end"] - allrohregions["start"]

    #rohfile, genesfile, vcffile, rohregions, sampleannot, force
    varcounts, classcounts, afcounts, vrohcounts = rohVarOverlap( rohfile, genesfile, vcffile, 
                                                     allrohregions, sampleannot, True ) 

    #plotVariantSlopes( varcounts, percentcov, dataset )

    #plotRohLengthBurden( classcounts, dataset )

    #plotRohSiteFreq( afcounts, dataset )

    vrohcounts["alleles"] = vrohcounts["In_hom"] + vrohcounts["Out_hom"]
    varrohcounts = vrohcounts[ vrohcounts["alleles"] > 3 ]
    varrohcounts["binom"] = [binom_test(s) for s in varrohcounts[["In","Out"]].values]
    print varrohcounts[ varrohcounts.In > 0 ].head()

    plotVarRohSiteFreq( vrohcounts, dataset )

# END MAIN
################################################################################
    #varcounts["delnonroh"] = varcounts["delvars"] - varcounts["delroh"]
    #varcounts["benignnonroh"] = varcounts["benignvars"] - varcounts["benignroh"]

