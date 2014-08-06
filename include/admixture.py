#!/usr/bin/env python

#-------------------------------------------------------------------#
#                       Import Packages                             #
#-------------------------------------------------------------------#
import os, sys
import getopt
import csv
from collections import Counter
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist, squareform, euclidean
import re
#from random import randint

from housekeeping import *
import patientInfo
import popgencommands
from localglobals import *

gridextra = importr('gridExtra')

######################################################################
# run_admixture
# cv - cross validation
# B - bootstrapping
######################################################################
def run_admixture( targetfile, K, bootstrap=False, forceFlag=False ) :
    write2log( " - Running:"+whoami(), True )
    targetdir, filename, suffix = getBasename(targetfile)

    #if bootstrap :
        #Qfile = targetdir+"/"+filename+"."+str(K)+".best.Q"
        #logfile = "%s/%s_admix.%d.best.log" % (LOGDIR,filename,K)
    #else :
    Qfile = targetdir+"/"+filename+"."+str(K)+".Q"
    logfile = "%s/%s_admix.%d.log" % (LOGDIR,filename,K)
    print "Qfile:",Qfile
    if os.path.exists( Qfile ) and forceFlag is False:
        write2log(" - File"+Qfile+" already exists. Skipping command!",False)
        with open(logfile, 'r') as f:
            runout = f.read() 
        return Qfile, runout

    if bootstrap :
        command = "admixture -B %s %s" % ( targetfile, K )
    else :
        command = "admixture --cv %s %s" % ( targetfile, K )

    runout = runCommand( command )
    OUT = open(logfile,"wb")
    OUT.write(runout)
    OUT.close()
    command = "mv %s*.Q %s ; mv %s*.P %s" % (filename, targetdir, filename, targetdir)
    runCommand( command )
    return Qfile, runout
# End run_admixture

######################################################################
# parseGroupsFile
######################################################################
def parseGroupsFile( groupsfile ):
    groups = {}
    breakdown = {}
    FILE = open(groupsfile,"r")
    reader = csv.reader(FILE,delimiter="\t")
    header = {}
    for row in reader :
        if len(row) == 0 :
            continue
        if len(header) == 0 :
            for i in range(len(row)) :
                header[row[i]] = i
            continue
        currgroup = row[header["groups"]].strip()
        currpop = row[header["ethnicity"]].strip()
        currid = row[header["Individual.ID"]].strip()
        if not groups.has_key(currgroup) :
            groups[currgroup] = []
        groups[currgroup].append(currid)
        if not breakdown.has_key(currgroup) :
            breakdown[currgroup] = {}
        if not breakdown[currgroup].has_key(currpop) :
            breakdown[currgroup][currpop] = 0
        breakdown[currgroup][currpop] += 1
    #with open(groupsfile, 'r') as f:
        #runout = f.read()
    targetdir, filename, suffix = getBasename(groupsfile)
    outfile = "%s/%s_grpbreakdown.txt" % (targetdir, filename)
    print "Breakdown file:",outfile
    OUT = open(outfile,"wb")
    for group in sorted(breakdown) :
        OUT.write( "Group: #"+str(group)+"\n" )
        for pop in breakdown[group] :
            OUT.write( " - "+pop+":"+str(breakdown[group][pop])+"\n" )
    return groups
# End parseGroupsFile

######################################################################
# runAdmixtureClustering
######################################################################
def runAdmixtureClustering(bedfile, currK, forceFlag=forceFlag,
                           annotfile="./resources/annotation/patientannotation.ped") :
    write2log( " - Running:"+whoami(), True )
    print "Target File:",bedfile , "Annot:",annotfile,"CurrK:",currK
    bedfile  = fixfilepath( bedfile  )
    annotfile = fixfilepath( annotfile )

    targetdir, filename, suffix = getBasename(bedfile )

    Qfile = targetdir + filename+"."+str(currK)+".Q"
    groupsfile = "%s/%s.%d.groups.txt" % (targetdir,filename,currK)
    logfile = "%s/%s_admixclustering.%d.log" % (LOGDIR,filename,currK)
    print groupsfile
    if os.path.exists( groupsfile ) and forceFlag is False:
        write2log(" - File"+groupsfile+" already exists. Skipping command!",False)
        groups = parseGroupsFile( groupsfile )
        return groups

    command = "./rscripts/admix.R %s %s %d" % (bedfile , annotfile, currK)
    print "#"*50
    print command
    runout = runCommand( command )
    print runout
    OUT = open(logfile,"wb")
    OUT.write(runout)
    OUT.close()
    groups = parseGroupsFile( groupsfile )
    return groups
# runAdmixtureClustering

######################################################################
# findBestK
######################################################################
def findBestK( targetfile, maxK ) :
    write2log( " - Running:"+whoami(), True )
    targetdir, filename, suffix = getBasename(targetfile)
    logfile = "%s/%s_findBestK.%d.log" % (LOGDIR,filename,maxK)
    bestk = 0
    lowesterr = 1.0
    errorrates = {}
    qfiles = {}
    for K in range(2,maxK+1) :
        qfile, runout = run_admixture( targetfile, K )
        qfiles[K] = qfile
        cverrorsearch = re.search( "CV error \WK=\d+\W+([\.\d]+)", runout )
        cverror = 0
        if cverrorsearch is not None :
            cverror = cverrorsearch.group(1)
        print "K:",K,"Error:",cverror
        errorrates[K] = cverror 
        if float(cverror) < lowesterr :
            lowesterr = float(cverror)
            bestk = K
    errlogfile = fixfilepath("%s/%s_error.log"%(LOGDIR, filename),shouldexist=False)
    OUT = open(errlogfile, "wb")
    OUT.write("K\tError\n")
    for K in sorted(errorrates) :
        OUT.write( "%s\t%s\n" % (K, errorrates[K]) )
    OUT.close()
    # MAKE ERROR RATE GRAPH HERE
    command = "./rscripts/errorrate.R %s" % errlogfile
    runout = runCommand(command)
    OUT = open(logfile,"wb")
    OUT.write(runout)
    OUT.close()
    return qfiles, bestk, errlogfile
# End findBestK

######################################################################
# mergeVcfs
######################################################################
def mergeVcfs( path, filelist, targetfile ) :
    write2log( " - Running:"+whoami(), True )
    print "filelist:",filelist

    if os.path.exists( targetfile ) and forceFlag is False:
        print " - File"+targetfile+" already exists. Skipping command!"
        return

    for eachfile in filelist :
        assert eachfile[-3:] == ".gz", "Error: file isnt bgzipped"
    #for vfile in allfiles :
        #newfile = compress(vfile)
        
    assert( existsInPath("vcf-merge") )
    command = "vcf-merge --ref-for-missing '0/0' -d %s > %s" \
             % (" ".join( filelist ), targetfile)
    print command
    runCommand( command, True )
    return
# END mergeVcfs

######################################################################
# parseSampleGroups
# annotfile="./resources/annotation/patientannotation.ped"
######################################################################
def parseSampleGroups( bedfile, patientdata ) :

    filepath, basename, suffix = getBasename( bedfile )
    famfile = "%s/%s.fam" % (filepath,basename)
    pats = patientInfo.currentFilePats( famfile )

    #print pats
    #patientdata = read_csv( annotfile, sep="\t" )
    #newids = [ x[:x.find('_')] if len(x) > 21 else x for x in patientdata["Individual.ID"].tolist()]
    #print newids
        #header = patientdata["header"]

    #print pats[:10]
    pats = DataFrame( pats, columns=["Individual.ID"])
    groups = merge(pats, patientdata[["Individual.ID","ethnicity"]], on="Individual.ID", how="left")
    #print groups.head(10)
    return groups
# END parseSampleGroups

def orderSamples( admixdf, K, levels, annotcol="GeographicRegions2", maxK=100 ):
    print "Running orderSamples, K:",K
    comporder = []
    for lvl in levels :
        subdf = admixdf[(admixdf[annotcol] == lvl)]
        if len(subdf) < 3 : print "Error: level -",lvl,"len",len(subdf); continue
        #rawdata = pivot_table(subdf, values="value", index="IID", columns="variable" )
        rawdata = subdf[[x for x in subdf.columns if re.match("K\d+", x) is not None]]
        rawdata.index = subdf.IID
        distxy = squareform(pdist(rawdata, metric='euclidean'))
        l = leaves_list(linkage( distxy ))
        #print "L",len(l),l[:10] 
        nodeorder =  DataFrame(index=l, data={"node_order":range(1,len(l)+1)})
        #print "Node order"
        #print nodeorder.head()
        rawdata.reset_index(inplace=True)
        rawdata["node_order"] = 0
        rawdata.update( nodeorder )
        print rawdata.head()
        rawdata["lvl"] = lvl
        comporder.append(rawdata[["IID","lvl","node_order"]])
        #meanvals = subdf.groupby("variable")["value"].mean().reset_index()
        #vorder = meanvals.sort("value",ascending=[False])["variable"]
        #print vorder
        #maxcols = subdf.ix[subdf.groupby( "IID" )["value"].idxmax()]
        #mostimportantcomp,count = Counter( maxcols["variable"]).most_common(1)[0]
        #print mostimportantcomp,count
        #taken.append(mostimportantcomp)
        #comporder.append( DataFrame({'lvl':lvl, 'component':vorder}))
                                  #'N':len(maxcols)}))

    comporder = concat(comporder).reset_index()
    #uniqcomps = comporder.groupby("component").size().reset_index()
    #uniqcomps["sort_idx"] = uniqcomps.index

    admixdf_melt = melt( admixdf, id_vars=["IID",annotcol], 
                        value_vars=["K%d"%x for x in range(1,K+1)])

    admixdf_melt = merge( admixdf_melt, comporder, left_on=[annotcol,"IID"], 
                    right_on=["lvl","IID"] )

    #print admixdf.head()
    admix_sort = admixdf_melt.sort( [annotcol,"node_order"] )

    return admix_sort.IID.tolist()
# END orderSamples

def reorderComponents(previousadmix, currentadmix, currentK ):
    print "Running reorderComponents", currentK
    if previousadmix is None : return {}

    #padmix = previousadmix[["K%d"%x for x in range(1,currentK)]].to_dict() 
    #cadmix = currentadmix[["K%d"%x for x in range(1,currentK+1)]].to_dict() 

    test = previousadmix["K2"]
    print type(test)
    print test.head()

    targetcolumns = ["K%d"%x for x in range(1,currentK+1)]
    padmix = {K:previousadmix[K].map(float).tolist() 
              for K in ["K%d"%x for x in range(1,currentK)]}
    cadmix = {K:currentadmix[K].tolist() 
              for K in ["K%d"%x for x in range(1,currentK+1)]}

    matching = []
    for prevK in padmix : 
        for currK in cadmix :
            #print "prevK",prevK,"currK",currK
            #dst = euclidean(previousadmix[prevK].tolist(),
                            #currentadmix[currK].tolist())
            assert len(padmix[prevK]) == len(cadmix[currK])
            dst = euclidean(padmix[prevK], cadmix[currK])
            matching.append(Series({'prevK':prevK,'currK':currK,'dist':dst}))

    matching = DataFrame( matching )
    bestmatches = DataFrame(columns=["prevK","currK","dist"])
    for idx, row in matching.sort("dist").iterrows() :
        if any( (bestmatches.prevK == row["prevK"]) | 
               (bestmatches.currK == row["currK"])  ) : continue
        bestmatches.loc[len(bestmatches)+1] = row

    #bestmatches = matching.ix[matching.groupby('prevK')["dist"].idxmin()]
    assert len(bestmatches.currK.unique()) == (currentK-1)
    newlbls = {y:x for x,y in bestmatches[['prevK','currK']].values}
    for currK in ["K%d"%x for x in range(1,currentK+1)] : 
        if not newlbls.has_key(currK) : 
            newlbls[currK] = 'K%d'%currentK
            break
    assert all([newlbls.has_key(x) for x in targetcolumns])
    return newlbls
    #if annotcol == "GeographicRegions2" : 
        #levels = target_geographic_regions2 
    #elif annotcol == "Continent" :
        #levels = target_continents
    #taken = []
    #comporder = []
    #for lvl in levels :
        #subdf = admixdf[(admixdf[annotcol] == lvl) & ~(admixdf["variable"].isin(taken))]
        #if len(subdf) == 0 : continue
        #maxcols = subdf.ix[subdf.groupby( "IID" )["value"].idxmax()]
        #mostimportantcomp,count = Counter( maxcols["variable"]).most_common(1)[0]
        #print mostimportantcomp,count
        #taken.append(mostimportantcomp)
        #comporder.append( Series({'lvl':lvl, 'component':mostimportantcomp,
                                  #'N':len(maxcols)}))
    #comporder = DataFrame(comporder)
    #print comporder.head()
    #uniqcomps = comporder.groupby("component").size().reset_index()
    #uniqcomps["sort_idx"] = uniqcomps.index
    #admixdf = merge( admixdf, uniqcomps, left_on="variable", right_on="component" )
    #print admixdf.head()
    #admix_sort = admixdf.sort( [annotcol,"sort_idx","value"] )
    #return admix_sort.IID.tolist(), uniqcomps.component
# END reorderComponents

def readQfile( qfile, annotcol=None ) :

    filepath, fullbasename, suffix = getBasename(qfile) 
    basename = fullbasename[:fullbasename.rfind(".")]
    K = int(fullbasename[fullbasename.rfind(".")+1:])
    famfile = os.path.join(filepath,basename+".fam")
    famdf = read_csv(famfile,delim_whitespace=True, 
                     header=None, names=["FID","IID","PID","MID","Aff","Gender"])

    #print famdf.head()
    datadf = read_csv( qfile, delim_whitespace=True, 
                      header=None, names=["K"+str(x) for x in range(1,K+1)])
    #print datadf.head()

    alldat = concat([famdf,datadf],axis=1)
    #alldat_melt = melt( alldat, id_vars=["FID","IID","PID","MID","Aff","Gender"])
    alldat_annot = addSampleAnnotation(alldat, mergecol="IID")
    if annotcol is not None :
        alldat_annot = alldat_annot[(alldat_annot[annotcol].notnull()) & 
                                    (alldat_annot[annotcol] != "Unknown")]
        alldat_annot[annotcol] = [x.strip().replace(" ",".") for x in alldat_annot[annotcol]]
    #finallevels = [x for x in levels if x in alldat_annot[annotcol].unique()]
    #alldat_annot = alldat_annot[alldat_annot[annotcol].isin(finallevels)]
    alldat_annot.to_csv(os.path.join(filepath,fullbasename+"_annot.tsv"),sep="\t",index=False)
    return alldat_annot, K
# END readQfile

################################################################################
def plotAdmixture( admixdf, K, figbasename, sampleorder=None, 
                  levels=None, annotcol="GeographicRegions2" ):
    print "Running plotAdmixture"
    
    if sampleorder is None :
        sampleorder = orderSamples( admixdf, K, levels, annotcol )

    if any( [x not in levels for x in admixdf[annotcol].unique()] ):
        print "Levels :",levels
        print "Unique atts",admixdf[annotcol].unique()
        sys.exit(1)

    admixdf = admixdf[admixdf.IID.isin(sampleorder)]
    admixdf_melt = melt( admixdf, id_vars=["IID",annotcol], 
                        value_vars=["K%d"%x for x in range(1,K+1)])

    admixdf_melt["prop"] = admixdf_melt.value.astype(float)
    admixdf_melt["currK"] = "K = "+str(K)
    s = admixdf_melt.groupby("IID")["prop"].sum().reset_index()

    #print [x for x in admixdf.IID if x not in sampleorder ]
    #print "$" *60
    r_dataframe = com.convert_to_r_dataframe(admixdf_melt)
    r_dataframe = fixRLevels( r_dataframe, "IID", sampleorder )
    r_dataframe = fixRLevels( r_dataframe, annotcol, levels )
    r_dataframe = fixRLevels( r_dataframe, "variable", ["K%d"%x for x in range(1,K+1)])
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(IID)", y="value", fill="factor(variable)") +
                ggplot2.geom_bar(stat="identity") +
                #ggplot2.ggtitle("K = "+str(K)) +
                ggplot2.scale_y_continuous("Components") +
                ggplot2.facet_grid(robjects.Formula('currK ~ '+annotcol),scale="free_x", space="free") +

                #ggplot2.scale_fill_brewer(palette="Set1") +
                ggplot2.scale_fill_manual(values=makeLargePalette(K)) +
                ggplot2.coord_fixed() +
                ggplot2.theme(**admixtheme) )

    #figurename = os.path.join(filepath,basename+"_"+str(K)+".png")
    figurename = figbasename+"_"+str(K)+".png"
    print "Making figure:",figurename
    grdevices.png(figurename, width=8, height=3,units="in",res=300)
    p.plot()
    grdevices.dev_off()
    return p
# END plotAdmixture

################################################################################
def plotAllAdmixture( admixdf, maxK, figbasename, sampleorder=None, 
                  levels=None, annotcol="GeographicRegions2" ):
    print "Running plotAllAdmixture"
    
    admixdf = admixdf[admixdf.IID.isin(sampleorder)]
    #admixdf_melt = melt( admixdf, id_vars=["IID",annotcol], 
                        #value_vars=["K%d"%x for x in range(1,K+1)])

    print admixdf.head()
    #admixdf["prop"] = admixdf.value.astype(float)
    #s = admixdf.groupby("IID")["prop"].sum().reset_index()

    #print [x for x in admixdf.IID if x not in sampleorder ]
    #print "$" *60
    r_dataframe = com.convert_to_r_dataframe(admixdf)
    r_dataframe = fixRLevels( r_dataframe, "IID", sampleorder )
    r_dataframe = fixRLevels( r_dataframe, annotcol, levels )
    r_dataframe = fixRLevels( r_dataframe, "variable", ["K%d"%x for x in range(1,maxK+1)])
    r_dataframe = fixRLevels( r_dataframe, "currK", ["K = %d"%x for x in range(1,maxK+1)])
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(IID)", y="value", fill="factor(variable)") +
                ggplot2.geom_bar(stat="identity") +
                #ggplot2.ggtitle("K = "+str(K)) +
                ggplot2.scale_y_continuous("Components") +
                ggplot2.facet_grid( robjects.Formula('currK ~ '+annotcol), scale="free_x", space="free" ) +
                #ggplot2.scale_fill_brewer(palette="Set1") +
                #ggplot2.scale_fill_discrete() +
                ggplot2.scale_fill_manual(values=makeLargePalette(maxK)) +
                ggplot2.coord_fixed() +
                ggplot2.theme(**admixtheme) )
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +

    #figurename = os.path.join(filepath,basename+"_"+str(K)+".png")
    figurename = figbasename+"_all2.png"
    print "Making figure:",figurename
    grdevices.png(figurename, width=8, height=1*maxK,units="in",res=300)
    p.plot()
    grdevices.dev_off()
    return p
# END plotAllAdmixture

def plotErrorRate( errorfile, figbasename ) :
    perror = read_csv(errorfile,sep="\t")
    r_dataframe = com.convert_to_r_dataframe(perror)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="K", y="Error") +
                ggplot2.geom_point(colour="red",size=3) +
                ggplot2.scale_y_continuous("Cross-validation Error") +
                ggplot2.scale_x_continuous("K runs", 
                                           breaks=robjects.IntVector(range(2,perror.K.max()+1))) +
                #ggplot2.ggtitle("Predicted Error") +
                ggplot2.theme(**pointtheme) )

    figurename = figbasename+"_error.png"
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# END plotErrorRate

######################################################################
# admixtureAnalysis
######################################################################
def admixtureAnalysis( bedfile, sampleannot, tlevels, maxK=None, 
                      force=False, annotcol="GeographicRegions" ) :

    #print "Working Dir",os.getcwd()
    print "Bedfile:",bedfile
    figbasename = bedfile[:-4]

    #print "Annotfile",annotfile
    #assert os.path.exists( annotfile )
    if maxK is None :
        groups = parseSampleGroups( bedfile, sampleannot )
        #maxK = len(uniqify(groups))
        maxK = len(groups["ethnicity"].unique().tolist())
        #print "Groups",len(groups),"Max K:",maxK

    print "MaxK:",maxK
    if maxK < 2: return maxK
    maxK = min(maxK,14)
    #maxK = min(maxK,4)
    qfiles,bestK,errorfile = findBestK( bedfile, maxK )

    plotErrorRate( errorfile, figbasename )
    print qfiles

    # Identify sample order using the "BestK"
    admixdf, K = readQfile( qfiles[bestK], annotcol )
    print hk.dfTable(admixdf[annotcol].tolist())
    sampleorder = orderSamples( admixdf, K, tlevels, annotcol )

    #allplots = []
    alldata = []
    prevdata = None
    for K in qfiles : 
        print K, qfiles[K]
        admixdf, K = readQfile( qfiles[K], annotcol )
        print admixdf[annotcol].unique()
        #sys.exit(1)
        comporder = reorderComponents( prevdata, admixdf, K )
        admixdf.rename(columns=comporder,inplace=True)
        admixdf_melt = melt( admixdf, id_vars=["IID",annotcol], 
                        value_vars=["K%d"%x for x in range(1,K+1)])
        admixdf_melt["currK"] = "K = "+str(K)
        alldata.append( admixdf_melt )
        plotAdmixture( admixdf, K, figbasename, sampleorder, tlevels, annotcol )
        #allplots.append(p)
        prevdata = admixdf

    alldata = concat(alldata).reset_index(drop=True)

    plotAllAdmixture( alldata, maxK, figbasename, sampleorder, tlevels, annotcol )
    #figurename = bedfile[:-4]+"_all.png"
    #print "Making figure:",figurename
    #grdevices.png(figurename, width=8, height=2*len(allplots),
                  #units="in",res=300)
    #gridextra.grid_arrange( *allplots )
    #grdevices.dev_off()
    sys.exit(1)

    qfile,runout = run_admixture( bedfile, bestK, True, forceFlag )

    runAdmixtureClustering(bedfile,bestK,True)
    return bestK
# End admixtureAnalysis

def seriousClean( vcffile, keepfile=None, rerun=False ):
    print "Running seriousClean"
    filepath, basename, suffix = hk.getBasename(vcffile)
    cleanped = "%s/%s"%(filepath,basename)
    command = ["vcftools","--gzvcf", vcffile,
               "--remove-filtered-all",
               "--remove-indels",
               "--maf",".005",
               "--min-alleles","2",
               "--max-alleles","2",
               "--plink-tped","--out",cleanped]

    if keepfile is not None : 
        assert os.path.exists(keepfile)
        command = command+["--keep",keepfile]

    print " ".join(command)
    tped = cleanped+".tped"
    if not os.path.exists( tped ) or rerun:
        out = subprocess.check_output( command )

    bedfile = popgencommands.run_plink_convert(tped, force=rerun)
    return bedfile
# END seriousClean
    #patientdata = sampleAnnotation()
    #targetdir = "%s/pca/" % filepath
    #hk.makeDir(targetdir)
    #cptargetfile = targetdir+basename+suffix+".gz"

    #if not os.path.exists(cptargetfile) :
        #out = hk.runCommand( "cp %s/%s%s* %s" % (filepath,basename,suffix, targetdir) )
        #print "Copy targetfile:",cptargetfile
    #recodefile = plink_recode12( tped, force=rerun )
    #print "#"*50
    #print "Made recodefile:",recodefile
    # Modify ped
    #modped = popgencommands.modify_ped( recodefile, patientdata, calcdistances=True, shorten=True, force=rerun )
    #frqfile = plink_makefrqfile( modped, force=rerun )
    #newfrq, excludebase = popgencommands.excludeSnps( modped, force=rerun )

def getLevels( annotcol, targetvcf ):
    print "Running getLevels"
    if annotcol == "GeographicRegions2" : 
        levels = target_geographic_regions2 
    elif annotcol == "Continent" :
        levels = target_continents
    elif annotcol == "GeographicRegions" :
        levels = target_geographic_regions
        if "Oceania" in levels : levels.remove("Oceania") 
        if "America" in levels : levels.remove("America")
        print levels
    else :
        print "Error: unknown levels for column -",annotcol
        sys.exit(1)

    targetpats = patientInfo.getPats( targetvcf )
    sampleannot = sampleAnnotation(targetpats)

    sampleannot = sampleannot[(sampleannot[annotcol].notnull()) & 
                                (sampleannot[annotcol] != "Unknown")]
    sampleannot[annotcol] = [x.strip().replace(" ",".") for x in sampleannot[annotcol]]

    finallevels = [x for x in levels if x in sampleannot[annotcol].unique()]
    sampleannot = sampleannot[sampleannot[annotcol].isin(finallevels)]
    #keepfiles = sampleannot["Individual.ID"].tolist()
    return finallevels, sampleannot#, keepfiles
# END getLevels

def subsectionVCFbyPop( vcffile, sampleannot, targetdir, basename="test", 
                       force=False, annotcol="Continent" ) :
    allpopfiles = {}
    for group, data in sampleannot.groupby("Continent") :
        pop = group.replace(" ","_")
        popdir = os.path.join(targetdir,pop)
        makeDir( popdir )
        #popfile = "%s_%s.txt" %(popsprefix, pop )
        popfile = os.path.join(popdir, "%s_%s.txt" %(basename,pop))
        data[["Individual.ID"]].to_csv(popfile,sep="\t",index=False,header=False)
        allpopfiles[pop] = popfile
    print allpopfiles

    popvcfs = {}
    for pop in allpopfiles :
        popdir = os.path.join(targetdir,pop)
        newvcffile = os.path.join(popdir, "%s_%s" %(basename,pop))
        finalvcf = newvcffile+".vcf.gz"
        popvcfs[pop] = finalvcf
        # Chekc if file already exists
        if os.path.exists( finalvcf ) and not force : continue

        print "Targetvcf:",newvcffile+".recode.vcf"
        print "Pop:",pop, "File:",allpopfiles[pop]
        command = "vcftools --gzvcf %s" \
            " --min-alleles 2 --max-alleles 2" \
            " --recode" \
            " --keep %s" \
            " --out %s" \
            % ( vcffile, allpopfiles[pop], newvcffile )

        out = runCommand(command)
        assert os.path.exists( newvcffile+".recode.vcf" ), out
        subprocess.call(['mv', newvcffile+".recode.vcf", newvcffile+".vcf"])
        subprocess.call(["bgzip",newvcffile+".vcf"])
        subprocess.call(["tabix","-p","vcf",finalvcf])

    return popvcfs
# subsectionVCFbyPop

def copyToSubDir1( tfile, foldername, force=False ) :
    #filepath, filename, suffix = getBasename(vcffile)
    filepath, filename = os.path.split( tfile )
    tdir = os.path.join(filepath, foldername)
    #if outdir is None : outdir = targetdir
    makeDir( tdir )

    targetfile = os.path.join(tdir,filename)
    if not os.path.exists(targetfile) or force: 
        command = ["cp",os.path.join(filepath,filename), tdir]
        retval = subprocess.call( command )
    assert os.path.exists(targetfile)
    return targetfile
# End copyToSubDir1
    #out = runCommand("cp %s/%s* %s" % (filepath,filename,tdir))

######################################################################
# Main
######################################################################
if __name__ == "__main__":

    optlist, args = getopt.getopt( sys.argv[1:], "bk")
    optlist = dict(optlist)

    # Change running directory
    os.chdir("..")

    #vcffile = "./rawdata/test2/main/test2.clean.vcf.gz"
    #vcffile = "./rawdata/daily/main/daily.clean.vcf.gz"
    #vcffile = "./rawdata/mevariome/main/variome.clean.vcf.gz"
    vcffile = "./rawdata/mergedaly/main/meceu.clean.vcf.gz"
    #vcffile = "./rawdata/merge1kg/main/me1000G.clean.vcf.gz"
    print "Using Vcffile:",vcffile

    #maxK = 4
    #maxK = optlist.get("-k",maxK)
    #bedfile = optlist.get("-b",bedfile)
    #assert os.path.exists(bedfile)
    #assert is_int( maxK )

    annotcol="GeographicRegions"

    targetvcf = copyToSubDir( vcffile, "admixture" )
    filepath, filename, suffix = getBasename(targetvcf)

    # Run for all data together

    levels, sampleannot = getLevels( annotcol, targetvcf )
    keepfile = os.path.join(filepath, filename+".keeppats")
    sampleannot[["Individual.ID"]].to_csv(keepfile, header=None,index=None)
    bedfile = seriousClean( targetvcf, keepfile, False )
    admixtureAnalysis( bedfile, sampleannot, levels, annotcol=annotcol )

    # Run for all continents separately
    popvcfs = subsectionVCFbyPop( targetvcf, sampleannot, filepath, filename )
    for pop in popvcfs :
        #targetdir, filename, suffix = getBasename(popvcfs[pop])
        #targetpats = patientInfo.getPats( popvcfs[pop] )
        #sampleannot = sampleAnnotation(targetpats)
        levels, sampleannot = getLevels( annotcol, popvcfs[pop] )
        #print "Targetpats:",targetpats
        if len(sampleannot) < 15 : continue
        keepfile = os.path.join(filepath, filename+".keeppats")
        sampleannot[["Individual.ID"]].to_csv(keepfile, header=None,index=None)
        bedfile = seriousClean( popvcfs[pop], keepfile, False )
        print "Running pop:",pop,"vcf:",popvcfs[pop],"bed:",popvcfs[pop]
        admixtureAnalysis( bedfile, sampleannot, levels, annotcol=annotcol )



# END Main
################################################################################
    #admixtureAnalysis( bedfile, annotfile, maxK )


