#!/usr/bin/python

import os, sys
import string
from pandas import *
import pandas.rpy.common as com
from collections import Counter
from scipy import stats
import math
from localglobals import *

#import rpy2.robjects as robjects
#from rpy2.robjects.packages import importr
#from rpy2.robjects.lib import grid
#from rpy2.robjects.lib import ggplot2
#r = robjects.r
#rprint = robjects.globalenv.get("print")
#rstats = importr('stats')
#grdevices = importr('grDevices')
#base = importr('base')

forceFlag = False

#--------------------------------------------------------------------#
#                               Annotation                           #
#--------------------------------------------------------------------#

######################################################################
# parseRegionsFile
######################################################################
def parseRegionsFile():
    #refseqgenesbed = "./bedfiles/refseq_genes_simple.bed"
    bedfile = "./resources/regionbeds/collab_percentcoverage.txt"

    regions = {}
    for row in csv.reader(open(bedfile), delimiter="\t") :
        if not regions.has_key(row[0]) : regions[row[0]] = []
        regions[row[0]].append([int(row[1]),int(row[2])])
        #print [int(row[1]),int(row[2]),row[3]]
    return regions
# End parseRegionsFile

######################################################################
# intersectRegions
######################################################################
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

######################################################################
# burdenAnalysis
######################################################################
def burdenAnalysis( tfile ) :
    head = {}
    patburden = {}
    for row in csv.reader(open(tfile),delimiter="\t") :
        if len(head) == 0 :
            for i in range(len(row)) :
                head[row[i]] = i
            continue
        sample = row[head["Sample"]]
        genotype = row[head["Genotype"]]
        if not patburden.has_key(sample): patburden[sample] = {"het":0,"hom":0}
        assert genotype in ["het","hom"], "Genotype:"+genotype
        patburden[sample][genotype] += 1
    return patburden
# End burdenAnalysis

######################################################################
# readAnnotData
######################################################################
def readAnnotData( ):
    pfile = "./patientannot/annotation/patientannotation_laser.txt"
    hd = {}
    annottable = {}
    targetannot = ["Gender","Phenotype","Source","ethnicity","continent"]
    for row in csv.reader(open(pfile),delimiter="\t") :
        if len(hd) == 0 :
            for i in range(len(row)) :
                hd[row[i]] = i
            #print hd
            #continue
        annottable[row[hd["Normalized"]]] = [row[hd[x]] for x in targetannot]
    return annottable, targetannot
# End readAnnotData

######################################################################
# allPatVarStats
######################################################################
def allPatVarStats( allclasses, outdir="./classes" ) :
    alldata = {}
    for vclass in allclasses :
        cfile = allclasses[vclass]
        patburdens = burdenAnalysis( cfile )
        for pat in patburdens :
            if not alldata.has_key(pat) : alldata[pat]={}
            alldata[pat][vclass] = patburdens[pat]

    patannot, header = readAnnotData()
    samplesfile = "%s/allpatsvarstats.txt" % outdir
    OUT = open(samplesfile, "wb")
    print allclasses, header
    OUT.write("Sample"+"\t".join([x+"(hom|het)" for x in allclasses])+
              "hets\thoms"+"\t".join(header)+"\n")
    for pat in alldata :
        entry = [pat]
        totals = [0,0]
        for vclass in allclasses :
            if alldata[pat].has_key(vclass) : 
                if isinstance(alldata[pat][vclass], basestring) :
                    entry.append(alldata[pat][vclass])
                else :
                    values = alldata[pat][vclass].values()
                    entry.append( "|".join([str(x) for x in values]) )
                    totals[0] += values[0]
                    totals[1] += values[1]
            else :
                entry.append("0|0")
        entry += [str(x) for x in totals]
        if patannot.has_key(pat): entry += patannot[pat]
        OUT.write( "\t".join(entry) + "\n" )
        print entry
    OUT.close()
    return
# END allPatVarStats

######################################################################
# calculateAF
######################################################################
def calculateAF( vdata ) :
    print "Starting to calculate AF!"
    print vdata.head(3)
    uniqsamples = len( vdata.Sample.unique().tolist() )
    vdata["key"] = vdata.chrom.map(str)+":"+vdata.pos.map(str)
    varcounts = vdata.groupby("key")
    print varcounts.count()
    sys.exit(1)

######################################################################
# calculateRVIS
######################################################################
def calculateRVIS( genecounts, prefix=None ):

    #xdata_uniq = xdata[["chrom","pos","gene"]].drop_duplicates()
    #Xgenecounts = Counter(xdata_uniq["gene"].tolist())
    #ydata_uniq = ydata[["chrom","pos","gene"]].drop_duplicates()
    #Ygenecounts = Counter(ydata_uniq["gene"].tolist())
    #ycounts = DataFrame([[x,Ygenecounts[x]] for x in Ygenecounts],columns=["gene","Ycount"])
    #xcounts = DataFrame([[x,Xgenecounts[x]] for x in Xgenecounts],columns=["gene","Xcount"])

    robjects.globalenv["Ycount"] = robjects.FloatVector(genecounts["Ycount"])
    robjects.globalenv["Xcount"] = robjects.FloatVector(genecounts["Xcount"])
    model_x = rstats.lm("Ycount ~ Xcount")
    resid = list(rstats.rstudent(model_x))
    genecounts["StdResid"] = resid
    genecounts["Colour"] = "Innerquartile"
    lowquantile = genecounts.StdResid.quantile(.05)
    highquantile = genecounts.StdResid.quantile(.95)
    genecounts["Colour"][genecounts.StdResid <= lowquantile] = "2% most intolerant"
    genecounts["Colour"][genecounts.StdResid >= highquantile] = "2% least intolerant"
    #print "Allcounts head:",genecounts.head(20)
    #print "Allcounts shape:",genecounts.shape
    #print genecounts[genecounts.Colour != "Innerquartile"].head(20)

    r_dataframe = com.convert_to_r_dataframe(genecounts)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x = "Xcount",y="Ycount" ) + \
                ggplot2.geom_jitter(ggplot2.aes_string(colour="factor(Colour)")) + \
                ggplot2.scale_colour_manual(values = robjects.StrVector(("blue", "red", "grey"))) + \
                ggplot2.ggtitle("Determination of RVIS") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.scale_y_continuous("Sum of all Common (MAF >.01) functional variants in a gene", limits=robjects.IntVector((0,120)))+ \
                ggplot2.scale_x_continuous("Sum of all variant sites in a gene", limits=robjects.IntVector((0,600))) + \
                ggplot2.stat_smooth(method="lm", se=False)
                #ggplot2.geom_abline(intercept=rstats.coef(model_x)[1], slope=rstats.coef(model_x)[2])
    if prefix is None : prefix = "test"
    if prefix.find(" ")>=0 : prefix = prefix.replace(" ","_")
    figname = "results/figures/rvis/%s_rvis_lm.png" % prefix
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
    return genecounts
# END calculateRVIS

######################################################################
# readVariantSet_old
######################################################################
def readVariantSet_old( prefix, vclasses=["HIGH","MODERATE","LOW"] ):
    allvars = []
    for vclass in vclasses :
        #targetfile = "./classes/ciliopathies.names.chimp.regions.annot_"+vclass+".tsv"
        targetfile = "./results/classes/"+prefix+".annot_"+vclass+".tsv"
        #print targetfile
        assert os.path.exists(targetfile) 
        data = read_csv(targetfile,sep="\t")
        #data = data[data["Sample"].isin(unrelated)]
        allvars.append(data)
    return concat(allvars)
# END readVariantSet_old
######################################################################
# CommentStripper
######################################################################
def CommentStripper2 (iterator):
    for line in iterator:
        if line [:1] == '#':
            continue
        if not line.strip ():
            continue
        if len(line) > 200 : 
            continue
        yield line 

################################################################################
# makeGeneCounts
#Sample  Genotype        chrom   pos     ref     mut     gene    AF      vclass  polyphen2
#HG00736 het     1       866422  C       T       SAMD11  0.0035  SYNONYMOUS_CODING       .
################################################################################
def makeGeneCounts( prefix, regions, regionlookup, vclasses=["HIGH"], maf=0.0 ) :
    allvars = {}
    hdr = {}
    #vclasses=["HIGH"]
    for vclass in vclasses :
        targetfile = "./results/classes/"+prefix+"_"+vclass+".tsv"
        assert os.path.exists(targetfile), "Targetfile:"+targetfile
        varcount = 0
        for row in csv.reader(CommentStripper2(open(targetfile)),delimiter="\t") :
            if len(hdr)==0 : 
                hdr = dict([[row[i],i] for i in range(len(row))])
                continue
            if not hk.is_float(row[hdr["AF"]]) : print row
            if regionlookup.has_key( row[hdr["Sample"]] ) and float(row[hdr["AF"]]) > maf:
                vkey = row[hdr["chrom"]]+":"+row[hdr["pos"]]+":"+row[hdr["gene"]]
                if not allvars.has_key(vkey) : 
                    allvars[vkey] = [0 for i in range(len(regions))]
                allvars[vkey][regionlookup[row[hdr["Sample"]]]] += 1
                varcount += 1
        print "Class:",vclass,"#vars:",varcount

    allvars = DataFrame([x.split(":")+allvars[x] for x in allvars], 
                        columns=["chrom","pos","gene"]+regions)
    return allvars
# END makeGeneCounts

################################################################################
# makeGeneCounts1
################################################################################
def makeGeneCounts1( prefix, regions, vclasses=[1], minmaf=0.0 ) :
    varfile = prefix+"_genes.tsv"
    csvread = csv.reader(open(varfile),delimiter="\t")
    header = csvread.next()
    foundregions = [x for x in regions if x in header]
    allvars = []
    for row in csvread :
        s = Series(data=row,index=header)
        if s["vclass"] not in vclasses : continue
        if s["AF"] < minmaf : continue
        for region in foundregions : 
            #ref,het,hom = [int(x) for x in s[region].split(":")]
            regionaf = (2*hom+het) / (2.*(ref+het+hom))
            s[region] = regionaf
            #s[region] = sum([int(x) for x in s[region].split(":")[1:]])
            #if regionaf < minmaf : continue
            #regseries = Series({"chrom":s["chrom"],"pos":s["pos"],"Gene":s["Gene"],
                        #"AF":regionaf,"region":region,"samps":(hom+het)})
            #allvars.append(regseries)
        allvars.append(s[["chrom","pos","Gene","vclass"]+foundregions])
        #print s[["chrom","pos","Gene"]+foundregions]

    allvars = DataFrame(allvars)
    #print allvars.head(10)
    return allvars
# END makeGeneCounts1

################################################################################
# makeGeneCounts2
################################################################################
def makeGeneCounts2( prefix, region, vclasses=[1], minmaf=0.0 ) :
    print "Running makeGeneCounts2 - region:",region
    varfile = prefix+"_genes.tsv"
    csvread = csv.reader(open(varfile),delimiter="\t")
    header = csvread.next()
    foundregions = [x for x in regions if x in header]
    allvars = []
    #limit = 100
    for row in csvread :
        s = Series(data=row,index=header)
        if s["vclass"] not in vclasses : continue
        #if s["AF"] < minmaf : continue
        #for region in foundregions : 
        ref,het,hom = [int(x) for x in s[region].split(":")]
        #print s[region]
        #print hom, het, ref
        #print (2*hom+het),(2.*(ref+het+hom))
        regionaf = (2*hom+het) / (2.*(ref+het+hom))
        #print regionaf
        if regionaf == 0. or regionaf < minmaf : continue
        #s[region] = sum([int(x) for x in s[region].split(":")[1:]])
        regseries = Series({"chrom":s["chrom"],"pos":s["pos"],"Gene":s["Gene"],
                    "AF":regionaf,"region":region,"nsamp":(hom+het)})
        allvars.append(regseries)
        #limit -= 1
        #if limit <= 0 : sys.exit(1)
        #print s[["chrom","pos","Gene"]+foundregions]

    allvars = DataFrame(allvars)
    #print allvars.head(10)
    return allvars
# END makeGeneCounts2

######################################################################
# Main
######################################################################
if __name__ == "__main__" :

    os.chdir("..")
    #Y - total number of common missense and LoF SNVs
    #X - total number of protein coding variants (including synonymous)
    #Regress Y on X then take studentized residual
    #sampleannot = read_csv("resources/annotation/patientannotation.ped", sep="\t")
    sampleannot = sampleAnnotation()

    regions = sorted(sampleannot[sampleannot["Continent"].notnull()]["Continent"].unique().tolist())
    #regionlookup = dict([[samp,regions.index(cont)] for samp,cont in sampleannot[["Individual.ID","Continent"]].values if cont in regions])

    # Make X
    #prefix = "merged.clean"
    #prefix = "onekg.clean"
    #prefix = "daily.clean"
    #prefix = "everything_set1.chr1.snp.clean"
    prefix = "variome"
    
    #prefix = "rawdata/daily/daily.clean"
    #prefix = "rawdata/merged/main/merged.clean"
    #prefix = "rawdata/mevariome/main/variome.clean"
    #prefix = "rawdata/mergedaly/main/meceu.clean"
    prefix = "rawdata/daily/main/daily.clean"
    tregion = "Europe"
    bigX = makeGeneCounts2( prefix, tregion, vclasses=["4","3","2","1"] )
    bigY = makeGeneCounts2( prefix, tregion, vclasses=["4","3","2"], minmaf=.01 )

    print bigX.head(10)
    print "X",bigX.shape,"Y",bigY.shape

    #print "#"*70,"\n", "Running :",cont
    subX = DataFrame({'Xcount' : bigX.groupby( "Gene" ).size()}).reset_index()
    #if len(subX) == 0 : print "Region:",cont,"not found"; continue
    subY = DataFrame({'Ycount' : bigY.groupby( "Gene" ).size()}).reset_index()
    genecounts = merge(subX,subY,on="Gene")
    print genecounts.head(10)
    crvis = calculateRVIS( genecounts, prefix=tregion ) 
    print crvis.head(3)
    sys.exit(1)

    #bigX = makeGeneCounts1( prefix, regions, vclasses=["5","4","3","2","1"] )
    #bigY = makeGeneCounts1( prefix, regions, vclasses=["5","4","3"], maf=.1 )
    #bigX = makeGeneCounts( prefix, regions, regionlookup, vclasses=["5","4","3","2","1"] )
    #bigY = makeGeneCounts( prefix, regions, regionlookup, vclasses=["5","4","3","2"], maf=.05 )

    #print bigX.head(10)
    #print bigY.head(10)
    #print bigX.shape, bigY.shape
    #sys.exit(1)

    foundregions = [x for x in regions if x in bigX]
    for cont in foundregions :
        print "#"*70,"\n", "Running :",cont
        subX = DataFrame({'Xcount' : bigX[bigX[cont]>0].groupby( ["Gene"] ).size()}).reset_index()
        if len(subX) == 0 : print "Region:",cont,"not found"; continue
        subY = DataFrame({'Ycount' : bigY[bigY[cont]>0].groupby( ["Gene"] ).size()}).reset_index()
        genecounts = merge(subX,subY,on="Gene")
        print genecounts.head(10)
        crvis = calculateRVIS( genecounts, prefix=cont ) 
        print crvis.head(3)

    bigX["total"] = bigX[foundregions].apply(sum,axis=1)
    bigY["total"] = bigY[foundregions].apply(sum,axis=1)
    subX = DataFrame({'Xcount' : bigX.groupby( ["Gene"] ).size()}).reset_index()
    subY = DataFrame({'Ycount' : bigY.groupby( ["Gene"] ).size()}).reset_index()
    genecounts = merge(subX,subY,on="Gene")
    print genecounts.head(10)
    crvis = calculateRVIS( genecounts, prefix="total" ) 
 
    sys.exit(1)

################################################################################
# END MAIN

    #bigX = makeGeneCounts( prefix, regions, regionlookup, vclasses=["HIGH","MODERATE","LOW"] )
    #bigY = makeGeneCounts( prefix, regions, regionlookup, vclasses=["HIGH","MODERATE"], maf=.1 )
 
        #DataFrame({'count' : df1.groupby( [ "Name", "City"] ).size()}).reset_index()
        #subY = bigY[["chrom","pos","gene",cont]]
        #subY.columns = ["chrom","pos","gene","Ycount"]
        #subX = bigX[["chrom","pos","gene",cont]]
        #subX.columns = ["chrom","pos","gene","Xcount"]
        #if subX.Xcount.sum() == 0 : print "Region:",cont,"not found"; continue

    bigX = readVariantSet(statsfile, prefix, sampleannot, ["HIGH","MODERATE","LOW"])
    bigY = readVariantSet(statsfile, prefix, sampleannot, ["HIGH","MODERATE"], force=False)

    
    bigX_annot = merge(bigX, sampleannot[["Individual.ID","ethnicity","continent","country"]], left_on="Sample",right_on="Individual.ID",how="left").groupby("continent")
    bigY_annot = merge(bigY, sampleannot[["Individual.ID","ethnicity","continent","country"]], left_on="Sample",right_on="Individual.ID",how="left").groupby("continent")

    print bigY_annot.indices.keys()
    print bigX_annot.indices.keys()

    for region in bigX_annot.indices.keys() :
        print "Region:",region
        if region in ["Oceania","East Asia","Africa"] : continue
        cbigx = bigX_annot.get_group(region)
        cbigy = bigY_annot.get_group(region)
        crvis = calculateRVIS( cbigx, cbigy,prefix=region )
        crvis["Region"] = region
        print crvis.head(3)

    crvis = calculateRVIS( cbigx, cbigy )
    sys.exit(1)
# END 
    #print rstats.coef(model_x)

    # eruption.lm = lm(eruptions ~ waiting, data=faithful)
    # eruption.res = resid(eruption.lm)
    # eruption.qtl = quantile(eruption.res,c(.02,.98))
    # faithful$Colour = "black"
    # faithful$Colour[eruption.res <= eruption.qtl[1]]="red"
    # faithful$Colour[eruption.res >= eruption.qtl[2]]="blue"
    # plot( faithful$waiting, faithful$eruptions, col=faithful$Colour )
    # abline( eruption.lm )

    #print "Xcount null?:",allcounts[allcounts.Xcount.isnull()].shape
    #print "Ycount null?:",allcounts[allcounts.Ycount.isnull()].shape

    rvisdata = read_csv("resources/rvis_data.txt", sep="\t")
    print rvisdata.head(3)
    merged = merge( allcounts, rvisdata, on="gene", how="left" )
    print "Failed:",merged[merged.RVIS.isnull()].shape
    print "Success!:",merged[merged.RVIS.notnull()].shape
    print merged[merged.RVIS.isnull()].head(4)
    print merged.shape
    merged = merge( allcounts, rvisdata, on="gene" )
    print merged.shape

    print merged['StdResid'].corr(merged['RVIS'])
    print merged['StdResid'].corr(merged['RVIS'], method='spearman')

    merged["Diff"] = merged['StdResid'] - merged['RVIS']
    mostpos = merged.Diff.quantile(.95)
    mostneg = merged.Diff.quantile(.05)

    r_dataframe = com.convert_to_r_dataframe(merged)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x = "RVIS",y="StdResid" ) + \
                ggplot2.geom_jitter() + \
                ggplot2.ggtitle("Comparison of RVIS scores") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.scale_x_continuous("External RVIS")+ \
                ggplot2.scale_y_continuous("Internal RVIS") + \
                ggplot2.stat_smooth(method="lm", se=False)

    print "Writing file test1.png"
    grdevices.png("test1.png")
    p.plot()
    grdevices.dev_off()

    #q5 <- quantile(x,.05)
    #q95 <- quantile(x,.95)
    #medx <- median(x)
    #x.dens <- density(x)
    #df.dens <- data.frame(x = x.dens$x, y = x.dens$y)
    #p + geom_area(data = subset(df.dens, x >= q5 & x <= q95), 
                            #aes(x=x,y=y), fill = 'blue') +
    #geom_vline(xintercept = medx)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x = "Diff" ) + \
                ggplot2.geom_density() + \
                ggplot2.ggtitle("Density of RVIS difference\n(Internal-External)") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.geom_vline(xintercept=mostpos,colour="blue") + \
                ggplot2.geom_vline(xintercept=mostneg,colour="red")
                #ggplot2.scale_x_continuous("External RVIS")+ \
                #ggplot2.scale_y_continuous("Internal RVIS") + \
                #ggplot2.stat_smooth(method="lm", se=False)

    grdevices.png("test2.png")
    p.plot()
    grdevices.dev_off()

    negdata = merged[merged.Diff <= mostneg]
    negdata.to_csv("results/negativegenes.txt",sep="\t",index=False)
    sub1 = negdata[['gene','RVIS']]
    sub1.columns = ['gene','score']
    sub1['Type'] = 'RVIS'
    sub2 = negdata[['gene','StdResid']]
    sub2.columns = ['gene','score']
    sub2['Type'] = 'MEvariome'
    neggenes = concat( [sub1, sub2], ignore_index=True )

    print "$"*70
    print sub1.head(3)
    print sub2.head(3)
    print neggenes.head(10)

    posdata = merged[merged.Diff >= mostpos]
    posdata.to_csv("results/positivegenes.txt",sep="\t",index=False)
    sub1 = posdata[['gene','RVIS']]
    sub1.columns = ['gene','score']
    sub1['Type'] = 'RVIS'
    sub2 = posdata[['gene','StdResid']]
    sub2.columns = ['gene','score']
    sub2['Type'] = 'MEvariome'
    posgenes = concat( [sub1, sub2], ignore_index=True )

    print "$"*70
    print "Neg:",neggenes.shape
    print "Pos:",posgenes.shape
    r_dataframe = com.convert_to_r_dataframe(neggenes)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="factor(gene)",y="score",fill="factor(Type)" ) + \
                ggplot2.geom_bar(stat = "identity",position="dodge") + \
                ggplot2.ggtitle("Comparison of Internal and External RVIS") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                ggplot2.scale_y_continuous("RVIS")

    print "Writing file test3.png"
    grdevices.png("test3.png")
    p.plot()
    grdevices.dev_off()

    r_dataframe = com.convert_to_r_dataframe(posgenes)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="factor(gene)",y="score",fill="factor(Type)" ) + \
                ggplot2.geom_bar(stat = "identity",position="dodge") + \
                ggplot2.ggtitle("Comparison of Internal and External RVIS") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                ggplot2.scale_y_continuous("RVIS")

    grdevices.png("test4.png")
    p.plot()
    grdevices.dev_off()


    knowngenes = read_csv("/home/escott/nextjensystem/data/omim/genelists/omimhits.txt",sep="\t",names=["gene","entrez","MIM","Class"])
    geneclasses = read_csv("resources/allgene_classes.tsv",sep="\t",names=["gene","entrez","Class"])
    geneclasses = concat([geneclasses, knowngenes])
    results = []
    for tcolumn in ["StdResid","RVIS"] :
        negset = merged[merged[tcolumn] < 0]
        posset = merged[merged[tcolumn] >= 0]
        for gclass, gdata in geneclasses.groupby('Class') :
            print gclass
            neginclass = sum(negset.gene.isin(gdata.gene))
            posinclass = sum(posset.gene.isin(gdata.gene))
            oddsratio,pvalue = stats.fisher_exact( [[neginclass, posinclass],[len(negset) - neginclass, len(posset) - posinclass]])
            results.append([tcolumn,gclass,oddsratio,pvalue])

    classsig = DataFrame( results, columns=['RVIS_type','Type','Oddsratio','Pvalue'] ).sort("Pvalue")
    classsig["Pval_cor"] = list(rstats.p_adjust(robjects.FloatVector(classsig['Pvalue'].tolist()), method = 'BH'))
    classsig["-log(Pvalue)"] = classsig['Pval_cor'].map(lambda x: -math.log(x))
    print classsig.head(3)

    classsig.to_csv("results/classpvals.txt",sep="\t",index=False)
    r_dataframe = com.convert_to_r_dataframe(classsig)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="factor(Type)",y="-log(Pvalue)" ) + \
                ggplot2.geom_bar(stat = "identity",position="dodge") + \
                ggplot2.ggtitle("Comparison of Internal and External RVIS") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                ggplot2.facet_grid( robjects.Formula('RVIS_type ~ .') )
                #ggplot2.scale_y_continuous("RVIS") + \

    print "Writing file test5.png"
    grdevices.png("test5.png")
    p.plot()
    grdevices.dev_off()


    # group by gene to get sizes.

    #vgroups = highaf.groupby("gene")
    #genesizes = vgroups.size()
    #genes = genesizes[ genesizes > 1 ] # Select genes with more than 1 var
    #print genesizes
    #print genes

    #allPatVarStats( outfiles, outdir="./classes" )

# END Main
