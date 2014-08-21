#!/usr/bin/python
# coding: utf-8

import os, sys
import csv
import commands
import getopt
import math
import gzip
import string
import re
import subprocess
import numpy as np
from pandas import *
from scipy import stats
#from ggplot import *

#from housekeeping import *
import housekeeping as hk
import patientInfo

forceFlag=False

import rpy2.robjects as robjects
import pandas.rpy.common as com
from rpy2.robjects.packages import importr
from rpy2.robjects.lib import grid
from rpy2.robjects.lib import ggplot2


r = robjects.r
rprint = robjects.globalenv.get("print")
rstats = importr('stats')
grdevices = importr('grDevices')
gridextra = importr('gridExtra')
rcolorbrewer = importr('RColorBrewer')
#base = importr('base')
#head = robjects.r['head']
#summary = robjects.r['summary']

utils = importr('utils')
#ggplot2 = importr('ggplot2')

def makeLargePalette(ncols=12) :
    set1cols = list(rcolorbrewer.brewer_pal(9,"Set1"))
    set3cols = list(rcolorbrewer.brewer_pal(12,"Set3"))
    set2cols = list(rcolorbrewer.brewer_pal(8,"Set2"))
    allcols = set1cols+set2cols+set3cols
    return robjects.StrVector(allcols[:ncols])
# END makeLargePallete

#print robjects.r('packageVersion("ggplot2")')

#--------------------------------------------------------------------#
#                               Annotation                           #
#--------------------------------------------------------------------#
mytheme = {
            'panel.background':ggplot2.element_rect(fill='white',colour='white'),
            'axis.text':ggplot2.element_text(colour="black",size=15),
            'axis.title':ggplot2.element_text(colour="black",size=15),
            'plot.title':ggplot2.element_text(face="bold", size=20,colour="black"),
            'panel.grid.minor':ggplot2.element_blank(),
            'panel.grid.major':ggplot2.element_blank(),
            'strip.text.y':ggplot2.element_text(colour="black",face="bold",size=15),
            'strip.text.x':ggplot2.element_text(colour="black",face="bold",size=15)
            }
            #'axis.line':ggplot2.theme_line(size = 1.2, colour="black"),
            #'panel.grid.major':ggplot2.theme_line(colour = "grey90"),
pointtheme = {
            'panel.background':ggplot2.element_rect(fill='white',colour='black',size=2),
            'axis.text':ggplot2.element_text(colour="black",size=15),
            'axis.title':ggplot2.element_text(colour="black",size=15),
            'plot.title':ggplot2.element_text(face="bold", size=20,colour="black"),
            'panel.grid.major':ggplot2.element_blank(),
            'panel.grid.minor':ggplot2.element_blank(),
            'axis.ticks':ggplot2.ggplot2.element_line(colour="black"),
            'strip.text.y':ggplot2.element_text(colour="black",face="bold",size=15),
            'strip.text.x':ggplot2.element_text(colour="black",face="bold",size=15)
            }
admixtheme = {
            'panel.background':ggplot2.element_rect(fill='white',colour='white'),
            'axis.text.x': ggplot2.element_blank(),
            'axis.text.y': ggplot2.element_blank(),
            'axis.title.x':ggplot2.element_blank(),
            'axis.title.y':ggplot2.element_blank(),
            'legend.position':"none",
            'strip.text.x':ggplot2.element_text(size=6, colour="black",angle=0),
            'strip.text.y':ggplot2.element_text(size=12, face="bold", colour="black"),
            'strip.background':ggplot2.element_rect(colour="white", fill="white"),
            'axis.ticks':ggplot2.element_blank()
            }
            #'strip.text':ggplot2.element_text(size=8, colour="blue",angle=90),

pointtheme_nolegend = {
            'panel.background':ggplot2.element_blank(),
            'axis.text':ggplot2.element_text(colour="black",size=15),
            'axis.title':ggplot2.element_text(colour="black",size=15),
            'plot.title':ggplot2.element_text(face="bold", size=20,colour="black"),
            'panel.grid.minor':ggplot2.element_blank(),
            'panel.grid.major':ggplot2.element_blank(),
            'legend.position':"none",
            'axis.text.x': ggplot2.element_text(angle = 45, hjust=1, vjust=1),
            'strip.text.y':ggplot2.element_text(colour="black",face="bold",size=15,angle=-90),
            'strip.text.x':ggplot2.element_text(colour="black",face="bold",size=15),
            }
            #'strip.background':ggplot2.element_rect(colour="white", fill="white")
            #'axis.title':ggplot2.element_blank(),

#'panel.background':ggplot2.element_rect(colour = "black"),
            #'panel.grid.minor':ggplot2.element_blank(),


#qq = function(pvector, ptitle="Quantile-quantile plot of p-values", spartan=F){
    #o <- -log10( sort(pvector,decreasing=F) )
    #e <- -log10( 1:length(o)/length(o) )
    #plot=(qplot(e,o, main=ptitle) + 
            #stat_abline(intercept=0,slope=1, col="red")+
            #scale_y_continuous(name=expression(Observed~~-log[10](italic(p))),limits=c(0,max(o)))+
            #scale_x_continuous(name="something",limits=c(0,max(e))))
    #if (spartan) 
        #plot=plot+opts(panel.backgorund=theme_rect(col="grey50"),panel.grid.minor=theme_blank())
    #plot 
#}

ANNOTATIONFILE = "/home/escott/workspace/variome/resources/annotation/patientannotation.ped"

target_geographic_regions_me = (["Northwest Africa","Northeast Africa","Arabian Peninsula",
                                        "Syrian Desert","Turkish Peninsula","Central Asia"])

target_geographic_regions = ([ 'Africa','Northwest.Africa','Northeast.Africa',
                                 'Arabian.Peninsula','Syrian.Desert','Turkish.Peninsula',
                                  'Central.Asia','Europe','East.Asia','America','Oceania' ])

target_me_countries = (["Morocco","Algeria","Tunisia","Libya","Egypt",
                        "Saudi Arabia","Oman","Qatar","UAE","Yemen",
                        "Jordan","Palestine","Lebanon","Syria","Kuwait",
                        "Iraq","Turkey"])

target_ca_countries = (["Iran","Pakistan","Afganistan"])

target_continent_regions = (["Middle.East","South.Asia","Europe","Africa","East.Asia"])

target_continents = (["Middle.East","South.Asia","Europe","Africa","East.Asia"])
#target_geographic_regions2 = (["Northwest Africa","Northeast Africa","Arabian Peninsula",
                                    #"Syrian Desert","Turkish Peninsula","Central Asia"]+
                             #onekg_ethnicities)

target_geographic_regions2 = (["YRI","LWK", "Northwest Africa","Northeast Africa",
                               "Arabian Peninsula", "Syrian Desert","Turkish Peninsula",
                               "Central Asia","TSI","IBS","GBR","CEU","FIN",
                               "CHB","CHS","JPT"])#,"CLM","ASW","MXL","PUR",

target_ethnicities_all = ([ 'MbutiPygmy', 'BiakaPygmy', 'Makrani', 'Mozabite', 'Bedouin',
                 'Palestinian', 'Druze', 'Adygei', 'Uygur', 'Pathan', 'Tuscan',
                 'Brahui', 'Balochi', 'Burusho', 'Sindhi','Sardinian', 'Italian',
                 'Basque','French',  'Russian','Orcadian', 'Hazara',  'Maya'  'Miao', 
                 'Mongola', 'Tujia', 'Melanesian', 'Cambodian' ])

target_me_ethnicities = (["Bedouin","Adygei","Druze","Mozabite","Palestinian"])

onekg_ethnicities = (["ASW","CEU","CHB","CHS","CLM","FIN","GBR","IBS",
                      "JPT","LWK","MXL","PUR","TSI","YRI"])


def geographicRegions( sampleannot ):
    fixed = []
    blacklist = read_csv("resources/blacklist.txt",header=None,names=["toremove"])
    for idx, row in sampleannot.iterrows() :
        country = row["Origin"]
        ethnicity = row["ethnicity"]
        continent = row["Continent"]
        Ethnicity = row["Ethnicity"]
        if row["Individual.ID"] in blacklist.toremove.tolist() :
            fixed.append("Unknown")
        elif ((country in ["Morocco","Algeria","Tunisia"] or ethnicity == "Mozabite") and 
              (country not in ["Senegal","Togo","Cameroun","Angola"])):
            fixed.append("Northwest Africa")
        elif country in ["Egypt","Libya"] : 
            fixed.append("Northeast Africa")
        elif (country in ["Saudi Arabia","Yemen","UAE","Oman","Qatar","Bahrain","Kuwait"]): 
            fixed.append("Arabian Peninsula")
        elif ((country in ["Iran","Pakistan","Afganistan"] or
              ethnicity in ["Balochi", "Brahui", "Burusho", "Hazara", "Kalash", 
                             "Makrani", "Pashtun", "Pathan", "Sindhi"]) and
             (Ethnicity not in ["Hispanic","Dominican","Asian-Indian"])):
            fixed.append( "Central Asia")
        elif (country in ["Iraq","Syria","Lebanon","Palestine","Jordan"] or
                ethnicity in ["Druze"]):
            fixed.append( "Syrian Desert")
        elif country in ["Turkey"] or ethnicity in ["Adygei","Yakut"]:
            fixed.append("Turkish Peninsula") 
        elif (continent == "Middle East" and ethnicity == "Bedouin" and 
              row["Source"] in ["BROAD","Casanova","Fowzan","Frazer","Murat"]) :
            fixed.append( "Northeast Africa" )
        else : 
            if continent == "Middle East" : 
                #print "Err: unknown:", row
                fixed.append("Unknown")
            else :fixed.append(continent)
    return fixed
# geographicRegions

def geographicRegions2( sampleannot ):
    fixed = []
    blacklist = read_csv("resources/blacklist.txt",header=None,names=["toremove"])
    for idx, row in sampleannot.iterrows() :
        country = row["Origin"]
        ethnicity = row["ethnicity"]
        Ethnicity = row["Ethnicity"]
        continent = row["Continent"]
        if Ethnicity in onekg_ethnicities :
            fixed.append( Ethnicity )
        elif row["Individual.ID"] in blacklist.toremove.tolist() :
            fixed.append("Unknown")
        elif ((country in ["Morocco","Algeria","Tunisia"] or ethnicity == "Mozabite") and 
              (country not in ["Senegal","Togo","Cameroun","Angola"])):
            fixed.append("Northwest Africa")
        elif country in ["Egypt","Libya"] : 
            fixed.append("Northeast Africa")
        elif (country in ["Saudi Arabia","Yemen","UAE","Oman","Qatar","Bahrain","Kuwait"]) :
            #or ethnicity in ["Palestinian"]):
            fixed.append("Arabian Peninsula")
        elif ((country in ["Iran","Pakistan","Afganistan"] or
              ethnicity in ["Balochi", "Brahui", "Burusho", "Hazara", "Kalash", 
                             "Makrani", "Pashtun", "Pathan", "Sindhi"]) and
             (Ethnicity not in ["Hispanic","Dominican","Asian-Indian"])):
            fixed.append( "Central Asia")
        elif (country in ["Iraq","Syria","Lebanon","Palestine","Jordan"] or
              ethnicity == "Druze" ):
            fixed.append( "Syrian Desert")
        elif country in ["Turkey"] or ethnicity in ["Adygei","Yakut"]:
            fixed.append("Turkish Peninsula") 
        elif (continent == "Middle East" and ethnicity == "Bedouin" and 
              row["Source"] in ["BROAD","Casanova","Fowzan","Frazer","Murat"]) :
            fixed.append( "Northeast Africa" )
        else : 
            if continent == "Middle East" : 
                #print "Err: unknown:", row
                fixed.append("Unknown")
            else :fixed.append(continent)

    #print fixed
    return fixed
# geographicRegions2

def fixRLevels( r_dataframe, column, tlevels ):
    replace = robjects.FactorVector( r_dataframe.rx2(column), 
                                    levels=robjects.StrVector(tlevels) )
    allcolumns = list(r_dataframe.colnames)
    allcolumns = [ x+"_old" if x == column else x for x in allcolumns]
    new_r_df = r_dataframe.cbind(replace)
    new_r_df.colnames = robjects.StrVector(allcolumns+[column])
    return new_r_df
# END fixRLevels

################################################################################
# qqplot
################################################################################
def qqplot( pvector, ptitle="Quantile-quantile plot of p-values", spartan=False, 
           figname="./results/test.png", pfactor=None):
    print "Running qqplot"

    #o = robjects.FloatVector([-math.log10(x) for x in sorted(pvector)])
    #veclen = float(len(pvector))
    #e = robjects.FloatVector([-math.log10(i/veclen) for i in range(1,len(pvector))])
    if pfactor is not None : 
        assert len(pvector) == len(pfactor)
        pvaldf = DataFrame({'pvec':pvector,'Group':pfactor})
        print pvaldf.Group.unique()
        obsexp = []
        for grp, data in pvaldf.groupby("Group") :
            veclen = float(len(data.pvec))
            tmp = DataFrame({
                'obs':[-math.log10(x) for x in sorted(data.pvec)],
                'exp':[-math.log10(i/veclen) for i in range(1,len(data.pvec)+1)],
                'Group':grp})
            obsexp.append(tmp)
        obsexp = concat( obsexp ).reset_index(drop=True)
        print obsexp.Group.unique()
    else :
        veclen = float(len(pvector))
        obsexp = DataFrame({'obs':[-math.log10(x) for x in sorted(pvector)],
                        'exp':[-math.log10(i/veclen) for i in range(1,len(pvector)+1)]})

    print "Current DF:"
    print obsexp.head(10)
    r_dataframe = com.convert_to_r_dataframe(obsexp)
    if pfactor is not None :
        p = (ggplot2.ggplot(r_dataframe) +
             ggplot2.aes_string(x="exp",y="obs",colour="factor(Group)"))
    else :
        p = (ggplot2.ggplot(r_dataframe) +
             ggplot2.aes_string(x="exp",y="obs"))
    p = (p + 
         ggplot2.geom_point() +
         ggplot2.stat_abline(intercept=0,slope=1, col="red") +
         ggplot2.ggtitle(ptitle) +
         ggplot2.scale_y_continuous(name="Observed -log[10]",
                                    limits=robjects.IntVector([0,max(obsexp.obs)])) +
         ggplot2.scale_x_continuous(name="Expected -log[10]",
                                    limits=robjects.IntVector([0,max(obsexp.exp)])) +
         ggplot2.theme(**mytheme) )

    print "Writing figure:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END qqplot
#plot=(ggplot2.qplot(e,o,)
      #ggplot2.ggtitle( ptitle ) +
#expression(Observed~~-log[10](italic(p)))
#if (spartan) 
    #plot=plot+opts(panel.backgorund=theme_rect(col="grey50"),panel.grid.minor=theme_blank())

######################################################################
# parseClinVar
######################################################################
def parseClinVar( vfilter=5 ):
    clinvarvcf = "/home/escott/snpEff/clinvar_20140303.vcf.gz"
    clindata = []
    for row in csv.reader( hk.CommentStripper(hk.conditionalOpen(clinvarvcf)), delimiter="\t" ) :
        #if row[0].startswith("#") or len(row) == 0 : continue
        chrom,pos,vid,ref,alt,qual,vfilter,info = row
        infodict = dict([x.split("=") if x.find("=") >0 else [x,''] for x in info.split(";")])
        clinsig = infodict["CLNSIG"]
        if clinsig.find("|") > 0 or clinsig.find(",")  > 0:
            tmp = [ int(x) for x in re.split("\W+", infodict["CLNSIG"]) if int(x) != 255]
            clinsig = max(tmp) if len(tmp) > 0 else 255
        newrow = [chrom,pos,vid,ref,alt,clinsig]
        clindata.append(newrow)
    clindata = DataFrame(clindata, columns=["chrom","pos","id","ref","alt","clinsig"])
    clindata.index = "chr"+clindata.chrom+":"+clindata.pos#+":"+clindata.ref+":"+clindata.alt
    if vfilter != None : clindata = clindata[clindata.clinsig == vfilter]
    return clindata
# END parseClinVar

def sampleAnnotation(targetpats=None, annotationfile=ANNOTATIONFILE ) :
    sampleannot = read_csv(annotationfile, sep="\t")
    hgdp = read_csv("/home/escott/workspace/variome/resources/HGDP_ethnicgroups.txt",sep="\t")
    hgdp.columns = ["ethnicity","Count","Continent","Country"]
    keepcols =  ["Family.ID","Individual.ID","Paternal.ID","Maternal.ID",
                 "Phenotype","Gender", "Consang", "Plate","Source","Ethnicity",
                 "Origin","ethnicity","continent","country"]
    sampleannot = merge(sampleannot[keepcols], hgdp[["ethnicity","Continent","Country"]], 
                        on="ethnicity",how="left")
    print "annot len:",sampleannot.shape
    sampleannot.index = sampleannot["Individual.ID"].tolist()
    sampleannot["GeographicRegions"] = geographicRegions(sampleannot)
    sampleannot["GeographicRegions2"] = geographicRegions2(sampleannot)
    sampleannot["Continent2"] = ["Middle East" if x in target_geographic_regions_me 
                                 else x for x in sampleannot["GeographicRegions"]]
    #( sampleannot.Origin.tolist(), #sampleannot.ethnicity.tolist(), #sampleannot.continent)
    #print sampleannot[sampleannot.Continent2 == "South Asia"]
    if targetpats is not None :
        print "Targetpats",len(targetpats),targetpats[:10]
        sampleannot = sampleannot[sampleannot["Individual.ID"].isin( targetpats )]
    return sampleannot
# END sampleAnnotation

def addSampleAnnotation( targetdf, mergecol=None, update=False ):
    patientdata = sampleAnnotation()

    #print targetdf.head(10)
    tids = patientdata["Individual.ID"].tolist()
    patientdata[mergecol] = [x[:x.find(".combined")] 
                          if x.find(".combined") > 0 
                          else x
                          for x in patientdata["Individual.ID"].tolist()]
                          #and x[:x.find(".combined")] in tids

    if update : tids = targetdf.index.tolist()
    else : tids = targetdf[mergecol].tolist() 
    patientdata["fixed_id"] = [x[:x.find("_")] 
                               if x[:x.find("_")] in tids and x.find("JL") == 0
                               else x
                               for x in patientdata["Individual.ID"]]
    keepcols = ["fixed_id", "Family.ID", "Gender", "Phenotype", "Origin",
                "Continent","GeographicRegions", "GeographicRegions2", "Continent2",
                "Ethnicity", "Source", "ethnicity", "continent", "country"] 
                #"Individual.ID", 
                #"Paternal.ID", "Maternal.ID",
                #"Consang", #"Plate", "continent_prob", , "country_prob" "ethnicity_prob",
    if update : 
        if mergecol is not None : targetdf.index = targetdf[mergecol]
        patientdata.index = patientdata.fixed_id
        targetdf.update( patientdata )
        return targetdf

    assert mergecol is not None
    alldata = merge(targetdf, patientdata[keepcols], left_on=mergecol, 
                    right_on="fixed_id", how="left")
    return alldata
# ENd addSampleAnnotation

######################################################################
# readVariantSet
######################################################################
def readVariantSet( prefix, vclasses=["HIGH","MODERATE","LOW"], force=False ):
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
# END readVariantSet

######################################################################
# readVariantSet1
######################################################################
def readVariantSet1( statsfile, prefix, sampleannot, minaf=1, maxaf=0, vclasses=["HIGH"], force=False ):
    print "Running readVariantSet1"
    limit = 10
    if os.path.exists(statsfile) and not force: return read_csv(statsfile, sep="\t")

    knownsamples = sampleannot.index.tolist()
    allvars = []
    for vclass in vclasses :
        regioncounts = {}
        #targetfile = "./results/classes/"+prefix+".annot_"+vclass+".tsv"
        targetfile = "./results/classes/"+prefix+"_"+str(vclass)+".tsv"
        #print "Reading:",targetfile
        FM = open(targetfile,"r")
        assert os.path.exists(targetfile)

        for line in FM :
            if len(line) > 200 : continue
            var = line.rstrip().split("\t")
            if var[0] == "Sample" or var[0] not in knownsamples : continue
            if len(var) == 1 : print limit,var; sys.exit(1)
            region = sampleannot.loc[var[0]]["Continent"]
            sample = var[0]
            af = var[7]
            if float(af) >= minaf : continue
            if float(af) <= maxaf : continue
            geno = 0 if var[1] == "het" else 1
            if not regioncounts.has_key(region) : regioncounts[region] = {}
            if not regioncounts[region].has_key(sample) : regioncounts[region][sample] = [0,0]

            regioncounts[region][sample][geno] += 1

        classvars = []
        for region in regioncounts :
            for samp in regioncounts[region] :
                #print [samp,vclass,region,regioncounts[region][samp][0],regioncounts[region][samp][1],sum(regioncounts[region][samp])]
                classvars.append([samp,vclass,region,regioncounts[region][samp][0],regioncounts[region][samp][1],sum(regioncounts[region][samp])])
        allvars.append( DataFrame( classvars, columns=["Individual.ID","class","region","het","hom","total"]) )

    print "Length:",len(allvars)
    if len(allvars) == 1 : allvars = allvars[0]
    else : 
        allvars = concat( allvars )
        allvars = allvars.groupby(["Individual.ID","region"],as_index=False).aggregate(np.sum)
        print "Length:",len(allvars)
        print allvars.head(10)
    varsannot = merge(allvars, sampleannot[["Individual.ID","Origin","ethnicity","Continent"]], on="Individual.ID")

    print "Writing:",statsfile
    varsannot.to_csv(statsfile, sep="\t", index=False)
    return varsannot
# END readVariantSet1

######################################################################
# readVariantSet2
######################################################################
def readVariantSet2( statsfile, prefix, sampleannot, vclasses=["HIGH"], force=False ):
    print "Running readVariantSet1"
    limit = 10
    if os.path.exists(statsfile) and not force: return read_csv(statsfile, sep="\t")
    #for var in csv.reader(open(targetfile),delimiter="\t"):

    knownsamples = sampleannot.index.tolist()
    regioncounts = {}
    seen = []
    for vclass in vclasses :
        targetfile = "./results/classes/"+prefix+".annot_"+vclass+".tsv"
        print "Reading:",targetfile
        FM = open(targetfile,"r")
        assert os.path.exists(targetfile)

        for line in FM :
            if len(line) > 200 : continue
            var = line.rstrip().split("\t")
            if var[0] == "Sample" or var[0] not in knownsamples : continue
            if len(var) == 1 : print limit,var; sys.exit(1)
            region = sampleannot.loc[var[0]].continent
            af = float(var[7])
            key = var[2]+":"+var[3]
            print "key",key
            if key in seen : continue
            if not regioncounts.has_key(region) : regioncounts[region] = {}
            if not regioncounts[region].has_key(af) : regioncounts[region][af] = 0
            regioncounts[region][af] += 1

    alldata = []
    for region in regioncounts:
        for af in sorted(regioncounts[region]) :
            alldata.append([region,af,regioncounts[region][af]])
    alldata = DataFrame( alldata, columns=["Region","AF","count"] )
    return alldata
# END readVariantSet2


if __name__ == "__main__" :

    patientlistfile = "/home/escott/workspace/variome/rawdata/merge1kg/allsamples.txt"
    parseClinVar( vfilter=5 )   
    samples = read_csv( patientlistfile, sep="\t", header=None, names=["sample"] )
    annotated = addSampleAnnotation( samples, mergecol="sample" )
    print annotated.head()

    print annotated.Ethnicity.unique()
    annotated[annotated.GeographicRegions2 == "Northwest Africa"][["sample","Origin","GeographicRegions","Ethnicity","Source","ethnicity"]].to_csv("ethnicitytest.txt",sep="\t")

# END MAIN




