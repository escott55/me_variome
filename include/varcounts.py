#!/usr/bin/env python

#-------------------------------------------------------------------#
#                       Import Packages                             #
#-------------------------------------------------------------------#
import os, sys
import getopt
import csv
import re
import subprocess

import housekeeping as hk
import popgencommands as pop
from localglobals import *

def plotHomozygosity( samplecounts, outdir, basename ): 
    homcols = [vc+"_hom" for vc in vclasses]
    hetcols = [vc+"_het" for vc in vclasses]
    samplecounts["homtotal"] = samplecounts[homcols].sum(axis=1)
    samplecounts["hettotal"] = samplecounts[hetcols].sum(axis=1)
    samplecounts["allvars"] = samplecounts[homcols+hetcols].sum(axis=1)
    samplecounts["homozygosity"] = samplecounts["homtotal"] / samplecounts["allvars"].map(float)

    r_dataframe = com.convert_to_r_dataframe(samplecounts)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(Region)", y="homozygosity", fill="factor(Region)") +
                ggplot2.geom_boxplot() +
                ggplot2.ggtitle("Boxplot of Homozygosity") +
                ggplot2.theme(**mytheme) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.scale_y_continuous("Homozygosity") +
                ggplot2.scale_x_discrete("Geographic Region") )

    figurename = "%s/%s_homozygosity.png" % (outdir,basename)
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
    print samplecounts[["Sample","Region","Country","GeographicRegions","homtotal","homozygosity"]].head(10)
    scounts = melt(samplecounts,id_vars=["Sample","Region","GeographicRegions","homozygosity","Source"],
                  value_vars=vclasses)

    print "Basename:", basename
    if basename.find( "meceu" ) >= 0:
        #scounts = scounts[scounts.Region.isin(["Europe","Middle East"])] #,"South Asia"])]
        #scounts_sub = scounts[scounts.Source.isin(["BROAD","Daly"])]
        scounts_sub = scounts[((scounts.Region=="Europe") & (scounts.Source=="Daly")) |
                          ((scounts.Region=="Middle East") & (scounts.Source=="BROAD"))] #,"South Asia"
    elif basename.find( "1000G" ) >= 0 :
        print "Using 1000 Genomes pops!"
        scounts_sub = scounts[((scounts.Region=="Europe") & (scounts.Source=="1KG")) |
                              ((scounts.Region=="Africa") & (scounts.Source=="1KG")) |
                          ((scounts.Region=="Middle East") & (scounts.Source=="BROAD"))] #,"South Asia"
    elif basename.find( "variome" ) >= 0 :
        scounts_sub = scounts[scounts.Region.isin(["Europe","Africa","Middle East"])] #,"South Asia"
    elif basename.find( "test2" ) >= 0 :
        scounts_sub = scounts[scounts.Region.isin(["Europe","Africa","Middle East"])] #,"South Asia"
    else :
        print "Error: Unknown basename",basename; sys.exit(1)

    r_dataframe = com.convert_to_r_dataframe(scounts_sub)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="homozygosity", y="value" ) +
                ggplot2.geom_point(ggplot2.aes_string(colour="factor(Source)")) +
                ggplot2.ggtitle("Boxplot of Variants by Variant Class") +
                ggplot2.theme(**mytheme) +
                ggplot2.scale_y_continuous("Class Burden") +
                ggplot2.scale_x_continuous("Sample Homozygosity") +
                ggplot2.facet_grid( robjects.Formula('variable ~ Region') , scale="free") )
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                #ggplot2.scale_x_discrete("Geographic Region") +
                #ggplot2.scale_y_log10() + ggplot2.xlim(0., 1.) + \
                #ggplot2.geom_histogram(breaks=robjects.FloatVector([0,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])) + \
                #ggplot2.aes_string(fill="factor(continent)")
                #fill="factor(group)"

    figurename = "%s/%s_homozburden.png" % (outdir,basename)
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# END plotHomozygosity

def plotVariantClassCounts( scounts, outdir, basename ) :
    r_dataframe = com.convert_to_r_dataframe(scounts)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(Region)", y="value", fill="factor(Region)") +
                ggplot2.geom_boxplot() +
                ggplot2.ggtitle("Boxplot of Variants by Variant Class") +
                ggplot2.theme(**mytheme) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.scale_y_continuous("Variant Burden") +
                ggplot2.scale_x_discrete("Geographic Region") +
                ggplot2.facet_grid( robjects.Formula('variable ~ .') , scale="free") )
                #ggplot2.scale_y_log10() + ggplot2.xlim(0., 1.) + \
                #ggplot2.geom_histogram(breaks=robjects.FloatVector([0,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])) + \
                #ggplot2.aes_string(fill="factor(continent)")
                #fill="factor(group)"

    figurename = "%s/%s_sampcounts.png" % (outdir,basename)
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# ENd plotVariantClassCounts

def plotKnockOutVarCounts( samplecounts, outdir, basename ) :
    vclasses = ["ClassNMD","ClassLOF"]#,"Class5"
    for vclass in vclasses :
        samplecounts[vclass] = (samplecounts[vclass+"_het"] +
                                samplecounts[vclass+"_hom"])
    scounts = melt(samplecounts,id_vars=["Sample","Region","GeographicRegions"],
                  value_vars=vclasses)
    scount = scounts[~scounts.Region.isin(["NA","Oceania"])]
    if basename.find( "meceu" ) >= 0:
        scounts = scounts[scounts.Region.isin(["Europe","Middle East","South Asia"])]
    elif basename.find( "1000G" ) >= 0 :
        scounts = scounts[scounts.Region.isin(["Europe","Africa","Middle East"])] #,"South Asia"
    #print scounts.head(10)

    scounts = scounts[scounts.value < 50]

    r_dataframe = com.convert_to_r_dataframe(scounts)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(Region)", y="value", fill="factor(Region)") +
                ggplot2.geom_boxplot() +
                ggplot2.ggtitle("High Impact Variants by Region") +
                ggplot2.theme(**mytheme) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.scale_y_continuous("Variant Burden" ) +#, limits=robjects.IntVector((0,50))) +
                ggplot2.scale_x_discrete("Geographic Region") +
                ggplot2.facet_grid( robjects.Formula('variable ~ .'), scale="free") )
    figurename = "%s/%s_lofcounts.png" % (outdir,basename)
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# ENd plotKnockOutVarCounts

################################################################################
# plotSampleCounts
################################################################################
def plotSampleCounts( samplecounts, outdir="./results/figures/varcounts", force=False) :

    plotHomozygosity( samplecounts, outdir, basename )

    vclasses = ["Class0","Class1","Class2","Class3","Class4"]#,"Class5"
    scounts = melt(samplecounts,id_vars=["Sample","Region","GeographicRegions","Source"],
                  value_vars=vclasses)

    plotVariantClassCounts( scounts, outdir, basename ) 

    plotKnockOutVarCounts( samplecounts, outdir, basename )
# END plotSampleCounts

def isOutlier( x, mean, std ):
    return abs((x-mean)/std) > 3

def lrVarClassPlots( varcounts, dataset, tcol="Continent3" ):
    #vcounts = melt( varcounts_merge, 
                   #id_vars=["Individual.ID","percent"],
                   #value_vars=["benignroh","benignnonroh","delroh","delnonroh"])
    #vcounts["rohclass"] = ["Outside" if x.find("non") >= 0 else "Within" for x in vcounts.variable]
    #vcounts["vtype"] = ["Benign" if x.find("benign") == 0 else "Deleterious" for x in vcounts.variable]

    #vcounts_annot = addSampleAnnotation( vcounts, mergecol="Individual.ID" )
    #vcounts = varcounts[ varcounts[tcol].isin(["Europe","Middle East","Anglo-American", "Africa"]) ]
    #vcounts = varcounts[ varcounts[tcol].isin(["Middle East","Anglo-American","Europe","Africa"]) ]
    vcounts = varcounts[ varcounts[tcol].isin(["Middle East","Anglo-American"]) ]
    #vcounts = varcounts[ varcounts[tcol].isin(["TP","SD","NWA","NEA","AP","CA"]) ]
    #print vcounts_annot.head()"Africa","Europe","Middle East"

    #mean = vcounts["Class1"].mean()
    #std = vcounts["Class1"].std()
    #sdout = [isOutlier(x,mean,std) for x in vcounts.Class1]
    #outliers = vcounts[["Sample","Class1"]][sdout]
    #vcounts = vcounts[[not x for x in sdout]]

    lmdata = []
    for idx, grp in vcounts.groupby([tcol,"tclass"]) :
        lm = stats.linregress( grp["Class0"].tolist(), grp["Vcount"].tolist() )
        lmdata.append( Series(data=idx+lm, index=[tcol,"tclass","slope",
                                                 "intercept", "rval", "pval", "stderr"]) )
    
    lmdata = DataFrame(lmdata)

    lmdata["xpos"] = vcounts.Class0.min()

    for idx, subdata in lmdata.groupby("tclass"): 
        maxY = vcounts[vcounts.tclass==idx]["Vcount"].max()
        lmdata.loc[lmdata.tclass==idx,"ypos"] = [maxY*(1.-.04*(i))
                                                 for i in range(1,len(subdata[tcol])+1)]

    lmdata["Slope"] = ["%s Slope: %.3f" % (x,y) 
                       for x,y in lmdata[[tcol,"slope"]].values]

    print lmdata
    r_dataframe = com.convert_to_r_dataframe(vcounts)
    r_lm = com.convert_to_r_dataframe(lmdata)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="Class0", y="Vcount", colour="factor("+tcol+")")+
                                   #group="factor(tclass)", colour="factor(tclass)") +
                ggplot2.geom_point(alpha=.6) +
                ggplot2.geom_abline(ggplot2.aes_string(intercept="intercept",slope="slope", 
                                                       colour=tcol),
                                     data=r_lm) +
                ggplot2.geom_text(ggplot2.aes_string(label="Slope",x="xpos",y="ypos"),
                                  color="black", size=4, hjust=0, vjust=0, data=r_lm) +
                ggplot2.scale_y_continuous("Burden of Damaging Variants") +
                ggplot2.scale_x_continuous("Benign variant burden") +
                ggplot2.scale_colour_brewer("Region",palette="Set1") +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45),
                                                         'legend.position':'top'}) +
                #ggplot2.scale_colour_manual(name="Variant Location",
                                            #values=robjects.StrVector(["red","blue"]),
                                            #breaks=robjects.StrVector(["benignroh","benignnonroh"]),
                                            #labels=robjects.StrVector(["Within RoH","Outside RoH"]) ) +
                ggplot2.facet_grid(robjects.Formula('tclass ~ .'), scale="free_y") +
                #ggplot2.stat_smooth(method="lm", se=False) +
                ggplot2.theme(**pointtheme) )

                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 90)}) +
                #ggplot2.ggtitle("Benign Variant counts within\nRuns of Homozygosity") +

    figname = "results/figures/varcounts/%s_%s_lr_varcounts.pdf" % (dataset,tcol)
    print "Writing file:",figname
    grdevices.pdf(figname)
    p.plot()
    grdevices.dev_off()
# END lrVarClassPlots

######################################################################
if __name__ == "__main__":

    #optlist, args = getopt.getopt( sys.argv[1:], "bk")
    #optlist = dict(optlist)
    #bedfile = optlist.get("-b",bedfile)

    # Change running directory
    os.chdir("..")

    #vcffile = "./rawdata/merge1kg/main/me1000G.clean.vcf.gz"
    #vcffile = "./rawdata/mergedaly/meceu.X.vcf.gz"
    #vcffile = "./rawdata/merge1kg/me1000G.X.vcf.gz"
    #samplecountsfile = "./rawdata/mevariome/main/variome.clean_indiv.tsv"
    samplecountsfile = "./rawdata/mevariome/classify/variome.sfilt_indiv.tsv"
    clustfile = "./rawdata/mevariome/main/clust/variome.clean.annot"
    samplecountsfile = "./rawdata/merge1kg/main/classify/me1000G.clean.sfilt_indiv.tsv"
    clustfile = "./rawdata/merge1kg/main/clust/me1000G.clean.annot"

    samplecountsfile = "./results/other/merged_indiv.tsv"
    clustfile = "./results/other/merged.annot"
    #samplecountsfile = "./rawdata/mergedaly/main/meceu.clean_indiv.tsv"

    filepath, basename, ext = hk.getBasename( samplecountsfile )
    dataset = basename[:basename.find(".")]
    samplecounts = read_csv(samplecountsfile, sep="\t")

    if os.path.exists(clustfile) :
        sampannot = read_csv(clustfile, sep='\t')
        samplecounts = merge(samplecounts, sampannot, 
                             left_on="Sample", right_on="Individual.ID")
        mapping_regions = {"NWA":"Middle East", "NEA":"Middle East", "AP":"Middle East",
                           "SD":"Middle East", "TP":"Middle East", "CA":"Middle East",
                           "LWK":"Africa", "YRI":"Africa", "IBS":"Europe",
                           "CEU":"Europe", "TSI":"Europe", "FIN":"Europe",
                           "GBR":"Europe", "CHB":"East Asia", "CHS":"East Asia",
                           "JPT":"East Asia", "Anglo-American":"Anglo-American"}
        samplecounts["Continent3"] = [mapping_regions[x] if mapping_regions.has_key(x) 
                                      else "Unknown"
                                      for x in samplecounts.GeographicRegions3]
    else :
        samplecounts = addSampleAnnotation( samplecounts, "Sample" )

    #print samplecounts.head()

    vclasses = ["Class0","Class1","Class2","Class3","Class4"]#,"Class5"
    for vclass in vclasses :
        samplecounts[vclass] = (samplecounts[vclass+"_het"] +
                                samplecounts[vclass+"_hom"])

    # Remove outliers
    outliers = []
    for cont,sc_sub in samplecounts.groupby("Continent2") :
        for vclass in vclasses :
            mean = sc_sub[vclass].mean()
            std = sc_sub[vclass].std()
            sdout = [isOutlier(x,mean,std) for x in sc_sub[vclass]]
            outliers.append(sc_sub[["Sample","Continent2"]][sdout])
 
    outliers = concat(outliers)
    toremove = outliers.groupby("Sample").size().reset_index()
    samplecounts = samplecounts[~samplecounts.Sample.isin(toremove.Sample)]
    #plotSampleCounts( samplecounts )
    print hk.dfTable( samplecounts.Continent2 )

    print samplecounts.head()
    vclasses = ["Class1","Class2","Class3","Class4"]
    varcounts = melt(samplecounts,id_vars=["Sample","Continent3","GeographicRegions3","Source","Class0"],
                  value_vars=vclasses)
    varcounts.rename(columns={"variable":"tclass", "value":"Vcount"},inplace=True)
    #print hk.dfTable( varcounts.tclass )

    lrVarClassPlots( varcounts, dataset, tcol="GeographicRegions3" )
    lrVarClassPlots( varcounts, dataset, tcol="Continent3" )

# END MAIN
