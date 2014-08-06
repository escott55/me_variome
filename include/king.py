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


######################################################################
# make_kinshipplots
######################################################################
def make_kinshipplots( kinshipfile ) :
    targetdir, filename, suffix = hk.getBasename(kinshipfile)
    # MAKE RELATEDNESS PLOTS
    command = ["./rscripts/kinship.R",kinshipfile]
    runout = subprocess.check_output(command)
    print runout
    return
# End make_kinshipplots

################################################################################
# run_king
# Input - bed
# Process - Run king kinship analysis and king unrelated
# Output - samplelist
################################################################################
def run_king( targetfile, force=False ) :
    targetdir, filename, suffix = hk.getBasename(targetfile)

    basename = "%s/%s" % (targetdir, filename )
    unrelfile = basename+"_unrelated.txt"
    relfile = basename+"_rel.kin"
    otherrelfile = basename+"_rel.kin0"
    if os.path.exists( unrelfile ) and not force:
        print "File "+unrelfile+" already exists!"
        return unrelfile, relfile, otherrelfile

    command = ["king",
               "-b",basename+".bed",
               "--kinship", 
               "--prefix",basename+"_rel"]

    subprocess.call( command )

    command = ["king",
               "-b",basename+".bed",
               "--unrelated",
               "--degree","2",
               "--prefix",basename+"_" ]

    subprocess.call( command )

    return unrelfile, relfile, otherrelfile
# End run_king
    #patientlist = []
    #FILE = open( unrelfile )
    #reader = csv.reader(FILE,delimiter="\t")
    #for row in reader :
        #if len(row) > 0 :
            #patientlist.append(row[1])
    #return patientlist, outfile, basename+"_rel.kin"
        #patientlist = []
        #FILE = open( outfile )
        #reader = csv.reader(FILE,delimiter="\t")
        #for row in reader :
            #if len(row) > 0 :
                #patientlist.append(row[1])

def isOutlier( x, mean, std ):
    return abs((x-mean)/std) > 3

################################################################################
def popOutliers( vcffile, force=False ) : #, sampleannot
    print "Running popOutliers"
    vtargetdir, vfilename, vsuffix = hk.getBasename(vcffile)
    #vtargetdir = outprefix if outprefix is not None else originaldir
    #print vtargetdir, outprefix, vfilename
    pseqfile = os.path.join(vtargetdir, "%s_istats.txt"%vfilename)
    pseqerr = os.path.join(vtargetdir, "%s_istats.err"%vfilename)
    outlierfile = os.path.join(vtargetdir, "%s.qcoutliers"%vfilename)
    if not os.path.exists( pseqfile ) or hk.fileIsEmpty( pseqfile ) or force:
        with open(pseqfile,"wb") as out, open(pseqerr,"wb") as err:
            print " ".join(["pseq",vcffile,"i-stats"])
            subprocess.call(["pseq",vcffile,"i-stats"],stdout=out,stderr=err)
            out.close(); err.close()
    
    #print pseqerr, pseqfile
    istats = read_csv( pseqfile, sep="\t" )
    alloutliers = []
    for col in istats.columns : 
        if col == "ID" : continue
        if col in ["DP"] : continue
        mean = istats[col].mean()
        std = istats[col].std()
        outliers = istats[["ID",col]][[isOutlier(x,mean,std) 
                                       for x in istats[col].tolist()]]
        for out,val in outliers[["ID",col]].values :
            alloutliers.append([out,val,col,mean,std])
    if len(alloutliers) == 0 :
        alloutliers = DataFrame( columns=["ID","Value","Col","Mean","STD"])
    else :
        alloutliers = DataFrame( alloutliers, columns=["ID","Value","Col","Mean","STD"])

    #print "Shape:",alloutliers.shape
    print "Going to remove",len(alloutliers.ID.unique()),"samples"
    print "Writing outlier file:",outlierfile
    alloutliers.to_csv(outlierfile,sep="\t")
    #print alloutliers[alloutliers.ID == "JL0688_ACAGTG_BD2GJBACXX_L006_001"]
    return alloutliers.ID.unique().tolist()
# END popOutliers

################################################################################
def seriousClean( vcffile, rerun=False ):
    print "Running seriousClean"
    patientdata = sampleAnnotation()
    filepath, basename, suffix = hk.getBasename(vcffile)

    cleanped = "%s/%s"%(filepath,basename)
    tped = cleanped+".tped"
    if os.path.exists( cleanped+".bed" ) and not rerun :
        return cleanped+".bed"

    # Annotation Outliers
    vcfpats = patientInfo.getPats( vcffile )
    vcfpats = DataFrame({'Individual.ID':vcfpats})
    vcfpats = addSampleAnnotation( vcfpats )
    sampleannot = vcfpats[(vcfpats.Continent2.notnull() 
                              & (vcfpats.Continent2 != "Unknown"))]

    #print sampleannot.Continent2.unique()
    
    #keepfile = os.path.join(filepath, basename+".annotkeeppats")
    #sampleannot[["Individual.ID"]].to_csv(keepfile, header=None,index=None)
    annotpats = sampleannot["Individual.ID"].tolist()

    toremove="./toremove.txt"
    autoremove = []
    if os.path.exists(toremove):
        autoremove = ([x.strip() for x in open(toremove).readlines()
                       if len(x.strip()) > 0])
    
    qcoutliers = popOutliers( vcffile )
    
    filepats = patientInfo.currentFilePats( vcffile )
    keepfile = os.path.join(filepath,basename+".keep")
    OUT = open(keepfile,"wb")
    for samp in filepats :
        if samp in autoremove : print "autoremove",samp; continue
        if samp in qcoutliers : print "qcoutlier",samp; continue # skip outliers
        if samp not in annotpats : print "no annotation",samp; continue 
        OUT.write(samp+"\n")
    OUT.close()

    command = ["vcftools","--gzvcf", vcffile,
               "--remove-filtered-all",
               "--remove-indels",
               "--keep",keepfile, 
               "--maf",".005",
               "--min-alleles","2",
               "--max-alleles","2",
               "--plink-tped","--out",cleanped]

    if not os.path.exists( tped ) or rerun:
        print " ".join(command)
        out = subprocess.check_output( command )

    print patientdata.head()
    pop.modify_tped( tped, patientdata, force=rerun )
                             #calcdistances=True, shorten=True, 

    bedfile = pop.run_plink_convert(tped, force=rerun)
    return bedfile
# END seriousClean

def dfTable( targetlist ):
    from sets import Set
    return DataFrame( [[x,targetlist.count(x)] for x in Set(targetlist)], 
                     columns=["Variable","Count"]).sort("Count",ascending=0)
# END dfTable

def dfProportion( targetlist ):
    from sets import Set
    datapts = float(len(targetlist))
    return DataFrame( [[x,targetlist.count(x)/datapts] for x in Set(targetlist)], 
                     columns=["Variable","Count"]).sort("Count",ascending=0)
# END dfProportion

def plotSecondDegreeBar( interfam ) :
    seconddeg = dfTable( interfam[interfam.IBS0 > .0884].ID1.tolist()+
                   interfam[interfam.IBS0 > .0884].ID2.tolist() )

    seconddeg = addSampleAnnotation( seconddeg, mergecol="Variable" )
    seconddegfilt = seconddeg[seconddeg.Count >1]
    r_dataframe = com.convert_to_r_dataframe(seconddegfilt)
    r_dataframe = fixRLevels( r_dataframe, "Variable", seconddegfilt.Variable.tolist() )
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(Variable)",
                                   y="Count",
                                   fill="factor(Source)") +
                ggplot2.geom_bar(stat="identity") +
                ggplot2.ggtitle("Count of Interfamily 2nd Degree relationship") +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.theme(**mytheme) )
                #ggplot2.stat_smooth(method="lm", se=False) +
                ##ggplot2.facet_grid(robjects.Formula('. ~ variable')) +

    figname = "results/figures/king/%s_2nddegree_bar.png" % (basename)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
    seconddegremove = seconddegfilt.Variable.tolist()
# ENd plotSecondDegreeBar

def king_workflow( vcffile, force=False ) :
    print "Using Vcffile:",vcffile
    assert os.path.exists(vcffile)

    targetvcf = hk.copyToSubDir( vcffile, "king" )
    #filepath, filename, suffix = hk.getBasename(targetvcf)
    # Run for all data together

    bedfile = seriousClean( targetvcf, force )
    unrelfile, famrelfile, interrelfile = run_king( bedfile, force )

    return unrelfile, famrelfile, interrelfile
# END king_workflow

def plotWithinBetweenRelatedness( interfam, basename ):
    print interfam.head()
    within_sub = interfam[interfam["GeoRegion1"] == interfam["GeoRegion2"]]
    r_dataframe = com.convert_to_r_dataframe(within_sub)
    #r_dataframe = fixRLevels( r_dataframe, "Individual.ID", classcounts["Individual.ID"].tolist() )
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="IBS0",
                                   y="Kinship",
                                   colour="factor(Comparison)") +
                ggplot2.geom_point() +
                ggplot2.ggtitle("Relatedness Coefficient by Population") +
                ggplot2.scale_y_continuous("Kinship Coefficient",
                                           limits=robjects.FloatVector((-.5,.2))) +
                ggplot2.scale_x_continuous("Pr(IBS=0)",
                                           limits=robjects.FloatVector((0,.06))) +
                #ggplot2.ylim(robjects.IntVector((0,.3))) +
                ggplot2.theme(**{'legend.position':"none"}) +
                ggplot2.facet_grid(robjects.Formula('. ~ GeoRegion1')) + #, scale="free") +
                ggplot2.theme(**mytheme) )
                #ggplot2.coord_flip() +
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                #ggplot2.theme(**{'axis.text.y': ggplot2.element_blank()}) +

    figname = "results/figures/king/%s_within_points.png" % (basename)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()

    #interfam["GeographicRegions"] = "Unknown"
    #interfam = addSampleAnnotation(in, "ID1", True)
    within_sub = interfam[interfam["GeoRegion1"] == interfam["GeoRegion2"]]
    r_dataframe = com.convert_to_r_dataframe(within_sub)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(group="factor(GeoRegion1)",
                                   x="Kinship",
                                   colour="factor(GeoRegion1)") +
                ggplot2.geom_density() +
                ggplot2.ggtitle("Kinship Coefficient Density") +
                #ggplot2.scale_y_continuous("Kinship Coefficient",
                                           #limits=robjects.FloatVector((-.5,.2))) +
                #ggplot2.scale_x_continuous("Pr(IBS=0)",
                                           #limits=robjects.FloatVector((0,.06))) +
                #ggplot2.ylim(robjects.IntVector((0,.3))) +
                #ggplot2.theme(**{'legend.position':"none"}) +
                #ggplot2.facet_grid(robjects.Formula('. ~ GeoRegion1')) + #, scale="free") +
                ggplot2.theme(**mytheme) )
                #ggplot2.coord_flip() +
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                #ggplot2.theme(**{'axis.text.y': ggplot2.element_blank()}) +

    figname = "results/figures/king/%s_within_density.png" % (basename)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()


    print interfam.head()
    between_sub = interfam[interfam["GeoRegion1"] != interfam["GeoRegion2"]]
    r_dataframe = com.convert_to_r_dataframe(between_sub)
    #r_dataframe = fixRLevels( r_dataframe, "Individual.ID", classcounts["Individual.ID"].tolist() )
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="IBS0",
                                   y="Kinship",
                                   colour="factor(Comparison)") +
                ggplot2.geom_point() +
                ggplot2.ggtitle("Something") +
                ggplot2.scale_y_continuous("Kinship Coefficient",
                                           limits=robjects.FloatVector((0,.2))) +
                ggplot2.scale_x_continuous("Pr(IBS=0)",
                                           limits=robjects.FloatVector((0,.06))) +
                #ggplot2.ylim(robjects.IntVector((0,.3))) +
                #ggplot2.theme(**{'legend.position':"none"}) +
                #ggplot2.facet_grid(robjects.Formula('. ~ GeoRegion1')) + #, scale="free") +
                ggplot2.theme(**mytheme) )
                #ggplot2.coord_flip() +
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                #ggplot2.theme(**{'axis.text.y': ggplot2.element_blank()}) +

    figname = "results/figures/king/%s_betweenrel.png" % (basename)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotWithinBetweenRelatedness

######################################################################
if __name__ == "__main__":

    #optlist, args = getopt.getopt( sys.argv[1:], "bk")
    #optlist = dict(optlist)
    #bedfile = optlist.get("-b",bedfile)

    # Change running directory
    os.chdir("..")

    vcffile = "./rawdata/mevariome/variome.vcf.gz"
    basename = "variome"
    #vcffile = "./rawdata/mevariome/main/variome.clean.vcf.gz"
    #vcffile = "./rawdata/merge1kg/main/me1000G.clean.vcf.gz"

    
    unrelfile, famrelfile, interrelfile = king_workflow( vcffile, force=False )

    #unrelated = [x.strip() for x in open(unrelfile).readlines() if len(x) > 0]
    unrelated = read_csv( unrelfile, sep="\t", header=None, names=["FID","Individual.ID"] )

    #>0.354, [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884] 
    # duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree
    #annotcol = "Continent2"
    interfam = read_csv( interrelfile, sep="\t" )
    fambins = [-1, 0.0442, 0.0884, 0.177, 0.354, 1]
    famlabels = ["4th & Greater","3rd degree","2nd degree","1st degree","duplicate" ]
    # duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree
    interfam["Fclass"] = cut(interfam.IBS0, bins=fambins, retbins=True, labels=famlabels)[0]

    # remove related samples
    interfam_sub = interfam[(interfam.ID1.isin(unrelated["Individual.ID"])) &
                              (interfam.ID2.isin(unrelated["Individual.ID"]))]

    idset1 = DataFrame({'ID1':interfam_sub["ID1"].unique().tolist(),
                        'Continent2':"Unknown", 'GeographicRegions':"Unknown"})
    idset1 = addSampleAnnotation( idset1, "ID1", True )
    idset1.columns = ["Cont1","GeoRegion1","ID1"]
    interfam_sub = merge(interfam_sub, idset1, on="ID1" )

    idset2 = DataFrame({'ID2':interfam_sub["ID2"].unique().tolist(), 
                        'Continent2':"Unknown", 'GeographicRegions':"Unknown"})
    idset2 = addSampleAnnotation( idset2, "ID2", True )
    idset2.columns = ["Cont2","GeoRegion2","ID2"]
    interfam_sub = merge(interfam_sub, idset2, on="ID2" )
    interfam_sub["Comparison"] = ["%s : %s" % tuple( sorted([x,y]) )
                              for x,y in interfam_sub[["GeoRegion1","GeoRegion2"]].values]

    plotWithinBetweenRelatedness( interfam_sub, basename )
    sys.exit(1)
    #allids = (interfam[(interfam.IBS0 > .0442)].ID1.tolist()+ 
              #interfam[(interfam.IBS0 > .0442)].ID2.tolist() )
    #thirddeg = dfTable(allids)
    #thirddegremove = thirddeg[thirddeg.Count > 300]
    #thirddeg = addSampleAnnotation( thirddeg, mergecol="Variable" )
    #thirdnew = third[~third.Variable.isin(seconddegfilt.Variable.tolist())]
    #thirdfilt = thirddeg[(thirddeg.Count >1) & (thirddeg.Count < 500 )]

    interfam_unrel = interfam[(interfam.ID1.isin(unrelated["Individual.ID"])) &
                              (interfam.ID2.isin(unrelated["Individual.ID"]))]

    fambins = [-1, 0.0442, 0.0884, 0.177, 0.354, 1]
    famlabels = ["4th & Greater","3rd degree","2nd degree","1st degree","duplicate" ]

    # duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree
    interfam_unrel["Fclass"] = cut(interfam_unrel.IBS0, bins=fambins, retbins=True, labels=famlabels)[0]

    id1annot = addSampleAnnotation( interfam_unrel, mergecol="ID1" )
    id2annot = addSampleAnnotation( interfam_unrel, mergecol="ID2" )

    interfam_unrel["continent1"] = id1annot["continent"]
    interfam_unrel["continent2"] = id2annot["continent"]

    classcounts = []
    for idx, data in interfam_unrel.groupby(["ID1","continent1","continent2"]) : #"Fclass",
        sid, continent1, continent2 = idx
        tmptbl = dfProportion(data.Fclass.tolist()) #data.ID1.tolist()+data.ID2.tolist()
        #tmptbl["Fclass"] = fclass
        tmptbl["Individual.ID"] = sid
        tmptbl["continent1"] = continent1
        tmptbl["continent2"] = continent2
        classcounts.append( tmptbl )
        
    classcounts = concat( classcounts ).reset_index()
    #classcounts = addSampleAnnotation( classcounts, mergecol="Variable" )
    print classcounts.head()
    classcounts.sort( ["Variable","Count"], inplace=True )

    #sub = classcounts[classcounts.continent1 == "Middle East"]
    #invcounts = sub.copy()
    #invcounts["continent1"] = sub["continent2"]
    #invcounts["continent2"] = sub["continent1"]
    #sub = concat( [sub, invcounts] ).reset_index(drop=True)
    #alldata = dfTable(interfam.ID1.tolist()+interfam.ID2.tolist())

    r_dataframe = com.convert_to_r_dataframe(classcounts)
    r_dataframe = fixRLevels( r_dataframe, "Individual.ID", classcounts["Individual.ID"].tolist() )
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(Individual.ID)",
                                   y="Count",
                                   fill="factor(Variable)") +
                ggplot2.geom_bar(stat="identity") +
                ggplot2.ggtitle("Count of Interfamily 3nd Degree relationship") +
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.theme(**{'axis.text.y': ggplot2.element_blank()}) +
                ggplot2.facet_grid(robjects.Formula('continent1 ~ continent2'), scale="free") +
                #ggplot2.coord_flip() +
                ggplot2.theme(**mytheme) )
    p.plot()
#stat="identity"
                #ggplot2.ggtitle("Variant counts within\nRuns of Homozygosity") +

    r_dataframe = com.convert_to_r_dataframe(classcounts)
    r_dataframe = fixRLevels( r_dataframe, "Individual.ID", classcounts["Individual.ID"].tolist() )
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(Variable)",
                                   y="Count",
                                   fill="factor(Variable)") +
                ggplot2.geom_boxplot() +
                #ggplot2.ggtitle("Count of Interfamily 3nd Degree relationship") +
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.theme(**{'axis.text.y': ggplot2.element_blank()}) +
                ggplot2.facet_grid(robjects.Formula('continent1 ~ continent2'), scale="free") +
                #ggplot2.coord_flip() +
                ggplot2.theme(**mytheme) )
    p.plot()

    figname = "results/figures/king/%s_3nddegree_bar.png" % (basename)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()


# END MAIN



