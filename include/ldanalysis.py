#!/usr/bin/env python

#-------------------------------------------------------------------#
#                       Import Packages                             #
#-------------------------------------------------------------------#
import os, sys
import getopt
import csv
import re
import subprocess
from random import randint

from housekeeping import *
from localglobals import *
import popgencommands as popc

def getLevels( annotcol, targetvcf, clustannotfile ):
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
    elif annotcol == "GeographicRegions3" :
        levels = ["Africa","Europe","East.Asia","NWA","NEA","AP","SD","TP","CA","Anglo-American"]
    else :
        print "Error: unknown levels for column -",annotcol
        sys.exit(1)

    targetpats = patientInfo.getPats( targetvcf )

    filepath,fbase,suffix = hk.getBasename(targetvcf)
    filepats = patientInfo.currentFilePats( targetvcf )

    #clustannotfile = os.path.join(filepath[:filepath.find("main")+4],"clust",fbase+".annot")
    if annotcol =="GeographicRegions3":
        sampleannot = read_csv( clustannotfile, sep="\t" )
        sampleannot.index = sampleannot["Individual.ID"]
        sampleannot = sampleannot.ix[filepats,:]

        assert os.path.exists(clustannotfile), ("Error: clust annot file not found:"+clustannotfile)

        mapping_regions = {"NWA":"Middle East", "NEA":"Middle East", "AP":"Middle East",
                           "SD":"Middle East", "TP":"Middle East", "CA":"Middle East",
                           "LWK":"Africa", "YRI":"Africa", "IBS":"Europe", "CEU":"Europe",
                           "TSI":"Europe", "FIN":"Europe", "GBR":"Europe",
                           "CHB":"East Asia", "CHS":"East Asia", "JPT":"East Asia",
                           "Anglo-American":"Anglo-American"
                          }

        sampleannot["Continent3"] = [mapping_regions[x] if mapping_regions.has_key(x)
                                      else "Unknown"
                                      for x in sampleannot.GeographicRegions3]

        mapping_regions = { "NWA":"NWA", "NEA":"NEA", "AP":"AP",
                            "SD":"SD", "TP":"TP", "CA":"CA",
                            "LWK":"Africa", "YRI":"Africa", "IBS":"Europe", "CEU":"Europe",
                            "TSI":"Europe", "FIN":"Europe", "GBR":"Europe",
                            "CHB":"East.Asia", "CHS":"East.Asia", "JPT":"East.Asia",
                            "Anglo-American":"Anglo-American"
                          }
        sampleannot["GeographicRegions3"] = [mapping_regions[x] if mapping_regions.has_key(x)
                                      else "Unknown"
                                      for x in sampleannot.GeographicRegions3]
        sampleannot = sampleannot[(sampleannot[annotcol].notnull()) &
                                    ~(sampleannot[annotcol].isin(["Unknown","Unknown1","Unknown2"]))]

    else :
        sampleannot = sampleAnnotation(targetpats)
        sampleannot = sampleannot[(sampleannot[annotcol].notnull()) &
                                    (sampleannot[annotcol] != "Unknown")]
        sampleannot[annotcol] = [x.strip().replace(" ",".") for x in sampleannot[annotcol]]

    print "Before levels:",levels
    finallevels = [x for x in levels if x in sampleannot[annotcol].unique()]
    sampleannot = sampleannot[sampleannot[annotcol].isin(finallevels)]

    #keepfiles = sampleannot["Individual.ID"].tolist()
    return finallevels, sampleannot #, keepfiles
# END getLevels

################################################################################
def seriousClean( vcffile, keepfile=None, rerun=False ):
    print "Running seriousClean"
    filepath, basename, suffix = hk.getBasename(vcffile)

    cleanped = "%s/%s"%(filepath,basename)
    #bedfile = cleanped+".bed"
    bedfile = cleanped+".filt.bed"
    if os.path.exists(bedfile) and not rerun :
        return bedfile

    command = ["vcftools","--gzvcf", vcffile,
               "--remove-filtered-all",
               "--remove-indels",
               "--hwe",".001",
               "--maf",".1",
               "--min-alleles","2",
               "--max-alleles","2",
               "--plink-tped","--out",cleanped]

    #if keepfile is not None :
        #assert os.path.exists(keepfile)
        #command = command+["--keep",keepfile]

    print " ".join(command)
    tped = cleanped+".tped"
    if not os.path.exists( tped ) or rerun:
        out = subprocess.check_output( command )

    bedfile = popc.run_plink_convert(tped, force=rerun)
    #bedfile = popc.run_plink_filter(tped, force=rerun)
    return bedfile
# END seriousClean

def splitRange( x, firstname="start", secondname="end" ) :
    return Series({firstname:int(x[1:x.find(",")]),
                   secondname:int(x[x.find(",")+2:len(x)-1])})
# ENd splitRange

################################################################################
def subsectionVCFbyPop( vcffile, sampleannot, targetdir, basename="test",
                       force=False, annotcol="Continent3" ) :
    allpopfiles = {}
    for group, data in sampleannot.groupby(annotcol) :
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
        command = [ "vcftools","--gzvcf",vcffile,
                    "--min-alleles","2","--max-alleles","2",
                    "--maf","0.0001","--max-maf","1.0",
                    "--recode","--keep",allpopfiles[pop],"--out",newvcffile]

        out = subprocess.check_output(command)

        assert os.path.exists( newvcffile+".recode.vcf" ), out
        subprocess.call(['mv', newvcffile+".recode.vcf", newvcffile+".vcf"])
        subprocess.call(["bgzip",newvcffile+".vcf"])
        subprocess.call(["tabix","-p","vcf",finalvcf])

    return popvcfs
# subsectionVCFbyPop

def runPlinkLD( bedfile, force=False ) :
    filepath, basename, suffix = hk.getBasename(bedfile)

    ldfile = filepath+"/"+basename+".ld"
    if not os.path.exists( ldfile ) or force :
        print "Bedfile:",bedfile
        command = ["plink","--noweb","--bfile",bedfile[:-4],"--r2","--ld-window-r2",
                   "0","--ld-window","99999","--ld-window-kb","70",
                   "--out",os.path.join(filepath,basename)]

        print "cmd:"," ".join(command)
        out = subprocess.check_output( command )
        if out.find("ERROR")  > 0: print out

    return ldfile
# END runPlinkLD

#LD <- read.table("your_plinkresult.ld", header = T) ##read your .ld data
#LD$distancekb <- with (LD, LD$BP_B-LD$BP_A)/1000 ## the distance between snp1 and snp2 in kb
#LD$grp <- cut(LD$distancekb, 0:70) ## bin 70kb
#r2means <- with (LD, tapply(LD$R2, LD$grp, FUN = mean)) ##r2 mean every 1kb
def calculateLDmeans( ldfile, classname="All" ) :
    print "Running calculateLDmeans"

    lddata = read_csv( ldfile, delim_whitespace=True )
    lddata["distancekb"] = (lddata.BP_B - lddata.BP_A) / 1000

    lddata["grp"] = cut( lddata.distancekb, range(0,71) )
    #print lddata.head()

    r2means = lddata.groupby("grp")["R2"].mean().reset_index()
    r2means["dist"] = [int(x[1:x.find(",")])*1000 for x in r2means.grp]
    print r2means.head()
    r2means["Class"] = classname
    return r2means
# END calculateLDmeans

def populationLdMeans( targetvcf, filepath, basename, sampleannot, annotcol, force=False ) :
    r2meansfile = os.path.join(filepath,basename+"_r2means.txt")
    if os.path.exists( r2meansfile ) and not force : 
        allr2 = read_csv( r2meansfile, sep="\t" )
        return allr2

    print sampleannot.head()

    allr2 = []
    popvcfs = subsectionVCFbyPop( targetvcf, sampleannot, filepath, 
                                 filename, force=False, annotcol=annotcol )
    for pop in popvcfs : 
        bedfile = seriousClean( popvcfs[pop], rerun=force )
        ldfile = runPlinkLD( bedfile, force=force )
        r2res = calculateLDmeans( ldfile, pop )
        allr2.append(r2res)

    allr2 = concat( allr2 ).reset_index()
    allr2.to_csv( r2meansfile, sep="\t", index=False )
    return allr2
# END populationLdMeans

def plotLDdecay( r2means, targetdir, dataset, annotcol, levels ) :

    print levels
    custpalette = makeLargePalette2(len(levels))
    print custpalette

    r_dataframe = com.convert_to_r_dataframe(r2means)
    #r_lm = com.convert_to_r_dataframe(lmdata)
    r_dataframe = fixRLevels( r_dataframe, "Class", levels )
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="dist", y="R2", colour="Class") +
                #ggplot2.geom_point() +
                ggplot2.geom_line(size=2) +
                #ggplot2.stat_smooth(method="loess", se=False) +
                ggplot2.scale_y_continuous("Average Correlation", 
                                           expand=robjects.IntVector((0,0))) +
                ggplot2.scale_x_continuous("Distance between SNPs (bp)", 
                                           expand=robjects.IntVector((0,0))) +
                ggplot2.scale_colour_manual("Region",values=custpalette) +
                #ggplot2.scale_colour_brewer("Region",palette="Set1") +
                ggplot2.theme(**pointtheme) )
                #ggplot2.scale_x_continuous("Percent of Genome in RoH") +
                #ggplot2.scale_colour_manual(name="Variant Location",
                #ggplot2.facet_grid(robjects.Formula('vtype ~ Continent'), scale="free_y") +
                                            #values=robjects.StrVector(["red","blue"]),
                                            #breaks=robjects.StrVector(["benignroh","benignnonroh"]),
                                            #labels=robjects.StrVector(["Within RoH","Outside RoH"]) ) +
                                   #group="factor(rohclass)",
                                   #colour="factor(rohclass)") +
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 90)}) +
                #ggplot2.geom_abline(ggplot2.aes_string(intercept="intercept",slope="slope"),
                                    #color="black", data=r_lm) +
                #ggplot2.geom_text(ggplot2.aes_string(label="Slope",x="xpos",y="ypos"),
                                  #size=4, hjust=0, data=r_lm) +
                #ggplot2.ggtitle("Benign Variant counts within\nRuns of Homozygosity") +


    figname = os.path.join(targetdir,"%s_%s_ld.pdf"%(dataset,annotcol))
    print "Writing file:",figname
    grdevices.pdf(figname, width=7, height=5)
    p.plot()
    grdevices.dev_off()
# END plotLDdecay

################################################################################
if __name__ == "__main__" :

    os.chdir("..")

    optlist, args = getopt.getopt( sys.argv[1:], "r:ot")
    optlist = dict(optlist)

    dataset = optlist.get("-r","test2")
    if dataset == "test2" :
        vcffile = "./rawdata/test2/main/test2.clean.vcf.gz"
        clustfile = "./rawdata/mevariome/main/clust/variome.clean.annot"
    elif dataset == "daily" :
        vcffile = "./rawdata/daily/main/daily.clean.vcf.gz"
        clustfile = "./rawdata/daily/main/clust/daily.clean.annot"
    elif dataset == "merge1kg" :
        vcffile = "./rawdata/merge1kg/main/me1000G.clean.vcf.gz"
        clustfile = "./rawdata/merge1kg/main/clust/me1000G.clean.annot"
    elif dataset == "onekg" :
        vcffile = "./rawdata/onekg/main/onekg.clean.vcf.gz"
        clustfile = "./rawdata/onekg/main/clust/onekg.clean.annot"
    elif dataset == "mevariome" :
        vcffile = "./rawdata/mevariome/main/variome.clean.vcf.gz"
        clustfile = "./rawdata/mevariome/main/clust/variome.clean.annot"
    else :
        print "Error: dataset not recognized -",dataset
        sys.exit(1)

    annotcol="GeographicRegions3"
    
    #filepath, filename = os.path.split(vcffile)
    #targetdir = os.path.join(filepath,"ld")
    #hk.makeDir(targetdir)

    #annotcol = "GeographicRegions3"
    targetvcf = hk.copyToSubDir( vcffile, "ld" )
    filepath, filename, suffix = hk.getBasename(targetvcf)

    basename = filename[:filename.find(".")] if filename.find(".") > 0 else filename

    levels, sampleannot = getLevels( annotcol, targetvcf, clustfile )

    allr2 = populationLdMeans( targetvcf, filepath, basename, sampleannot, annotcol, True )

    plotLDdecay( allr2, filepath, basename, annotcol, levels )

# END MAIN
################################################################################

    #levels, sampleannot = getLevels( annotcol, targetvcf )
    #print levels
    #print hk.dfTable(sampleannot[annotcol])


