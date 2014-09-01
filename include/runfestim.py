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
from popgencommands import *
from localglobals import *
import patientInfo

################################################################################
# vcftools_splitbypop
################################################################################
def vcftools_splitbypop( vcffile, patientannot, force=forceFlag ):
    write2log(" - Running "+whoami(), True)

    filepath, basename, suffix = getBasename(vcffile)

    targetdir = os.path.join(filepath,"pops")
    makeDir(targetdir)
    patientannot.index = patientannot["Individual.ID"].tolist()
    samples = patientInfo.getPats( vcffile ) 
    assert len(samples) > 0
    samples = [ x for x in samples if x in patientannot["Individual.ID"].tolist() ]
    annot = patientannot.loc[samples,["Individual.ID","Continent2"]]
    #print annot.head(3)

    patfiles = {}
    for region, data in annot.groupby("Continent2") : 
        if len(data) < 10 : continue
        patfiles[region.replace(" ","_")] = "%s/%s_pats.txt" % (targetdir, region.replace(" ","_"))
        patsfile = os.path.join(targetdir,region.replace(" ","_")+"_pats.txt")
        data["Individual.ID"].to_csv( patsfile, index=False, header=False)

    newfiles = {}
    for pop, popfile in patfiles.iteritems() :
        plinkfile = "%s/%s.%s" % (targetdir, basename, pop.replace(" ","_"))
        newfiles[pop] = plinkfile+".tped"
        if os.path.exists( plinkfile+".tped" ) and not force : print "Already exists"; continue
        print plinkfile
        command = [ "vcftools","--gzvcf",vcffile,
                    "--keep",popfile,
                    "--remove-filtered-all",
                    "--maf",".005",
                    "--remove-indels","--min-alleles","2",
                    "--max-alleles","2",
                    "--plink-tped","--out",plinkfile]
                #% ( vcffile, popfile, plinkfile )
                #" --max-maf .99" \
                #"--plink-tped" \
                #" --thin 1000" \
                #" --FILTER-summary %s" \

        print " ".join(command)
        out = subprocess.check_output( command )
        print out
    return newfiles
# END vcftools_splitbypop

######################################################################
# vcftools_makeplink
######################################################################
def vcftools_makeplink( vcffile, force=forceFlag ):
    write2log(" - Running "+whoami(), True)
    targetdir, filename, suffix = getBasename(vcffile)
    plinkfile = "%s/%s.plink" % (targetdir, filename)
    print plinkfile

    if os.path.exists( plinkfile+".tped" ) and not force :
        print " - File"+plinkfile+".tped already exists. Skipping command!"
        return plinkfile+".tped"

    command = [ "vcftools","--gzvcf",vcffile,
                "--remove-filtered-all",
                "--maf",".005",
                "--remove-indels","--min-alleles","2",
                "--max-alleles","2",
                "--plink-tped","--out",plinkfile ]
        #% ( vcffile, plinkfile )
        #"--plink-tped" \
        #" --thin 1000" \
        #" --FILTER-summary %s" \

    print " ".join(command)
    out = subprocess.check_output( command )
    #print out
    return plinkfile+".tped"
# END vcftools_makeplink

######################################################################
# excludeSnps
# Goal - remove variants from dataset that fall 
# within the same centimorgan
################################################################################
################################################################################
def excludeSnps_old( pedfile, force=False ) :
    write2log(" - Running "+whoami(), True)
    filepath, basename, suffix = getBasename(pedfile)
    frqfile = filepath+"/"+basename+".frq"
    mapfile = filepath+"/"+basename+".map"
    excludefile = filepath+"/"+basename+".exclude.txt"
    pfile = plink_fileinput( pedfile )
    print "Fixing mapfile:",mapfile

    newfrq = filepath+"/"+basename+".exclude.frq"
    #newmap = filepath+"/"+basename+".exclude.map"
    thinped = filepath+"/"+basename+".exclude"

    if os.path.exists(thinped+".ped") and os.path.exists(newfrq) and not force : 
        return newfrq, thinped

    print "Reading:",frqfile
    # VID - variant id
    freqdata = read_csv(frqfile,sep=" ",header=None, 
                        names=["chrom","vid","ref","alt","freq","count"], 
                        dtype={"vid":str})
    freqdata["pos"] = [int(x.split(".")[1]) for x in freqdata.vid]
    zerofreq = freqdata.vid[freqdata.freq.isin([0,1])].tolist()
    print "Zerofreq:",zerofreq

    #evars = {}
    #for row in csv.reader(open(mapfile),delimiter=" ") :
    # http://www.nature.com/nature/journal/v409/n6822/fig_tab/409951a0_T1.html
    # calculate Centimorgans!!!!! 1.30 cM/Mb
    print "Reading mapfile:",mapfile
    mapdata = read_csv( mapfile, header=None, delim_whitespace=True, 
                        names=["chrom","vid","centimorgans","pos"],
                        dtype={"vid":str})

    mapdata["centimorgans"] = ["%.1f" % (float(pos) / 1000000 * 1.3) 
                               for pos in mapdata.pos]
    alldata = merge(mapdata, freqdata, on="vid" )
    print "mapdata:",mapdata.shape
    print "Alldata:",alldata.shape
    for i in range(1,1000) :
        if freqdata.ix[i,"vid"]  != mapdata.ix[i,"vid"] :
            print "Error",i,freqdata.ix[i,"vid"], mapdata.ix[i,"vid"]

    maxfreqbycm = alldata[alldata.freq > 0.0].groupby("centimorgans")["freq"].idxmax()
    newfrqdata = freqdata.ix[maxfreqbycm,:].sort(["chrom","pos"])
    newfrqdata[["chrom","vid","ref","alt","freq","count"]].to_csv(
        newfrq,sep=" ",index=False,header=False)

    # Write out exclude variants
    newfrqdata["vid"].to_csv(excludefile, index=False,header=False)
    filetype = "--file" if pedfile[-3:] == "ped" else "--bfile"
    command = ["plink","--noweb",filetype,pedfile[:-4],"--allow-no-sex",
               "--extract",excludefile,"--recode","--out",thinped]
    print " ".join(command)
    try :
        out = subprocess.check_output( command )
    except subprocess.CalledProcessError as e :
        print "Error:",e
        sys.exit(1)

    print "Thin ped!:",thinped
    assert os.path.exists(thinped+".ped")
    return newfrq, thinped
# END excludeSnps_old


################################################################################
# plink_recode12
################################################################################
def plink_recode12( pedfile, force=forceFlag ) :
    write2log(" - Running "+whoami(), True)
    targetdir,filename,suffix = getBasename(pedfile)

    basename = "%s/%s" % (targetdir, filename)
    recodefile = basename+".recode12"
    # --1
    if not os.path.exists( recodefile+".ped" ) or force:
        command = "plink --noweb --tfile %s" \
                " --recode12 --out %s" \
                % ( basename, recodefile )
                #" --geno 0.1" \
                #" --maf 0.01" \
        print command
        out = runCommand( command )
        print out
    return recodefile+".ped"
# END plink_recode12

################################################################################
# plink_runhaplotyping
################################################################################
def plink_runhaplotyping( pedfile ) :
    write2log(" - Running "+whoami(), True)
    # Check if outfile already exists
    targetdir, filename, suffix = getBasename(pedfile)
    basename = "%s/%s" % (targetdir, filename)

    hapfile = basename + ".hap"
    print "#" *70
    print "Running haplotyping!!", hapfile
    print "t",targetdir,"f",filename,"s",suffix,"b",basename,hapfile
    # --1
    if not os.path.exists( hapfile ) :
        command = "plink --noweb --file %s --blocks --out %s" % (basename, basename)
        print command
        out = runCommand( command )
        print out

        blockfile = basename +".blocks"
        command = "plink --noweb --file %s --hap %s --hap-freq --out %s" % (basename, blockfile, basename)
        print command
        out = runCommand( command )
        print out

    print "#" *70
    return hapfile
# END plink_runhaplotyping

################################################################################
# plink_makefrqfile
################################################################################
def plink_makefrqfile( pedfile, force=forceFlag ) :
    write2log(" - Running "+whoami(), True)
    filepath,basename,suffix = getBasename(pedfile)
    # Make frqfile
    freqfile = filepath+"/"+basename+".frq"
    tmpfile = filepath+"/"+basename+".1.frq"
    pfile = plink_fileinput( pedfile )
    #print "Checking frq file:",freqfile, os.path.exists( freqfile )
    # --1
    if not os.path.exists( freqfile ) or force :
        command = "plink --noweb %s" \
                " --freq --maf 0.005 --out %s" \
                % ( pfile, filepath+"/"+basename+".1" )
        # --filter-controls 
        #--recode12 
        print command
        out = runCommand( command )
        print out
        command = "awk '$1 !~ /CHR/ {print $1,$2,$3,$4,$5,$6}' %s > %s" % (tmpfile, freqfile )
        print command
        out = runCommand( command )
        print out
    return freqfile
# END plink_makefrqfile

################################################################################
# plink_outliers
################################################################################
def plink_outliers( pedfile ) :
    write2log(" - Running "+whoami(), True)
    # Run outlier analysis
    filepath, basename, suffix = getBasename( pedfile)
    pfile = plink_fileinput( pedfile )
    clusterfile = basename+".cluster0"
    # --1
    if not os.path.exists( clusterfile ) or forceFlag is True:
        command = "plink --noweb %s" \
                " --cluster --neighbour 1 5 --out %s" \
                % ( pfile, filepath+"/"+basename )
        print command
        out = runCommand( command )
        print out
    return clusterfile
# END plink_outliers

################################################################################
# run_festim
################################################################################
def run_festim( pedfile, force=False, iterations="2000" ) :
    write2log(" - Running "+whoami(), True)
    filepath, basename, suffix = getBasename( pedfile )
    festimfile = "%s/%s_festim_snp.out" % (filepath,basename)
    ibcfile = "%s/%s_inbreedings.out" % (filepath,basename)
    if os.path.exists(ibcfile) and not force : return ibcfile
    print "FEstimfile:",festimfile,"IBCfile:",ibcfile
    #command = "FEstim --mplink %s/%s" \
            #" --SE --iterMC 2000" \
            #" --ped %s > %s" \
            #% ( filepath, basename, pedfile, festimfile )

    command = ["FEstim", "--mplink", "%s/%s" % (filepath,basename),"--SE", "--iterMC",iterations, "--ped", pedfile]
            #% ( filepath, basename, pedfile, festimfile )
    #print "#"*70
    print " ".join(command)
    #write2log( command, True )
    #out = runCommand( command )
    try:
        stdout = subprocess.check_output(command, stderr=subprocess.STDOUT)
        print stdout
    except subprocess.CalledProcessError as e:
        print "Caught error:",e

    #print out
    runCommand("mv lodscore.out %s/%s_lodscore.out" % (filepath,basename) )
    runCommand("mv locusIBD.out %s/%s_locusIBD.out" % (filepath,basename) )
    runCommand("mv family_lodscore.out %s/%s_family_lodscore.out" % (filepath,basename) )
    runCommand("mv inbreedings.out %s/%s_inbreedings.out" % (filepath,basename) )
    return ibcfile
# END run_festim

# findMid
# Find middle occurence
################################################################################
def findMid(s, ch):
    allindexes= [i for i, ltr in enumerate(s) if ltr == ch]
    return allindexes[len(allindexes)/2]

################################################################################
# plotIBC
################################################################################
def plotIBC( ibcfile, patientdata, factor="continent" ):
    print "Running plotIBC:",ibcfile

    filepath, basename, suffix = getBasename( ibcfile )
    #patientannotation = "./resources/annotation/patientannotation.ped"
    #patientdata = read_csv( patientannotation, sep="\t" )
    ibcdata = read_csv( ibcfile, sep="\t" )

    #ibcdata["Family"] = [x[:findMid(x,"_")] for x in ibcdata["Individual"].tolist()]
    #ibcdata["Individual.ID"] = [x[findMid(x,"_")+1:] for x in ibcdata["Individual"].tolist()]
    alldata = merge(ibcdata[["Individual.ID","F","StdE","A"]],
                    patientdata[["Individual.ID",factor]],how="left",on="Individual.ID")
    filtdata = alldata[alldata[factor].notnull()]
    filtdata.F = filtdata.F.map(float)
    print filtdata.head(10)

    r_dataframe = com.convert_to_r_dataframe(filtdata)
    #print summary(r_dataframe)
    #print r_dataframe.rx2("F")
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor("+factor+")",y="F" ) +
                ggplot2.geom_boxplot(notch=True) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.ggtitle("IBC by "+factor.capitalize()) +
                ggplot2.scale_y_continuous("Inbreeding Coefficient (festim)") +
                ggplot2.theme(**mytheme) )
                #ggplot2.aes_string(x = "factor(continent)",y="F" ) +

                #ggplot2.scale_colour_manual(values = robjects.StrVector(("blue", "red", "grey"))) + \
                #ggplot2.scale_x_continuous("Sum of all variant sites in a gene", limits=robjects.IntVector((0,600))) + \
                #ggplot2.stat_smooth(method="lm", se=False)
    if basename is None : basename = "test"
    if basename.find(" ")>=0 : basename = basename.replace(" ","_")
    figname = "results/figures/festim/%s_%s_festim_ibc.png" % (basename,factor)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotIBC

################################################################################
# plotPlinkIBC
# FID       Family ID
# IID       Individual ID
# O(HOM)    Observed number of homozygotes
# E(HOM)    Expected number of homozygotes
# N(NM)     Number of non-missing genotypes
# F         F inbreeding coefficient estimate
################################################################################
def plotPlinkIBC( ibcfile, patientdata, factor="continent" ):
    print "Running plotPlinkIBC:",ibcfile

    filepath, basename, suffix = getBasename( ibcfile )

    #patientannotation = "./resources/annotation/patientannotation.ped"
    #patientdata = read_csv( patientannotation, sep="\t" )
    ibcdata = read_csv( ibcfile, sep="\t" )
    #print ibcdata.head(10)
    #print ibcdata.columns
    #ibcdata.columns = ["blank","Family","Individual.ID","Homozygotes","Expected","Nonmissing","F"]
    alldata = merge(ibcdata[["Individual.ID","F","Homozygotes","Expected"]],patientdata[["Individual.ID",factor]],how="left",on="Individual.ID")
    alldata = addSampleAnnotation( ibcdata, "Individual.ID" )
    #print alldata.head(3)
    #print alldata[alldata[factor].isnull()].shape
    filtdata = alldata[alldata[factor].notnull()]

    r_dataframe = com.convert_to_r_dataframe(filtdata)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor("+factor+")",y="F" ) +
                ggplot2.geom_boxplot(notch=True) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.ggtitle("IBC by "+factor.capitalize()) +
                ggplot2.scale_y_continuous("Inbreeding Coefficient (plink)") +
                ggplot2.theme(**mytheme) )
                #ggplot2.scale_colour_manual(values = robjects.StrVector(("blue", "red", "grey"))) + \
                #ggplot2.scale_x_continuous("Sum of all variant sites in a gene", limits=robjects.IntVector((0,600))) + \
                #ggplot2.stat_smooth(method="lm", se=False)
                #ggplot2.geom_abline(intercept=rstats.coef(model_x)[1], slope=rstats.coef(model_x)[2])
    if basename is None : basename = "test"
    if basename.find(" ")>=0 : basename = basename.replace(" ","_")
    figname = "results/figures/festim/%s_%s_plink_ibc.png" % (basename, factor)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotPlinkIBC

################################################################################
# ibcCorrelation
################################################################################
def ibcCorrelation( plinkibc, festimibc, patientdata, factor="continent" ):
    print "Running ibcCorrelation:",plinkibc

    filepath, basename, suffix = getBasename( plinkibc )

    #patientannotation = "./resources/annotation/patientannotation.ped"
    #patientdata = read_csv( patientannotation, sep="\t" )

    plinkibcdata = read_csv( plinkibc, sep="\t" )
    #print plinkibcdata.head(3)
    #plinkibcdata.columns = ["blank","Family","Individual.ID","Homozygotes","Expected","Nonmissing","F"]
    plinkall = merge(plinkibcdata[["Individual.ID","F"]],patientdata[["Individual.ID",factor]],how="left",on="Individual.ID")
    plinkall = plinkall[plinkall[factor].notnull()]
    plinkall["Type"] = "Plink"

    feibcdata = read_csv( festimibc, sep="\t" )
    #feibcdata["Family"] = [x[:x.rfind("_")] for x in feibcdata["Individual"].tolist()]
    #feibcdata["Individual.ID"] = [x[x.rfind("_")+1:] for x in feibcdata["Individual"].tolist()]
    festimall = merge(feibcdata[["Individual.ID","F"]],patientdata[["Individual.ID",factor]],how="left",on="Individual.ID")
    festimall = festimall[festimall[factor].notnull()]
    festimall["Type"] = "FEstim"
    #print alldata[alldata[factor].isnull()].shape

    alldata = concat( [plinkall,festimall] ).reset_index()
    alldata.F = alldata.F.map(float)

    r_dataframe = com.convert_to_r_dataframe(alldata)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor("+factor+")",y="F" ) +
                ggplot2.geom_boxplot(ggplot2.aes_string(fill="factor(Type)"),notch=True) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.ggtitle("IBC by Regions") +
                ggplot2.scale_y_continuous("Inbreeding Coefficient") +
                ggplot2.theme(**mytheme) )
                #ggplot2.scale_colour_manual(values = robjects.StrVector(("blue", "red", "grey"))) + \
                #ggplot2.scale_x_continuous("Sum of all variant sites in a gene", limits=robjects.IntVector((0,600))) + \
                #ggplot2.stat_smooth(method="lm", se=False)
                #ggplot2.geom_abline(intercept=rstats.coef(model_x)[1], slope=rstats.coef(model_x)[2])
    if basename is None : basename = "test"
    if basename.find(" ")>=0 : basename = basename.replace(" ","_")
    figname = "results/figures/festim/%s_%s_ibc.png" % (basename,factor)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()

    allibcdata = merge(plinkibcdata[["Individual.ID","F"]],feibcdata[["Individual.ID","F"]], on="Individual.ID")
    print "Allibcdata"
    print allibcdata.head(30)

    alldata = merge(allibcdata,patientdata[["Individual.ID",factor]], on="Individual.ID")
    if len(alldata) == 0 :
        shortids = [x[:x.find('_')] if len(x) > 21 and x.find('_') > 0 else x for x in patientdata["Individual.ID"].tolist()]
        patientdata["Individual.ID.short"] = shortids
        alldata = merge(allibcdata,patientdata[["Individual.ID.short",factor]], right_on="Individual.ID.short",left_on="Individual.ID")

    print "Alldata for IBC correlation"
    print alldata.head(3)
    alldata.F_x = alldata.F_x.map(float)
    alldata.F_y = alldata.F_y.map(float)
    r_dataframe = com.convert_to_r_dataframe(alldata)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x = "F_x",y="F_y") + \
                ggplot2.geom_point(ggplot2.aes_string(colour="factor("+factor+")")) + \
                ggplot2.theme(**mytheme) + \
                ggplot2.ggtitle("Correlation between Inbreeding Coefficients") + \
                ggplot2.scale_x_continuous("F (plink)") + \
                ggplot2.scale_y_continuous("F (festim)") + \
                ggplot2.stat_smooth(method="lm", se=False)
                #ggplot2.scale_colour_manual(values = robjects.StrVector(("blue", "red", "grey"))) + \
                #ggplot2.scale_x_continuous("Sum of all variant sites in a gene", limits=robjects.IntVector((0,600))) + \
                #ggplot2.geom_abline(intercept=rstats.coef(model_x)[1], slope=rstats.coef(model_x)[2])
    if basename is None : basename = "test"
    if basename.find(" ")>=0 : basename = basename.replace(" ","_")
    figname = "results/figures/festim/%s_%s_ibccorr.png" % (basename,factor)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END ibcCorrelation

def wavg(group):
    d = group['ConsangPercent']
    w = group['N']
    return (d * w).sum() / w.sum()
# END wavg

def getMarriageRateAnnot( annotcol="Continent2" ):
    print "Running getMarriageRateAnnot",annotcol
    countryrates = read_csv("./resources/consangrates_norm.txt",sep="\t")
    countryrates = countryrates[~countryrates.Country.isin(["Sudan","Israel"])]
    region2countriesdf = DataFrame([ {'Region':x,'Country':y}
                                      for x in region2countries 
                                      for y in region2countries[x] ])

    regionmerge = merge(countryrates, region2countriesdf, on="Country")
    if annotcol == "ethnicity" :
        hgdp = read_csv("resources/HGDP_ethnicgroups.txt",sep="\t")
        ethnicities = merge(hgdp, countryrates, on="Country",how="left")
        s = countryrates.groupby("ethnicity",as_index=False).apply(wavg)
        s.rename(columns={'ethnicity':annotcol,0:'ConsangPercent'}, inplace=True)
    elif annotcol == "country" or annotcol == "Origin" :
        s = countryrates[['Country','ConsangPercent']]
        s.rename(columns={'Country':annotcol}, inplace=True)
    elif annotcol == "Continent2" :
        subset = regionmerge[regionmerge.Region.isin( [x.replace("."," ") 
                                                       for x in target_continents] )] 
        s = subset.groupby("Region",as_index=True).apply(wavg).reset_index()
        s.rename(columns={'Region':annotcol,0:'ConsangPercent'}, inplace=True)
    elif annotcol == "Continent" :
        s = countryrates.groupby("Continent",as_index=True).apply(wavg).reset_index()
        s.rename(columns={'Continent':annotcol,0:'ConsangPercent'}, inplace=True)
    elif annotcol == "GeographicRegions" :
        subset = regionmerge[regionmerge.Region.isin( [x.replace("."," ") 
                                                       for x in target_geographic_regions] )] 
        s = subset.groupby("Region",as_index=True).apply(wavg).reset_index()
        s.rename(columns={'Region':annotcol,0:'ConsangPercent'}, inplace=True)
    elif annotcol == "GeographicRegions2" :
        subset = regionmerge[regionmerge.Region.isin( [x.replace("."," ") 
                                                       for x in target_geographic_regions2] )] 
        s = subset.groupby("Region",as_index=True).apply(wavg).reset_index()
        s.rename(columns={'Region':annotcol,0:'ConsangPercent'}, inplace=True)
    else :
        print "Error: Unknown annotation column:", annotcol
        sys.exit(1)
    
    print s.head()
    s.sort('ConsangPercent', inplace=True)
    return s
    #countryrates = merge(hgdp[["Continent","Country"]].drop_duplicates(),consangrates,on="Country",how="left")
    #countryrates.to_csv("countryrates.tsv",sep="\t",index=False)

    #print countryrates.head(3)
    #s = countryrates.groupby("Continent",as_index=False).apply(wavg)
    #regionrates = countryrates.groupby("Continent",as_index=False).mean()
    #regionrates.index = regionrates.Continent
    #regionrates.update(DataFrame(s,columns=["ConsangPercent"]))
    #print regionrates.head(3)
    #regionrates = regionrates[regionrates["ConsangPercent"].notnull()].sort("ConsangPercent")
    #countryrates.rename(columns={'ConsangPercent':'CountryConsangPercent'}, inplace=True)
    #regionrates.rename(columns={'ConsangPercent':'ContConsangPercent'}, inplace=True)
# END getMarriageRateAnnot

################################################################################
# ibcCountryCorrelation
################################################################################
def ibcCountryCorrelation( ibcdata, basename, algorithm="plink", annotcol="Continent" ):
    print "Running ibcCountryCorrelation:",basename

    #sampleannot = read_csv( "./resources/annotation/patientannotation.ped", sep="\t" )
    #sampleannot = merge(sampleannot[["Individual.ID","continent","Origin","ethnicity"]], hgdp[["Ethnicity","Continent","Country"]], left_on="ethnicity",right_on="Ethnicity",how="left")

    #print plinkibcdata.head(3)
    #plinkibcdata.columns = ["blank","Family","Individual.ID","Homozygotes","Expected","Nonmissing","F"]
    print ibcdata.shape
    #ibcall = merge(ibcdata[["Individual.ID","F"]],sampleannot[["Individual.ID","Continent","Origin","ethnicity"]],how="left",on="Individual.ID")
    ibcall = addSampleAnnotation(ibcdata[["Individual.ID","F"]],"Individual.ID")

    print ibcall.shape
    ibcall = ibcall[(ibcall[annotcol].notnull()) & (ibcall[annotcol] != "Unknown")]
    ibcall["Type"] = "Plink"

    regionrates = getMarriageRateAnnot( annotcol )
    #print "Region"
    print regionrates.head(10)
    #print "Country"
    #countrylist = countryrates.Country.unique().tolist()
    #ibclist = ibcall.Origin.unique().tolist()
    #print [x for x in ibclist if x not in countrylist]
    withrates = merge(ibcall, regionrates, on=annotcol) #, how="left"
    #withrates = withrates[(withrates["ConsangPercent"].notnull())  &
                          #(withrates["ConsangPercent"] != "NA")]
    #withrates = merge(withrates, countryrates[["Country","CountryConsangPercent"]], left_on="Origin", right_on="Country", how="inner")

    print "Testing without internal Europeans"
    withrates = withrates[~((withrates.GeographicRegions == "Europe") & 
                          (withrates.GeographicRegions2 == "Europe"))]

    withrates = withrates[~((withrates.GeographicRegions == "Africa") & 
                          (withrates.GeographicRegions2 == "Africa"))]
    withrates = withrates[~((withrates.GeographicRegions == "East Asia") & 
                          (withrates.GeographicRegions2 == "East Asia"))]
    print "finished"
    withrates.to_csv("results/ibc/%s_%s_%s_ibcrates.txt" %(basename,algorithm,annotcol),
                     sep="\t",index=False)
    #print "Making continent ibcrates file"
    #print withrates.shape
    withrates.F = withrates.F.astype(float)
    #withrates.ContConsangPercent = withrates.ContConsangPercent.astype(float)

    colmeans = withrates.groupby(annotcol)["F"].mean().reset_index()
    #colmeans.rename(columns={0:"Colmean"},inplace=True)
    colmeans["Label"] = ["%.3f"%x for x in colmeans.F]
    print colmeans.head()
    r_dataframe = com.convert_to_r_dataframe(withrates)
    rmeans = com.convert_to_r_dataframe(colmeans)
    r_dataframe = fixRLevels( r_dataframe, annotcol, 
                             regionrates[annotcol].tolist() )
    rmeans = fixRLevels( rmeans, annotcol, 
                             regionrates[annotcol].tolist() )
    for col in r_dataframe:
        col.rclass = None
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor("+annotcol+")",y="F" ) +
                ggplot2.geom_boxplot(ggplot2.aes_string(fill="ConsangPercent"),notch=True) +
                ggplot2.geom_text(ggplot2.aes_string(label="Label", x="F"),y=0, size=.6, data=rmeans ) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45, hjust=1,vjust=1)}) +
                #ggplot2.ggtitle("IBC by %s"%annotcol) +
                ggplot2.scale_y_continuous("Inbreeding Coefficient") +
                ggplot2.scale_x_discrete("Geographic Region") +
                ggplot2.scale_fill_continuous("Percent of\nConsanguineous\nMarriages") +
                ggplot2.theme(**mytheme) )
                #ggplot2.scale_colour_manual(values = robjects.StrVector(("blue", "red", "grey"))) + \
    #if basename is None : basename = "test"
    if basename.find(" ")>=0 : basename = basename.replace(" ","_")
    figname = "results/figures/festim/%s_%s_%s_percent_ibc.png" % (basename,algorithm,annotcol)
    print "Writing file:",figname
    grdevices.png(figname, width=5, height=3, units="in",res=300)
    p.plot()
    grdevices.dev_off()

    minF = withrates.F.min()
    maxF = withrates.F.max()
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string( x="F" ) +
                ggplot2.stat_density(ggplot2.aes_string(ymax="..density..", ymin="-..density..",
                                                        fill="ConsangPercent"),
                                     geom="ribbon", position="identity") + #,notch=True 
                ggplot2.geom_text(ggplot2.aes_string(label="Label", x="F"),y=0, size=3.5, data=rmeans ) +
                #ggplot2.ggtitle("HWE Probability by AF") +
                ggplot2.facet_grid(robjects.Formula('. ~ '+annotcol), scales="free") +
                ggplot2.scale_x_continuous("Inbreeding Coefficient", 
                                           limits=robjects.FloatVector((minF,maxF))) +
                ggplot2.scale_y_discrete("Geographic Region") +
                ggplot2.scale_fill_continuous("Percent of\nConsanguineous\nMarriages") +
                ggplot2.coord_flip() +
                ggplot2.theme(**ribbontheme))

    
    #figurename = "%s/%s_hweribbon.png" % (outdir, basename)
    figname = "results/figures/festim/%s_%s_%s_percent_density_ibc.pdf" % (basename,algorithm,annotcol)
    print "Making figure:",figname
    #grdevices.png(figname, width=7, height=4, units="in",res=300)
    grdevices.pdf(figname, width=7.5, height=4)
    p.plot()
    grdevices.dev_off()

    return
# END ibcCountryCorrelation

def ibcCountryCorrelationold( ibcdata, basename, annotcol="Continent" ):
    if sum(ibcall.Origin.isin(countryrates.Country.tolist())) == 0 :
        print "No countries with known consang rates! Exiting"
        return

    withrates = merge(ibcall, countryrates, left_on="Origin", right_on="Country", how="inner")
    #withrates = withrates[withrates["ConsangPercent"].notnull()]
    withrates["CountryConsangPercent"] = withrates["CountryConsangPercent"].fillna(0).astype(int)
    withrates.to_csv("results/ibc/%s_country_ibcrates.txt" %basename,sep="\t",index=False)
    print "Making country ibcrates file"
    print withrates.head(10)
    withrates.F = withrates.F.map(float)

    r_dataframe = com.convert_to_r_dataframe(withrates)
    r_dataframe = fixRLevels( r_dataframe, "Origin", countryrates.Country.tolist() )
    for col in r_dataframe:
        col.rclass = None
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor(Origin)",y="F" ) +
                ggplot2.geom_boxplot(ggplot2.aes_string(fill="CountryConsangPercent"),notch=True) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.ggtitle("IBC by Origin") +
                ggplot2.scale_y_continuous("Inbreeding Coefficient") +
                ggplot2.scale_fill_continuous("Percent of\nConsanguineous\nMarriage") +
                ggplot2.theme(**mytheme) )
                #ggplot2.scale_colour_manual(values = robjects.StrVector(("blue", "red", "grey"))) + \
                #ggplot2.scale_x_continuous("Sum of all variant sites in a gene", limits=robjects.IntVector((0,600))) + \
    if basename is None : basename = "test"
    if basename.find(" ")>=0 : basename = basename.replace(" ","_")
    figname = "results/figures/festim/%s_%s_percent_ibc.png" % (basename,"Origin")
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()

    out = runCommand("cd rscripts/; ./ibcplots.R %s" %basename)
    print out
# END ibcCountryCorrelationold

################################################################################
# plotPlinkIBCold
################################################################################
#def plotPlinkIBCold( bedfile ) :
    #write2log( " - Running:"+whoami(), True )
    #print "Target File:",bedfile 
    #targetdir, filename, suffix = getBasename(bedfile )
    #command = "./rscripts/inbreedCoeff.R %s %s" % ( targetdir, filename )
    #print "Command",command
    #output = runCommand( command )
    #return
# plotPlinkIBCold

################################################################################
# mergeFEstimFiles
################################################################################
def mergeFEstimFiles( festimfiles, targetfile, patientdata ) :
    print "Running mergeFEstimFiles", festimfiles
    shortids = [x[:x.find('_')] if len(x) > 21 and x.find('_') > 0 else x for x in patientdata["Individual.ID"].tolist()]
    patientdata["Individual.ID.short"] = shortids

    alldata = []
    for festimibc in festimfiles :
        feibcdata = read_csv( festimibc, sep=r"\s+", dtype={'Individual':str, 'F':float} )
        #feibcdata = read_csv( festimibc, header=None, sep=r"\s+", names=["Individual","F","StdE","A","StdE2"], dtype={'Individual':str, 'F':float} )
        #feibcdata["Family"] = [x.split("_",1)[0] for x in feibcdata["Individual"].tolist()]
        #feibcdata["Individual.ID"] = [x.split("_",1)[1] for x in feibcdata["Individual"].tolist()]
        print "Filename:",festimibc
        print feibcdata.head(3)
        print [x for x in feibcdata.Individual]
        feibcdata["Family"] = [x[:x.rfind("_")] if x.find("_") > 0 else x 
                               for x in feibcdata["Individual"].tolist() ]
        feibcdata["Individual.ID.short"] = [x[x.rfind("_")+1:] for x in feibcdata["Individual"].tolist()]
        alldata.append(feibcdata)

    #"Individual","F","StdE","A","StdE.1","Family","Individual.ID
    alldata = concat( alldata )
    alldata_merge = merge( alldata, patientdata, on="Individual.ID.short", how="left" )
    alldata_filt = alldata_merge[["Family.ID","Individual.ID","F","StdE","A","StdE.1","continent","Origin","Ethnicity","ethnicity"]]
    print alldata_filt[alldata_filt["Family.ID"].isnull()].head(3)
    alldata_filt.to_csv( targetfile, sep="\t", index=False )
# END mergeFEstimFiles 

################################################################################
# mergehetfiles
################################################################################
def mergehetfiles( hetfiles, targetfile, patientdata ) :

    shortids = [x[:x.find('_')] if len(x) > 21 and x.find('_') > 0 else x for x in patientdata["Individual.ID"].tolist()]
    patientdata["Individual.ID.short"] = shortids

    alldata =[]
    for plinkibc in hetfiles :
        plinkibcdata = read_csv( plinkibc, delim_whitespace=True )#, columns=None
        plinkibcdata.columns = ["Family","Individual.ID.short","Homozygotes","Expected","Nonmissing","F"]
        print plinkibcdata.head(4)
        alldata.append(plinkibcdata)

    alldata = concat(alldata)
    alldata_merge = merge( alldata, patientdata, on="Individual.ID.short", how="left" )
    alldata_filt = alldata_merge[["Family.ID","Individual.ID","Homozygotes","Expected","Nonmissing","F","continent","Origin","Ethnicity","ethnicity"]]
    print alldata_filt[alldata_filt["Family.ID"].isnull()].head(3)
    alldata_filt.to_csv( targetfile, sep="\t", index=False )
# END mergehetfiles

################################################################################
# calculateIBC
################################################################################
def calculateIBC( targetfile, rerun=False ) :
    write2log(" - Running "+whoami(), True)
    print "FEstim Targetfile:",targetfile

    # Parse Patient data
    #patientannotation = "./resources/annotation/patientannotation.ped"
    #patientdata = read_csv( patientannotation, sep="\t" )
    targetpats = patientInfo.getPats( targetfile )
    patientdata = sampleAnnotation(targetpats)

    filepath, basename, suffix = getBasename( targetfile )
    print filepath, basename, suffix
    targetdir = "%s/ibc/" % filepath
    makeDir(targetdir)
    filepath, basename, suffix = getBasename( targetfile )
    cptargetfile = targetdir+basename+suffix+".gz"

    if not os.path.exists(cptargetfile) :
        out = runCommand( "cp %s/%s%s* %s" % (filepath,basename,suffix, targetdir) )
        print "Copy targetfile:",cptargetfile
    if suffix == ".vcf" :
        vcffile = cptargetfile
        # Convert vcffile into ped
        pedfiles = vcftools_splitbypop( vcffile, patientdata, force=rerun )
        #pedfile = vcftools_makeplink(vcffile, force=rerun )
        #print pedfile
    elif suffix == ".ped" : 
        pedfiles = {"all":targetfile}
    else :
        print "Error: Unknown File type:", basename, suffix
        sys.exit(1)

    hetfiles = []
    festimfiles = []
    for region,ped in pedfiles.iteritems() :
        print region,ped
        recodefile = plink_recode12( ped, force=rerun )
        print "#"*50
        print "Made recodefile:",recodefile
        # Modify ped
        modped = modify_ped( recodefile, patientdata, calcdistances=True, shorten=True, force=rerun )
        bedfile = run_plink_convert(modped, force=rerun)
        #hetfile = run_plink_het( bedfile, force=rerun )
        frqfile = plink_makefrqfile( modped, force=rerun )
        #hapfile = plink_runhaplotyping( modped )
        newfrq, excludebase = excludeSnps( modped, force=rerun )
        bedfile = run_plink_convert(excludebase+".ped", force=rerun)
        hetfile = run_plink_het( bedfile, force=rerun )
        ibcfestim = run_festim( excludebase+".ped", force=rerun )
        hetfiles.append(hetfile)
        festimfiles.append(ibcfestim)
    
    hetfile = "%s/%s.het" % (targetdir, basename)
    mergehetfiles( hetfiles, hetfile, patientdata )
    ibcfestim = "%s/%s.ibc" % (targetdir, basename)
    mergeFEstimFiles( festimfiles, ibcfestim, patientdata )
    print "Hetfile:",hetfile
    print "Festim:",ibcfestim

    plinkibcdata = read_csv( hetfile, sep="\t" )
    festimibcdata = read_csv( ibcfestim, sep="\t" )
    filepath, basename, suffix = getBasename( hetfile )

    #targetfactors = ["Continent2","country", "GeographicRegions", "GeographicRegions2"]
    targetfactors = ["GeographicRegions"]
    for factor in targetfactors :
        #plotPlinkIBC( hetfile, patientdata, factor )
        #plotIBC( ibcfestim, patientdata, factor )
        #ibcCorrelation( hetfile, ibcfestim, patientdata, factor )
        ibcCountryCorrelation( plinkibcdata, basename, "plink", factor )
        #ibcCountryCorrelation( festimibcdata, basename, "festim", factor )

    return recodefile
# End calculateIBC

################################################################################
# calculateIBC_all
################################################################################
def calculateIBC_all( targetfile, rerun=False ) :
    write2log(" - Running "+whoami(), True)
    print "FEstim Targetfile:",targetfile

    # Parse Patient data
    #patientannotation = "./resources/annotation/patientannotation.ped"
    #patientdata = read_csv( patientannotation, sep="\t" )
    targetpats = patientInfo.getPats( targetfile )
    patientdata = sampleAnnotation(targetpats)

    filepath, basename, suffix = getBasename( targetfile )
    print filepath, basename, suffix
    targetdir = "%s/ibc/" % filepath
    makeDir(targetdir)
    filepath, basename, suffix = getBasename( targetfile )
    cptargetfile = targetdir+basename+suffix+".gz"

    if not os.path.exists(cptargetfile) :
        print "Copy targetfile:",cptargetfile
        out = runCommand( "cp %s/%s%s* %s" % (filepath,basename,suffix, targetdir) )
    #if suffix != ".vcf" :
        #print "Error: Unknown File type:", basename, suffix
        #sys.exit(1)

    vcffile = cptargetfile
    pedfile = vcftools_makeplink(vcffile, force=rerun )

    hetfiles = []
    festimfiles = []

    recodefile = plink_recode12( pedfile, force=rerun )
    print "#"*50
    print "Made recodefile:",recodefile
    # Modify ped
    modped = modify_ped( recodefile, patientdata, calcdistances=True, shorten=True, force=rerun )
    bedfile = run_plink_convert(modped, force=rerun)
    #hetfile = run_plink_het( bedfile, force=rerun )
    frqfile = plink_makefrqfile( modped, force=rerun )
    #hapfile = plink_runhaplotyping( modped )
    newfrq, excludebase = excludeSnps( modped, force=rerun )
    hetfile = run_plink_het( excludebase+".ped", force=rerun )
    ibcfestim = run_festim( excludebase+".ped", force=rerun )
    hetfiles.append(hetfile)
    festimfiles.append(ibcfestim)
    
    hetfile = "%s/%s.het" % (targetdir, basename)
    mergehetfiles( hetfiles, hetfile, patientdata )
    ibcfestim = "%s/%s.ibc" % (targetdir, basename)
    mergeFEstimFiles( festimfiles, ibcfestim, patientdata )
    print "Hetfile:",hetfile
    print "Festim:",ibcfestim

    plinkibcdata = read_csv( hetfile, sep="\t" )
    festimibcdata = read_csv( ibcfestim, sep="\t" )
    filepath, basename, suffix = getBasename( hetfile )

    targetfactors = ["Continent2", "country", "GeographicRegions", "GeographicRegions2"]
    for factor in targetfactors :
        plotPlinkIBC( hetfile, patientdata, factor )
        plotIBC( ibcfestim, patientdata, factor )
        ibcCorrelation( hetfile, ibcfestim, patientdata, factor )
        ibcCountryCorrelation( plinkibcdata, basename, "plink", factor )
        ibcCountryCorrelation( festimibcdata, basename, "festim", factor )

    return recodefile
# End calculateIBC_all

################################################################################
# Main
################################################################################
if __name__ == "__main__" :
    os.chdir("..")
    optlist, args = getopt.getopt( sys.argv[1:], "r:ot")
    optlist = dict(optlist)
    filename = None

    dataset = optlist.get("-r",None)
    if dataset is None and not optlist.has_key("-t"):
        print "Error: no dataset provided"
        print " -r <string> - the name of the dataset"
        print " -o flag - overwrite everything"
        print " -t flag - run a test"
        sys.exit(1)

    print "Using dataset:", dataset
    path = os.path.abspath("./rawdata/")
    #vcffile = path+"ciliopathies/ciliopathies.unfilt.vcf.gz"
    #if dataset == "ciliopathies" :
        #filename = "/home/escott/workspace/inbreed/rawdata/ciliopathies/ciliopathies.chimp.recode.vcf.gz"
    if dataset == "daily" :
        filename = "/home/escott/workspace/variome/rawdata/daily/daily.clean.vcf.gz"
    elif dataset == "onekg" :
        filename = "/home/escott/workspace/variome/rawdata/onekg/main/onekg.clean.vcf.gz"
    elif dataset == "mevariome" :
        filename = "/home/escott/workspace/variome/rawdata/mevariome/main/variome.clean.vcf.gz"
    elif dataset == "merge1kg" :
        filename = "/home/escott/workspace/variome/rawdata/merge1kg/main/me1000G.clean.vcf.gz"
    elif dataset == "mergedaly" :
        filename = "/home/escott/workspace/variome/rawdata/mergedaly/main/meceu.clean.vcf.gz"
    elif dataset == "test2" :
        filename = "/home/escott/workspace/variome/rawdata/test2/main/test2.clean.vcf.gz"
    elif dataset == "casanova" :
        filename = "/home/escott/workspace/variome/rawdata/casanova/casanova.snp.recal.chimp.regions.filt.samp.vcf.gz"
    else : print "Error: no dataset provided:",dataset

    outdir = "./results/liftover/tmp"

    #plinkout = "results/plink_ibc.het"
    #plinkout = "rawdata/merged/ibc/merged.chimp.regions.filt.samp.plink.recode12.mod.het"

    changeLogFile( LOGDIR+"/runfestim_test.log" )

    print optlist,"Dataset:",dataset,"Filename:",filename
    if optlist.has_key("-t") : 
        festimout = "results/test_inbreedings.out"
        plinkout = "rawdata/onekg/ibc/onekg.chimp.regions.filt.samp.samp.plink.recode12.mod.het"
        ibcCountryCorrelation( plinkout, festimout )
    elif filename is not None : 
        #calculateIBC_all( filename, rerun=optlist.has_key("-o") )
        calculateIBC( filename, rerun=optlist.has_key("-o") )

    #plotIBC() 
    #ibcCorrelation( plinkout, festimout )
# END Main
    #if len(sys.argv) > 1 :
        #filename = sys.argv[1]
    #else :
        #write2log( "Error: no filename submitted" )

    #elif dataset == "fowzan" :
        #filename = "/home/escott/workspace/variome/rawdata/fowzan/fowzan.snp.recal.chimp.regions.recode.recode.recode.vcf.gz"
    #elif dataset == "test" :
        #filename = "/home/escott/workspace/variome/rawdata/test/everything_set1.chr1.snp.chimp.regions.filt.samp.samp.vcf.gz"
    #elif dataset == "variome" :
        #filename = "/home/escott/workspace/variome/rawdata/variome/variome_snps.regions.filt.samp.samp.vcf.gz"
    #elif dataset == "variome1" :
        #filename = "/home/escott/workspace/variome/rawdata/variome1/variome.regions.filt.samp.samp.vcf.gz"
    #elif dataset == "merged" :
        #filename = "/home/escott/workspace/variome/rawdata/merged/merged.chimp.regions.filt.samp.samp.vcf.gz"
    #elif dataset == "hgdp" :
        #filename = "/home/escott/workspace/variome/rawdata/hgdp/HGDP_938.filt.samp.samp.vcf.gz"


