#!/usr/bin/env python

#-------------------------------------------------------------------#
#                       Import Packages                             #
#-------------------------------------------------------------------#
import os, sys
import getopt
import csv
import re
import gzip
from collections import Counter
from pandas import *
#from random import randint

from housekeeping import *
from localglobals import *
import popgencommands as pop
import patientInfo

######################################################################
# compress
######################################################################
def compress( tfile ) :
    print " - Running compress"
    if tfile[-3:] == ".gz" :
        return tfile

    outfile = tfile+".gz"
    if os.path.exists( outfile ) and forceFlag is False:
        return outfile
    command = "bgzip %s; tabix -p vcf %s" % (tfile,outfile)
    runCommand( command, True )
    return outfile
# END compress

#def parseClinVar( ) :
    #clinvarvcf = "/home/escott/snpEff/clinvar_20140303.vcf.gz"
    #clindata = []
    #for row in csv.reader( CommentStripper(conditionalOpen(clinvarvcf)), delimiter="\t" ) :
        ##if row[0].startswith("#") or len(row) == 0 : continue
        #chrom,pos,vid,ref,alt,qual,vfilter,info = row
        #infodict = dict([x.split("=") if x.find("=") >0 else [x,''] for x in info.split(";")])
        #clinsig = infodict["CLNSIG"]
        #if clinsig.find("|") > 0 or clinsig.find(",")  > 0: 
            #tmp = [ int(x) for x in re.split("\W+", infodict["CLNSIG"]) if int(x) != 255]
            #clinsig = max(tmp) if len(tmp) > 0 else 255
        #newrow = [chrom,pos,vid,ref,alt,clinsig]
        #clindata.append(newrow)
    #clindata = DataFrame(clindata, columns=["chrom","pos","id","ref","alt","clinsig"])
    #return clindata

def parseAnnotation(targetpats, annotationfile="./resources/annotation/patientannotation.ped" ) : 
    sampleannot = read_csv(annotationfile, sep="\t")
    hgdp = read_csv("resources/HGDP_ethnicgroups.txt",sep="\t")
    sampleannot = merge(sampleannot, hgdp[["Ethnicity","Continent","Country"]], left_on="ethnicity",right_on="Ethnicity",how="left")
    print "Targetpats",len(targetpats),targetpats[:10]
    print "annot len:",sampleannot.shape

    sampleannot = sampleannot[sampleannot["Individual.ID"].isin( targetpats )]
    return sampleannot

def subsectionVCFbyPop( vcffile, sampleannot, outprefix, popsprefix ) :
    allpopfiles = {}
    for group, data in sampleannot.groupby("Continent2") :
        pop = group.replace(" ","_")
        popfile = "%s_%s.txt" %(popsprefix, pop )
        data[["Individual.ID"]].to_csv(popfile,sep="\t",index=False,header=False)
        allpopfiles[pop] = popfile
    print allpopfiles

    popvcfs = {}
    for pop in allpopfiles : 
        newvcffile = "%s.%s" % (outprefix, pop)
        popvcfs[pop] = newvcffile+".vcf.gz"
        if os.path.exists(newvcffile+".vcf.gz" ) : continue
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
        runCommand( 'mv %s.recode.vcf %s.vcf' % (newvcffile,newvcffile) )
        cfile = compress( newvcffile+".vcf" )
        print "New compressedfile:", cfile

    return popvcfs
# END subsectionVCFbyPop

def readDbSnp( vcffile ):
    dbsnp = []
    for line in conditionalOpen( vcffile ) :
        if line[0] == "#" : continue
        else : 
            row = line.strip().split("\t",8)
            dbsnp.append(row[2])
    return dbsnp
# END readDbSnp

def vcftools_freq( vcffile, outprefix ) :
    print "Calculate Variant Frequency"
    if os.path.exists( outprefix+".frq" ) : return outprefix+".frq"
    command = "vcftools --gzvcf %s" \
            " --freq" \
            " --out %s" \
            % (vcffile, outprefix)
            #--derived
    runout = runCommand(command)
    if os.path.exists( outprefix+".frq" ) : return outprefix+".frq"
    print "Command:",command
    print "Runout:",runout
# END vcftools_freq

def vcftools_counts( vcffile, outprefix ) :
    print "Calculate Variant Frequency"
    if os.path.exists( outprefix+".frq.count" ) : return outprefix+".frq.count"
    command = "vcftools --gzvcf %s" \
            " --counts2" \
            " --out %s" \
            % (vcffile, outprefix)
            #--derived
    runout = runCommand(command)
    if os.path.exists( outprefix+".frq.count" ) : return outprefix+".frq.count"
    print "Command:",command
    print "Runout:",runout
# END vcftools_counts

def plotFreqCounts( countfile, sampleannot, outprefix="test", factor="ethnicity" ) :
    print "plotFreqCounts",countfile
    chist = {}
    for row in csv.reader(open(countfile),delimiter="\t") :
        if row[0] == "CHROM" : continue
        vcount = int(row[5])
        if not chist.has_key(vcount) : chist[vcount] = 0
        chist[vcount] += 1

    histdf = DataFrame([[x,chist[x]] for x in sorted(chist) if int(x) < 16 and int(x) != 0],columns=["Vcount","Freq"])
    print histdf.head(4)
    total = histdf.Freq.map(int).sum()
    print "Total",total
    histdf["Prop"] = histdf.Freq.map(int) / total
    print histdf.head(3)
    r_dataframe = com.convert_to_r_dataframe(histdf)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x = "factor(Vcount)",y="Prop" ) + \
                ggplot2.geom_bar(stat="identity") + \
                ggplot2.theme(**mytheme)# + \
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})  + \
                #ggplot2.ggtitle("Singletons by "+factor.capitalize()) #+ \
                #ggplot2.scale_y_continuous("Singletons")
                #ggplot2.geom_boxplot(ggplot2.aes_string(fill="Continent")) \
                 #, notch=True ) + \
    figname = "%s_sitefreq.png" % (outprefix)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotFreqCounts

def vcftools_singletons( vcffile, outprefix ) :
    print "Calculate singletons"
    if os.path.exists( outprefix+".singletons" ) : return outprefix+".singletons"
    command = "vcftools --gzvcf %s" \
            " --singletons" \
            " --out %s" \
            % (vcffile, outprefix)
    runout = runCommand(command)
    if os.path.exists( outprefix+".singletons" ) : return outprefix+".singletons"
    print "Command:",command
    print "Runout:",runout
# END vcftools_singletons

def plotSingletons( singletonfile, sampleannot, outprefix="test", factor="GeographicRegions" ) :
    print "plotSingletons",singletonfile
    singletons = read_csv(singletonfile, sep="\t")
    singletons.columns =  ["CHROM","POS","Vtype","ALLELE","INDV"]
    # DataFrame({'count' : df1.groupby( [ "Name", "City"] ).size()}).reset_index()
    singleton_counts = DataFrame({'count': singletons.groupby(["INDV","Vtype"]).size()}).reset_index()
    print singleton_counts.head(3)
    singleton_counts["Vtype"] = ["Singleton" if x == "S" else "Doubleton" 
                                 for x in singleton_counts.Vtype.tolist()]

    #print singletons.head(10)
    singletons_merge = merge(singleton_counts,
                             sampleannot[["Individual.ID",factor,"Source","Continent2"]], 
                             left_on="INDV", right_on="Individual.ID", how="inner")
    #singletons_merge = singletons_merge[singletons_merge[factor].notnull()]
    singletons_merge["count"] = singletons_merge["count"].map(int)
    r_dataframe = com.convert_to_r_dataframe(singletons_merge)

    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor("+factor+")",y="count",fill="factor(Continent2)" ) +
                ggplot2.geom_boxplot() +
                ggplot2.ggtitle("Singletons by "+factor.capitalize()) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})  +
                ggplot2.scale_y_continuous("Variant count") +
                ggplot2.scale_x_discrete(string.capitalize(factor)) +
                ggplot2.theme(**mytheme) +
                ggplot2.facet_grid( robjects.Formula('Vtype ~ .'), scale="free") )
                #ggplot2.geom_boxplot(ggplot2.aes_string(fill="Continent")) \
                 #, notch=True ) + \
    figname = "%s_%s_singletons.png" % (outprefix,factor)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotSingletons

def plotDbSnp( countfile, dbsnpids, header="", outprefix="test", climit=1 ) :
    countdat = read_csv( countfile, delim_whitespace=True, header=None, skiprows=[0], names=["chrom","pos","alleles","N","Ref","Alt"] )
    print countdat.head(3)
    print countdat.columns
    print countdat.index[:3]

    #countdat.columns = ["chrom","pos","alleles","N","Ref","Alt"]
    countdat["dbsnp"] = ["Novel" if x is "." else "DBSNP" for x in dbsnpids]
    countdat = countdat[countdat["Alt"] > climit]
    print countdat.head(3)

    #pie <- ggplot(mtcars, aes(x = factor(1), fill = factor(cyl))) +
    #geom_bar(width = 1)
    #pie + coord_polar(theta = "y")
    r_dataframe = com.convert_to_r_dataframe(countdat)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x = "factor(1)", fill = "factor(dbsnp)") + \
                ggplot2.geom_bar(width=1) + \
                ggplot2.coord_polar(theta="y") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.ggtitle(header.replace("_"," ")+" DBSNP Proportion") #+ \

    figname = "%s_dbsnp.png" % (outprefix)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotDbSnp

def siteFreqSpectrum( frqfile, dbsnpids, header="", outprefix="results/figures/test" ) :
    countdat = read_csv( frqfile, delim_whitespace=True, header=None, skiprows=[0], names=["chrom","pos","alleles","N","Ref","Alt"] )
    countdat["dbsnp"] = [ "no" if x == "." else "yes" for x in dbsnpids ]
    countdat["AF"] = [float(x.split(":")[1]) for x in countdat.Alt.tolist()]
    bins=[0.,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]
    #afbinned = DataFrame({'vcount':countdat.groupby(["dbsnp",np.digitize(countdat.AF, bins)]).size()}).reset_index()
    afbinned = DataFrame({'vcount':countdat.groupby(["dbsnp",cut(countdat.AF, bins)]).size()}).reset_index()
    afbinned["proportions"] = afbinned["vcount"] / float(len(countdat))

    r_dataframe = com.convert_to_r_dataframe(afbinned)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="factor(AF)", y="proportions", fill="factor(dbsnp)") + \
                ggplot2.geom_bar(stat="identity",position="dodge") + \
                ggplot2.ggtitle("Proportions of DBSNP by AF") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                ggplot2.scale_y_continuous("Proportion")
                #ggplot2.facet_grid( robjects.Formula('group ~ .') )

    figurename = "%s_sf_dbsnp.png" % outprefix
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# END siteFreqSpectrum

def plotDBSnpBar( frqfile, dbsnpids, header="", outprefix="results/figures/test" ) :
    countdat = read_csv( frqfile, delim_whitespace=True, header=None, skiprows=[0], names=["chrom","pos","alleles","N","Ref","Alt"] )
    countdat["dbsnp"] = [ "no" if x == "." else "yes" for x in dbsnpids ]
    countdat["AF"] = [float(x.split(":")[1]) for x in countdat.Alt.tolist()]
    print countdat.head(10)
    bins=[0.,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]
    #afbinned = DataFrame({'vcount':countdat.groupby(["dbsnp",np.digitize(countdat.AF, bins)]).size()}).reset_index()
    afbinned = DataFrame({'vcount':countdat.groupby(["dbsnp",cut(countdat.AF, bins)]).size()}).reset_index()
    vcounttotals = afbinned[["AF","vcount"]].groupby("AF",as_index=False).sum() 
    vcounttotals.columns = ["AF","vtotals"]
    afbinned = merge(afbinned,vcounttotals,on="AF")
    print vcounttotals.head(10)
    #afbinned["proportions"] = afbinned["vcount"] / float(len(countdat))
    afbinned["proportions"] = afbinned["vcount"] / afbinned["vtotals"].map(float)
    print afbinned.head(10)

    ptitle = "Proportions of DBSNP by AF"
    if len(header) > 0 : ptitle = "Proportions of DBsnp by AF\n(%s)" %header
    r_dataframe = com.convert_to_r_dataframe(afbinned)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="factor(AF)", y="proportions", fill="factor(dbsnp)") + \
                ggplot2.geom_bar(stat="identity") + \
                ggplot2.ggtitle(ptitle) + \
                ggplot2.theme(**mytheme) + \
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                ggplot2.scale_y_continuous("Proportion") + \
                ggplot2.scale_x_discrete("AF binned")
                #ggplot2.facet_grid( robjects.Formula('group ~ .') )

    figurename = "%s_sf_dbsnp.png" % outprefix
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# END plotDBSnpBar

def vcftools_het( vcffile, outprefix) :
    print "Calculate heterozygosity"
    if os.path.exists( outprefix+".het" ) : return outprefix+".het"
    command = "vcftools --gzvcf %s" \
            " --het" \
            " --out %s" \
            % (vcffile, outprefix)
    runout = runCommand(command)
    if os.path.exists( outprefix+".het" ) : return outprefix+".het"
    print "Command:",command
    print "Runout:",runout
# END vcftools_het

def vcftools_relatedness( vcffile, outprefix ) :
    print "Calculate relatedness"
    if os.path.exists( outprefix+".relatedness" ) : return outprefix+".relatedness"
    command = "vcftools --gzvcf %s" \
            " --relatedness" \
            " --out %s" \
            % (vcffile, outprefix)
    runout = runCommand(command)
    if os.path.exists( outprefix+".relatedness" ) : return outprefix+".relatedness"
    #print "Command:",command
    #print "Runout:",runout
# END vcftools_relatedness
 
def vcftools_relatedness2( vcffile, outprefix ) :
    print "Calculate relatedness"
    if os.path.exists( outprefix+".relatedness2" ) : return outprefix+".relatedness2"
    command = "vcftools --gzvcf %s" \
            " --relatedness2" \
            " --out %s" \
            % (vcffile, outprefix)
    runout = runCommand(command)
    if os.path.exists( outprefix+".relatedness2" ) : return outprefix+".relatedness2"
    #print "Runout:",runout
# END vcftools_relatedness2
 
def vcftools_indelhist( vcffile, outprefix ) :
    print "Calculate indel hist"
    if os.path.exists( outprefix+".indel.hist" ) : return outprefix+".indel.hist"
    command = "vcftools --gzvcf %s" \
            " --hist-indel-len" \
            " --out %s" \
            % (vcffile, outprefix)
    runout = runCommand(command)
    if os.path.exists( outprefix+".indel.hist" ) : return outprefix+".indel.hist"
    print "Command:",command
    print "Runout:",runout
# END vcftools_indelhist

def vcftools_genor2( vcffile, outprefix ) :
    print "Calculate geno r^2"
    if os.path.exists( outprefix+".geno" ) : return outprefix+".geno"
    command = "vcftools --gzvcf %s" \
            " --geno-r2" \
            " --out %s" \
            % (vcffile, outprefix)
    runout = runCommand(command)
    if os.path.exists( outprefix+".geno" ) : return outprefix+".geno"
    print "Command:",command
    print "Runout:",runout
# END vcftools_genor2

def vcftools_titv( vcffile, outprefix ) :
    print "Calculate tstv summary"
    if os.path.exists( outprefix+".titv" ) : return outprefix+".titv"
    command = "vcftools --gzvcf %s" \
            " --TsTv-summary" \
            " --out %s" \
            % (vcffile, outprefix)
    runout = runCommand(command)
    if os.path.exists( outprefix+".titv" ) : return outprefix+".titv"
    print "Command:",command
    print "Runout:",runout
# END vcftools_titv

def vcftools_hardy( vcffile, outprefix, force=False ) :
    print "Calculate hardy wineburg per site"
    if os.path.exists( outprefix+".hwe" ) and not force : return outprefix+".hwe"
    command = "vcftools --gzvcf %s" \
            " --hardy" \
            " --out %s" \
            % (vcffile, outprefix)
    runout = runCommand(command)
    if os.path.exists( outprefix+".hwe" ) : return outprefix+".hwe"
    print "Command:",command
    print "Runout:",runout
# END vcftools_hardy

def vcftools_sitespi( vcffile, outprefix ) :
    print "Calculate Nucleotide divergency"
    if os.path.exists( outprefix+".sites.pi" ) : return outprefix+".sites.pi"
    command = "vcftools --gzvcf %s" \
            " --site-pi" \
            " --out %s" \
            % (vcffile, outprefix)
    runout = runCommand(command)
    if os.path.exists( outprefix+".sites.pi" ) : return outprefix+".sites.pi"
    print "Command:",command
    print "Runout:",runout
# END vcftools_sitespi

def vcftools_fst( vcffile, sampleannot, outprefix, popsprefix ) :
    allpopfiles = {}
    for group, data in sampleannot.groupby("ethnicity") :
        popfile = "%s_%s.txt" %(popsprefix, group)
        data[["Individual.ID"]].to_csv(popfile,sep="\t",index=False,header=False)
        allpopfiles[group] = popfile
    print allpopfiles
    command = "vcftools --gzvcf %s" \
            " --weir-fst-pop %s" \
            " --out %s" \
            % (vcffile, " --weir-fst-pop ".join(allpopfiles.values()), outprefix)
    runout = runCommand(command)
    if os.path.exists( outprefix+".weir.fst" ) : return outprefix+".weir.fst"
    print "command:",command
    print runout

    #singletonsfile =  "%s/%s.singletons" % (targetdir, filename)
# END vcftools_fst

def CNV_plot(cnv_file, subdir, png):
    print """Make CNV plots using R cnv library"""

    os.system('mkdir -p working/cnv_seq/CNV')
    os.system('mkdir -p working/cnv_seq/CNV/' 
              + subdir + '/')
    output_cnv_file = 'working/cnv_seq/CNV/' + subdir + '/' + cnv_file.split('.hits')[0].rstrip('T')  + '.cnvs'
    os.system('rm ' + output_cnv_file)
    rtmp = 'rtmp' + str(random.randint(0,1000))
    with open(rtmp, 'w') as f:
        f.write("source('funcs.R')\n")
        f.write('library(cnv)\n')
        f.write("data<-read.delim('"
                + cnv_file + "')\n")
        f.write("png('" + png + "')\n")
        f.write("plot.cnv.all.perry(data,colour=9)\n")
        f.write('dev.off()\n')
        f.write("cnv.print(data, file='" + output_cnv_file + "')\n")
        f.write('q()\n')
    if utils.check_input(cnv_file):
        os.system('R CMD BATCH --vanilla ' + rtmp + ' tmpLog')
        os.system('rm ' + rtmp + ' tmpLog')
    os.system('mv ' + cnv_file + ' working/cnv_seq/CNV/' + subdir + '/')
    os.system('mv ' + cnv_file.replace('cnv', 'count') 
              + ' working/cnv_seq/CNV/' + subdir + '/')

def run_vcftools_custom( vcffile, keepfile, outprefix, force=False ) :
    assert os.path.exists(vcffile)

    if os.path.exists( outprefix+".clean.vcf.gz" ) and not force :
        print " - File"+outprefix+".clean.vcf.gz : already exists. Skipping command!"
        return outprefix+".clean.vcf.gz"

    # remove-indels
    # mix/max-alleles - Preven n allelic sites
    # snps <filename> include snps based on list of ids
    # Geno - removes genotypes without sufficient coverage
    command = "vcftools --gzvcf %s" \
        " --min-alleles 2 --max-alleles 2" \
        " --max-missing 0.95" \
        " --remove-filtered-all" \
        " --recode" \
        " --out %s" \
        % ( vcffile, outprefix )
        #" --min-meanDP 4 --minQ 20" \

    print command
    out = runCommand( command )

    assert os.path.exists( outprefix+".recode.vcf" ), out
    runCommand( 'mv %s.recode.vcf %s.clean.vcf' % (outprefix,outprefix) )
    cfile = compress( outprefix+".clean.vcf" )
    print " - End run_vcftools"
    return cfile
# END run_vcftools_custom

def plotHardy( filedf, prefix ) :
    allpvals = []
    for index, row in filedf.iterrows() :
        print row
        print row["pop"], row["hardyfile"]
        hardydata = read_csv(row["hardyfile"], sep="\t" )
        print hardydata.head(10)
        hardyfilt = hardydata[hardydata.P > 0]
        print "Before:",len(hardydata), "After:",len(hardyfilt)
        allpvals.append(DataFrame({'pval':hardyfilt.P, 'grp':row['pop']}))
        qqplot( hardyfilt.P, figname=prefix+"_"+row["pop"]+".png", 
               ptitle="QQ of HWE in "+row["pop"]+" variants" )
    allpvals = concat( allpvals ).reset_index(drop=True)
    print allpvals.head(10)
    qqplot( allpvals.pval, pfactor=allpvals.grp, figname=prefix+"_qq.png", 
               ptitle="QQ of HWE for all pops" )
# END plotHardy

################################################################################
# run_vcftools_analyze
################################################################################
def run_vcftools_analyze( vcffile, keepfile=None, outdir=None, force=False ) :
    print " - Running: run_vcftools_analyze"

    targetdir, filename, suffix = getBasename(vcffile)
    if outdir is None : outdir = targetdir
    makeDir("%s/vstats/" % outdir )
    makeDir("%s/vstats/pops/" % outdir )
    outprefix = "%s/vstats/%s" % (outdir, filename)
    popsprefix = "%s/vstats/pops/%s" % (outdir, filename)

    outfilefile = "%s.vcftoolsfiles" % outprefix
    if os.path.exists(outfilefile) and not force : 
        outfiles = read_csv( outfilefile, sep="\t" )
        return outfiles

    #if outfilefile is not None : 
        #print "Fixing vcffile:",vcffile
        #vcffile = pop.run_vcftools_clean( vcffile, keepfile, outprefix )
        #print "New vcffile:", vcffile
        
    targetpats = patientInfo.getPats( vcffile )
    #sampleannot = parseAnnotation(targetpats)
    sampleannot = sampleAnnotation( targetpats )

    outfiles = []

    singlefile = vcftools_singletons( vcffile, outprefix )
    plotSingletons( singlefile, sampleannot, outprefix, factor="GeographicRegions2")

    popvcfs = subsectionVCFbyPop( vcffile, sampleannot, outprefix, popsprefix )
    for pop in popvcfs :
        targetdir, filename, suffix = getBasename(popvcfs[pop])
        outprefix = "%s/vstats/%s" % (outdir, filename)
        dbsnpinfo = readDbSnp( popvcfs[pop] )
        vcountfile= vcftools_counts( popvcfs[pop], outprefix )
        plotDbSnp( vcountfile, dbsnpinfo, pop, outprefix )
        plotFreqCounts( vcountfile, sampleannot, outprefix )
        vfreqfile= vcftools_freq( popvcfs[pop], outprefix )
        siteFreqSpectrum( vfreqfile, dbsnpinfo, outprefix=outprefix )
        plotDBSnpBar( vfreqfile, dbsnpinfo, header=pop, outprefix=outprefix )
        singlefile= vcftools_singletons( popvcfs[pop], outprefix )
        #plotSingletons( singlefile, sampleannot, outprefix )
        hardyfile= vcftools_hardy( popvcfs[pop], outprefix, force=force )
        outfiles.append([pop,popvcfs[pop],vcountfile,vfreqfile,singlefile,hardyfile])

    outfiles=DataFrame(outfiles,columns=
                       ["pop","vcffile","vcountfile","vfreqfile","singlefile","hardyfile"])
    
    outfiles.to_csv( outfilefile, sep="\t", index=False )
    print outfiles.head(10)
    plotHardy( outfiles, popsprefix )
    return outfiles
# END run_vcftools_analyze
    #s= vcftools_relatedness( vcffile, outprefix )
    #s= vcftools_fst( vcffile, sampleannot, outprefix, popsprefix )
    #s= vcftools_counts( vcffile, outprefix )
    #s= vcftools_het( vcffile, outprefix )
    #s= vcftools_genor2( vcffile, outprefix )
    #s= vcftools_titv( vcffile, outprefix )
    #s= vcftools_sitespi( vcffile, outprefix )

def clinvarAnalysis( countfile ) :
    countdat = read_csv( countfile, delim_whitespace=True, header=None, skiprows=[0], names=["chrom","pos","alleles","N","Ref","Alt"] )
    clinvar["key"] = "chr"+clinvar.chrom.map(str)+":"+clinvar.pos.map(str)
    print clinvar.head(3)
    countdat["key"] = "chr"+countdat.chrom.map(str)+":"+countdat.pos.map(str)
    print countdat.head(3)

    clincounts = merge( clinvar, countdat, on="key", how="left" )
    print "Missing:",sum(clincounts.Alt.isnull())
    print "Non-Missing:",sum(clincounts.Alt.notnull())
    #clincounts["Alt"] = clincounts["Alt"].fillna(0)
    clincounts_filt = clincounts[clincounts.Alt.notnull()]
    r_dataframe = com.convert_to_r_dataframe(clincounts_filt)
    p = ggplot2.ggplot(r_dataframe) + \
            ggplot2.aes_string(x = "factor(1)", fill = "factor(clinsig)") + \
            ggplot2.geom_bar(width=1) + \
            ggplot2.coord_polar(theta="y") + \
            ggplot2.theme(**mytheme) + \
            ggplot2.ggtitle("Clinvar overlap") #+ \

    #figname = "%s_dbsnp.png" % (outprefix)
    figname = "results/clinvar.png"
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
    clincounts["Found"] = ["no" if math.isnan(found) else clin for clin, found in clincounts[["clinsig","Alt"]].values ]
    print clincounts.head(10)
    r_dataframe = com.convert_to_r_dataframe(clincounts)
    p = ggplot2.ggplot(r_dataframe) + \
            ggplot2.aes_string(x = "factor(1)", fill = "factor(Found)") + \
            ggplot2.geom_bar(width=1) + \
            ggplot2.coord_polar(theta="y") + \
            ggplot2.theme(**mytheme) + \
            ggplot2.ggtitle("Clinvar overlap") #+ \

    #figname = "%s_dbsnp.png" % (outprefix)
    figname = "results/clinvar_found.png"
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END clinvarAnalysis

################################################################################
# identifyUnrelatedSamples
################################################################################
def identifyUnrelatedSamples( vcffile, targetdir=None, force=False ): 

    filepath, filename, suffix = getBasename(vcffile)
    if targetdir is not None :
        assert os.path.exists(targetdir)
        outprefix = "%s/%s" % (targetdir, filename) 
    else : outprefix = "%s/%s" % (filepath, filename)
    filepats = patientInfo.currentFilePats(vcffile)
    sampleannot = sampleAnnotation(filepats)

    unrelatedfile = "%s.unrelated" % outprefix
    if os.path.exists(unrelatedfile) and not force :
        return unrelatedfile

    print "Starting Unrelated:",len(filepats)
    #relatednessfile = vcftools_relatedness( vcffile, outprefix )
    relatednessfile = vcftools_relatedness2( vcffile, outprefix )
    maxrelatedness = 0.025 # Degree = 5
    maxrelatedness = 0.05 # Degree = 4
    maxrelatedness = 0.13 # Degree = 3
    reldf = read_csv( relatednessfile, delim_whitespace=True, header=None, skiprows=1,
                     names=["INDV1","INDV2","RELATEDNESS_AJK"])
    relationships = reldf[reldf.INDV1 != reldf.INDV2] # remove exact comparisons
    relationships.RELATEDNESS_AJK = relationships.RELATEDNESS_AJK.map(float)

    #print "Could prioritize sample selection here"
    are_related = relationships.RELATEDNESS_AJK > maxrelatedness

    fambins = [0, 0.05, 0.13, 0.25, .55, 1]
    famlabels = ["4th & Greater","3rd degree","2nd degree","1st degree","duplicate" ]

    # duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree
    relationships["Fclass"] = cut(relationships.RELATEDNESS_AJK, 
                                  bins=fambins, retbins=True, labels=famlabels)[0]
    print relationships.head()
    
    while sum(are_related) > 1 :
        #relatedsub = relationships[are_related]
        currsamps = (relationships[are_related].INDV1.tolist() + 
                     relationships[are_related].INDV2.tolist() )
        relcounts = Counter(currsamps).most_common(1)
        print "Remove", relcounts
        tsamp, count = relcounts[0]
        relationships = (relationships[(relationships.INDV1 != tsamp) & 
                                      (relationships.INDV2 != tsamp)])
        targetindexes = relationships.index[(relationships.INDV1 == tsamp) & 
                                      (relationships.INDV2 == tsamp)].tolist()
        print "Length of targetindexes:",len(targetindexes)
        relationships.drop(targetindexes, inplace=True)
        #print "Related",relationships.shape
        #unrelatedset = hk.uniqify(relationships.INDV1.tolist() + relationships.INDV2.tolist())
        #print "Remaining",len(unrelatedset)
        are_related = relationships.RELATEDNESS_AJK > maxrelatedness
        #print "are relationships count:", sum(are_related)

    unrelatedset = hk.uniqify(relationships.INDV1.tolist() + relationships.INDV2.tolist())
    print "Finalcount:",len(unrelatedset)
    #print unrelatedset[:10]

    OUT = open(unrelatedfile,"wb")
    for samp in unrelatedset :
        OUT.write(samp+"\n")
    OUT.close()

    return unrelatedfile
# end identifyUnrealtedSamples

################################################################################
# Main
################################################################################
if __name__ == "__main__" :
    os.chdir("..")

    optlist, args = getopt.getopt( sys.argv[1:], "r:ot")
    optlist = dict(optlist)

    dataset = optlist.get("-r",None)
    if dataset is None and not optlist.has_key("-t"):
        print "Error: no dataset provided"
        print " -r <string> with the name of the dataset"
        print " -o flag to overwrite everything"
        sys.exit(1)

    #filename = "/home/escott/workspace/variome/rawdata/onekg/onekg.chimp.regions.filt.samp.samp.vcf.gz"
    #filename = "/home/escott/workspace/variome/rawdata/test/everything_set1.chr1.snp.chimp.regions.filt.samp.samp.vcf.gz"
    #filename = "/home/escott/workspace/variome/rawdata/merged/merged.chimp.regions.filt.samp.samp.vcf.gz"
    #filename = "/home/escott/workspace/variome/rawdata/daily/daily.chimp.regions.filt.samp.samp.vcf.gz"
    #filename = "/home/escott/workspace/variome/rawdata/variome1/variome.chimp.regions.filt.samp.samp.vcf.gz"
    path = os.path.abspath("./rawdata/")
    if dataset == "test" :
        vcffile = path+"/test/everything_set1.chr1.snp.clean.vcf.gz"
    #elif dataset == "onekg" :
        #vcffile = path+"/onekg/onekg.clean.vcf.gz"
    elif dataset == "daily" :
        vcffile = path+"/daily/main/daily.clean.vcf.gz"
    elif dataset == "onekg" :
        vcffile = path+"/onekg/main/onekg.clean.vcf.gz"
    #elif dataset == "eichler" :
        #vcffile = path+"/eichler/eichler.clean.vcf.gz"
    #elif dataset == "turks" :
        #vcffile = path+"/turks/turks.clean.vcf.gz"
    #elif dataset == "variome1" :
        #vcffile = path+"/variome1/variome.clean.vcf.gz"
    #elif dataset == "variome" :
        #vcffile = path+"/variome/main/variome.clean.vcf.gz"
    #elif dataset == "hgdp" :
        #vcffile = path+"/hgdp/HGDP_938.clean.vcf.gz"
    #elif dataset == "merged" :
        #vcffile = path+"/merged/main/merged.clean.vcf.gz"
    elif dataset == "casanova" :
        vcffile = path+"/casanova/casanova.snp.recal.clean.vcf.gz"
    elif dataset == "meceu" :
        vcffile = path+"/mergedaly/main/meceu.clean.vcf.gz"
    elif dataset == "merge1kg" :
        vcffile = path+"/merge1kg/main/me1000G.clean.vcf.gz"
    elif dataset == "mevariome" :
        vcffile = path+"/mevariome/main/variome.clean.vcf.gz"


    if optlist.has_key("-t") :
        clinvar = parseClinVar()
        print clinvar.head(10)
        #vcffile = "/home/escott/workspace/variome/rawdata/test/vstats/everything_set1.chr1.snp.Middle_East.vcf.gz"
        #countfile = "/home/escott/workspace/variome/rawdata/test/vstats/everything_set1.chr1.snp.Middle_East.frq.count"
        #countfile = "/home/escott/workspace/variome/rawdata/daily/vstats/daily.Europe.frq.count"
        #clinvarAnalysis( countfile )

        #vcffile = "/home/escott/workspace/variome/rawdata/daily/vstats/daily.Middle_East.vcf.gz"
        #frqfile = "/home/escott/workspace/variome/rawdata/daily/vstats/daily.Middle_East.frq"
        #vcffile = "/home/escott/workspace/variome/rawdata/daily/daily.clean.vcf.gz"
        vcffile = "/home/escott/workspace/variome/rawdata/test/everything_set1.chr1.snp.vcf.gz"
        vcffile = "/media/data/workspace/variome/rawdata/variome2/qc/variome.recal.clean.vcf.gz"

        #dbsnpinfo = readDbSnp( vcffile )
        #siteFreqSpectrum( frqfile, dbsnpinfo )
        #plotDBSnpBar( frqfile, dbsnpinfo, header="Middle East" )
        unrelated = identifyUnrelatedSamples( vcffile, force=True )
        print "Unrelated:",len(unrelated)
    else :
        print "Using dataset:", dataset
        #keepfile = "%s/%s/tokeep.txt" % (path, dataset)
        #assert os.path.exists( keepfile )
        run_vcftools_analyze( vcffile, force=True )

# END Main

