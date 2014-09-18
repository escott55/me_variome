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
gtable = importr('gtable',on_conflict="warn")

#genotypename: input genotype file (in any format: see ../CONVERTF/README)
#snpname:      input snp file      (in any format: see ../CONVERTF/README)
#indivname:    input indiv file    (in any format: see ../CONVERTF/README)
#evecoutname:  output file of eigenvectors.  See numoutevec parameter below.
#evaloutname:  output file of all eigenvalues
#numoutevec:     number of eigenvectors to output.  Default is 10.
#numoutlieriter: maximum number of outlier removal iterations.
  #Default is 5.  To turn off outlier removal, set this parameter to 0.
#numoutlierevec: number of principal components along which to 
  #remove outliers during each outlier removal iteration.  Default is 10.
#outliersigmathresh: number of standard deviations which an individual must 
  #exceed, along one of the top (numoutlierevec) principal components, in
  #order for that individual to be removed as an outlier.  Default is 6.0.

################################################################################
# modifyFamFile
################################################################################
def modifyFamFile( famfile, targetfile, sampleannot ): 
    print "Making file:",targetfile
    famdata = read_csv(famfile, delim_whitespace=True, header=None,
                               names=["Fam","IID","MID","PID","Gender","AFF"]) 
    famdata = famdata.fillna(0)
    famdata.index = famdata.IID
    famdata["Region"] = "Unk"
    famdata.update(sampleannot)
    #famdata = addSampleAnnotation( famdata, update=True )
    #famdata[annotcol] = famdata[annotcol].apply(lambda x: 
                        #x.strip().replace(' ','.'))
    famdata["Region"] = famdata["Region"].apply(lambda x:
                        x.strip().replace(' ','.'))
           
    famdata["Gender"] = famdata.Gender.astype(int)
    #famdata["Gender"] = [1 if x == "M" else 2 for x in famdata.Gender]

    famdata.to_csv(targetfile, sep=" ", header=False, index=False, 
                   columns=["Fam","IID","MID","PID","Gender","Region"])
    assert os.path.exists(targetfile)

    return famdata['Region'].unique().tolist()
# END modifyFamFile

################################################################################
# makeParFile
################################################################################
def makeParFile( bedfile, sampleannot, targetdir=None, 
                numvec=20, fstonly=False ):
    filepath,basename,ext = hk.getBasename(bedfile)
    hk.makeDir(targetdir)
    bedfile = "%s/%s.bed" % (filepath,basename)
    bimfile = "%s/%s.bim" % (filepath,basename)
    famfile = "%s/%s.fam" % (filepath,basename)

    mapfile = "%s/%s.map" % (targetdir,basename)
    pedind = "%s/%s.pedind" % (targetdir,basename)
    #subprocess.call(["cp", famfile, pedind] )
    #if annotcol in ["GRsub"] : annotcol = "GeographicRegions"
    regions = modifyFamFile( famfile, pedind, sampleannot )
    subprocess.call(["cp", bimfile, mapfile] )

    eigenvec = "%s/%s.evec" % (targetdir, basename)
    eigenval = "%s/%s.eval" % (targetdir, basename)
    grmout = "%s/%s.grm" % (targetdir, basename)
    outliers = "%s/%s.outliers" % (targetdir, basename)

    parfile = "%s/%s.par" % (targetdir, basename)
    OUT = open(parfile,"wb")
    OUT.write(  ("genotypename:\t%s\n"%bedfile +
                "snpname:\t%s\n" %mapfile +
                "indivname:\t%s\n"% pedind +
                "evecoutname:\t%s\n"%eigenvec + 
                "evaloutname:\t%s\n"%eigenval +
                "outlieroutname:\t%s\n"%outliers +
                "altnormstyle:\tNO\n"+
                "fstonly:\t%s\n"%("YES" if fstonly else "NO") +
                "outliersigmathresh:\t6.0\n"+
                "numoutevec:\t%d\n"%numvec +
                "numoutlierevec:\t5\n" +
                "familynames:\tNO\n" +
                "grmoutname:\t%s\n"%grmout ))
    OUT.close()

    return parfile,outliers,eigenvec,eigenval,regions
# END makeParFile

################################################################################
# parseOutliers
################################################################################
def parseOutliers(outliers) :
    samps = []
    if os.path.exists(outliers) :
        outliers = open(outliers,"r").readlines()
        for line in outliers :
            row = line.split(" ")
            samps.append( row[2] )
    print samps
    print "# outliers:",len(samps)
    return samps
# END parseOutliers

def getTargetRegions( regions, annotcol ):
    print "Regions:",regions
    if annotcol=="Continent" :
        targetregions = ["Middle.East","South.Asia","Europe","Africa","East.Asia"]
        regions = [x for x in regions if x in targetregions]
    if annotcol=="Continent2" :
        #targetregions = ["Middle.East","South.Asia","Europe","Africa","East.Asia"]
        excluderegions = ["Unknown","America","Oceania"]
        regions = sorted([x for x in regions if x not in excluderegions])
    elif annotcol=="ethnicity" :
        targetregions = ["Bedouin","Adygei","Druze","Mozabite","Palestinian"]
        regions = [x for x in regions if x in targetregions]
    elif annotcol=="Origin" :
        targetregions = ["Morocco","Algeria","Tunisia","Libya","Egypt","Saudi Arabia","Oman","Qatar","UAE","Yemen","Jordan","Palestine","Lebanon","Syria","Kuwait","Iraq","Turkey"] + ["Iran","Pakistan","Afganistan"]
        regions = [x for x in regions if x in targetregions]
    elif annotcol=="GeographicRegions" :
        exclude = ["Oceania","America","Unknown"]
        targetregions = [x.replace(" ",".") 
                         for x in target_geographic_regions if x not in exclude]
        regions = [x for x in targetregions if x in regions]
    elif annotcol=="GeographicRegions2" :
        targetregions = [x.replace(" ",".") for x in target_geographic_regions2]
        print targetregions
        regions = [x for x in regions if x in targetregions if x != "Unk"]
    elif annotcol=="GRsub" :
        #targetregions = [ 'Northwest.Africa', 'Northeast.Africa', 'Arabian.Peninsula', 'Syrian.Desert',
                        #'Turkish.Peninsula', 'Central.Asia' ]
        targetregions = [x.replace(" ",".") for x in target_geographic_regions_me]
        regions = [x for x in regions if x in targetregions]
    elif annotcol=="Source" :
        regions = regions
    else :
        print "Error: Unknown annotation:",annotcol

    print "Regions:",regions
    return regions
# END getTargetRegions

################################################################################
# pcaOutlierAnalysis
################################################################################
def pcaOutlierAnalysis( eigfiles, force=False, 
                       annotcol ="Continent", tcols="1:2" ):
    print "Running pcaOutlierAnalysis" 

    parfile,outliers,eigenvec,eigenval,regions = eigfiles
    targetdir,basename,ext = hk.getBasename(parfile)
    #targetdir = targetdir if targetdir is not None else filepath
    #print "Using targetdir:",targetdir
    #hk.makeDir(targetdir)
    outbase = "%s_%s" %(basename, tcols.replace(":","_"))

    if not os.path.exists(outliers) or force :
        print "Parfile:",parfile
        try :
            print "Running smartpca"
            subout = subprocess.check_output( ["smartpca", "-p",parfile] )
        except subprocess.CalledProcessError as e :
            print "Caught error:",e

    tregions = getTargetRegions( regions, annotcol )

    assert len(tregions) > 0
    plotfile = "%s/%s.xtxt" % (targetdir, outbase)
    command = ["ploteig","-i",eigenvec,"-c",tcols,
                     "-p",":".join(tregions),
                     "-x","-o",plotfile]

    print " ".join(command)
    print "Trying subprocess"
    print "Making file:",outbase,annotcol,".pdf"
    retval = subprocess.call(command)
    retval = subprocess.call(["mv","%s.pdf"%outbase, 
                              "results/pcaplots/%s_%s.pdf" 
                              %(outbase,annotcol)])
    #print retval
    #print "Trying system!"
    #retval = os.system( " ".join(command) )
    #print retval
    twtable = "/home/escott/Packages/EIG5.0.2/POPGEN/twtable"
    twout = "%s/%s.twout" % (targetdir, outbase)
    subprocess.call(["twstats","-t",twtable,"-i",eigenval,"-o",twout])
    outliersamps = parseOutliers(outliers)

    return outliersamps
# END pcaOutlierAnalysis

def makeFstPlot( fstdata, annotcol, basename, targetdir ):
    newdata = []
    allregions = hk.uniqify(fstdata.Region2.tolist()+fstdata.Region1.tolist())
    print "All Regions",allregions
    for region in allregions :
        newdata.append(
            Series({'Region1':region, 'Region2':region,'Fst':0.0,'stderr':0.0})
        )

    invfst = fstdata.copy()
    invfst["Region1"] = fstdata["Region2"]
    invfst["Region2"] = fstdata["Region1"]

    fstdata = concat([fstdata, DataFrame(newdata), invfst]).reset_index()

    #print [fstdata.Region2.tolist().count(x) for x in fstdata.Region2.unique()]
    fstdata = fstdata[(fstdata.Region1 != "Unk") & (fstdata.Region2 != "Unk")]
    maxfst = fstdata.Fst.max()
    midfst = maxfst / 2
    #print fstdata.head()
    #fstdata = concat( [fstdata, invfst] ).reset_index(drop=True)

    targetlevels = None
    if annotcol == "GeographicRegions" :
        targetlevels = target_geographic_regions

    elif annotcol == "GRsub" :
        targetlevels = target_geographic_regions_me

    elif annotcol == "GeographicRegions2" :
        targetlevels = target_geographic_regions2
        targetlevels = ["FIN","CEU","GBR","IBS","TSI", 
                        "Turkish.Peninsula","Syrian.Desert",
                        "Central.Asia","Arabian.Peninsula",
                        "Northeast.Africa","Northwest.Africa"]

    elif annotcol == "Continent" :
        targetlevels = [ 'Africa',
                         'Middle.East',
                         'Central.Asia',
                         'Europe',
                         'East.Asia',
                         'America',
                         'Oceania' ]
    elif annotcol == "Continent2" :
        targetlevels = [ 'Africa',
                         'Middle.East',
                         'Central.Asia',
                         'Europe',
                         'East.Asia']
                         #'America',
                         #'Oceania' ]

    elif annotcol == "Country" :
        targetlevels = target_me_countries

    elif annotcol == "ethnicity" :
        targetlevels = target_ethnicities_all 

    if targetlevels is not None :
        targetlevels = [x.strip().replace(" ",".") for x in targetlevels]
        newfst = []
        for i in range(0,len(targetlevels)) :
            subset = targetlevels[:i+1]
            newfst.append(fstdata[(fstdata.Region1==targetlevels[i]) 
                                  & (fstdata.Region2.isin(subset))])
        newfst = concat( newfst )
        maxfst = newfst.Fst.max()
        midfst = maxfst / 2
        r_dataframe = com.convert_to_r_dataframe(newfst)
        r_dataframe = fixRLevels( r_dataframe, "Region1", targetlevels )
        r_dataframe = fixRLevels( r_dataframe, "Region2", targetlevels )
    else : 
        r_dataframe = com.convert_to_r_dataframe(fstdata)

    print "Tlevels:",targetlevels
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(Region2)",y="factor(Region1)") +
                ggplot2.geom_tile(ggplot2.aes_string(fill="unclass(Fst)")) +
                ggplot2.scale_fill_gradient2(name="Fst", low="blue",high="red",
                                            mid="yellow",midpoint=midfst) +
                ggplot2.theme(**{
                    'axis.text.x': ggplot2.element_text(angle=-90, size=15, hjust=0),
                    'axis.text.y': ggplot2.element_text(angle=0, size=15, hjust=1),
                    'legend.title' :ggplot2.element_text(size=15),
                    'axis.title.x': ggplot2.element_blank(),
                    'axis.title.y': ggplot2.element_blank()
                }) +
                ggplot2.coord_fixed() +
                ggplot2.theme(**mytheme) )
                #ggplot2.ggtitle("Fst") +
                #ggplot2.scale_x_discrete(position="top") +

    lofhetfile = "%s/%s_%s_heat.png"%(targetdir,basename,annotcol)
    print "Writing file %s" % lofhetfile
    grdevices.png(lofhetfile)
    p.plot()
    grdevices.dev_off()
# END makeFstPlot

def calcFst( eigenfiles, force=True, annotcol="Continent" ) :
    print "Running calcFst" 

    parfile,outliers,eigenvec,eigenval,regions = eigfiles
    filepath,basename,ext = hk.getBasename(parfile)
    targetdir = targetdir if targetdir is not None else filepath
    print "Using targetdir:",targetdir
    hk.makeDir(targetdir)
    #eigfiles = makeParFile(filepath,basename,
                             #targetdir=targetdir,
                             #annotcol=annotcol,
                             #fstonly=True ) 

    fstfile = os.path.join( targetdir, basename+".fst")
    print "Fst file:", fstfile
    outfile = os.path.join( targetdir, basename+".runout")

    if not os.path.exists(outliers) or force :
        print "Parfile:",parfile
        print "Running smartpca"
        subout = subprocess.check_output( ["smartpca", "-p",parfile] )
        ALLOUT = open(outfile,"w")
        FSTOUT = open(fstfile,"w")
        print >> FSTOUT, "Region1\tRegion2\tFst\tstderr"
        print >> ALLOUT, "Smartpca out"
        infst = False
        fstdata = []
        for line in subout.split("\n") :
            print >> ALLOUT, line
            if line.find("## Fst statistics") > -1 : infst = True
            elif len(line.strip()) == 0 : infst = False
            elif infst : print >> FSTOUT, "\t".join( line.strip().split() )
        ALLOUT.close()
        FSTOUT.close()
    
    assert os.path.exists(fstfile)
    fstdata = read_csv( fstfile, sep="\t" )
    if len(fstdata) > 0 :
        makeFstPlot( fstdata, annotcol, basename, targetdir )
    else :
        print "Error! - no fstdata?"

    return fstdata
# END calcFst

######################################################################
# plink_recode12
######################################################################
def plink_recode12( pedfile, force=forceFlag ) :
    print "Running plink_recode12"
    targetdir,filename,suffix = hk.getBasename(pedfile)

    basename = "%s/%s" % (targetdir, filename)
    recodefile = basename+".recode12"
    # --1
    if not os.path.exists( recodefile+".ped" ) or force:
        command = ["plink","--noweb","--tfile", basename, \
                "--maf","0.01","--recode12","--out", recodefile ] \
                #" --geno 0.1" \
                #% ( basename, recodefile )
        print " ".join(command)
        #out = hk.runCommand( command )
        out = subprocess.check_output( command )
        print out
    return recodefile+".ped"
# END plink_recode12

################################################################################
# plink_makefrqfile
################################################################################
def plink_makefrqfile( pedfile, force=forceFlag ) :
    print "Running plink_makefrqfile"
    filepath,basename,suffix = hk.getBasename(pedfile)
    print "filepath",filepath,"basename",basename,"suffix",suffix
    if suffix == ".tped" : filetype = "--tfile"
    elif suffix == ".bed" : filetype = "--bfile"
    elif suffix == ".ped" : filetype = "--file"
    else : print "Error! Unknown filetype", suffix; sys.exit(1)

    # Make frqfile
    freqfile = filepath+"/"+basename+".frq"
    tmpfile = filepath+"/"+basename+".1.frq"
    #print "Checking frq file:",freqfile, os.path.exists( freqfile )
    # --1
    if not os.path.exists( freqfile ) or force :
        command = "plink --noweb %s %s" \
                " --freq --out %s" \
                % ( filetype, os.path.join(filepath,basename), 
                   os.path.join(filepath,basename+".1") )
        # --filter-controls 
        #--recode12 
        print command
        out = hk.runCommand( command )
        print out
        command = "awk '$1 !~ /CHR/ {print $1,$2,$3,$4,$5,$6}' %s > %s" % (tmpfile, freqfile )
        print command
        out = hk.runCommand( command )
        print out
    return freqfile
# END plink_makefrqfile

def seriousClean( vcffile, rerun=False ): 
    print "Running seriousClean"
    patientdata = sampleAnnotation()
    filepath, basename, suffix = hk.getBasename(vcffile)
    targetdir = "%s/pca/" % filepath
    hk.makeDir(targetdir)
    cptargetfile = targetdir+basename+suffix+".gz"
    cleanped = "%s/%s"%(targetdir,basename)

    if not os.path.exists(cptargetfile) :
        out = hk.runCommand( "cp %s/%s%s* %s" % (filepath,basename,suffix, targetdir) )
        print "Copy targetfile:",cptargetfile

    print "Making sample missingness estimates", cleanped
    #command = ["vcftools","--gzvcf", cptargetfile,
               #"--missing-indv","--out",cleanped]
    #out = subprocess.check_output( command )

    command = ["vcftools","--gzvcf", cptargetfile,
               "--remove-filtered-all",
               "--remove-indels",
               "--maf",".001",
               "--min-alleles","2",
               "--max-alleles","2",
               "--plink-tped","--out",cleanped]

    print " ".join(command)
    tped = cleanped+".tped"
    if not os.path.exists( tped ) or rerun:
        out = subprocess.check_output( command )

    recodefile = plink_recode12( tped, force=rerun )
    # Modify ped
    modped = pop.modify_ped( recodefile, patientdata, calcdistances=True, shorten=True, force=rerun )
    #pop.modify_tped( tped, patientdata, force=rerun )
    #print modped
    frqfile = plink_makefrqfile( recodefile, force=rerun )
    newfrq, excludebase = pop.excludeSnps( recodefile, force=rerun )
    bedfile = pop.run_plink_convert(excludebase+".ped", force=rerun)
    #bedfile = pop.run_plink_convert(tped, force=rerun)
    return bedfile
    #print vcffile
# END seriousClean

def seriousClean2( vcffile, rerun=False ): 
    print "Running seriousClean2"
    patientdata = sampleAnnotation()
    filepath, basename, suffix = hk.getBasename(vcffile)
    targetdir = "%s/pca/" % filepath
    hk.makeDir(targetdir)
    cptargetfile = targetdir+basename+suffix+".gz"
    cleanped = "%s/%s"%(targetdir,basename)

    if not os.path.exists(cptargetfile) :
        out = hk.runCommand( "cp %s/%s%s* %s" % (filepath,basename,suffix, targetdir) )
        print "Copy targetfile:",cptargetfile

    print "Making sample missingness estimates", cleanped
    #command = ["vcftools","--gzvcf", cptargetfile,
               #"--missing-indv","--out",cleanped]
    #out = subprocess.check_output( command )

    command = ["vcftools","--gzvcf", cptargetfile,
               "--remove-filtered-all",
               "--remove-indels",
               "--maf",".001",
               "--min-alleles","2",
               "--max-alleles","2",
               "--tped-plink","--out",cleanped]

    print " ".join(command)
    ped = cleanped+".ped"
    #bedfile = pop.run_plink_convert(tped, force=rerun)
    if not os.path.exists( ped ) or rerun:
        tped = cleanped+".tped"
        out = subprocess.check_output( command )
        subprocess.call( ["plink","--noweb","--recode","--tfile",cleanped,"--out",cleanped] )

    #recodefile = plink_recode12( ped, force=rerun )
    # Modify ped
    modped = pop.modify_ped( ped, patientdata, calcdistances=True, shorten=True, force=rerun )
    #pop.modify_tped( tped, patientdata, force=rerun )
    #print modped
    frqfile = plink_makefrqfile( ped, force=rerun )
    newfrq, excludebase = pop.excludeSnps( ped, force=rerun )
    bedfile = pop.run_plink_convert(excludebase+".ped", force=rerun)
    return bedfile
    #print vcffile
# END seriousClean2


def plotPCA( tdf, basename, pc1, pc2, smallplot=False ):
    tregions = tdf.Region.unique()
    r_dataframe = com.convert_to_r_dataframe(tdf)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x=pc1, y=pc2) +
                ggplot2.geom_point(ggplot2.aes_string(
                    shape="factor(Region)",colour="factor(Region)"), size=1) +
                #ggplot2.ggtitle("AF comparison between CEU and ME") +
                #ggplot2.stat_smooth(method="lm", se=False)+
                ggplot2.scale_y_continuous(pc2)+
                ggplot2.scale_x_continuous(pc1)+
                ggplot2.scale_shape_manual("Region",values=robjects.IntVector(range(1,len(tregions)+1)))+
                ggplot2.scale_colour_brewer("Region",palette="Set1") +
                ggplot2.theme(**pcatheme) )

    psmall = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x=pc1, y=pc2) +
                ggplot2.geom_point(ggplot2.aes_string(
                    shape="factor(Region)",colour="factor(Region)"), size=.8) +
                #ggplot2.ggtitle("AF comparison between CEU and ME") +
                #ggplot2.stat_smooth(method="lm", se=False)+
                ggplot2.scale_y_continuous(pc2)+
                ggplot2.scale_x_continuous(pc1)+
                ggplot2.scale_shape_manual("Region",values=robjects.IntVector(range(1,len(tregions)+1)))+
                ggplot2.scale_colour_brewer("Region",palette="Set1") +
                ggplot2.theme(**pcathemesmall) )

                #ggplot2.scale_shape_discrete(solid=False)+

    #psmall = p + ggplot2.theme(**pcathemesmall)
    #plarge = p + ggplot2.theme(**pcatheme)
    #if ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})
    figname = "results/figures/pca/%s_%s%s_pcs.pdf" % (basename, pc1, pc2)
    print "Writing file:",figname
    grdevices.pdf(figname, width=5, height=4)
    p.plot()
    grdevices.dev_off()
    #p = p  + ggplot2.theme(**{'legend.position':"none"}) 
    return psmall
# ENd plotPCA

def plotPCAGrid( eigfiles ) :
    parfile,outliers,eigenvec,eigenval,regions = eigfiles

    targetdir,basename,ext = hk.getBasename(parfile)

    tregions = getTargetRegions( regions, annot )

    evec = read_csv(eigenvec, delim_whitespace=True, header=None, skiprows=[0])
    evec.index = evec[0].tolist()
    del evec[0]
    evec.columns = ["V%d"%x for x in range(1,21)] + ["Region"]

    evecfilt = evec[evec.Region.isin(tregions)]

    #print evec.head()
    plotlist = []
    ng = grid.grid.nullGrob()
    legend = None
    for pc2 in [2,3,4] :
        #plotlist.append([])
        for pc1 in [1,2,3] :
            if pc1 >= pc2 : plotlist.append(ng); continue
            p = plotPCA( evecfilt, basename, "V%d"%pc1, "V%d"%pc2 )
            legend = gtable.gtable_filter(
                ggplot2.ggplot2.ggplot_gtable(
                    ggplot2.ggplot2.ggplot_build(p)),"guide-box")
            plotlist.append(p+ ggplot2.theme(**{'legend.position':"none"}) )

    print plotlist
    figname = "results/figures/pca/%s_grid_pcs.pdf" % (basename)
    print "Making figure:",figname
    ph = 5
    grdevices.pdf(figname, width=ph+1.5, height=ph)#res=300)
    #gridextra.grid_arrange( *plotlist , ncol=3 )
    gridextra.grid_arrange( gridextra.arrangeGrob(*plotlist , ncol=3), legend, nrow=1,
                                widths=grid.grid.unit_c( grid.unit(ph, "inches"), 
                                            grid.unit(1.5,"inches")))

    grdevices.dev_off()
# END plotPCAGrid

######################################################################
if __name__ == "__main__":

    #optlist, args = getopt.getopt( sys.argv[1:], "bk")
    #optlist = dict(optlist)
    #bedfile = optlist.get("-b",bedfile)

    # Change running directory
    os.chdir("..")

    #bedfile = "./rawdata/variome1/variome.chimp.regions.filt.samp.plink.mod.filt.uniq.bed"
    #bedfile = "./rawdata/onekg/onekg.chimp.regions.filt.samp.plink.filt.bed"
    #bedfile = "./rawdata/test/everything_set1.chr1.snp.clean.plink.filt.bed"
    #bedfile = "./rawdata/mevariome/main/ibc/variome.clean.Middle_East.recode12.mod.bed"
    #bedfile = "./rawdata/mevariome/main/variome.clean.plink.filt.bed"
    #bedfile = "./rawdata/mevariome/main/variome.clean.plink.filt.bed"
    #vcffile = "./rawdata/merge1kg/main/me1000G.clean.vcf.gz"
    #vcffile = "./rawdata/merge1kg/main/me1000G.clean.vcf.gz"
    vcffile = "./rawdata/mergedaly/meceu.X.vcf.gz"
    vcffile = "./rawdata/merge1kg/me1000G.X.vcf.gz"
    vcffile = "./rawdata/mevariome/variome.X.vcf.gz"
    #vcffile = "./rawdata/mevariome/main/variome.clean.vcf.gz"
    vcffile = "./rawdata/merge1kg/main/me1000G.clean.vcf.gz"
    #vcffile = "./rawdata/onekg/main/onekg.clean.vcf.gz"
    #vcffile = "./rawdata/test2/main/test2.clean.vcf.gz"

    #bedfile = seriousClean2( vcffile, True )
    bedfile = seriousClean( vcffile, False )

    print "Bedfile:", bedfile
    assert os.path.exists(bedfile)

    filepath, filename = os.path.split(bedfile)
    targetdir = os.path.join(filepath,"pca")

    clustfile = "./rawdata/merge1kg/main/clust/me1000G.clean.annot"
    filepats = patientInfo.currentFilePats( vcffile )
    if os.path.exists(clustfile) :
        sampleannot = read_csv( clustfile, sep="\t" )
        sampleannot.index = sampleannot["Individual.ID"]
        sampleannot = sampleannot.ix[filepats,:]
        mapping_regions = {
                   #"NWA":"Middle East", "NEA":"Middle East", "AP":"Middle East",
                   #"SD":"Middle East", "TP":"Middle East", "CA":"Middle East",
                   "LWK":"Africa", "YRI":"Africa", "IBS":"Europe", "CEU":"Europe",
                   "TSI":"Europe", "FIN":"Europe", "GBR":"Europe",
                   "Europe":"Unknown", "Africa":"Unknown", "East Asia":"East Asia",
                   "CHB":"East Asia", "CHS":"East Asia", "JPT":"East Asia"}

        sampleannot["Continent2"] = [mapping_regions[x] if mapping_regions.has_key(x)
                                      else x
                                      for x in sampleannot.GeographicRegions3]
    else :
        sampleannot = sampleAnnotation( filepats )

    annot = "Continent2"
    sampleannot["Region"] = sampleannot[annot]
    targetdir = os.path.join(filepath,"pca/"+annot)
    #fstdata = calcFst( bedfile, targetdir, force=False, annotcol=annot )
    #fixedannotname = ("GeographicRegions2" if annot == "GRsub" else annot)

    eigfiles = makeParFile(bedfile, sampleannot,targetdir=targetdir)
                           #annotcol=fixedannotname) 
    plotPCAGrid( eigfiles )
    
    #outliers = pcaOutlierAnalysis( eigfiles, force=True, 
                                  #annotcol=annot, tcols="1:2" )
    #print outliers

    

    #gtable_filter(ggplot_gtable(ggplot_build(p1)), "guide-box")
    #grid.arrange(plots[[1]], ng,         ng,
                 #plots[[2]], plots[[3]], ng,
                 #plots[[4]], plots[[5]], plots[[6]])

    #annotationcolumns = ["GeographicRegions","GeographicRegions2","Country","Continent2","ethnicity","Source"]
    #annotationcolumns = ["GeographicRegions2","GRsub"]
    #annotationcolumns = ["GeographicRegions","GeographicRegions2"]
        #for annot in annotationcolumns :


    #outliers = pcaOutlierAnalysis( bedfile, targetdir, force=False, annotcol="GeographicRegions" )
    #outliers = pcaOutlierAnalysis( bedfile, targetdir, force=False, annotcol="GeographicRegions" )
    #outliers = pcaOutlierAnalysis( bedfile, targetdir, force=True, annotcol="ethnicity" )
# END Main
     #cmd = ('unit.c(unit(2, "lines"), '+
                        #'unit(1, "npc") - unit(2, "lines") - %s["width"], '+
                        #'%s["width"]') % (legend.r_repr(), legend.r_repr() )
#
    #cmd = u"gtable_filter(ggplot_gtable(ggplot_build(%s)))" % p.r_repr()
#
    #cmd = ('unit.c(unit(2, "lines"), unit(1, "npc") - unit(3, "lines"), unit(3, "lines")')
    #lwidths = robjects.r( cmd )
    #lwidths = robjects.r( '%s["width"]'% (legend.r_repr()) )
#
                 
