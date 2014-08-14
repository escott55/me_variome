#!/usr/bin/python

import os, sys
from pandas import *
from collections import Counter
import subprocess

from localglobals import *
import popgencommands as pop
forceFlag = False

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
    #command = "bgzip %s; tabix -p vcf %s" % (tfile,outfile)
    command = ["bgzip",tfile]
    subprocess.call(command)
    command = ["tabix","-p","vcf",outfile]
    subprocess.call(command)
    return outfile
# END compress

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
def calculateRVIS( xdata, ydata, prefix=None ):

    xdata_uniq = xdata[["chrom","pos","gene"]].drop_duplicates()
    Xgenecounts = Counter(xdata_uniq["gene"].tolist())
    #ydata["newAF"] = [ydata.x/totalsamp for 
    #calculateAF( ydata )
    ydata_uniq = ydata[["chrom","pos","gene"]].drop_duplicates()
    Ygenecounts = Counter(ydata_uniq["gene"].tolist())
    ycounts = DataFrame([[x,Ygenecounts[x]] for x in Ygenecounts],columns=["gene","Ycount"])
    xcounts = DataFrame([[x,Xgenecounts[x]] for x in Xgenecounts],columns=["gene","Xcount"])
    allcounts = merge(ycounts,xcounts,on="gene")
    print allcounts.head(3)
    #p = ggplot2.ggplot(faithful_data) + \
                #ggplot2.aes_string(x = "eruptions") + \
                #ggplot2.geom_histogram(fill = "lightblue") + \
                #ggplot2.geom_density(ggplot2.aes_string(y = '..count..'), colour = "orange") + \
                #ggplot2.geom_rug() + \
                #ggplot2.scale_x_continuous("Eruption duration (seconds)") + \
                #ggplot2.opts(title = "Old Faithful eruptions")
    #p.plot()

    #r_dataframe = com.convert_to_r_dataframe(allcounts)

    #print allcounts["Ycount"]
    robjects.globalenv["Ycount"] = robjects.FloatVector(allcounts["Ycount"])
    robjects.globalenv["Xcount"] = robjects.FloatVector(allcounts["Xcount"])
    model_x = rstats.lm("Ycount ~ Xcount")
    resid = list(rstats.rstudent(model_x))
    print resid
    allcounts["StdResid"] = resid
    allcounts["Colour"] = "Innerquartile"
    lowquantile = allcounts.StdResid.quantile(.02)
    highquantile = allcounts.StdResid.quantile(.98)
    #allcounts["Colour"][allcounts.StdResid <= lowquantile] = "2% most intolerant"
    #allcounts["Colour"][allcounts.StdResid >= highquantile] = "2% least intolerant"
    print "Allcounts head:",allcounts.head(20)
    print "Allcounts shape:",allcounts.shape
    print allcounts[allcounts.Colour != "Innerquartile"].head(20)

    r_dataframe = com.convert_to_r_dataframe(allcounts)
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
    figname = "results/figures/%s_rvis_lm.png" % prefix
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
    return allcounts
# END calculateRVIS

######################################################################
# makeROHfile
######################################################################
def makeROHfile( vcffile, sampleannot, outprefix, famfile=None, rerun=False ) :
    print "FEstim Targetfile:",vcffile
    assert os.path.exists(vcffile)

    tpedfile = pop.run_vcftools_plink( vcffile, force=rerun )
    print "Pedfile :",tpedfile

    # Modify Ped File annotation
    pop.modify_tped( tpedfile, sampleannot, force=rerun )

    # LD vars, MAF >=.05, Missingness >=.05
    bedfile_filt = pop.run_plink_filter( tpedfile, force=rerun )

    roh_calls = pop.run_plink_homozyg( bedfile_filt, force=rerun )
    return roh_calls
# END makeROHfile
    #filepath, basename, suffix = hk.getBasename( vcffile )
    #targetdir = "%s/roh/" % filepath
    #makeDir(targetdir)
    #if famfile is not None and os.path.exists(famfile ):
        #famdata = read_csv(famfile, delim_whitespace=True, header=None, names=["family","Individual.ID","Paternal","Maternal","Gender","Aff"] )
        #keeplist = "%s/%s_keeplist.txt" % (targetdir,basename)
        #print "Writing Keep file:",keeplist
        #famdata[["Individual.ID"]].to_csv(keeplist,header=False,index=False)

# Goal 1 plot Individual Num segments
################################################################################
def makeRohBoxPlots( rohfile, sampleannot, outprefix="./results/figures/rohbox" ) :
    indivstats = read_csv(rohfile+".indiv",delim_whitespace=True)
    print "indivstats:",indivstats.shape
    #print allintervals.head(10)

    merged = merge( sampleannot[["Individual.ID","Continent2","ethnicity"]], indivstats[["IID","KB","KBAVG","NSEG"]], right_on="IID", left_on="Individual.ID", how="right" )
    merged = merged[merged.Continent2.notnull()]

    levels = merged[["Continent2","ethnicity"]].drop_duplicates().sort('Continent2')
    levels = levels.ethnicity.unique().tolist()
    print "Levels",levels

    #print "Merged:",merged.shape
    #print "sampleannot:",sampleannot.shape
    print merged.columns

    r_dataframe = com.convert_to_r_dataframe(merged)
    s = robjects.FactorVector(r_dataframe.rx2("ethnicity"), levels=robjects.StrVector(levels))
    new_r_df = r_dataframe.cbind(s)
    #new_r_df.colnames = robjects.StrVector(["Continent2",feature,"ethnicity","Ethnicity"])
    new_r_df.colnames = robjects.StrVector(['Individual.ID', 'Continent2', 'ethnicity', 'IID', 'KB', 'KBAVG', 'NSEG','Ethnicity'])
    #print robjects.r.head(new_r_df)

    p = ggplot2.ggplot(new_r_df) + \
                ggplot2.aes_string(x = "Ethnicity",y="KB" ) + \
                ggplot2.geom_boxplot(ggplot2.aes_string(fill="factor(Continent2)")) + \
                ggplot2.ggtitle("Comparison of ROH") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                ggplot2.scale_y_continuous("Total ROH")# + \
                #ggplot2.scale_x_continuous("External RVIS")+ \
                #ggplot2.stat_smooth(method="lm", se=False)

    print "Writing file %s_eth_kb.png" % outprefix
    grdevices.png("%s_eth_kb.png" % outprefix)
    p.plot()
    grdevices.dev_off()


    p = ggplot2.ggplot(new_r_df) + \
                ggplot2.aes_string(x = "Ethnicity",y="KBAVG" ) + \
                ggplot2.geom_boxplot(ggplot2.aes_string(fill="factor(Continent2)")) + \
                ggplot2.ggtitle("ROH KBAVG") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                ggplot2.scale_y_continuous("Average ROH")# + \

    print "Writing file %s_eth_kbave.png" %outprefix
    grdevices.png("%s_eth_kbave.png" %outprefix)
    p.plot()
    grdevices.dev_off()

    p = ggplot2.ggplot(new_r_df) + \
                ggplot2.aes_string(x = "Ethnicity",y="NSEG" ) + \
                ggplot2.geom_boxplot(ggplot2.aes_string(fill="factor(Continent2)")) + \
                ggplot2.ggtitle("Number of Segments") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                ggplot2.scale_y_continuous("NSEG")# + \

    print "Writing file %s_eth_nseg.png" %outprefix
    grdevices.png("%s_eth_nseg.png" %outprefix)
    p.plot()
    grdevices.dev_off()
# END makeRohBoxPlots

################################################################################
################################################################################
def makeRohCumPlot( rohfile, sampleannot, outprefix="results/rohcumplot" ) :
    allintervals = read_csv(rohfile,delim_whitespace=True)
    #indivstats = read_csv(rohfile+".indiv",delim_whitespace=True)
    #print "indivstats:",indivstats.shape

    intdf = merge( sampleannot[["Individual.ID","Continent2","ethnicity"]], allintervals[["IID","KB","NSNP","DENSITY"]], right_on="IID", left_on="Individual.ID", how="right" )
    intdf = intdf[intdf.Continent2.notnull()]

    allsizedist = []
    binstart = 1000
    bins = np.linspace(binstart, 15000, 100)
    binsize = bins[1] - bins[0]
    binskb = [binstart+binsize*i for i in range(len(bins)) ]
    for region, rints in intdf.groupby("Continent2") :
        groups = rints.groupby(np.digitize(rints.KB, bins))
        sizedist = DataFrame([ [group,len(data)] for group,data in groups ], columns=["Bin","Freq"])
        allsizes = sizedist.Freq.tolist()
        lentotal = float(sum( allsizes ))
        sizedist["CumProp"] = [sum(allsizes[:i])/lentotal for i in range(1,len(allsizes)+1) ]
        sizedist["RohSize"] = [binskb[i-1] for i in sizedist.Bin.tolist()]
        sizedist["Region"] = region
        allsizedist.append(sizedist)

    allsizedist = concat( allsizedist ).reset_index(drop=True)
    #print allsizedist[allsizedist.Region == "Middle East"].tail()
    r_dataframe = com.convert_to_r_dataframe(allsizedist)

    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x = "RohSize",y="CumProp", colour="factor(Region)" ) + \
                ggplot2.geom_point() + ggplot2.geom_line() + \
                ggplot2.ggtitle("Cumulative Proportion of ROH by Length") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                ggplot2.scale_x_continuous("ROH Length")+ \
                ggplot2.scale_y_continuous("Cumulative Proportion") #+ \
                #ggplot2.scale_x_log10("ROH Length")+ \

    print "Writing file %s_rohcumul.png" % outprefix
    grdevices.png("%s_rohcumul.png"%outprefix)
    p.plot()
    grdevices.dev_off()
# END makeRohCumPlot

################################################################################
################################################################################
def run_vcftools_custom( vcffile, keepfile, outprefix, force=False ) :
    assert os.path.exists(vcffile)
    
    if os.path.exists( outprefix+".clean.vcf.gz" ) and not force :
        print " - File"+outprefix+".clean.vcf.gz : already exists. Skipping command!"
        return outprefix+".clean.vcf.gz"
    
    # remove-indels
    # mix/max-alleles - Preven n allelic sites
    # snps <filename> include snps based on list of ids
    # Geno - removes genotypes without sufficient coverage
    command = ["vcftools","--gzvcf",vcffile, 
               "--min-alleles","2",
               "--max-alleles","2", 
               "--max-missing","0.95",
               "--remove-filtered-all",
               "--recode",
               "--out",outprefix]
        #" --min-meanDP 4 --minQ 20" \
    
    print " ".join(command)
    out = subprocess.check_output( command )
        
    assert os.path.exists( outprefix+".recode.vcf" ), out
    subprocess.call( ["mv","%s.recode.vcf"%outprefix, "%s.clean.vcf"%outprefix] )
    cfile = pop.compress( outprefix+".clean.vcf" )
    print " - End run_vcftools"
    return cfile
# END run_vcftools_custom

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
        subprocess.call( ["sort","-k1,1","-k2,2n",bedfile], stdout=outfile )
    return sortedfile
# END sortBed

################################################################################
# exomeCoverage
################################################################################
def exomeCoverage( rohfile, sampleannot, outprefix="results/rohcov", figureprefix=None ):
    print "Running exomeCoverage"
    rohdf = read_csv(rohfile,delim_whitespace=True)
    print "rohdf:",rohdf.shape
    rohannot = merge( sampleannot[["Individual.ID","Continent2","ethnicity"]], 
                   rohdf[["IID","CHR","POS1","POS2","KB"]], right_on="IID", 
                   left_on="Individual.ID", how="right" )
    rohannot = rohannot[(rohannot.Continent2.notnull()) & (rohannot.Continent2 != "Unknown")]

    rohbedfile = "%s.bed" % outprefix
    print "Writing file:",rohbedfile
    rohannot.to_csv(rohbedfile, cols=["CHR","POS1","POS2","IID"], sep="\t",index=False,header=False)
    rohbedsort = sortBed( rohbedfile )

    rohbases = calculateRegionsCovered( rohbedsort ) 

    exomeregions = "resources/regionbeds/refseqgenes_b37_sorted.bed"  
    exometotals = calculateRegionsCovered( exomeregions ) 
    print exometotals

    #sort -k1,1 -k2,2n refseqgenes_b37.bed
    #sort -k1,1 -k2,2n in.bed
    #
    allbases = {}
    limit = 10
    for samp, data in rohannot.groupby( "IID" ) : 
        print "Processing:",samp
        tmprohbedfile = "%s/%s_sort.bed" % (outprefix,samp)
        tmproh = data.sort(["CHR","POS1","POS2"])
        tmproh.to_csv(tmprohbedfile, cols=["CHR","POS1","POS2","IID"], sep="\t",index=False,header=False)
        command = ["bedtools","intersect","-a",tmprohbedfile,"-b",exomeregions,"-sorted"]
        #print " ".join(command)
        stdout = subprocess.check_output( command )
        intersectdata = (DataFrame( [x.split("\t") for x in stdout.split("\n")], 
                                  columns=["chrom","start","end","Individual.ID"] )
                         .drop_duplicates())
        intersectbed = "%s/%s_exomeroh.bed" % (outprefix, samp)
        #print "Writing file:",intersectbed
        intersectdata.to_csv( intersectbed, header=False, index=False, sep="\t" )
        totalbases = calculateRegionsCovered( intersectbed ) 
        allbases[samp] = totalbases
        limit -= 1
        print limit
        if limit == 0 : break

    allbases = DataFrame( allbases ).transpose()

    overlapfile = "%s/chrom_overlap.txt" % (outprefix )
    print "Writing file:",overlapfile
    allbases.to_csv( overlapfile, sep="\t", index=False)
    percentAnalysis( allbases, figureprefix )
# END exomeCoverage

def percentAnalysis( allbases, figureprefix=None ) :
    percents = melt(DataFrame({col:(allbases[int(col)].fillna(0)/float(exometotals[col])*100) 
                          for col in exometotals.index if col in allbases.columns}))
    #allbases.divide( exometotals.map(float), fill_value=0 )
    print percents.head(10)
    r_dataframe = com.convert_to_r_dataframe(percents)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x = "factor(variable)",y="value" ) + \
                ggplot2.geom_boxplot() + \
                ggplot2.ggtitle("Chromosomal distribution of ROH") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.scale_y_continuous("RoH (%)" )+ \
                ggplot2.scale_x_discrete("Chromosomes" )# + \

    if figureprefix is None : figureprefix = "test"
    if figureprefix.find(" ")>=0 : figureprefix = figureprefix.replace(" ","_")
    figname = "results/figures/roh/%s_roh_chr.png" % figureprefix
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()

    sums = allbases.head(10).fillna(0).apply(sum,axis=1) 
    percents = (sums / float(exometotals.sum()) *100).reset_index()
    percents.columns = ["Individual.ID","percent"]
    percents.sort("percent", inplace=True)
    levels = percents["Individual.ID"].unique().tolist()
    print levels
    print percents.head(10)
    r_dataframe = com.convert_to_r_dataframe(percents)
    s = robjects.FactorVector(r_dataframe.rx2("Individual.ID"), levels=robjects.StrVector(levels))
    new_r_df = r_dataframe.cbind(s)
    new_r_df.colnames = robjects.StrVector(['Individual.ID', 'percent', 'Sample'])
    p = ggplot2.ggplot(new_r_df) + \
                ggplot2.aes_string(x = "Sample",y="percent" ) + \
                ggplot2.geom_bar(stat="identity") + \
                ggplot2.ggtitle("Exome coverage") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                ggplot2.scale_y_continuous("Exome coverage (%)" )+ \
                ggplot2.scale_x_discrete("Samples" )# + \

    figname = "results/figures/roh/%s_roh_samp.png" % figureprefix
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END percentAnalysis

################################################################################
# runRohAnalysis
################################################################################
def runRohAnalysis( vcffile, sampleannot, keepfile=None ) :

    filepath, filename, suffix = hk.getBasename(vcffile)
    hk.makeDir("%s/roh/" % filepath )
    targetdir = "%s/roh/" % filepath
    outprefix = "%s/roh/%s" % (filepath, filename)
    figprefix = "./results/figures/roh/%s" % (filename)

    if keepfile is not None : newvcf = pop.run_vcftools_clean( vcffile, keepfile, outprefix )
    else :
        print "Copy vcffile:",vcffile
        newvcf = targetdir+filename+suffix+".gz"
        if not os.path.exists(newvcf) :
            command = ["cp",vcffile,targetdir]
            print " ".join(command)
            out = subprocess.check_output(command)
            print out

    roh_calls = makeROHfile( newvcf, sampleannot, outprefix )

    print "Roh calls:",roh_calls

    makeRohBoxPlots( roh_calls, sampleannot, figprefix )
    makeRohCumPlot( roh_calls, sampleannot, figprefix ) 
    return roh_calls
# END runRohAnalysis


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
        print " -t run a test"
        sys.exit(1)

    hk.changeLogFile( hk.LOGDIR+"/roh_analysis_test.log", True )

    print "Using dataset:", dataset
    path = os.path.abspath("./rawdata/")
    if dataset == "test" :
        vcffile = path+"/test/everything_set1.chr1.snp.vcf.gz"
    elif dataset == "onekg" :
        vcffile = path+"/onekg/onekg.vcf.gz"
    elif dataset == "daily" :
        vcffile = path+"/daily/daily.vcf.gz"
    elif dataset == "eichler" :
        vcffile = path+"/eichler/eichler.vcf.gz"
    #elif dataset == "variome" :
        #vcffile = path+"/variome/variome.clean.vcf.gz"
    elif dataset == "merge1kg" :
        vcffile = path+"/merge1kg/main/me1000G.clean.vcf.gz"
    elif dataset == "mevariome" :
        vcffile = path+"/mevariome/main/variome.clean.vcf.gz"
    #elif dataset == "variome1" :
        #vcffile = path+"/variome1/variome.clean.vcf.gz"
    #elif dataset == "hgdp" :
        #vcffile = path+"/hgdp/HGDP_938.vcf.gz"
    #elif dataset == "merged" :
        #vcffile = path+"/merged/merged.vcf.gz"
    elif dataset == "casanova" :
        vcffile = path+"/casanova/casanova.snp.recal.vcf.gz"
    else :
        print "Error: Unknown dataset provided"
        sys.exit(1)

    #sampleannot = read_csv("resources/annotation/patientannotation.ped", sep="\t")
    #hgdp = read_csv("resources/HGDP_ethnicgroups.txt",sep="\t")
    #sampleannot = merge(sampleannot, hgdp[["Ethnicity","Continent","Country"]], left_on="ethnicity",right_on="Ethnicity",how="left")
    filepats = patientInfo.currentFilePats( vcffile )
    sampleannot = sampleAnnotation(filepats)

    if optlist.has_key("-t") :
        #rohfile = "rawdata/test/roh/everything_set1.chr1.snp.clean.plink.filt.hom"
        rohfile = "rawdata/daily/roh/daily.clean.plink.filt.hom"
        #makeRohCumPlot( rohfile, sampleannot )
        #makeRohBoxPlots( rohfile, sampleannot )
        
        exomeCoverage( rohfile, sampleannot )
    else :
        #keepfile = "%s/%s/tokeep.txt" % (path, dataset)
        #assert os.path.exists( keepfile )
        calls = runRohAnalysis( vcffile, sampleannot, None )
        print calls

# END MAIN

