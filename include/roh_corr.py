#!/usr/bin/python

import os, sys
from localglobals import *
import glob
from roh_figures import *

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

        if not breakpoints.has_key( chrom ) : breakpoints[chrom] = 0
        breakpoints[chrom] += abs(end-start)
    return Series(breakpoints)
# END calculateRegionsCovered

################################################################################
# makeIntersectBed
################################################################################
def makeIntersectBed2( sample, bedfile, targetdir, exomeregions, force=False ):
    print "Running makeIntersectBed2:",sample
    intersectbed = os.path.join(targetdir,"%s_exomeroh.bed"%(sample))
    print "Intersectbed:",intersectbed
    if os.path.exists(intersectbed) and not force:
        print "Warning: intersectbed already exists. Skipping..."
        return intersectbed

    data = read_csv( bedfile, sep="\t", header=None, 
                    names=["CHR","POS1","POS2","att","vars"])
                    #dtype={'POS1':int,'POS2':int} )

    data.CHR[ data.CHR == "X"] = 23
    data.CHR[ data.CHR == "Y"] = 24
    data[["CHR","POS1","POS2"]] = data[["CHR","POS1","POS2"]].astype(int)
    tmprohbedfile = "%s/%s_sort.bed" % (targetdir,sample)
    tmproh = data.sort(["CHR","POS1","POS2"])
    tmproh.to_csv(tmprohbedfile, cols=["CHR","POS1","POS2"], 
                  sep="\t",index=False,header=False)
    command = ["bedtools","intersect","-a",tmprohbedfile,"-b",exomeregions]#,"-sorted"
    print " ".join(command)
    try :
        stdout = subprocess.check_output( command )
    except subprocess.CalledProcessError as e:
        print "Caught error:",e
        sys.exit(1)

    if len(stdout.strip()) == 0 :
        print "Error: stdout - ",stdout
        print "Command -"," ".join(command)
        return "Failed"
    
    intersectdata = (DataFrame( [x.split("\t") for x in stdout.split("\n")],
                              columns=["chrom","start","end"] )#,"Individual.ID"
                     .drop_duplicates())
    intersectdata["Individual.ID"] = sample
    print "Writing file:",intersectbed
    intersectdata.to_csv( intersectbed, header=False, index=False, sep="\t" )
    return intersectbed
# END makeIntersectBed2

################################################################################
# exomeCov
################################################################################
def exomeCov( bedpath, targetdir, sampleannot, force=False, targetsamples=None, 
             outprefix="h3m2") :
    print "Running exomeCov:",bedpath
    rawrohbedfiles = glob.glob( os.path.join( bedpath, "*.bed" ) )
    overlapfile = "results/roh_h3m2/chrom_overlap.txt"# % (outprefix )

    # Read Annotation
    exomeregions = "/home/escott/workspace/variome/resources/regionbeds/refseq_condensed.bed"
    exometotals = calculateRegionsCovered( exomeregions )
    print exometotals

    if os.path.exists(overlapfile) and not force :
        print "Final file already exists:",overlapfile
        #rawrohbedfiles = retrieveBedFiles(outprefix)
        allbases = read_csv( overlapfile, sep="\t" )
        return allbases, exometotals

    allbases = {}
    bedfiles = {}
    limit = -1
    for bedfile in rawrohbedfiles : 
        filepath, filename = os.path.split(bedfile)
        sample = filename[:filename.find("_1e")]
        #print sample
        if targetsamples is not None :
            if sample not in targetsamples: 
                print "Error: skipping sample:",sample
                continue
        print "Sample passed!",sample
        intersectbed = makeIntersectBed2( sample, bedfile, targetdir, exomeregions, force )
        if intersectbed == "Failed" : print "Error!: intersectbed failed!", sample

        #print "Made intersectbed -", intersectbed
        totalbases = calculateRegionsCovered( intersectbed )
        #print "calc:",totalbases
        bedfiles[sample] = intersectbed
        allbases[sample] = totalbases.ix[range(1,23)]
        limit -= 1
        if limit == 0 : break

    print "Finished!:"
    #print "Allbases - ",allbases
    allbases = DataFrame( allbases ).transpose().astype(int).reset_index()
    allbases.rename(columns={'index':'Individual.ID'}, inplace=True)
    #print "After renaming head!:"
    #print allbases.head(10)[range(10,15)]

    print "Writing file:",overlapfile
    foundcols = [x for x in exometotals.index.tolist() if x in allbases.columns]
    print "Found columns!:",foundcols
    samplesums = allbases[foundcols].fillna(0).apply(sum,axis=1)

    print samplesums.head(10)[foundcols]
    exomesum = float(exometotals[foundcols].sum())
    print "Exomesum:",exomesum
    allbases["percent"] = (samplesums / exomesum *100)
    print allbases.head(10)

    allbases.to_csv( overlapfile, sep="\t", index=False)
    return allbases, exometotals
# END exomeCov

def collapseBedFiles( bedpath, targetdir, targetsamples=None, outprefix="h3m2", 
                     force=False ): 
    print "Running collapseBedFiles"
    targetbedfile = os.path.join( targetdir, outprefix+"_all.bed" )
    targetfile = os.path.join( targetdir, outprefix+"_all.hom" )

    if os.path.exists(targetfile) and not force:
        print "Warning: final file already exists:",targetfile,". skipping..."
        return targetbedfile, targetfile

    print "Making file:",targetfile
    rawrohbedfiles = glob.glob( os.path.join( bedpath, "*.bed" ) )
    alldata = []
    for bedfile in rawrohbedfiles : 
        filepath, filename = os.path.split(bedfile)
        sample = filename[:filename.find("_1e")]
        #print sample
        if targetsamples is not None :
            if sample not in targetsamples: 
                print "Error: skipping sample:",sample
                continue
        print "Sample passed!",sample
        regions = read_csv( bedfile, sep="\t", header=None, 
                           names=["chrom","POS1","POS2","cnt","NSNP"])
        regions["IID"] = sample
        regions["KB"] = [(end-start)/1000. 
                         for start,end in regions[["POS1","POS2"]].values]
        regions["DENSITY"] = regions["NSNP"] / regions["KB"]
        alldata.append( regions )

    alldata = concat( alldata )
    alldata[["chrom","POS1","POS2","IID"]].to_csv(targetbedfile, header=None, sep="\t" )
    alldata.to_csv(targetfile, sep="\t" )
    return targetbedfile,targetfile
# END collapseBedFiles

################################################################################
# Read all h3m2 files
# Calculate exome coverage for h3m2 files
# calculate overlap with other samples
# plot percentage overlap
# plot exome coverage corr
################################################################################
if __name__ == "__main__" :
    os.chdir("..")
    basedir = "/home/escott/workspace/variome/results/roh_h3m2/"
    bedpath = os.path.join(basedir,"raw/")
    exomepath = os.path.join(basedir,"exome/")
    figuredir = "results/figures/roh/h3m2/"
    
    outprefix = "h3m2"
    prefix = "test"
    prefix = "meceu"
    targetsamples= [x.strip() 
                    for x in open("./rawdata/mergedaly/qc/meceu.finalpats").readlines() 
                    if len(x) > 0]
    sampleannot = sampleAnnotation(targetsamples)
    allbases, exometotals = exomeCov( bedpath, exomepath, sampleannot, force=False, 
                                     targetsamples=targetsamples, outprefix=outprefix )

    percentDensity( allbases, sampleannot, outprefix )
    #percentAnalysis( allbases, exometotals, outprefix )
    percentAnalysisRegions( allbases, exometotals, sampleannot, outprefix )

    rohbedfile, allrohfile = collapseBedFiles(bedpath, basedir, targetsamples=targetsamples, 
                                  outprefix=outprefix, force=False)
    makeRohCumsumPlot( allrohfile, sampleannot, os.path.join(figuredir,prefix) )

    sys.exit(1)
    print allbases.head(10)

    #prefix = "variome1"
    otherrohfile = "./results/rohcov/%s_chrom_overlap.txt" % prefix
    print otherrohfile
    allbases2 = read_csv(otherrohfile, sep="\t" ) 
    print "Data1:",
    print allbases.head(10)
    print "Data2:"
    print allbases2.head(10)

    print allbases.shape
    allbases["Algorithm"] = "H3M2"
    allbases2["Algorithm"] = "Plink"
    print allbases2.shape

    keepcols = ["Individual.ID","percent"]
    mrg = merge( allbases[keepcols], allbases2[keepcols], on="Individual.ID" )
    mrg.percent_x = mrg.percent_x.map(float)
    mrg.percent_y = mrg.percent_y.map(float)
    percents_annot = merge(mrg, sampleannot[["Individual.ID","Continent2","Origin"]], on="Individual.ID")

    print hk.dfTable(percents_annot.Continent2)
    percents_annot = percents_annot[ percents_annot.Continent2.isin(
                                    ["East Asia","Europe","Middle East","Africa"]) ]

    print "Alldata for IBC correlation"
    r_dataframe = com.convert_to_r_dataframe(percents_annot)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "percent_x",y="percent_y", colour="factor(Continent2)") +
                ggplot2.geom_point() +
                ggplot2.ggtitle("Correlation between Inbreeding Coefficients") +
                ggplot2.scale_x_continuous("H3M2 Exome coverage (%)") +
                ggplot2.scale_y_continuous("Plink Exome coverage (%)") +
                ggplot2.stat_smooth(method="lm", se=False) +
                ggplot2.theme(**mytheme) )
                #ggplot2.geom_point(ggplot2.aes_string(colour="factor("+factor+")")) + \
                #ggplot2.scale_colour_manual(values = robjects.StrVector(("blue", "red", "grey"))) + \
                #ggplot2.scale_x_continuous("Sum of all variant sites in a gene", limits=robjects.IntVector((0,600))) + \
                #ggplot2.geom_abline(intercept=rstats.coef(model_x)[1], slope=rstats.coef(model_x)[2])
    figname = "%s/%s_ibccorr.png" % (figuredir,prefix)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()

    sys.exit(1)
    keepcols = ["Individual.ID","percent","Algorithm"]
    percents_annot = merge(concat([allbases[keepcols],allbases2[keepcols]]), 
                           sampleannot[["Individual.ID","Continent2","Origin"]], 
                           on="Individual.ID")
    percents_annot = percents_annot[percents_annot.Continent2.isin(
                                    ["Europe","Middle East","Africa"])]
    r_dataframe = com.convert_to_r_dataframe(percents_annot)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor(Continent2)",y="percent", colour="factor(Continent2)") +
                ggplot2.geom_boxplot() +
                ggplot2.ggtitle("Runs of Homozygosity Exome overlap\nBetween Algorithms") +
                ggplot2.facet_grid(robjects.Formula('. ~ Algorithm')) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.theme(**mytheme) )
                #ggplot2.scale_x_continuous("H3M2 Exome coverage (%)") +
                #ggplot2.scale_y_continuous("Plink Exome coverage (%)") +
                #ggplot2.stat_smooth(method="lm", se=False) +
                #ggplot2.geom_point(ggplot2.aes_string(colour="factor("+factor+")")) + \
                #ggplot2.scale_colour_manual(values = robjects.StrVector(("blue", "red", "grey"))) + \
                #ggplot2.scale_x_continuous("Sum of all variant sites in a gene", limits=robjects.IntVector((0,600))) + \
                #ggplot2.geom_abline(intercept=rstats.coef(model_x)[1], slope=rstats.coef(model_x)[2])
    figname = "%s/%s_ibccorr_box.png" % (figuredir,prefix)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()

    percents_annot = merge(concat([allbases[keepcols],allbases2[keepcols]]), 
                           sampleannot[["Individual.ID","GeographicRegions","Origin"]], 
                           on="Individual.ID")

    percents_annot = percents_annot[percents_annot["GeographicRegions"].isin(
                                    target_geographic_regions)]

    r_dataframe = com.convert_to_r_dataframe(percents_annot)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor(GeographicRegions)",y="percent", colour="factor(GeographicRegions)") +
                ggplot2.geom_boxplot() +
                ggplot2.ggtitle("Runs of Homozygosity Exome overlap\nBetween Algorithms") +
                ggplot2.facet_grid(robjects.Formula('. ~ Algorithm')) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.theme(**mytheme) )
    figname = "%s/%s_ibccorr_geobox.png" % (figuredir,prefix)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()


# END Main



