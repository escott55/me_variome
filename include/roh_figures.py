#!/usr/bin/python

from localglobals import *

################################################################################
# percentAnalysis
################################################################################
def percentAnalysis( allbases, exometotals, figureprefix=None ) :
    percents = melt(DataFrame(
                {col:(allbases[int(col)].fillna(0)/float(exometotals[col])*100)
                      for col in exometotals.index if col in allbases.columns}))
    #allbases.divide( exometotals.map(float), fill_value=0 )
    r_dataframe = com.convert_to_r_dataframe(percents)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x = "factor(variable)",y="value" ) + \
                ggplot2.geom_boxplot() + \
                ggplot2.geom_jitter(colour="light blue") + \
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

    #sums = allbases.fillna(0).apply(sum,axis=1)
    #percents = (sums / float(exometotals.sum()) *100).reset_index()
    #percents.columns = ["Individual.ID","percent"]
    percents = allbases[["Individual.ID","percent"]]
    percents.sort("percent", inplace=True)
    levels = percents["Individual.ID"].unique().tolist()
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
# percentAnalysisRegions
################################################################################
def percentAnalysisRegions( allbases, exometotals, sampleannot, figureprefix=None ) :

    allbases.columns = allbases.columns.map(str)
    exometotals.index = exometotals.index.map(str)
    percents_dict = {col:(allbases[col].fillna(0)/float(exometotals[col])*100)
                      for col in exometotals.index if col in allbases.columns.map(str)}
    percents = DataFrame( percents_dict )
    percents["Individual.ID"] = allbases["Individual.ID"]
    percents_melt = melt( percents, id_vars=["Individual.ID"] )
    percents_melt.columns = ["Individual.ID","chrom","percent"]
    percents_annot = merge(percents_melt, sampleannot[["Individual.ID","Continent"]], on="Individual.ID")
    #print percents_annot.head(10)
    r_dataframe = com.convert_to_r_dataframe(percents_annot)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x = "factor(chrom)",y="percent" ) + \
                ggplot2.geom_boxplot() + \
                ggplot2.ggtitle("Chromosomal distribution of ROH") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.scale_y_continuous("RoH (%)" )+ \
                ggplot2.scale_x_discrete("Chromosomes" ) + \
                ggplot2.facet_grid( robjects.Formula('Continent ~ .') )

    if figureprefix is None : figureprefix = "test"
    if figureprefix.find(" ")>=0 : figureprefix = figureprefix.replace(" ","_")
    figname = "results/figures/roh/%s_roh_chr_regions.png" % figureprefix
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()

    percents = allbases[["Individual.ID","percent"]]
    percents_annot = merge(percents, sampleannot[["Individual.ID","Continent"]], on="Individual.ID")
    percents_annot.sort(["Continent","percent"], inplace=True)
    levels = percents_annot["Individual.ID"].unique().tolist()
    #print levels
    #print percents_annot.head(10)
    r_dataframe = com.convert_to_r_dataframe(percents_annot)
    s = robjects.FactorVector(r_dataframe.rx2("Individual.ID"), levels=robjects.StrVector(levels))
    new_r_df = r_dataframe.cbind(s)
    new_r_df.colnames = robjects.StrVector(['Individual.ID', 'percent', 'Continent', 'Sample'])
    p = ggplot2.ggplot(new_r_df) + \
                ggplot2.aes_string(x = "Sample",y="percent",fill="factor(Continent)" ) + \
                ggplot2.geom_bar(stat="identity") + \
                ggplot2.ggtitle("Exome coverage") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.scale_y_continuous("Exome coverage (%)" )+ \
                ggplot2.theme(**{'axis.text.x': ggplot2.element_blank()}) + \
                ggplot2.scale_x_discrete("Samples" )# + \
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \

    figname = "results/figures/roh/%s_roh_samp_regions.png" % figureprefix
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()

    percents_annot = merge(percents, 
                           sampleannot[["Individual.ID","GeographicRegions"]], 
                           on="Individual.ID")
    percents_annot.sort(["GeographicRegions","percent"], inplace=True)
    levels = percents_annot["Individual.ID"].unique().tolist()
    #print levels
    #print percents_annot.head(10)
    r_dataframe = com.convert_to_r_dataframe(percents_annot)
    s = robjects.FactorVector(r_dataframe.rx2("Individual.ID"), levels=robjects.StrVector(levels))
    new_r_df = r_dataframe.cbind(s)
    new_r_df.colnames = robjects.StrVector(['Individual.ID', 'percent', 'GeographicRegions', 'Sample'])
    p = ggplot2.ggplot(new_r_df) + \
                ggplot2.aes_string(x = "Sample",y="percent",fill="factor(GeographicRegions)" ) + \
                ggplot2.geom_bar(stat="identity") + \
                ggplot2.ggtitle("Exome coverage") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.scale_y_continuous("Exome coverage (%)" )+ \
                ggplot2.theme(**{'axis.text.x': ggplot2.element_blank()}) + \
                ggplot2.scale_x_discrete("Samples" )# + \
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \

    figname = "results/figures/roh/%s_roh_samp_georegions.png" % figureprefix
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
    
# END percentAnalysisRegions

################################################################################
# percentDensity
################################################################################
def percentDensity( allbases, sampleannot, figureprefix=None ) :

    percents = allbases[["Individual.ID","percent"]]
    percents_annot = merge(percents, sampleannot[["Individual.ID","Continent",
                                                  "Origin","GeographicRegions"]], 
                           on="Individual.ID")
    percents_annot.sort(["Continent","percent"], inplace=True)
    levels = percents_annot["Individual.ID"].unique().tolist()
    #percents_annot = percents_annot[(
        #percents_annot.Continent.isin(["Europe","Middle East"]))]

    #targetregions = ["Morocco","Algeria","Tunisia","Libya","Egypt","Saudi Arabia","Oman","Qatar","UAE","Yemen","Jordan","Palestine","Lebanon","Syria","Kuwait","Iraq","Turkey"] + ["Iran","Pakistan","Afganistan"]
    #percents_annot = percents_annot[( percents_annot.Origin.isin(targetregions) )]
    percents_annot = percents_annot[( percents_annot["GeographicRegions"].isin(target_geographic_regions_me) )]

    #print levels
    #print percents_annot.head(10)
    r_dataframe = com.convert_to_r_dataframe(percents_annot)
    r_dataframe = fixRLevels( r_dataframe, "GeographicRegions", target_geographic_regions_me )
    #s = robjects.FactorVector(r_dataframe.rx2("Individual.ID"), levels=robjects.StrVector(levels))
    #new_r_df = r_dataframe.cbind(s)
    #new_r_df.colnames = robjects.StrVector(['Individual.ID', 'percent', 'Origin', 'Sample'])
    p = (ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="factor(GeographicRegions)", y="percent", 
                                   fill="factor(GeographicRegions)") + \
                ggplot2.geom_boxplot() + \
                ggplot2.ggtitle("Percent of Exome within\nRuns of Homozygosity") + \
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                ggplot2.theme(**mytheme))
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_blank()}) + \
                #ggplot2.scale_y_continuous("Exome coverage (%)" )+ \
                #ggplot2.scale_x_discrete("Samples" )# + \

    figname = "results/figures/roh/%s_box_roh_samp_regions.png" % figureprefix
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()

    p = (ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="percent" ) + \
                ggplot2.geom_density(ggplot2.aes_string(colour="factor(GeographicRegions)"), 
                                     fill="NA", size=2) + \
                ggplot2.ggtitle("Percent of Exome within\nRuns of Homozygosity") + \
                ggplot2.theme(**mytheme))
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_blank()}) + \
                #ggplot2.scale_y_continuous("Exome coverage (%)" )+ \
                #ggplot2.scale_x_discrete("Samples" )# + \

    figname = "results/figures/roh/%s_density_roh_samp_regions.png" % figureprefix
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END percentDensity

################################################################################
################################################################################
def makeRohCumsumPlot( rohfile, sampleannot, outprefix="results/rohcumplot" ) :
    allintervals = read_csv(rohfile,delim_whitespace=True)
    #indivstats = read_csv(rohfile+".indiv",delim_whitespace=True)
    #print "indivstats:",indivstats.shape

    intdf = merge( sampleannot[["Individual.ID","Continent2","ethnicity"]], 
                  allintervals[["IID","KB","NSNP","DENSITY"]], right_on="IID", 
                  left_on="Individual.ID", how="right" )
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

    allsizedist = concat( allsizedist ).reset_index(drop=True).sort("CumProp")
    print allsizedist[allsizedist.Region == "Oceania"]
    sys.exit(1)
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
    sys.exit(1)
# END makeRohCumsumPlot

def plotDelVars( rohvars ) :
    r_dataframe = com.convert_to_r_dataframe(rohvars)
    p = (ggplot2.ggplot(r_dataframe) +
            ggplot2.aes_string(x = "percent",y="delroh" ) +
            ggplot2.point() +
            ggplot2.ggtitle("Distribution of del vars in ROH") +
            ggplot2.scale_y_continuous("Burden of Variants" ) +
            ggplot2.scale_x_continuous("RoH (%)" ) +
            ggplot2.theme(**mytheme))
            #ggplot2.geom_boxplot() + \
            #ggplot2.geom_jitter(colour="light blue") + \

    #if figureprefix is None : figureprefix = "test"
    #if figureprefix.find(" ")>=0 : figureprefix = figureprefix.replace(" ","_")
    figname = "results/figures/roh/delvars_roh.png"
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotDelVars

################################################################################
# plotIdeogram
################################################################################
def plotIdeogram( population, popsize, bedfile ):
    print "Running plotIdeogram"
    print "Pop",population
    print "popsize",popsize
    print "bedfile",bedfile
    command = ["./rscripts/ideogram.R","single",
               os.path.abspath(bedfile),str(popsize),population.replace(" ","_")]
    print "Command:"," ".join(command)
    out = subprocess.check_output(command)
    print out
# END plot Ideogram

################################################################################
# plotIdeogramDiff
################################################################################
def plotIdeogramDiff( population, popsize, bedfile, population1, popsize1, bedfile1 ):
    print "Running plotIdeogramDiff"
    #print "Pop",population
    #print "popsize",popsize
    #print "bedfile",bedfile
    command = ["./rscripts/ideogram.R","difference",
               os.path.abspath(bedfile), os.path.abspath(bedfile1),
               str(popsize),str(popsize1),
               population +"-"+ population1 ]
    print "Command:"," ".join(command)
    out = subprocess.check_output(command)
    print out
# END plotIdeogramDiff
        
################################################################################
# makeIdeograms
################################################################################
def makeIdeograms( prefix, bedfiles, sampleannot ) :
    #bedfiles = retrieveBedFiles( prefix )
    print "Found bedfiles:",bedfiles
    popfiles = summarizeOverlap( prefix, sampleannot, bedfiles )
    popfiles.index = popfiles.Pop
    for index,pop in popfiles.iterrows() :
        plotIdeogram( pop["Pop"], pop["found"], pop["file"] )
    
    pop1 = popfiles.ix["Middle East"]
    pop2 = popfiles.ix["Europe"]

    plotIdeogramDiff( pop1["Pop"].replace(" ",'.'), pop1["found"], pop1["file"],
                      pop2["Pop"].replace(" ",'.'), pop2["found"], pop2["file"]
                    )
# END makeIdeograms

if __name__ == "__main__" :
    print "Should write main"
# END Main
