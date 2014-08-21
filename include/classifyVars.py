#!/usr/bin/python

import os, sys
import string
from scipy import stats
from scipy.stats import ttest_ind


from localglobals import *
import housekeeping as hk
import king
#
forceFlag = False

#--------------------------------------------------------------------#
#                               Annotation                           #
#--------------------------------------------------------------------#

######################################################################
# annotateVCF
# java -Xmx4G -jar snpEff.jar -i vcf -o vcf GRCh37.69 demo.1kg.vcf > demo.1kg_snpEff.vcf
######################################################################
def annotateVCF( vcffile, targetdir=None, force=False ):
    hk.write2log("#"*70+"\n# Running"+hk.whoami()+"\n"+"#"*70, True )
    print "annotateVCF:",vcffile
    assert os.path.exists( vcffile ), "No file found"

    filepath, basename, suffix = hk.getBasename( vcffile )
    if targetdir is None :
        targetdir = filepath
    effvcf = "%s/%s.eff.vcf" % (targetdir, basename)
    scorevcf = "%s/%s.annot.vcf" % (targetdir, basename)
    #knownvcf = "%s/%s.known.vcf" % (targetdir, basename)
    genevcf = "%s/%s.gene.vcf" % (targetdir, basename)

    allfiles={"eff":effvcf+".gz", "scores":scorevcf+".gz", "gene":genevcf+".gz"}#"known:":knownvcf+".gz", 
    #if os.path.exists( scorevcf+".gz" )  and os.path.exists( genevcf+".gz" ) and not force :
    if all([os.path.exists(x) for x in allfiles.values()]) and not force :
        hk.write2log(" - File"+scorevcf+" already exists. Skipping command!",True)
        return allfiles

    snpeffjar = "/home/escott/snpEff/snpEff.jar"
    snpeffconfig = "/home/escott/snpEff/snpEff.config"

    # Refs : GRCh37.69, HG19
    command = "java -Xmx8G -jar %s -lof -c %s -i vcf -o vcf GRCh37.72 %s > %s " \
            % ( snpeffjar, snpeffconfig, vcffile, effvcf )
    if not os.path.exists( effvcf+".gz" ) :
        print command
        out = hk.runCommand( command )
        #print out
        effvcf = hk.tabixbgzip( effvcf, force )

    # Prediction annotation
    snpsiftjar = "/home/escott/snpEff/SnpSift.jar"
    dblink = "/home/escott/snpEff/dbNSFP2.3.txt.gz"
    #java -jar SnpSift.jar dbnsfp [-q|-v] [-a] dbNSFP.txt.gz file.vcf > newFile.vcf
    command = "java -Xmx4G -jar %s dbnsfp -a %s %s > %s " \
            % ( snpsiftjar, dblink, effvcf, scorevcf )
    if not os.path.exists( scorevcf+".gz" ) :
        print command
        out = hk.runCommand( command )
        #print out
        scorevcf = hk.tabixbgzip( scorevcf, force )

    # DBsnp Annotation
    #snpsiftjar = "/home/escott/snpEff/SnpSift.jar"
    #dbsnp = "/home/escott/snpEff/dbsnp_138.vcf.gz"
    ##java -jar SnpSift.jar annotate dbSnp132.vcf variants.vcf > variants_annotated.vcf 
    #command = "java -Xmx4G -jar %s annotate %s %s > %s " \
            #% ( snpsiftjar, dbsnp, effvcf, knownvcf )
    #if not os.path.exists( knownvcf+".gz" ) :
        #print command
        #out = hk.runCommand( command )
        #print out
        #knownvcf = hk.tabixbgzip( knownvcf, force )

    # GeneSets
    dblink = "/home/escott/snpEff/msigdb.v4.0.symbols.gmt"
    #java -jar SnpSift.jar miSigDb.gmt file.vcf > file.geneSets.vcf
    command = "java -Xmx4G -jar %s geneSets %s %s > %s " \
            % ( snpsiftjar, dblink, effvcf, genevcf )
    if not os.path.exists( genevcf+".gz" ) :
        print command
        out = hk.runCommand( command )
        #print out
        genevcf = hk.tabixbgzip( genevcf, force )
    return allfiles
    #return {'effvcf':effvcf,'scorevcf':scorevcf,'genevcf':genevcf}
# End annotateVCF
#Score (dbtype)  # variants in LJB23 build hg19  Categorical Prediction
#SIFT (sift)     77593284        D: Deleterious (sift<=0.05); T: tolerated (sift>0.05)
#PolyPhen 2 HDIV (pp2_hdiv)      72533732        D: Probably damaging (>=0.957), P: possibly damaging (0.453<=pp2_hdiv<=0.956); B: benign (pp2_hdiv<=0.452)
#PolyPhen 2 HVar (pp2_hvar)      72533732        D: Probably damaging (>=0.909), P: possibly damaging (0.447<=pp2_hdiv<=0.909); B: benign (pp2_hdiv<=0.446)
#LRT (lrt)       68069321        D: Deleterious; N: Neutral; U: Unknown
#MutationTaster (mt)     88473874        A" ("disease_causing_automatic"); "D" ("disease_causing"); "N" ("polymorphism"); "P" ("polymorphism_automatic"
#MutationAssessor (ma)   74631375        H: high; M: medium; L: low; N: neutral. H/M means functional and L/N means non-functional
#FATHMM (fathmm) 70274896        D: Deleterious; T: Tolerated
#MetaSVM (metasvm)       82098217        D: Deleterious; T: Tolerated
#MetaLR (metalr) 82098217        D: Deleterious; T: Tolerated
#GERP++ (gerp++) 89076718        higher scores are more deleterious
#PhyloP (phylop) 89553090        higher scores are more deleterious
#SiPhy (siphy)   88269630        higher scores are more deleterious

#There are two databases for PolyPhen2: HVAR and HDIV. They are explained below:

#ljb2_pp2hvar should be used for diagnostics of Mendelian diseases, which requires distinguishing mutations with drastic effects from all the remaining human variation, including abundant mildly deleterious alleles.The authors recommend calling "probably damaging" if the score is between 0.909 and 1, and "possibly damaging" if the score is between 0.447 and 0.908, and "benign" is the score is between 0 and 0.446.

#ljb2_pp2hdiv should be used when evaluating rare alleles at loci potentially involved in complex phenotypes, dense mapping of regions identified by genome-wide association studies, and analysis of natural selection from sequence data. The authors recommend calling "probably damaging" if the score is between 0.957 and 1, and "possibly damaging" if the score is between 0.453 and 0.956, and "benign" is the score is between 0 and 0.452.

def retrievePolyphenScore( scoredata ):
    #D: Probably damaging (>=0.909), P: possibly damaging (0.447<=pp2_hdiv<=0.909); B: benign
    score = 'unk'
    scoretypes = {'D':'Probably damaging','P':'Possibly damaging','B':'Benign','unk':'.','.':'.'}
    if len(scoredata) == 1: score = scoretypes[scoredata[0]]
    elif "D" in scoredata : score = scoretypes["D"]
    elif "P" in scoredata : score = scoretypes["P"]
    elif "B" in scoredata : score = scoretypes["B"]
    else : score = '.'
    #print "Score:",score
    return score

################################################################################
# parseInfoColumn
# Effect_Impact|Codon_Change|Amino_Acid_change|Gene_Name|Gene_BioType|Coding|
# Transcript|Rank|ERRORS|WARNINGS
################################################################################
def parseInfoColumn( chrom, pos, infodata, key, clinvar ) :
    info = "."
    scoredata = {}
    #print infodata
    for annotsection in infodata.split(";") :
        #print annotsection[:annotsection.find("=")]
        secsplit = annotsection.find("=")
        if secsplit < 0 :
            #print "Annotsection",annotsection
            continue

        sec_name = annotsection[:secsplit]
        sec_value = annotsection[secsplit+1:]
        assert secsplit > 0, ("Error: Info field not equal to anything:")
        if sec_name == "EFF" :
            info = sec_value
        elif sec_name.find("dbNSFP") == 0 :
            scorename = sec_name[sec_name.find("dbNSFP")+7:]
            scoredata[scorename] = sec_value.split("|")
        elif sec_name in ["LOF","NMD"] :
            scoredata[sec_name] = [sec_value]
        #elif annotsection[:annotsection.find("=")] == "dbNSFP_Polyphen2_HVAR_pred" :
            #scoredata = annotsection[annotsection.find("=")+1:].split("|")
    if not scoredata.has_key("LOF") : scoredata["LOF"] = [""]
    if not scoredata.has_key("NMD") : scoredata["NMD"] = [""]
    scoredata = Series( scoredata  )
    
    #if infodata.find("LOF") > 0 : 
        #print infodata
        #print scoredata
    #varscore = retrievePolyphenScore(scoredata["Polyphen2_HVAR_pred"])
    if len(info) <2 : print "Failed info:",info; return None, None, None, '.'

    aheadlist = ["chrom","pos","FunctionGVS","Priority","Effect","Codon","AA","AAlength","Gene",
                 "Vartype","Coding","Transcript","Rank","Errors","Warnings"]
    #print [re.split("[\(\)\|]",x)[:13] for x in info.split(",") ]
    vardata = DataFrame( ([[chrom,pos] + re.split("[\(\)\|]",x)[:13] for x in info.split(",") ]), 
                        columns=aheadlist )
    vardata = vardata[vardata.Vartype == "protein_coding"]

    priorities = DataFrame({"Priority":["HIGH","MODERATE","LOW","MODIFIER"],"vprior":[1,2,3,4]})
    vardata = merge(vardata,priorities,how="left",on="Priority")
    varcounts = DataFrame({'penetrance' : 
                           vardata.groupby(["chrom","pos","FunctionGVS","Gene","Priority","vprior"])
                           .size()}).reset_index()
    tokeep = varcounts[["Gene","vprior"]].groupby("Gene",as_index=False).min()
    isocounts =  varcounts[["Gene","penetrance"]].groupby("Gene",as_index=False).sum()
    isocounts.columns = ["Gene","isoforms"]
    vcounts_filt = varcounts.ix[tokeep.index]
    vcounts_filt = merge(vcounts_filt, isocounts, on="Gene")
    for col in scoredata.index : 
        #print col,scoredata[col]
        if len(scoredata[col]) == len(vcounts_filt) : vcounts_filt[col] = scoredata[col]
        else : vcounts_filt[col] = scoredata[col][0]

    #for col in vcounts_filt.columns : 
        #print col, vcounts_filt[col].values
    vclass = []
    for index,iso in vcounts_filt.iterrows() :
        if "Polyphen2_HVAR_pred" not in iso.index : vclass.append("0"); continue
        vclass.append(getVarClass( chrom, pos, iso["Priority"], iso["Polyphen2_HVAR_pred"], clinvar ))

    vcounts_filt["vclass"] = vclass
    return vcounts_filt
# END parseInfoColumn

################################################################################
# parseRegionsFile
################################################################################
def parseRegionsFile():
    #refseqgenesbed = "./bedfiles/refseq_genes_simple.bed"
    #bedfile = "./resources/regionbeds/percentcoverage.bed"
    bedfile = "./resources/regionbeds/collab_percentcoverage.txt"

    regions = {}
    for row in csv.reader(open(bedfile), delimiter="\t") :
        if not regions.has_key(row[0]) : regions[row[0]] = []
        regions[row[0]].append([int(row[1]),int(row[2])])
        #print [int(row[1]),int(row[2]),row[3]]
    return regions
# End parseRegionsFile

################################################################################
# intersectRegions
################################################################################
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

################################################################################
# plotAFDistClass
################################################################################
def plotAFDistClass( alldata, figurename ) :
    #alldata["distance"] = alldata.het + alldata.hom
    print "Running AFDistClass"

    print alldata.head(10)
    bins=[0.,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]
    tbin = []
    for af in alldata.AF.tolist() : 
        found = None
        for s, e in zip(bins[:-1], bins[1:]):
            if af >= s and af <= e : 
                found = str(e)
        if found is not None : tbin.append(found)
        else : print "AF:",af

    alldata["Hist"] = tbin
    #chrom     pos     Vtype                 Vclass        AF  Hist
    #0     1  861329  MODERATE  NON_SYNONYMOUS_CODING  0.000865  0.01
    summarized = []
    #for group, cdata in alldata.groupby("Vtype") :
    alldata = alldata[alldata.vclass > -1]
    groupdf = alldata.groupby("vclass")
    #for group in ["HIGH","MODERATE","MODIFIER","LOW"] :
    for group in groupdf.groups :
        cdata = groupdf.get_group(group)
        #print "running group:",group
        dfcounts = cdata.groupby(["Hist"]).size()
        total = float(len(cdata))
        dfcounts = DataFrame({'Attr':dfcounts.index, 'Freq':dfcounts.tolist()})
        dfcounts["proportions"] = dfcounts["Freq"]/total
        dfcounts["group"] = group
        summarized.append(dfcounts)
    summarized = concat( summarized ).reset_index()

    #hist,bin_edges=np.histogram(alldata.,bins=bins)
    r_dataframe = com.convert_to_r_dataframe(summarized)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="factor(Attr)", y="proportions", fill="factor(group)") + \
                ggplot2.geom_bar(stat="identity",position="dodge") + \
                ggplot2.ggtitle("Proportions of Mutations by Variant Class") + \
                ggplot2.theme(**mytheme) #+ \
                #ggplot2.facet_grid( robjects.Formula('group ~ .') )
                #ggplot2.scale_y_log10() + ggplot2.xlim(0., 1.) + \
                #ggplot2.geom_histogram(breaks=robjects.FloatVector([0,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])) + \
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                #ggplot2.aes_string(fill="factor(continent)")

    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# END plotAFDistClass

def calcAF(genecounts):
    #print [x for x in genecounts if x.find(":") == -1]
    ref,het,hom = [int(x) for x in genecounts.split(":")]
    if ref+het+hom == 0 : print "Err:",genecounts; return 0
    return (1.*hom+het) / (1.*(ref+het+hom))
# END calcAF

################################################################################
# plotAFDistClassGroups
################################################################################
def plotAFDistClassGroups( distdf, figurename ) :
    print "Running plotAFDistClassGroups"

    distdf["vtotals"] = distdf.value.apply( lambda x: 
                                           sum([int(y) for y in x.split(":")]))

    distdf["AF"] = [calcAF(x) for x in distdf["value"].tolist()]
    distdf_filt = distdf[(distdf.AF != 0.0) &
                         (distdf.vtotals > 10) &
                         (distdf.vclass > 0)]
    print distdf_filt.head(30)

    bins=[0.,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]
    afbinned = []
    for index, groupdf in distdf_filt.groupby(["variable","vclass"]) :
        region, vclass = index
        afbinsub = (DataFrame({
                            'vcount':groupdf.groupby([np.digitize(groupdf.AF, bins)]).size()
                           }).reset_index())
        print afbinsub.head(10)
        afbinsub["afbin"] = [bins[int(x)] for x in afbinsub.index]
        afbinsub["proportions"] = afbinsub["vcount"].map(float) / len(groupdf)
        afbinsub["Region"] = region
        afbinsub["vclass"] = vclass
        print afbinsub.head(10)
        afbinned.append(afbinsub)
    afbinned = concat( afbinned ).reset_index()

    #afbinned = afbinned[afbinned.vclass.isin([1,2,3,4])]

    print afbinned.head(10)
    #hist,bin_edges=np.histogram(alldata.,bins=bins)
    r_dataframe = com.convert_to_r_dataframe(afbinned)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(afbin)", y="proportions", fill="factor(vclass)") +
                ggplot2.geom_bar(stat="identity",position="dodge") +
                ggplot2.ggtitle("Proportions of Mutations by Variant Class") +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.facet_grid( robjects.Formula('Region ~ .') ) +
                ggplot2.theme(**mytheme) )
                #ggplot2.scale_y_log10() + ggplot2.xlim(0., 1.) + \
                #ggplot2.aes_string(fill="factor(continent)")

    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# END plotAFDistClassGroups

################################################################################
# plotAFDistClassDensity
################################################################################
def plotAFDistClassDensity( distdf, figurename ) :
    print "Running plotAFDistClassDensity"
    print distdf.head(30)

    #distdf["vtotals"] = distdf.value.apply( lambda x: 
                                           #sum([int(y) for y in x.split(":")]))
    #distdf["AF"] = [calcAF(x) for x in distdf["value"].tolist()]
    #distdf_filt = distdf[(distdf.AF != 0.0) &
                         #(distdf.vtotals > 10) &
                         #(distdf.vclass > 0)]
    #print distdf_filt.head(30)

    #bins=[0.,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]
    #afbinned = []
    #for index, groupdf in distdf_filt.groupby(["variable","vclass"]) :
        #region, vclass = index
        #afbinsub = (DataFrame({
                            #'vcount':groupdf.groupby([np.digitize(groupdf.AF, bins)]).size()
                           #}).reset_index())
        #print afbinsub.head(10)
        #afbinsub["afbin"] = [bins[int(x)] for x in afbinsub.index]
        #afbinsub["proportions"] = afbinsub["vcount"].map(float) / len(groupdf)
        #afbinsub["Region"] = region
        #afbinsub["vclass"] = vclass
        #print afbinsub.head(10)
        #afbinned.append(afbinsub)
    #afbinned = concat( afbinned ).reset_index()
    print distdf[distdf.variable == "Africa"].head()
    sys.exit(1)
    distdf_filt = distdf[(distdf.vclass.isin([1,2,3,4])) & (distdf.AF <.05) & (distdf.AF >.0)]

    print distdf_filt.groupby("variable").size().head()

    #hist,bin_edges=np.histogram(alldata.,bins=bins)
    r_dataframe = com.convert_to_r_dataframe(distdf_filt)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="AF", group="factor(variable)", colour="factor(variable)") +
                ggplot2.geom_density() +
                ggplot2.ggtitle("AF Density by Variant Class") +
                ggplot2.facet_grid( robjects.Formula('vclass ~ .'), scale="free" ) +
                ggplot2.theme(**mytheme) )

                #ggplot2.scale_y_log10() + ggplot2.xlim(0., 1.) + \
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                #ggplot2.aes_string(fill="factor(continent)")

    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# END plotAFDistClassDensity

################################################################################
# getVarClass
################################################################################
def getVarClass( chrom, pos, vprior, score, clinvar ) :
    vclass = -1
    if "chr"+str(chrom)+":"+str(pos) in clinvar.index : vclass = 5
    elif vprior == "HIGH" : vclass = 4
    elif vprior == "MODERATE" and score == "D" : vclass =3
    elif (vprior == "MODERATE" and 
          (score == "P" or score is np.nan)) : vclass =2
    elif vprior == "MODERATE" and score == "B" : vclass =1
    elif vprior == "LOW" : vclass =1
    elif vprior == "MODIFIER" : vclass =0
    #else : print "Info:",vprior,":",score,"Class:",vclass
    return vclass
# END getVarClass


def calcRegionAF( regionlist, sampregions, genotypes ) :
    regionafs = {}
    for region in regionlist :
        print region, str(region)
        if str(region) == "nan" : continue
        subgeno = [ genotypes[i] for i in range(len(genotypes)) if sampregions[i] == region ]
        regionafs[region]= "%d:%d:%d" %(subgeno.count('0/0'),
                                      subgeno.count('0/1'),
                                      subgeno.count('1/1'))
    return Series(regionafs)
# END calcRegionAF

################################################################################
# classifyVars
################################################################################
def classifyVars( effvcf, callvcf, sampleannot, regionlist, filepats, 
                 outdir="./classes", limit=-1, force=False ) :
    hk.write2log("#"*70+"\n# Running "+hk.whoami()+"\n"+"#"*70, True )

    print "Callvcf:",callvcf
    print "Effvcf:",effvcf

    clinvar = parseClinVar(vfilter=5)
    print "Finished parsing clinvar"

    filepath, basename, suffix = hk.getBasename( callvcf )
    isoformfile = "%s/%s_genes.tsv" % (filepath, basename )
    if os.path.exists(isoformfile) and not force: 
        print "Warning: isoformfile already exists:",isoformfile
        return isoformfile

    regions = dict( [[samp,cont] 
                     for samp,cont in sampleannot[["Individual.ID","Continent2"]].values] )
    sampregions = [regions[x] if x in sorted(regions) else "None" for x in filepats]

    vcfcols = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Meta"]
    print "Attempting to read file",effvcf
    annotation = (read_csv(effvcf, sep="\t", comment='#', compression="gzip",
                          names=vcfcols,header=None,low_memory=False)
                  .dropna(how='all')
                  .reset_index(drop=True)) #skiprows=range(4),

    print "Finished annotation"
    print annotation.head()

    annotation["CHROM"] = [str(int(x)) if hk.is_int(x) else str(x) 
                           for x in annotation["CHROM"].tolist()]


    annotation.index = ("chr"+annotation["CHROM"]+":"+
                        annotation["POS"].map(int).map(str)+':'+
                        annotation["REF"]+':'+annotation["ALT"])

    targetcolumns = ['chrom', u'pos', u'FunctionGVS', u'Gene', u'Priority', 
                     u'vprior', u'penetrance', u'isoforms', u'1000Gp1_AF', 
                     u'1000Gp1_AFR_AF', u'1000Gp1_AMR_AF', u'1000Gp1_ASN_AF', 
                     u'1000Gp1_EUR_AF', u'29way_logOdds', u'ESP6500_AA_AF', 
                     u'ESP6500_EA_AF', u'GERP++_NR', u'GERP++_RS', 
                     u'Interpro_domain', u'LOF', u'LRT_pred', u'MutationTaster_pred', 
                     u'NMD', u'Polyphen2_HDIV_pred', u'Polyphen2_HVAR_pred', 
                     u'SIFT_pred', u'vclass', u'AF', 
                     u'Missing', u'Africa', u'America', u'East Asia', u'Europe', 
                     u'Middle East', u'South Asia']
                    #u'Uniprot_acc',
    #regions = parseRegionsFile()
    vartypes = {}
    varclasses = {}
    samples = []

    comments = []
    header = []
    if os.path.splitext( callvcf )[1] == ".gz" : FH = gzip.open( callvcf )
    else : FH = open(callvcf)

    print "Writing isoformfile", isoformfile
    iso_out = open(isoformfile, "wb")
    #for row in csv.reader( FH, delimiter="\t" ) :
    linecount = 0
    for line in FH :
        if line[:2] == "##" or len(line) < 1:
            comments.append( line )
            continue
        row = line.rstrip().split("\t")
        if row[0] == "#CHROM" :
            header = row
            samples = [x.strip() for x in row[9:]]
            continue

        linecount += 1
        chrom,pos,dbsnp,ref,mut,qual,passed = row[:7]
        pos,dbsnp,end = (int(pos),str(dbsnp!="."),(int(pos)+max(len(ref),len(mut))))
        if len(mut)+len(ref) > 100 : continue
        key = "chr%s:%s:%s:%s" % (chrom, pos, ref, mut )
        # retrieve var class
        if key not in annotation.index : continue #print "Key not found!"
        info, meta = annotation.loc[key][["INFO","Meta"]]
        meta = meta.split(",") if meta is not np.nan else ""
        for att in meta : 
            if att.find("=") <= 0 : continue
            key,value = att.split("=",1)
            if key == "AF" : AF = float(value)
            if key == "Missing": missingness = value
        if AF == 1. : print "Everybody has it"; continue

        genotypes = [ x[:3].replace("|","/") for x in row[9:] ]
        #isoformdata, scoredata = parseInfoColumn( info, key )
        isoformdata = parseInfoColumn( chrom, pos, info, key, clinvar )
        regionaf = calcRegionAF( regionlist, sampregions, genotypes )
        isoformdata["AF"] = AF
        isoformdata["Missing"] = missingness
        # what is going on here?
        for col in regionaf.index : isoformdata[col] = regionaf[col]

        for col in targetcolumns : 
            if col not in isoformdata.columns: isoformdata[col] = "" 
        #print isoformdata[ isoformdata.FunctionGVS == "SYNONYMOUS_CODING" ].values
        isoformdata.to_csv(iso_out, header=(linecount ==1), index=False, 
                           sep="\t", columns=targetcolumns)
        
        for index, iso in isoformdata.iterrows() :
            if not varclasses.has_key(iso["FunctionGVS"]) : 
                varclasses[iso["FunctionGVS"]] = 0
            varclasses[iso["FunctionGVS"]] += 1

        #if linecount > 40 :
            #print "Breaking due to hard linecount limit"
            #break

    iso_out.close()
    FH.close()

    return isoformfile
# END classifyVars

################################################################################
# plotClassDist
################################################################################
def plotClassDist( genefile, regionlist, outdir="./results/figures/classes" ) :
    print "Running plotClassDist"
    filepath, basename, ext = hk.getBasename( genefile )
    print "Reading file:",genefile
    genedata = read_csv(genefile, sep="\t")
    regionlist = [x for x in regionlist if x in genedata.columns]
    afdist = genedata[["chrom","pos","FunctionGVS","vclass","AF","Polyphen2_HVAR_pred"]+regionlist]
    figurename = "%s/%s_AFdistclasses.png" % (outdir,basename)
    #afdist = DataFrame( afdist, columns=["chrom","pos","FunctionGVS","vclass","AF"]+regionlist )
    keepcols = ["chrom","pos","AF","vclass"] + regionlist
    print afdist.shape
    afdist_filt = (afdist[keepcols]
                   .groupby(["chrom","pos","AF"])
                   .max().reset_index())
    print afdist_filt.shape
    print afdist_filt.head()

    plotAFDistClass( afdist_filt, figurename )

    afdistmelt = melt(afdist_filt,id_vars=["chrom","pos","vclass"],value_vars=regionlist)
    afdistmelt.to_csv("%s/%s_dist.tsv" %(outdir,basename), sep="\t", index=False)

    figurename = "%s/%s_AFdistclasses_grp.png" % (outdir,basename)
    plotAFDistClassGroups( afdistmelt, figurename )

    figurename = "%s/%s_AFdistclasses_density.png" % (outdir,basename)
    plotAFDistClassDensity( afdistmelt, figurename )

    varclasses = genedata.groupby("FunctionGVS").size().to_dict()
    statsfile = "%s/%s_stats.txt" % (outdir, basename)
    STATSOUT = open(statsfile,"wb")
    STATSOUT.write( "Varclasses:\n" )
    for vfunction in varclasses : 
        STATSOUT.write( "%s\t%s\n" % (vfunction,varclasses[vfunction]) )

    vartypes = genedata.groupby("vclass").size().to_dict()
    STATSOUT.write( "\nVartypes:\n" )
    for vclass in sorted(vartypes) :
        STATSOUT.write( "%s\t%s\n" % (vclass,vartypes[vclass]) )

    #STATSOUT.write( "\nVarclasses:\n" )
    #STATSOUT.write( "%s\tCounts:%d\n" % (vclass,len(classified[vclass]) ))

    STATSOUT.close()

    return 
# END plotClassDist

######################################################################
# plotAFDist
######################################################################
def plotAFDist( alldata, figurename ) :
    #alldata["distance"] = alldata.het + alldata.hom

    r_dataframe = com.convert_to_r_dataframe(alldata)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="AF") + \
                ggplot2.geom_histogram() + \
                ggplot2.ggtitle("AF distribution") + \
                ggplot2.scale_y_log10() + ggplot2.xlim(0., 1.) + \
                ggplot2.theme(**mytheme) #+ \
                #ggplot2.geom_histogram(breaks=robjects.FloatVector([0,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])) + \
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) + \
                #ggplot2.facet_grid( robjects.Formula('RVIS_type ~ .') )
                #ggplot2.aes_string(fill="factor(continent)")

    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# END plotAFDist

######################################################################
# makeSmallVCF
######################################################################
def makeSmallVCF( annotvcf, figuredir="./results/figures/classes", force=False ) :
    hk.write2log("#"*70+"\n# Running"+hk.whoami()+"\n"+"#"*70, True )
    print "Running makeSmallVCF",annotvcf

    filepath, basename, suffix = hk.getBasename( annotvcf )
    vcfoutfile = "%s/%s.small.vcf" % (filepath, basename)
    if os.path.exists(vcfoutfile+".gz") and not force : 
        return vcfoutfile+".gz"

    if os.path.splitext( annotvcf )[1] == ".gz" :
        FH = gzip.open( annotvcf )
    else :
        FH = open(annotvcf)

    regions = parseRegionsFile()

    print "Making file:",vcfoutfile
    OUT = open( vcfoutfile, "wb" )
    OUT.write("##fileformat=VCFv4.1\n")

    linecount = 0
    limit = -1
    afdist = []
    for line in FH :
        linecount += 1
        if line[:2] == "##" or len(line) < 1:
            #comments.append( line )
            continue
        row = line.split("\t")
        if row[0] == "#CHROM" :
            OUT.write('\t'.join(row[:9])+'\tMeta\n')
            continue
        chrom,pos,dbsnp,ref,mut,qual,passed = row[:7]
        if len(mut)+len(ref) > 100 : continue

        # Calculate meta info
        genotypes = [ x[:3].replace("|","/") for x in row[9:] ]
        missingness = float(genotypes.count("./.")) / len(genotypes)
        totalgeno = []
        for geno in genotypes :
            totalgeno += geno.split("/")

        if totalgeno.count("0") + totalgeno.count("1") == 0 : 
            print "Skipping 0 count"; continue

        AF = float(totalgeno.count("1")) / (totalgeno.count("0") + totalgeno.count("1"))
        if AF == 0.0 or AF == 1.0 : continue #print "Filter by AF"; continue
        afdist.append([chrom,pos,AF,missingness])
        filters = []
        if passed == False : filters.append("Passed")
        if AF > .1 : filters.append("AF")
        if missingness > .05 : filters.append("Missingness")
        if intersectRegions( chrom, pos, regions ) : filters.append("Regions")
        OUT.write( "\t".join(row[:9]) )
        if len(filters) > 0 : OUT.write( "\tMissing=%f,AF=%f,Filters=%s\n" %(missingness,AF,";".join(filters)) )
        else : OUT.write( "\tMissing=%f,AF=%f,Filters=\n" %(missingness,AF) )

        limit -= 1
        if limit == 0 : break

    OUT.close()
    figurename = "%s/%s_AFdist.png" % (figuredir, basename)
    afdist = DataFrame( afdist, columns=["chrom","pos","AF","missingness"] )
    plotAFDist( afdist, figurename )
    smallvcf = hk.tabixbgzip( vcfoutfile, force=True )
    return smallvcf
# END makeSmallVCF

######################################################################
# makeGeminiDB
######################################################################
def makeGeminiDB( ):
    return
# End makeGeminiDB

################################################################################
# burdenAnalysis
################################################################################
def burdenAnalysis( tfile ) :
    print "Running burdenAnalysis"
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

################################################################################
# readAnnotData
################################################################################
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

################################################################################
# allPatVarStats
################################################################################
def allPatVarStats( allclasses, outdir="./classes" ) :
    print "Running allPatVarStats"
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
    #print allclasses, header
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
        #print entry
    OUT.close()
    return
# END allPatVarStats

################################################################################
# individualVarStats
################################################################################
def individualVarStats( vcffile,geneannotfile,sampleannot,regionlist,filepats,force=False):
    print "Running individualVarStats"

    filepath, basename, suffix = hk.getBasename( vcffile )
    indivstatsfile = "%s/%s_indiv.tsv" % (filepath, basename )

    if os.path.exists(indivstatsfile) and not force: 
        print "Warning: File already exists:",indivstatsfile
        return indivstatsfile

    genedata = read_csv(geneannotfile, sep="\t")
    regionlist = [x for x in regionlist if x in genedata.columns]
    print "Regions included:",regionlist

    keepcols = ["chrom","pos","AF","vclass","NMD","LOF"] + regionlist
    keep_idx = genedata.groupby(["chrom","pos"])["vclass"].idxmax()
    varfilt = genedata.iloc[keep_idx][keepcols]
    varfilt.index = varfilt.chrom.map(str)+":"+varfilt.pos.map(str)
    vclasses = varfilt.vclass.unique().tolist()

    if os.path.splitext( vcffile )[1] == ".gz" :
        FH = gzip.open( vcffile )
    else : FH = open(vcffile)

    patcounts = {}
    for pat in filepats :
        if not patcounts.has_key(pat) : patcounts[pat] ={}
        for vclass in vclasses+["LOF","NMD"] :
            patcounts[pat][vclass] = [0,0]

    #for row in csv.reader( FH, delimiter="\t" ) :
    comments, header = [],[]
    linecount = 0
    for line in FH :
        linecount += 1
        if line[:2] == "##" or len(line) < 1:
            comments.append( line );continue
        row = line.rstrip().split("\t")
        if row[0] == "#CHROM" :
            header = row
            samples = [x.strip() for x in row[9:]]
            continue

        chrom,pos,dbsnp,ref,mut,qual,passed = row[:7]
        pos = int(pos)
        dbsnp = str(dbsnp != ".")
        end = pos + max(len(ref), len(mut))
        vid = chrom+':'+str(pos)
        if vid not in varfilt.index : 
            print "Error: cant find var:",vid;continue

        vclass = varfilt.ix[vid]["vclass"]
        islof = varfilt.ix[vid]["LOF"] is not np.nan
        isnmd = varfilt.ix[vid]["NMD"] is not np.nan
        if type(vclass) is Series : 
            print "Error: vclass isnt series!",vclass
            continue

        if vclass not in vclasses: 
            print "Error: unrecognized vclass",vclass
            continue

        for i in range(9,len(row)) : 
            genotype = sum([int(x) for x in row[i][:3].split("/") 
                            if hk.is_int(x)])
            if genotype == 0 : continue
            csamp = samples[i-9]
            #print "I:",i, csamp,"-",vclass,"-", genotype
            patcounts[csamp][vclass][genotype-1] += 1
            if islof : 
                patcounts[csamp]["LOF"][genotype-1] += 1
            if isnmd : 
                patcounts[csamp]["NMD"][genotype-1] += 1

        #if linecount == 500 : print "Breaking!"; break

    FH.close()
    allstats = []

    for pat in patcounts :
        if pat not in sampleannot.index :
            print patcounts[pat]; continue
        region = sampleannot.ix[pat]["Continent2"]
        country = sampleannot.ix[pat]["Origin"]
        ethnicity = sampleannot.ix[pat]["ethnicity"]
        sampdat = {'Sample':pat,'Region':region,
                   'Country':country,'Ethnicity':ethnicity}
        for vclass in patcounts[pat] :
            sampdat["Class%s_het"%vclass] = patcounts[pat][vclass][0]
            sampdat["Class%s_hom"%vclass] = patcounts[pat][vclass][1]
        sampdat = Series(sampdat)
        allstats.append(sampdat)

    allstats = DataFrame( allstats )
    #print "LOF count:",sum(allstats.ClassLOF_het > 0),sum(allstats.ClassLOF_hom > 0)
    #print "NMD count:",sum(allstats.ClassNMD_het > 0),sum(allstats.ClassNMD_hom > 0)

    cols = ["Sample","Region","Country","Ethnicity",
            "Class-1_het", "Class-1_hom","Class0_het","Class0_hom",
            "Class1_het","Class1_hom","Class2_het","Class2_hom",
            "Class3_het","Class3_hom","Class4_het","Class4_hom",
            "ClassLOF_het","ClassLOF_hom","ClassNMD_het","ClassNMD_hom" ]
    print "Writing file:",indivstatsfile
    allstats.to_csv(indivstatsfile,sep="\t",index=False,columns=cols)
    return indivstatsfile
# END individualVarStats

################################################################################
# plotSampleCounts
################################################################################
def plotSampleCounts( samplecountsfile, outdir="./results/figures/classes", force=False) :
    filepath, basename, ext = hk.getBasename( samplecountsfile )
    samplecounts = read_csv(samplecountsfile, sep="\t")
    #print "$"*70
    #print "Annotating!"
    samplecounts["Source"] = "Unk"
    samplecounts.index = samplecounts.Sample
    samplecounts = addSampleAnnotation( samplecounts, "Sample", True )
    #print samplecounts_annot.head()
    vclasses = ["Class1","Class2","Class3","Class4"]#,"Class5"
    for vclass in vclasses :
        samplecounts[vclass] = (samplecounts[vclass+"_het"] + 
                                samplecounts[vclass+"_hom"])
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

    print samplecounts[["Sample","Region","Country","Ethnicity","homtotal","homozygosity"]].head(10)
    scounts = melt(samplecounts,id_vars=["Sample","Region","Ethnicity","homozygosity","Source"],
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

    vclasses = ["ClassNMD","ClassLOF"]#,"Class5"
    for vclass in vclasses :
        samplecounts[vclass] = (samplecounts[vclass+"_het"] + 
                                samplecounts[vclass+"_hom"])
    scounts = melt(samplecounts,id_vars=["Sample","Region","Ethnicity"],
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
# END plotSampleCounts

def calcPairwisePvalues( tdf, groupcol, factcol, onetailed=False ) :
    pvals = []
    for grp in tdf[groupcol].unique():
        euro = tdf[(tdf.Continent2=="Europe") & (tdf[groupcol] ==grp)][factcol].tolist()
        me = tdf[(tdf.Continent2=="Middle East") & (tdf[groupcol] ==grp)][factcol].tolist()
        maxY = max(euro+me)
        ts,pv = ttest_ind(euro, me)
        #p/2 < alpha and t < 0
        #if onetailed : print "Todo"
        pvals.append( Series(data=[grp,pv,ts,"Europe",maxY], 
                             index=[groupcol,"pvalue","T-statistic","x","y"]))
    pvals = DataFrame(pvals)
    pvals["Pvalue"] = ["P-value: %.3g\nT-statistic: %.2f" %(x,y) 
                       for x,y in pvals[["pvalue","T-statistic"]].values]
    print pvals.head()
    return pvals
# END calcPairwisePvalues

def calcHomoz( varcounts ):
    ref,het,hom = [int(x) for x in varcounts.split(":")]
    if het+hom == 0 : return 0
    return hom/float(het+hom)
# END calcHomoz

################################################################################
# plotHomozygosity
################################################################################
def plotHomozygosity( geneannotfile, regionlist, outdir) :
    print "Running plotHomozygosity"
    vardf = read_csv( geneannotfile, sep="\t" )

    filepath, basename, suffix = hk.getBasename( geneannotfile )

    homzdf = melt(vardf,id_vars=["chrom","pos","vclass"],value_vars=regionlist)
    # Calculate AF

    homzdf = homzdf[homzdf["value"].notnull()]

    bins=[0.,.005,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]
    if basename.find("variome") >= 0 :
        levels = ["Europe","Middle East"]
    elif basename.find("meceu") >= 0 :
        levels = ["Europe","Middle East"]
    elif basename.find("me1000G") >= 0 :
        bins=[0.,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]
        levels = ["Africa","Europe","Middle East"]

    assert len(homzdf[homzdf["value"].isnull()]) == 0

    homzdf["AF"] = [calcAF(x) for x in homzdf["value"]]

    # Calculate Homozygosity
    homzdf["Homozygosity"] = [calcHomoz(x) for x in homzdf["value"]]

    #afbinned = DataFrame({'vcount':countdat.groupby(["dbsnp",np.digitize(countdat.AF, bins)]).size()}).reset_index()
    print homzdf.head(10)
    #print cut(homzdf.AF, bins)
    print "$"*60
    #print  (homzdf[homzdf.Homozygosity > 0][["variable","AF","Homozygosity"]]
                          #.groupby(["variable",cut(homzdf.AF, bins)])
                          #.mean())
    afbinned = homzdf[homzdf.Homozygosity > 0][["variable","AF","Homozygosity"]]
    afbinned["bins"] = cut(afbinned.AF, bins)

    print afbinned.head(10)
    r_dataframe = com.convert_to_r_dataframe(afbinned)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(bins)", y="Homozygosity", fill="factor(variable)") +
                ggplot2.geom_boxplot() + \
                ggplot2.ggtitle("Fraction of Homozygous Variants") +
                ggplot2.facet_grid( robjects.Formula('variable ~ .') ) +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45),
                                 'legend.position':"none"}) +
                ggplot2.scale_x_discrete("AF Binned") +
                ggplot2.scale_y_continuous("Homozygous Fraction") +
                ggplot2.theme(**mytheme))
                #ggplot2.scale_y_log10() + ggplot2.xlim(0., 1.) + \
                #ggplot2.geom_histogram(breaks=robjects.FloatVector([0,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])) + \
                #ggplot2.aes_string(fill="factor(continent)")

    figurename = "%s/%s_homozaf.png" % (outdir, basename)
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()
# END plotHomozygosity

def getLevels( annotcol, basename, samplecounts ):

    samplecounts.index = samplecounts.Sample
    if "Ethnicity" in samplecounts.columns : 
        samplecounts.drop("Ethnicity",1,inplace=True)

    if "Region" in samplecounts.columns : 
        samplecounts.drop("Region",1,inplace=True)

    scounts_annot = addSampleAnnotation( samplecounts, "Sample" )
    scounts = scounts_annot[(scounts_annot[annotcol].notnull()) & 
                                (scounts_annot[annotcol] != "Unknown")]
    scounts[annotcol] = [x.strip().replace(" ",".") for x in scounts[annotcol]]

    levels = None
    if annotcol == "GeographicRegions2" :
        meregions = target_geographic_regions_me
        levels = ["Africa","Europe"]+tregions
    elif annotcol == "Continent": 
        levels = target_continents
    elif annotcol == "Continent2": 
        levels = target_continents
    elif annotcol == "GeographicRegions" :
        levels = target_geographic_regions

    #print levels
    scounts_sub = scounts[scounts[annotcol].isin(levels)]

    if basename.find( "meceu" ) >= 0: 
        scounts_sub = scounts_sub[((scounts_sub.Continent2 =="Europe") &
                                   (scounts_sub.Source=="Daly")) |
                                  ((scounts_sub.Continent2 =="Middle East") &
                                   (scounts_sub.Source=="BROAD"))] 
    elif basename.find( "1000G" ) >= 0 :

        scounts_sub = scounts_sub[((scounts_sub[annotcol].isin(["Europe","Africa"]))) | 
                              ((scounts_sub.Continent2== "Middle East") &
                               (scounts_sub.Source=="BROAD"))] 
                                #"Middle East","South Asia" (scounts.Source=="1KG")) |

    elif basename.find( "variome" ) == 0 :
        print "nothing"
        scounts_sub = (scounts_sub[scounts_sub.Continent2.isin(
                                        ["Europe","Africa","Middle East"])])
    else : 
        print "Error: Unknown basename",basename; sys.exit(1)


    finallevels = [x for x in levels if x in scounts_sub[annotcol].unique()]
    return finallevels, scounts_sub
# END getLevels

################################################################################
# plotSampleCounts2
################################################################################
def plotSampleCounts2( samplecountsfile, outdir="./results/figures/classes", 
                      annotcol="GeographicRegions", force=False) :

    print "Running plotSampleCounts2 : annotcol-", annotcol
    filepath, basename, ext = hk.getBasename( samplecountsfile )
    samplecounts = read_csv(samplecountsfile, sep="\t")
    vclasses = ["Class1","Class2","Class3","Class4"]#,"Class5"
    lofclasses =["ClassNMD","ClassLOF"]
    #samplecounts["Source"] = "Unk"
    #for vclass in vclasses+lofclasses : 
        #samplecounts[vclass] = samplecounts[vclass+"_het"]+samplecounts[vclass+"_hom"]
    valcols = (  [x+"_het" for x in vclasses+lofclasses] 
               + [x+"_hom" for x in vclasses+lofclasses] )

    print "Before annotation"
    print samplecounts.head(10)
    tlevels, scounts_annot = getLevels(annotcol, basename, samplecounts )
    print "After annotation"
    print samplecounts.head(10)

    #print samplecounts.head()
    #print samplecounts[annotcol].unique()

    scounts = melt(scounts_annot,id_vars=["Sample",annotcol,"Continent2","Source"],
                   value_vars=valcols).reset_index()
    scounts["vclass"] = [x[:x.find("_")] for x in scounts["variable"].tolist()]
    scounts["vtype"] = [x[x.find("_")+1:] for x in scounts["variable"].tolist()]

    print scounts.head()


    scounts_sub = scounts[scounts.vclass.isin(["Class1","Class2","Class3","Class4"])]
    print "Scounts after filtering"
    scounts_sub["Burden"] = scounts_sub["value"].map(int)
    print scounts_sub.head(10)

    r_dataframe = com.convert_to_r_dataframe(scounts_sub)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor("+annotcol+")", y="Burden" )+ #, fill="factor(Region)") +
                ggplot2.geom_boxplot() +
                ggplot2.ggtitle("Distribution of Variants by Variant Class") +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.scale_y_continuous("Variant Burden") +
                ggplot2.scale_x_discrete(annotcol) +
                ggplot2.facet_grid( robjects.Formula('vclass ~ vtype') , scale="free") +
                ggplot2.theme(**mytheme) )

    figurename = "%s/%s_sampcounts2.png" % (outdir,basename)
    print "Making figure:",figurename
    grdevices.png(figurename)
    p.plot()
    grdevices.dev_off()

    hetcounts = scounts_sub[scounts_sub.vtype == "het"]
    pvals = calcPairwisePvalues( hetcounts, "vclass", "Burden" )
    pvals["x"] = tlevels[0]
    pvals["y"] = pvals["y"]*1.10
    print "$"*60
    print "Hetcounts"
    print hetcounts.head()

    assert tlevels is not None
    r_dataframe = com.convert_to_r_dataframe(hetcounts)
    r_dataframe = fixRLevels( r_dataframe, annotcol, tlevels )
    r_pvals = com.convert_to_r_dataframe(pvals)
    phet = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor("+annotcol+")", y="Burden")+#, color="factor(Region)")+ 
                ggplot2.geom_boxplot(ggplot2.aes_string(color="factor(Continent2)"), notch=True) +
                ggplot2.ggtitle("Heterozygous") +
                ggplot2.theme(**{ 'axis.title.x': ggplot2.element_blank(),
                                 'axis.title.y': ggplot2.element_blank()}) +
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", x="x", y="y"), 
                                  hjust=0, vjust=1, size=3, data = r_pvals ) +
                ggplot2.scale_y_continuous("Variant Burden") +
                ggplot2.scale_x_discrete(annotcol) +
                ggplot2.facet_grid( robjects.Formula('vclass ~ .') , scale="free") +
                ggplot2.theme(**pointtheme_nolegend) )

    figurename = "%s/%s_sampcounts2_het.png" % (outdir,basename)
    print "Making figure:",figurename
    grdevices.png(figurename)
    phet.plot()
    grdevices.dev_off()

    homcounts = scounts_sub[scounts_sub.vtype == "hom"]
    pvals = calcPairwisePvalues( homcounts, "vclass", "Burden" )
    pvals["x"] = tlevels[0]
    pvals["y"] = pvals["y"]*1.10

    r_dataframe = com.convert_to_r_dataframe(homcounts)
    r_dataframe = fixRLevels( r_dataframe, annotcol, tlevels )
    r_pvals = com.convert_to_r_dataframe(pvals)
    phom = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor("+annotcol+")", y="Burden") +
                ggplot2.geom_boxplot(ggplot2.aes_string(color="factor(Continent2)"), notch=True) +
                ggplot2.ggtitle("Homozygous") +
                ggplot2.theme(**{ 'axis.title.x': ggplot2.element_blank(),
                                 'axis.title.y': ggplot2.element_blank()}) +
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", x="x", y="y"), 
                                  hjust=0, vjust=1, size=3, data = r_pvals ) +
                #ggplot2.scale_y_continuous("Variant Burden") +
                #ggplot2.scale_x_discrete(annotcol) +
                ggplot2.facet_grid( robjects.Formula('vclass ~ .') , scale="free") +
                ggplot2.theme(**pointtheme_nolegend) )

    figurename = "%s/%s_sampcounts2_hom.png" % (outdir,basename)
    print "Making figure:",figurename
    grdevices.png(figurename)
    phom.plot()
    grdevices.dev_off()


    figurename = "%s/%s_sampcounts2_both.png" % (outdir,basename)
    print "Making figure:",figurename
    grdevices.png(figurename, width=6, height=9, units="in",res=300)
    gridextra.grid_arrange( *[phet,phom], ncol=2 )
    grdevices.dev_off()

    print "Starting LOF plots:"
    print scounts.head(10)
    #print scounts.vclass.unique()
    lof_het = scounts[(scounts.vclass.isin(["ClassNMD","ClassLOF"])) &
                      (scounts.vtype == "het")]
    assert len(lof_het) > 0
    print hk.dfTable(lof_het.vclass)
    lof_het["Burden"] = lof_het["value"].map(int)
    #print "$"*50
    #print lof_het.head()
    #print "$"*50

    #pvals = []
    #for vtype, subdat in lofcounts.groupby("vtype") :
        #tmpdf = calcPairwisePvalues( subdat, "vclass", "Burden" )
        #tmpdf["vtype"] = vtype
        #pvals.append( tmpdf )
    #pvals = concat( pvals ).reset_index(drop=True)
    pvals = calcPairwisePvalues( lof_het, "vclass", "Burden" )
    pvals["x"] = tlevels[0]
    pvals["y"] = pvals["y"]*1.10
    print pvals

    r_dataframe = com.convert_to_r_dataframe(lof_het)
    r_dataframe = fixRLevels( r_dataframe, annotcol, tlevels )
    r_pvals = com.convert_to_r_dataframe(pvals)
    phet = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor("+annotcol+")", y="Burden") +
                ggplot2.geom_boxplot(ggplot2.aes_string(color="factor(Continent2)"), 
                                     notch=True) +
                ggplot2.ggtitle("Heterozygous") +
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", x="x", y="y"), 
                                  hjust=0, vjust=1, size=3, data = r_pvals ) +
                ggplot2.theme(**{ 'axis.title.x': ggplot2.element_blank(),
                                 'axis.title.y': ggplot2.element_blank()}) +
                #ggplot2.scale_y_continuous("Variant Burden" ) +#, limits=robjects.IntVector((0,50))) +
                #ggplot2.scale_x_discrete(annotcol) +
                ggplot2.facet_grid( robjects.Formula('vclass ~ .'), scale="free") +
                ggplot2.theme(**pointtheme_nolegend) )
    figurename = "%s/%s_lofcounts_het.png" % (outdir,basename)
    print "Making figure:",figurename
    grdevices.png(figurename)
    phet.plot()
    grdevices.dev_off()


    lof_hom = scounts[(scounts.vclass.isin(["ClassNMD","ClassLOF"])) & 
                     (scounts.vtype == "hom")]
    print hk.dfTable(lof_hom.vclass)
    lof_hom["Burden"] = lof_hom["value"].map(int)

    pvals = calcPairwisePvalues( lof_hom, "vclass", "Burden" )
    pvals["x"] = tlevels[0]
    pvals["y"] = pvals["y"]*1.10
    print pvals

    r_dataframe = com.convert_to_r_dataframe(lof_hom)
    r_dataframe = fixRLevels( r_dataframe, annotcol, tlevels )
    r_pvals = com.convert_to_r_dataframe(pvals)
    phom = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor("+annotcol+")", y="Burden") +
                ggplot2.geom_boxplot(ggplot2.aes_string(color="factor(Continent2)"), 
                                     notch=True) +
                ggplot2.ggtitle("Homozygous") +
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", x="x", y="y"), 
                                  hjust=0, vjust=1, size=3, data = r_pvals ) +
                ggplot2.theme(**{ 'axis.title.x': ggplot2.element_blank(),
                                 'axis.title.y': ggplot2.element_blank()}) +
                ggplot2.facet_grid( robjects.Formula('vclass ~ .'), scale="free") +
                ggplot2.theme(**pointtheme_nolegend) )
    figurename = "%s/%s_lofcounts_hom.png" % (outdir,basename)
    print "Making figure:",figurename
    grdevices.png(figurename)
    phom.plot()
    grdevices.dev_off()

    figurename = "%s/%s_lofcounts.png" % (outdir,basename)
    print "Making figure:",figurename
    grdevices.png(figurename, width=6, height=5, units="in",res=300)
    gridextra.grid_arrange( *[phet,phom], ncol=2 )
    grdevices.dev_off()
    sys.exit(1)
# END plotSampleCounts2

################################################################################
def basicClean( vcffile, rerun=False ):
    print "Running basicClean"
    patientdata = sampleAnnotation()
    filepath, basename, suffix = hk.getBasename(vcffile)

    cleanbase = "%s/%s"%(filepath,basename)
    recodevcf = cleanbase+".recode.vcf"
    tvcf = cleanbase+".sfilt.vcf"
    tvcfgz = cleanbase+".sfilt.vcf.gz"
    if os.path.exists( tvcfgz ) and not rerun :
        return tvcfgz

    # Annotation Outliers
    vcfpats = patientInfo.getPats( vcffile )
    vcfpats = DataFrame({'Individual.ID':vcfpats})
    vcfpats = addSampleAnnotation( vcfpats, mergecol="Individual.ID" )
    sampleannot = vcfpats[(vcfpats.Continent2.notnull()
                              & (vcfpats.Continent2 != "Unknown"))]

    #keepfile = os.path.join(filepath, basename+".annotkeeppats")
    #sampleannot[["Individual.ID"]].to_csv(keepfile, header=None,index=None)
    annotpats = sampleannot["Individual.ID"].tolist()
    
    toremove="./toremove.txt"
    autoremove = []
    if os.path.exists(toremove):
        autoremove = ([x.strip() for x in open(toremove).readlines()
                       if len(x.strip()) > 0])
    
    qcoutliers = king.popOutliers( vcffile )

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
               "--keep",keepfile,
               "--min-alleles","2","--max-alleles","2",
               "--recode","--out",cleanbase]

    #print " ".join(command)
    out = subprocess.check_output( command )

    assert os.path.exists(recodevcf)

    #bedfile = pop.run_plink_convert(tped, force=rerun)
    #bedfile = pop.run_plink_filter(tped, force=rerun)
    subprocess.call(['mv', recodevcf, tvcf])
    subprocess.call(["bgzip",tvcf])
    subprocess.call(["tabix","-p","vcf",tvcfgz])

    return tvcfgz
# END seriousClean

################################################################################
# runClassifyWorkflow
################################################################################
def runClassifyWorkflow( vcffile, sampleannot, regionlist, figuredir ):
    print "Running runClassifyWorkflow"
    targetvcf = hk.copyToSubDir( vcffile, "classify" )
    # filter variants to remove poor quality
    cleanvcf = basicClean( targetvcf, rerun=False )

    filepats = patientInfo.currentFilePats( cleanvcf )
    smallvcf = makeSmallVCF(cleanvcf,figuredir=figuredir, force=False )
    allvcffiles = annotateVCF( smallvcf, force=False )
    print allvcffiles

    geneannotfile = classifyVars( allvcffiles["scores"], cleanvcf, 
                                 sampleannot, regionlist, filepats,
                                 "./results/classes", force=True )
    samplecounts = individualVarStats( cleanvcf, geneannotfile, 
                              sampleannot, regionlist, filepats, force=True)
    #plotSampleCounts( samplecounts, outdir=figuredir )
    plotSampleCounts2( samplecounts, outdir=figuredir )
    if cleanvcf.find("meceu") >= 0 : 
        regionlist = ["Europe","Middle East"]
    elif cleanvcf.find( "1000G" ) >= 0 :
        regionlist = ["Europe", "Middle East", "Africa"]
                   
    plotHomozygosity( geneannotfile, regionlist, outdir=figuredir )
    plotClassDist( geneannotfile, regionlist, outdir=figuredir )
# END runClassifyWorkflow

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

    hk.changeLogFile( hk.LOGDIR+"/classifyVars_test.log", True )

    print "Using dataset:", dataset
    path = os.path.abspath("./rawdata/")
    if dataset == "test" :
        vcffile = path+"/test/everything_set1.chr1.snp.clean.vcf.gz"
    elif dataset == "onekg" :
        vcffile = path+"/onekg/onekg.clean.vcf.gz"
    elif dataset == "test2" :
        vcffile = path+"/test2/main/test2.clean.vcf.gz"
        #vcffile = path+"/test2/test2.clean.vcf.gz"
    elif dataset == "daily" :
        vcffile = path+"/daily/daily.clean.vcf.gz"
    elif dataset == "merge1kg" :
        vcffile = path+"/merge1kg/main/me1000G.clean.vcf.gz"
    elif dataset == "mevariome" :
        vcffile = path+"/mevariome/variome.vcf.gz"
        #vcffile = path+"/mevariome/main/variome.clean.vcf.gz"
    #elif dataset == "variome1" :
        #vcffile = path+"/variome1/variome.clean.vcf.gz"
        #vcffile = path+"/variome1/variome.vcf.gz"
    #elif dataset == "variome" :
        #vcffile = path+"/variome/variome.clean.vcf.gz"
        #vcffile = path+"/variome1/variome.vcf.gz"
    elif dataset == "mergedaly" :
        #vcffile = path+"/mergedaly/main/meceu.clean.vcf.gz"
        vcffile = path+"/mergedaly/meceu.vcf.gz"
    elif dataset == "casanova" :
        vcffile = path+"/casanova/casanova.snp.recal.clean.vcf.gz"
    #elif dataset == "CEU" :
        #vcffile = path+"/fixed/daily.clean.Europe.vcf.gz"
    #elif dataset == "ME" :
        #vcffile = path+"/fixed/variome.clean.Middle_East.vcf.gz"
    #elif not optlist.has_key("-t") :
    else :
        print "Error: unknown dataset",dataset
        sys.exit(1)

    filepats = patientInfo.currentFilePats( vcffile )
    sampleannot = sampleAnnotation( filepats )
    regionlist = sorted(sampleannot[sampleannot.Continent2.notnull()]
                        .Continent2.unique().tolist())

    figuredir = "./results/figures/classes"

    if optlist.has_key("-t") : 
        print "Running test"
        #geneannotfile = "/media/data/workspace/variome/rawdata/test/everything_set1.chr1.snp_genes.tsv"
        geneannotfile = "/media/data/workspace/variome/rawdata/merged/merged.clean_genes.tsv"
        geneannotfile = "/media/data/workspace/variome/rawdata/mergedaly/main/meceu.clean_genes.tsv"
        #plotClassDist( geneannotfile, regionlist, outdir=figuredir )
    
        plotHomozygosity( geneannotfile, regionlist, outdir=figuredir )
    elif dataset is not None :
        runClassifyWorkflow( vcffile, sampleannot, regionlist, figuredir=figuredir)
    #allPatVarStats( outfiles, outdir="./classes" )

# END Main
################################################################################
    #targetfile = "./rawdata/onekg/onekg.chimp.regions.recode.vcf.gz"
    #targetfile = "./rawdata/variome1/variome.chimp.regions.filt.samp.samp.vcf.gz"
    #targetfile = "./rawdata/test/everything_set1.chr1.snp.chimp.regions.filt.samp.samp.vcf.gz"
    #targetfile = "./rawdata/merged/roh/merged.clean.vcf.gz"
    #targetfile = "./rawdata/daily/roh/daily.clean.vcf.gz"
    #targetfile = "./rawdata/casanova/roh/casanova.snp.recal.clean.vcf.gz"
    #targetfile = "./rawdata/variome1/roh/variome.clean.vcf.gz"
    #targetfile = "./rawdata/test/roh/everything_set1.chr1.snp.clean.vcf.gz"


