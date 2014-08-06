#!/usr/bin/python

import os, sys
import subprocess

from localglobals import *
import housekeeping as hk
import popgencommands as pop

def retrieveAdmixfiles( vcffile ) :
    filepath, mainbasename, suffix = hk.getBasename( vcffile )
    print "p",filepath, "b",mainbasename, "s",suffix

    allfiles = []
    admixfiles = hk.filesInDir( filepath, "Q" )
    print admixfiles
    for afile in admixfiles : 
        cpath, cbase, csuf = hk.getBasename( afile )
        kval = cbase[len(mainbasename)+1:]
        print kval
        allfiles.append(Series({"file":afile, "K":kval,
                                "path":cpath,"cbase":cbase}))

    allfiles = DataFrame(allfiles)
    return allfiles
# END retrieveAdmixfiles

if __name__ == "__main__" :

    originalfile = "/home/escott/workspace/variome/rawdata/mevariome/main/admixture/variome.clean.vcf.gz"
    admixfiles = retrieveAdmixfiles( originalfile ) 

    print admixfiles
    
    filepats = pop.currentFilePats1( originalfile )
    for idx, afile in admixfiles.iterrows() :
        print idx, afile

    colnames = ["a"+str(x) for x in range(1,int(afile["K"])+1)]
    print afile["file"]
    qdata = read_csv( afile["file"], delim_whitespace=True, header=None,
                     names=colnames )
    qdata["sid"] = filepats
    qdata_melt = melt( qdata, id_vars=["sid"] )
    qdata_annot = addSampleAnnotation( qdata_melt, "sid" )
    #qdata.reset_index()
    print qdata_annot.head()

    sys.exit(1)
    r_dataframe = com.convert_to_r_dataframe(qdata_annot)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x="factor(Individual.ID)", y="value", 
                                   fill="factor(variable)" ) + \
                ggplot2.geom_bar( stat="identity" ) + \
                ggplot2.ggtitle("Admixture Distribution") + \
                ggplot2.theme(**mytheme) #+ \
                #ggplot2.scale_y_continuous("Internal RVIS") + \
                #ggplot2.stat_smooth(method="lm", se=False)
                #ggplot2.scale_x_continuous("External RVIS")+ \
    #barplot ="results/figures/%s_pops.png" %(prefix)
    barplot ="%s_pops.png" %(figureprefix)
    print "Writing file %s"%barplot 
    grdevices.png(barplot) 
    p.plot()
    grdevices.dev_off()

# END MAIN




