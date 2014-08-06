#!/usr/bin/env python

#-------------------------------------------------------------------#
#                       Import Packages                             #
#-------------------------------------------------------------------#
import os, sys
#import getopt
#import csv
#import re
#import subprocess

#import housekeeping as hk
#import popgencommands as pop
from localglobals import *


if __name__ == "__main__" :
    #incidencefile = "/media/data/workspace/variome/resources/disincidence/noncomunicable_daly.txt"
    countryfile = "/media/data/workspace/variome/resources/disincidence/countrycodes.txt"
    countryinfo = read_csv( countryfile, sep="\t" )
    countryinfo["Population_int"] = [int(x.replace(",","")) 
                                     for x in countryinfo.Population.tolist()]

    incidencefile = "/media/data/workspace/variome/resources/disincidence/noncomunicable_yld.txt"
    yldincidence = read_csv(incidencefile, sep="\t" )
    yld_melt = melt(yldincidence, id_vars="Disease" ) 
    ylddata = merge( yld_melt, countryinfo, left_on="variable", right_on="Country" )
    ylddata["Stat"] = "YLD"

    incidencefile = "/media/data/workspace/variome/resources/disincidence/noncomunicable_yll.txt"
    yllincidence = read_csv(incidencefile, sep="\t" )
    yll_melt = melt(yllincidence, id_vars="Disease" ) 
    ylldata = merge( yll_melt, countryinfo, left_on="variable", right_on="Country" )
    ylldata["Stat"] = "YLL"

    alldat = concat([ylddata, ylldata])
    alldat_filt = alldat[(alldat.Region != "Unknown") & (alldat.value != ".") & (alldat.value != "..")]

    dalydat = []
    for idx, regiondata in  alldat_filt.groupby(["Disease","Region","Stat"]) :
        dis, region, stat = idx
        print dis, region
        print regiondata.head()
        #print regiondata.Population_int.sum()
        #print regiondata.value.sum()
        meandaly = regiondata.value.map(float).sum()/regiondata.Population_int.map(int).sum()
        dalydat.append( Series({'Disease':dis, 'Region':region, 'aveDaly':meandaly, 'Stat':stat}) )
    dalydat = DataFrame( dalydat )

    tdiseases = [ 'Congenital anomalies',
                 #'Neurological conditions',
                 #'Mental and behavioral disorders',
                 #'Skin diseases'
                 #'Respiratory diseases',
                 #'Diabetes mellitus', 
                 #'Digestive diseases',
                 #'Endocrine, blood, immune disorders', 
                 #'Malignant neoplasms', 
                 #'Musculoskeletal diseases', 
                 #'Genitourinary diseases',
                 #'Cardiovascular diseases', 
                 #'Sense organ diseases', 
                 #'Oral conditions', 
                ]
                 #'Other neoplasms', 

    tregions = [
        'Subsaharan Africa',
        'Northwest Africa',
        'Northeast Africa',
        'Syrian Desert',
        'Ariabian Peninsula',
        'Turkish Peninsula',
        'Central Asia',
        'Europe',
        'East Asia',
        'North America',
        'South America',
        'Oceania' ]
        #'South Asia',

    dalydat_filt = dalydat[ dalydat.Disease.isin( tdiseases ) ]

    r_dataframe = com.convert_to_r_dataframe(dalydat_filt)
    r_dataframe = fixRLevels( r_dataframe, "Region", tregions )
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="factor(Region)",
                                   y="aveDaly",
                                   fill="factor(Stat)") +
                ggplot2.geom_bar(stat="identity") +
                ggplot2.ggtitle("Disability-Adjusted Life Years (DALYs) Lost\nDue to Congenital Abnormalities") +
                ggplot2.scale_x_discrete("Geographic Region") +
                ggplot2.scale_y_continuous("DALYs Lost per 1000 people") +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 90)}) +
                ggplot2.coord_flip() +
                ggplot2.theme(**mytheme) )
                #ggplot2.stat_smooth(method="lm", se=False) +
                ##ggplot2.facet_grid(robjects.Formula('. ~ variable')) +

    figname = "../results/figures/disincidence/several_daly.png"
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()

# END MAIN

