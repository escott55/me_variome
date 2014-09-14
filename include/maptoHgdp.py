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
#import popgencommands as pop
from localglobals import *

vcffile = "./rawdata/mevariome/main/variome.clean.vcf.gz"

hgdpped = "/home/escott/workspace/variome/resources/hgdppanel4/panel4.ped"
hgdpmap = "/home/escott/workspace/variome/resources/hgdppanel4/panel4.map"
hgdpids = "/home/escott/workspace/variome/resources/hgdppanel4/sample_panel4.txt"

