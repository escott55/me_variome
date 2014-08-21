#!/usr/bin/env Rscript

library(RColorBrewer)
palette(brewer.pal(8, "Dark2"))
library(ggplot2)
library(reshape)
library(VennDiagram)
library(HardyWeinberg)

source("themes.R")

######################################################################
# Main
######################################################################
setwd( "/home/escott/workspace/variome/" )

varfile <- "rawdata/mevariome/main/variome.clean_genes.tsv"



