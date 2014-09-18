#!/usr/bin/env Rscript


library(RColorBrewer)
#palette(brewer.pal(8, "Dark2"))
library(ggplot2)
library(reshape)
library(VennDiagram)


file1 <- "/media/data/workspace/variome/rawdata/mevariome/main/lof/variome.clean_config1_LoF_genes.tsv"
file2 <- "/media/data/workspace/variome/rawdata/mevariome/main/lof/variome.clean_config2_LoF_genes.tsv"
file3 <- "/media/data/workspace/variome/rawdata/mevariome/main/lof/variome.clean_config3_LoF_genes.tsv"

ginfo1 <- read.csv(file1, sep="\t", comment.char="");
ginfo2 <- read.csv(file2, sep="\t", comment.char="");
ginfo3 <- read.csv(file3, sep="\t", comment.char="");

#genedata <- rbind( c())

venn.plot <- venn.diagram(
          list("All LoF" = ginfo1$Gene, "Hom LoF" = ginfo2$Gene, ">1 Hom LoF" = ginfo3$Gene),
          fill = brewer.pal(3,"Set1"),
          #main="Loss of Function Gene Proportions",
          #sub="Middle Easterns Versus Europeans",
          main.cex=3,sub.cex=2.3,cat.cex=2.1,
          alpha = c(0.5, 0.5, 0.5), 
          cex = 2, cat.fontface = 4,lty =2, fontfamily =3,
          filename=NULL)
          #filename = "figures/kos/chets.png");
png("figures/lof_class_proportions.png");
grid.draw(venn.plot);
dev.off();

