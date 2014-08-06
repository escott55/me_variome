#!/bin/bash
cd /home/escott/workspace/variome/rscripts/

#rawdata/test/everything_set1.chr1.snp.chimp.regions.filt.samp.vcf.gz
#rawdata/test/everything_set1.chr1.snp.chimp.regions.filt.samp.samp.vcf.gz

./admix.R rawdata/daily/admixture/daily.chimp.regions.filt.samp.plink.mod.filt.uniq.fam
./admix.R rawdata/merged/admixture/merged.chimp.regions.filt.samp.plink.mod.filt.uniq.fam
./admix.R rawdata/onekg/admixture/onekg.chimp.regions.filt.samp.plink.mod.filt.uniq.fam
./admix.R rawdata/test/admixture/everything_set1.chr1.snp.chimp.regions.filt.samp.plink.mod.filt.uniq.fam
./admix.R rawdata/casanova/admixture/casanova.snp.recal.chimp.regions.filt.samp.plink.mod.filt.uniq.fam
./admix.R rawdata/variome1/admixture/variome.chimp.regions.filt.samp.plink.mod.filt.uniq.fam

