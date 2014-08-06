#!/bin/bash
cd ..

#rawdata/test/everything_set1.chr1.snp.chimp.regions.filt.samp.vcf.gz
#rawdata/test/everything_set1.chr1.snp.chimp.regions.filt.samp.samp.vcf.gz

#./classifyVars.py rawdata/daily/daily.chimp.regions.filt.samp.vcf.gz
#./classifyVars.py rawdata/daily/daily.chimp.regions.filt.samp.samp.vcf.gz
#./classifyVars.py rawdata/daily/daily.regions.filt.samp.samp.vcf.gz
#./classifyVars.py rawdata/daily/daily.regions.filt.samp.vcf.gz
#./classifyVars.py rawdata/merged/merged.regions.filt.samp.samp.vcf.gz
#./classifyVars.py rawdata/merged/merged.chimp.regions.filt.samp.samp.vcf.gz
#./classifyVars.py rawdata/merged/merged.regions.filt.samp.vcf.gz
#./classifyVars.py rawdata/merged/merged.chimp.regions.filt.samp.vcf.gz
#./classifyVars.py rawdata/onekg/onekg.regions.filt.samp.vcf.gz
#./classifyVars.py rawdata/onekg/onekg.chimp.regions.filt.samp.samp.vcf.gz
#./classifyVars.py rawdata/onekg/onekg.chimp.regions.filt.samp.vcf.gz
#./classifyVars.py rawdata/variome1/variome.chimp.regions.filt.samp.samp.vcf.gz
#./classifyVars.py rawdata/variome1/variome.chimp.regions.filt.samp.vcf.gz
#./classifyVars.py rawdata/variome1/variome.regions.filt.samp.vcf.gz
#./classifyVars.py rawdata/variome1/variome.regions.filt.samp.samp.vcf.gz

./classifyVars.py ./rawdata/daily/roh/daily.clean.vcf.gz
./classifyVars.py ./rawdata/merged/roh/merged.clean.vcf.gz
./classifyVars.py ./rawdata/onekg/roh/onekg.clean.vcf.gz
./classifyVars.py ./rawdata/test/roh/everything_set1.chr1.snp.clean.vcf.gz
./classifyVars.py ./rawdata/casanova/roh/casanova.snp.recal.clean.vcf.gz
./classifyVars.py ./rawdata/variome1/roh/variome.clean.vcf.gz

