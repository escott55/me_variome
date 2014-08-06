#!/bin/bash
cd ..

#rawdata/test/everything_set1.chr1.snp.chimp.regions.filt.samp.vcf.gz
#rawdata/test/everything_set1.chr1.snp.chimp.regions.filt.samp.samp.vcf.gz

./distances.py rawdata/daily/daily.chimp.regions.filt.samp.vcf.gz
./distances.py rawdata/daily/daily.chimp.regions.filt.samp.samp.vcf.gz
./distances.py rawdata/daily/daily.regions.filt.samp.samp.vcf.gz
./distances.py rawdata/daily/daily.regions.filt.samp.vcf.gz
./distances.py rawdata/merged/merged.regions.filt.samp.samp.vcf.gz
./distances.py rawdata/merged/merged.chimp.regions.filt.samp.samp.vcf.gz
./distances.py rawdata/merged/merged.regions.filt.samp.vcf.gz
./distances.py rawdata/merged/merged.chimp.regions.filt.samp.vcf.gz
./distances.py rawdata/onekg/onekg.regions.filt.samp.vcf.gz
./distances.py rawdata/onekg/onekg.chimp.regions.filt.samp.samp.vcf.gz
./distances.py rawdata/onekg/onekg.chimp.regions.filt.samp.vcf.gz
#./distances.py rawdata/variome1/variome.chimp.regions.filt.samp.samp.vcf.gz
#./distances.py rawdata/variome1/variome.chimp.regions.filt.samp.vcf.gz
#./distances.py rawdata/variome1/variome.regions.filt.samp.vcf.gz
#./distances.py rawdata/variome1/variome.regions.filt.samp.samp.vcf.gz

./distances.py rawdata/hgdp/HGDP_938.filt.samp.samp.vcf.gz


