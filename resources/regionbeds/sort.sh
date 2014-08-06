#!/bin/bash

awk '$1 !~ /_/ {sub("X","23",$1) ;sub("Y","24",$1); print}' refseqgenes_b37.bed | sort -k1,1n -k2,2n > refseqgenes_b37_sorted.bed

