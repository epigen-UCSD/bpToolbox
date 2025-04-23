#!/usr/bin/env bash

GENOME="hg38"

BL="/tscc/projects/ps-epigen/GENOME/hg38/hg38.blacklist.bed.gz"
#BL="/tscc/projects/ps-epigen/GENOME/mm10/mm10.blacklist.bed.gz"

CS="/tscc/projects/ps-epigen/GENOME/hg38/hg38.chrom.sizes"
#CS="/tscc/projects/ps-epigen/GENOME/mm10/mm10.chrom.sizes"


Rscript iterative_overlap_peak_merging.R -i summits.txt \
        -g $GENOME\
        --blacklist $BL\
        --chromSize $CS \
        -d merge_out/ -o LFNIH
