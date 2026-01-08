#!/usr/bin/env bash

GENOME="hg38"

BL="/tscc/projects/ps-epigen/GENOME/hg38/hg38.blacklist.bed.gz"
#BL="/tscc/projects/ps-epigen/GENOME/mm10/mm10.blacklist.bed.gz"

CS="/tscc/projects/ps-epigen/GENOME/hg38/hg38.chrom.sizes"
#CS="/tscc/projects/ps-epigen/GENOME/mm10/mm10.chrom.sizes"

OUT='/tscc/lustre/ddn/scratch/rlancione/mergepeak_out/'

NAME='multiome_unionNarrow'

Rscript /tscc/projects/ps-epigen/bpToolbox/snATAC/iterative_overlap_peak_merging.R -i my_summits.txt \
        -g $GENOME\
        --blacklist $BL\
        --chromSize $CS \
        -d $OUT -o $NAME
