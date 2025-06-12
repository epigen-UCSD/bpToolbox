#!/bin/bash

# Run with direct resources 5cpu, 100G 
# activate cond env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate scenicplus

# Template replacements
	# 1. PROJECT
	# 2. BED_PATH
	# 3. GENOME

REGION_BED="[BED_PATH]"
ORG="[GENOME]"

if [ "$ORG" == "hg38" ]
then
	GENOME_FASTA="/tscc/projects/ps-epigen/GENOME/hg38/genome.fa"
	CHROMSIZES="/tscc/projects/ps-epigen/GENOME/hg38/hg38.chrom.sizes"
elif [ "$ORG" == "mm10" ]
then
	GENOME_FASTA="/tscc/projects/ps-epigen/GENOME/mm10/refdata-cellranger-atac-mm10-1.2.0/fasta/genome.fa"
	CHROMSIZES="/tscc/projects/ps-epigen/GENOME/mm10/mm10.chrom.sizes"
fi


SCRIPT_DIR="/tscc/nfs/home/rlancione/software/create_cisTarget_databases"
${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        consensus.fa\
        1000 \
        yes



