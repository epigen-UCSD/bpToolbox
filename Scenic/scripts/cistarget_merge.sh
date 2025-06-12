#!/bin/bash
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -J cistarget_merge
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem 100G
#SBATCH -t 12:00:00
#SBATCH -o /tscc/nfs/home/rlancione/jeo//ct_merge.out
#SBATCH -e /tscc/nfs/home/rlancione/jeo//ct_merge.err
#SBATCH --mail-user rlancione@ucsd.edu
#SBATCH --mail-type FAIL
#SBATCH -A csd772

set -e
source ~/.bashrc
conda activate scenicplus

IDIR="[DATABASE_PATH]"
ODIR="[OUT_PATH]"
NAME="[PROJECT]"


SCRIPT_DIR="/tscc/nfs/home/rlancione/software/create_cisTarget_databases/"
${SCRIPT_DIR}/combine_partial_motifs_or_tracks_vs_regions_or_genes_scores_cistarget_dbs.py \
	-i $IDIR \
	-o $ODIR

${SCRIPT_DIR}//convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py \
	-i ${ODIR}${NAME}.motifs_vs_regions.scores.feather \
