#!/bin/bash
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -J cistarget_array
#SBATCH -N 3
#SBATCH -c 8
#SBATCH --mem 400G
#SBATCH -t 12:00:00
#SBATCH -o /tscc/nfs/home/rlancione/jeo//ct_array.out_%a
#SBATCH -e /tscc/nfs/home/rlancione/jeo//ct_array.err_%a
#SBATCH --mail-user rlancione@ucsd.edu
#SBATCH --mail-type FAIL
#SBATCH -A csd772

set -e
source ~/.bashrc
conda activate scenicplus

FASTA_FILE="[FASTA_PATH]"
OUT_DIR="[OUT_PATH]"
DATABASE_PREFIX="[PROJECT]"

CBDIR="/tscc/nfs/home/rlancione/misc/aertslab_motif_colleciton/v10nr_clust_public/singletons"
MOTIF_LIST="/tscc/nfs/home/rlancione/ps-epigen/users/rlan/Preissel_Multiome/Scenic/V1/cistarget_db/motifs.txt"


SCRIPT_DIR="/tscc/nfs/home/rlancione/software/create_cisTarget_databases/"
${SCRIPT_DIR}/create_cistarget_motif_databases.py \
	-f ${FASTA_FILE} \
	-M ${CBDIR} \
	-m ${MOTIF_LIST} \
	-o ${OUT_DIR}${DATABASE_PREFIX} \
	-p ${SLURM_ARRAY_TASK_ID} 8 \
	--bgpadding 1000 \
	-t 20
