#!/bin/bash
#SBATCH -J splitfrags 
#SBATCH -N 1                            #Nodes
#SBATCH -c [CORES]                           #Cores
#SBATCH --mem=[MEM]G                       #Memory
#SBATCH -t 48:00:00
#SBATCH -o [WDIR]logs/splitfrags.o
#SBATCH -e [WDIR]logs/splitfrags.e
#SBATCH -p platinum
#SBATCH -q hcp-csd772
#SBATCH -A csd772
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=[USR]@ucsd.edu

# Navigate to oasis scratch
cd /tscc/lustre/ddn/scratch/rlancione/

# activate cond env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate renv43

FRAG_PATHS='[WDIR]/frags.txt'
CT_COL='[CTYPE_COL]'
META='[META_PATH]'
OUT='[OUTDIR]'


### Run python splitfrags support 
python [WDIR]support/splitfrags_support.py \
	--paths_txt $FRAG_PATHS \
	--groups $CT_COL \
	--metadata $META \
	--out $OUT 	



