#!/bin/bash
#SBATCH -J splitfrags 
#SBATCH -N 1                            #Nodes
#SBATCH -c 3                           #Cores
#SBATCH --mem=30G                       #Memory
#SBATCH -t 48:00:00
#SBATCH -o /tscc/projects/ps-epigen/bpToolbox/Callpeak_to_Matrix/logs/splitfrags.o
#SBATCH -e /tscc/projects/ps-epigen/bpToolbox/Callpeak_to_Matrix/logs/splitfrags.e
#SBATCH -p platinum
#SBATCH -q hcp-csd772
#SBATCH -A csd772
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rlancione@ucsd.edu

# Navigate to oasis scratch
cd /tscc/lustre/ddn/scratch/rlancione/

# activate cond env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate renv43

FRAG_PATHS='/tscc/projects/ps-epigen/bpToolbox/Callpeak_to_Matrix//frags.txt'
CT_COL='celltype'
META='/tscc/projects/ps-epigen/users/rlan/LFNIH_Pt2/Subclustering/SCVI/multi_metadata.csv'
OUT='/tscc/lustre/ddn/scratch/rlancione/splitfrags_out/'


### Run python splitfrags support 
python /tscc/projects/ps-epigen/bpToolbox/Callpeak_to_Matrix/support/splitfrags_support.py \
	--paths_txt $FRAG_PATHS \
	--groups $CT_COL \
	--metadata $META \
	--out $OUT 	



