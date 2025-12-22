#!/bin/bash
#SBATCH -J splitfrags 
#SBATCH -N [NODES]                            #Nodes
#SBATCH -c [CORES]                           #Cores
#SBATCH --mem=[MEM]                       #Memory
#SBATCH -t 48:00:00
#SBATCH -o [WDIR]/logs/splitfrags.o
#SBATCH -e [WDIR]/logs/splitfrags.e
#SBATCH -p platinum
#SBATCH -q hcp-csd772
#SBATCH -A csd772
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=[USR]@ucsd.edu

# Navigate to oasis scratch
cd /tscc/lustre/ddn/scratch/rlancione/

# activate cond env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate renv4


### args.. in-frags, ct col, meta 
### still need argparse 
python ../support/splitfrags_support.py 


