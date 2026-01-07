#!/bin/bash
#SBATCH -J bmat_frag
#SBATCH -N 1                            #Nodes
#SBATCH -c 2                            #Cores
#SBATCH --mem=30G                      #Memory
#SBATCH -t 24:00:00
#SBATCH -o /tscc/nfs/home/rlancione/jeo/bmat_frags_ac.sh.o-%a
#SBATCH -e /tscc/nfs/home/rlancione/jeo/bmat_frags_ac.sh.e-%a
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd772
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rlancione@ucsd.edu

# Navigate to oasis scratch
cd /tscc/lustre/ddn/scratch/rlancione/

# activate cond env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate epscan

### User Vars 
MANIFEST='/tscc/nfs/home/rlancione/ps-epigen/users/rlan/LFNIH_Pt2/Pairtag/Build_ac_mat/manifest.txt'
REGIONS='/tscc/projects/ps-epigen/users/rlan/LFNIH_Pt2/Pairtag/Callpeak_ac/ac_bedmerge_out/ac_unionNarrow_merged.bed'
OUTDIR='/tscc/nfs/home/rlancione/ps-epigen/users/rlan/LFNIH_Pt2/Pairtag/Build_ac_mat/sample_out/'
METADATA='/tscc/projects/ps-epigen/users/rlan/LFNIH_Pt2/Pairtag/Build_ac_mat/corr_ac_meta.csv'

mkdir -p $OUTDIR

# extract path from manifest 
FRAG=`cat $MANIFEST | sed -n ${SLURM_ARRAY_TASK_ID}p | cut -f1 -d' '`

### Run bld_mat.py for each frag file
PYSCRIPT='/tscc/projects/ps-epigen/bpToolbox/snATAC/bld_mat.py'
python3 $PYSCRIPT -f $FRAG -r $REGIONS -o $OUTDIR -m $METADATA 



