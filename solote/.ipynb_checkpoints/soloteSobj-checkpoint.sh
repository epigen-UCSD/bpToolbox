#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --account=csd772
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH -e /tscc/projects/ps-epigen/users/kdang/solote/script_eo/aging_solote2.e
#SBATCH -o /tscc/projects/ps-epigen/users/kdang/solote/script_eo/aging_solote2.o
#SBATCH --mail-type END
#SBATCH --mail-user k8dang@health.ucsd.edu

# include these 2 lines since bash commands aren't recognized outside my base environment (subshells)
source /tscc/nfs/home/k1dang/miniconda3/etc/profile.d/conda.sh
conda activate seuratv5

echo "Starting SoloTE Report:"

### USER SPECIFY ###
OUTPUT_DIRECTORY_PATH="/tscc/projects/ps-epigen/users/kdang/solote/Dp16_aging/" # directory where all solote output directories are; each solote odir should have dset as title
echo "OUTPUT_DIRECTORY_PATH: ${OUTPUT_DIRECTORY_PATH}"
SOLOTE_MTX="subfamilytes" # options: classtes, familytes, subfamilytes, legacytes, locustes ; BUT recommend legacy or locus
echo "SOLOTE_MTX: ${SOLOTE_MTX}"
SOBJ="/tscc/projects/ps-epigen/users/kdang/ASAP/rds_objs/aging_merged_processed_v3_wwn_nodoub_v2_batch_cpeaks_batch_finalishWNN_frozenVersion.rds"
echo "SOBJ: ${SOBJ}"

# please note - this can either one or multiple metacols
# ex: ADDMETACOL=("celltype") or ADDMETACOL=("celltype" "conditions" "age")
# if list; must sep by spaces, NO COMMAS
ADDMETACOL=("sub_celltypes")
echo "ADDMETACOL: ${ADDMETACOL[@]}"

PROJECTNAME="Aging_SoloTE_lognorm_subfamilytes" # also used as prefix for solote rds name
echo "PROJECTNAME: ${PROJECTNAME}"
UMAPREDUCTION="umap.hrm.wnn"
echo "UMAPREDUCTION: ${UMAPREDUCTION}"
echo "-----------------------"
### DONE ###

cd $OUTPUT_DIRECTORY_PATH

# this will merge all solote dsets into one solote sobj
# also matches barcodes to reference sobj and copies over umap coords and metacol(s)
Rscript /tscc/projects/ps-epigen/users/kdang/solote/solote/solote_to_seurat.R -mtx $SOLOTE_MTX -obj $SOBJ -col "${ADDMETACOL[@]}" -proj $PROJECTNAME -odir $OUTPUT_DIRECTORY_PATH -umap $UMAPREDUCTION
