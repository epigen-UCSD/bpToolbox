#!/bin/sh
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --account=csd772
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH -e /tscc/projects/ps-epigen/users/kdang/solote/script_eo/260324_subfamBulkMtx.e
#SBATCH -o /tscc/projects/ps-epigen/users/kdang/solote/script_eo/260324_subfamBulkMtx.o
#SBATCH --mail-type END
#SBATCH --mail-user k8dang@health.ucsd.edu

source /tscc/nfs/home/k1dang/miniconda3/etc/profile.d/conda.sh
conda activate seuratv5

### USER SPECIFY ###
### shared and unique vars for bulk mtxs and deseq scripts
# solote sobj path
SEU_OBJ="/tscc/projects/ps-epigen/users/kdang/solote/Dp16_aging/Aging_SoloTE_lognorm_subfamilytes.rds"
# Dataset column name in object metadata. E.G. "orig.ident"
DSET_COL="orig.ident"
COND_COL="Ref_conditions"
CT_COL="Ref_sub_celltypes_v2"
#Output path - bulk mtx directory
BM_ODIR="/tscc/projects/ps-epigen/users/kdang/solote/Dp16_aging/DESeq_OUT/bulkMtx/"
# Project: To name files
PROJ="Aging_SoloTE"
# Extra covariate columns for meta data. (optional)
#CV1="genotype"
#CV2="age"
#CV3="TOD"

# path to DESeq2 output directory
DE_ODIR='/tscc/projects/ps-epigen/users/kdang/solote/Dp16_aging/DESeq_SoloTE.sobj/' # usually DESeq_out
# path to tsv containing comparisons between conditions to run DESeq on. If left null (commed out) all conditions will be run pairwise.
COMP_PATH="/tscc/projects/ps-epigen/users/kdang/ASAP/DESeq/comparisons.tsv"
# formula used for DESeq - used to be manually written in Rscript command below - but since i use it in both Script and Report make a var for it
FORM="~ condition"
# adjusted p-value cutoff for countsummary, reports, and hmaps
PADJ_CUTOFF=0.05
# 4. log2FC threshold cutoff for countsummary, reports, and hmaps
# .58=log2(1.5)  1=log2(2)
LOG2FC_CUTOFF=0.5849625
### DONE ###

echo "Printing variables..."
echo "SEU_OBJ: $SEU_OBJ"
echo "DSET_COL: $DSET_COL"
echo "COND_COL: $COND_COL"
echo "CT_COL: $CT_COL"
echo "BULK.MTX.ODIR: $BM_ODIR"
echo "PROJ: $PROJ"
echo "DESEQ.ODIR: $DE_ODIR"
echo "COMP_PATH: $COMP_PATH"
echo "FORM: $FORM"
echo "-----------------------"

# making bulk mtxs - has both genes and TEs
echo "Making Bulk Mtxs..."
BM_Script='/tscc/projects/ps-epigen/users/kdang/solote/solote/solote_bulk_mat_create.R' 
#Rscript $Script -sobj $SEU_OBJ -dcol $DSET_COL -ccol $COND_COL -ctcol $CT_COL --project_name $PROJ -o $BM_ODIR \
#	-cov1 $CV1 -cov2 $CV2 -cov3 $CV3 
Rscript $BM_Script -sobj $SEU_OBJ -dcol $DSET_COL -ccol $COND_COL -ctcol $CT_COL --project_name $PROJ -o $BM_ODIR

echo "Running DESeq2 Pipeline..."
DE_Script='/tscc/projects/ps-epigen/users/kdang/solote/solote/solote_DEseq2_i10.R'
Rscript $DE_Script -sobj $SEU_OBJ -ctcol $CT_COL -mpath $BM_ODIR -o $DE_ODIR -proj $PROJ -comp $COMP_PATH -nf 10 -form $FORM

echo "Creating DE.reports..."
RE_Script='/tscc/projects/ps-epigen/users/kdang/solote/solote/solote_DE_report_12.R'
Rscript $RE_Script -sp $SEU_OBJ -rp $DE_ODIR -form $FORM -mp $BM_ODIR -ap $PADJ_CUTOFF -lfc $LOG2FC_CUTOFF \
    -pn $PROJ -ct $CT_COL

