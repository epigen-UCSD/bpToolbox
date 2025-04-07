#!/bin/bash
#SBATCH -J LTB_DA
#SBATCH -N 1                            #Nodes
#SBATCH -c 2                            #Cores
#SBATCH --mem=50G                       #Memory
#SBATCH -t 08:00:00
#SBATCH -o /tscc/nfs/home/rlancione/jeo/LTB_deseq.sh.o
#SBATCH -e /tscc/nfs/home/rlancione/jeo/LTB_deseq.sh.e
#SBATCH -p platinum
#SBATCH -q hcp-csd772
#SBATCH -A csd772
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rlancione@ucsd.edu

# Navigate to oasis scratch
#cd /tscc/lustre/scratch/rlancione/
# activate cond env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate renv43

# 1. Input merged object path
SEU_OBJ="/tscc/projects/ps-epigen/users/rlan/Liver_TB/PeakCalling/Objects/sobj_nfrags.rds"

# Record Input
# 2. pathway to matrices and metadata created from MatrixCreation.R script
MATRIX_PATH="/tscc/projects/ps-epigen/users/rlan/Liver_TB/DASeq/Matrices/"
# 3. adjusted p-value cutoff for EnrichR
PADJ_CUTOFF=0.05
# 4. log2FC threshold cutoff for EnrichR
# 0.5849625=log2(1.5)  1=log2(2)
LOG2FC_CUTOFF=0.5849625
# 5. path to output directory
OUTDIR='/tscc/projects/ps-epigen/users/rlan/Liver_TB/DASeq/Results/Res_batch/'
# 6. project name (must be the same from matrix creation script)
PROJ='LTB'
#7. Cell type column. Default seurat_clusters
CT_COL="CellType"
#8. Path to tsv containing comparisons between conditions to run DESeq on. If left null (commed out) all conditions will be run pairwise.
COMP_PATH="/tscc/projects/ps-epigen/users/rlan/Liver_TB/DASeq/comparisons.txt"
#9. Feature filter
FF=10
#10. Assay RNA or ATAC. Determines whether to run CP or Great 
ASSAY="ATAC"

# Path to script
Script='/tscc/projects/ps-epigen/bpToolbox/DESeq2-GO/DESeq2.R'
######## Run R Script #########
Rscript $Script -mpath $MATRIX_PATH -ap $PADJ_CUTOFF -lfc $LOG2FC_CUTOFF \
        -o $OUTDIR -proj $PROJ -comp $COMP_PATH -nf $FF


conda deactivate 
conda activate CProfiler
#### Now call report script 
Report='/tscc/projects/ps-epigen/bpToolbox/DESeq2-GO/DE_report.R'
######## Run R Script #########
Rscript $Report -rp $OUTDIR -mp $MATRIX_PATH -ap $PADJ_CUTOFF -lfc $LOG2FC_CUTOFF -org "human" \
	        -pn $PROJ -sp $SEU_OBJ -ct $CT_COL -a $ASSAY







