#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00 # USER might need to play around with; i had one dset run for ~14hrs; avg usually 3-4
#SBATCH --account=csd772
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH -e /tscc/projects/ps-epigen/users/kdang/solote/script_eo_ALD/240417_ALD.e-%a
#SBATCH -o /tscc/projects/ps-epigen/users/kdang/solote/script_eo_ALD/240417_ALD.o-%a
#SBATCH --mail-type END
#SBATCH --mail-user k8dang@health.ucsd.edu

# include these 2 lines since bash commands aren't recognized outside my base environment (subshells)
source /tscc/nfs/home/k1dang/miniconda3/etc/profile.d/conda.sh
conda activate solote

echo "Starting SoloTE Pipeline:"

# 1. generating the TE annotation BED file for specific genome (mm10)
# note: already done for mm10 and hg38; no need to repeat until you want another genome
# note: run python /tscc/projects/ps-epigen/users/kdang/solote/SoloTE/SoloTE_RepeatMasker_to_BED.py -l to get all genomes available
#cd /tscc/projects/ps-epigen/users/kdang/solote/SoloTE
#python /tscc/projects/ps-epigen/users/kdang/solote/SoloTE/SoloTE_RepeatMasker_to_BED.py -g mm10

# 2. running SoloTE pipeline on each BAM file from CellRanger
### USER SPECIFY ###
INPUT_FILE="/tscc/projects/ps-epigen/users/kdang/solote/lung_ALD/inSoloTE_array.tsv"
echo "INPUT_FILE: ${INPUT_FILE}"
OUTPUT_DIRECTORY_PATH="/tscc/projects/ps-epigen/users/kdang/solote/lung_ALD/"
echo "OUTPUT_DIRECTORY_PATH: ${OUTPUT_DIRECTORY_PATH}"
NUM_THREAD=24
echo "NUM_THREAD: ${NUM_THREAD}"
GENOME="hg38" # hg38 or mm10
echo "GENOME: ${GENOME}"
echo "-----------------------"
### DONE ###

# grabbing matching genome RepeatMasker
GENOME_BED=""
if [ $GENOME = "mm10" ]
then
    GENOME_BED="/tscc/projects/ps-epigen/users/kdang/solote/RepeatMasker_BEDfiles/mm10_rmsk.bed"
else
    GENOME_BED="/tscc/projects/ps-epigen/users/kdang/solote/RepeatMasker_BEDfiles/hg38_rmsk.bed"
fi

# lets do batch arrays instead
SAMPLE_ROW=`cat $INPUT_FILE | sed -n ${SLURM_ARRAY_TASK_ID}p`
echo "Running..."
BAM_FILE=$(echo $SAMPLE_ROW | cut -f1 -d' ')
OUTPUT_PREFIX=$(echo $SAMPLE_ROW | cut -f2 -d' ')
echo $OUTPUT_PREFIX
OUTPUT_DIRECTORY="${OUTPUT_DIRECTORY_PATH}${OUTPUT_PREFIX}/"
echo $OUTPUT_DIRECTORY
mkdir $OUTPUT_DIRECTORY
cd $OUTPUT_DIRECTORY

python /tscc/projects/ps-epigen/users/kdang/solote/SoloTE/SoloTE_pipeline.py -t $NUM_THREAD -b $BAM_FILE -a $GENOME_BED -o $OUTPUT_PREFIX -d $OUTPUT_DIRECTORY

cd ..
echo "Done!"
