### To bedtools merge
# 1. concat bed files
# 2. sort bed filed
# 3. bedtools merge 


IDIR='ac_callpeak_out/'
ODIR='ac_bedmerge_out/'
NAME='ac_unionNarrow'

mkdir -p $ODIR

FILES=$(find "$IDIR" -type f -name "*narrow*" -readable)
echo $FILES


cat $FILES > ${ODIR}${NAME}_concat.bed
sort -k1,1 -k2,2n ${ODIR}${NAME}_concat.bed > ${ODIR}${NAME}_concat_sorted.bed

bedtools merge -i ${ODIR}${NAME}_concat_sorted.bed > ${ODIR}${NAME}_merged.bed


