SoloTE Pipeline README
https://github.com/bvaldebenitom/SoloTE/tree/main

General steps:
1. Generate SoloTE RepeatMasker BED files
    * I already generated them for mm10 and hg38, but there are more genomes available
    * located in "/tscc/projects/ps-epigen/users/kdang/solote/RepeatMasker_BEDfiles" directory
    * modules required: argparse, pandas, requests, os, colorama, & tqdm (just in case you don't already have them installed)

2. Run SoloTE pipeline
    * script takes in tsv file and runs SoloTE pipeline on each line
    * tsv input format: <full pathway to CellRanger BAM file>\t<OUTPUT_PREFIX>
    * OUTPUT_PREFIX = prefix for all SoloTE output files and name of output directory (flags -o/-d respectively); should be dataset/seq lib name

### Steps 1 and 2 completed in solotePipeline.sh (note: this is the longest step) ###

3. Merge SoloTe mtxs and make Seurat object
    * matches cells of solote sobj to reference sobj
    * adds in umap reduction and metacols from reference
    * saves final solote sobj

### Step 3 completed in soloteSobj.sh ###

4. Run DESeq2 pipeline on SoloTE sobj
    * pairwise comparison DESeq2 on solote
    * originally DESeq2 output will have both genes and TEs - ran together for transcriptome wide variability
    * later will filter DESeq2 output files to have ONLY TEs - run DE reports on these
    * output for these are filtered features bulk mtxs, bulk mtxs, deseq2 csv output files (with both genes and TEs - will be filtered i next step)

5. Create SoloTE version of DE.report
    * TE.reports, summary barplots, heatmaps, solote count summary
    * takes all DESeq_out csv files and filters for only TEs

### Step 4 and 5 completed in soloteBulkMtxDeseq2.sh ###

DONE :)
