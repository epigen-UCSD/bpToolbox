import pandas as pd
import numpy as np
from datetime import datetime
import re
import sys
import glob
import subprocess
import os.path
import os

## Primary vars 
# 1. Working dir.. to derive Seurat_Data
# 2. genome hg38/mm10
# 3. Project "ZM"


def write_pfasta(wdir, genome):
    # bed path
    bed = f"{wdir}Seurat_Data/consensus.bed"    # Format project path 
    
    # Replace strings in file with extracted idir & pdir
    with open("/tscc/projects/ps-epigen/users/rlan/sandbox/scenic/templates/pfasta.sh", 'rt') as fh:
        shell = fh.read()
        #replace all occurrences of the required string
        shell = shell.replace('[BED_PATH]', bed)
        shell = shell.replace('[GENOME]', genome)

    
    if not os.path.exists(f"{wdir}prepare_fasta/"): os.mkdir(f"{wdir}prepare_fasta/")     
    with open(f"{wdir}prepare_fasta/pfasta.sh", "w") as fh:
        fh.write(shell)


    # subprocess.run(f"./{wdir}prepare_fasta/pfasta.sh")



def cistarget_build(wdir, project):
    # fasta path
    fasta = f"{wdir}prepare_fasta/consensus.fa"    # Format project path 
    
    # Replace strings in file with extracted idir & pdir
    with open("/tscc/projects/ps-epigen/users/rlan/sandbox/scenic/templates/cistarget_array.sh", 'rt') as fh:
        shell = fh.read()
        #replace all occurrences of the required string
        shell = shell.replace('[FASTA_PATH]', fasta)
        shell = shell.replace('[PROJECT]', project)
        shell = shell.replace('[OUT_DIR]', f"{wdir}/cistarget_db/out/partial_db/")

    
    if not os.path.exists(f"{wdir}cistarget_db/out/partial_db/"): os.mkdirs"{wdir}cistarget_db/out/partial_db/")     
    with open(f"{wdir}cistartet_db/cistarget_array.sh", "w") as fh:
        fh.write(shell)


    #sbatch_cmd = f"sbatch --array=1-{ndataset} {proj_path}bATAC_pipe.sh"
    #print(sbatch_cmd)
    #subprocess.run(sbatch_cmd, shell = True)



### During setup intitally these two scripts are built and launched
## 




