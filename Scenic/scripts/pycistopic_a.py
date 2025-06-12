import os
import pycisTopic
from pycisTopic.cistopic_class import create_cistopic_object
pycisTopic.__version__

import pandas as pd
import numpy as np
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import pickle

from scipy import io
import scipy.sparse as sp
from scipy import io
import subprocess


### Define project and output vars 
wdir = sys.argv[1]
project = sys.argv[2]
genome = sys.argv3[3]

outdir = f"{wdir}/pycistopic_outputs/"
if not os.path.exists(outdir): os.mkdir(f"pycistopic_outputs")

if genome == "hg38":
    path_to_blacklist = "/tscc/projects/ps-epigen/GENOME/hg38/hg38.blacklist.bed.gz"
elif genome == "mm10":
    path_to_blacklist = "/tscc/projects/ps-epigen/GENOME/mm10/mm10.blacklist.bed.gz"
else
    print(f"{genome} not supported")

### Build object
meta_df = pd.read_csv("Seurat_data/meta_data.tsv", header=0, sep="\t", index_col=0)

cell_names_tb = pd.read_csv("Seurat_data/barcodes.tsv", header=None, sep="\t")
cell_names = cell_names_tb[0].tolist()
region_names_tb = pd.read_csv("Seurat_data/region_names.tsv", header=None, sep="\t")
region_names = region_names_tb[0].tolist()
region_names = ["{}:{}-{}".format(*x.split("-")) for x in region_names]

atac_count_matrix = io.mmread("Seurat_data/atac_matrix.mtx")

cistopic_object = create_cistopic_object(
    fragment_matrix=atac_count_matrix.tocsr(),
    cell_names=cell_names,
    region_names=region_names,
    path_to_blacklist=path_to_blacklist,
    project=project,
    tag_cells=False
)

pickle.dump(
    cistopic_object,
    open("pycistopic_outputs/cistopic_object.pkl", "wb")
)

if not os.path.exists(f"{outdir}pycistopic_slurm/" ): os.mkdir(f"{outdir}pycistopic_slurm/")
if not os.path.exists(f"{outdir}model_out/" ): os.mkdir(f"{outdir}model_out/")
if not os.path.exists(f"{outdir}model_logs/" ): os.mkdir(f"{outdir}model_logs/")



# Use combination of f'strings' and .format to build template 
slurm_script_template = (
    "#!/bin/bash\n"
    "#SBATCH -p platinum\n"
    "#SBATCH -q hcp-csd772\n"
    "#SBATCH -J {{0}}\n"
    "#SBATCH -N 1\n"
    "#SBATCH -c 8\n"
    "#SBATCH --mem 100G\n"
    "#SBATCH -t 16:00:00\n"
    f"#SBATCH -o {outdir}/model_logs/sys_modeling_{{0}}.out\n"
    "#SBATCH -e /tscc/nfs/home/rlancione/jeo/psys_modeling_{0}.err\n"
    "#SBATCH --mail-user rlancione@ucsd.edu\n"
    "#SBATCH --mail-type FAIL\n"
    "#SBATCH -A csd772\n"
    "\n"
    "set -e\n"
    "source ~/.bashrc\n"
    "conda activate scenicplus\n"
    "\n"
    f"python3 {wdir}/mallet_script.py {{0}} {wdir}\n"
)


for n_topic in [2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]:
    with open(f"{outdir}pycistopic_slurm/model_{n_topic}_slurm.sh", "w") as f:
        f.write(slurm_script_template.format(n_topic))


#for n_topic in [2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]:
#    sbatch_cmd = f"sbatch {outdir}pycistopic_slurm/model_{n_topic}_slurm.sh"
#    print(sbatch_cmd)
#    subprocess.run(sbatch_cmd, shell = True)




