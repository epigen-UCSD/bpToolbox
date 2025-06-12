import scanpy as sc
import pandas as pd 
import os
import sys
from pathlib import Path

### Export RNA to h5ad

wdir = sys.argv[1]
rna_counts = sc.read_mtx(f"{wdir}/Seurat_Data/rna_matrix.mtx")
rna_counts = rna_counts.T.copy()

rna_cell_meta = pd.read_csv(f"{wdir}/Seurat_Data/meta_data.tsv", sep="\t", header=0, index_col=0)
rna_gene_names = pd.read_csv(f"{wdir}/Seurat_Data/gene_names.tsv", sep="\t", header=None, index_col=0)
rna_gene_names.index.name = None

rna_counts.raw = rna_counts.copy()
sc.pp.normalize_total(rna_counts, target_sum=1e4)
sc.pp.log1p(rna_counts)
sc.pp.highly_variable_genes(rna_counts, min_mean=0.0125, max_mean=3, min_disp=0.5)
rna_counts = rna_counts[:, rna_counts.var.highly_variable]
sc.pp.scale(rna_counts, max_value=10)
sc.tl.pca(rna_counts)
sc.pp.neighbors(rna_counts)
sc.tl.umap(rna_counts)

os.mkdir("Scanpy_out")
rna_counts.write("Scanpy_out/rna_counts.h5ad")


