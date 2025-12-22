
import sys
sys.path.append("/tscc/projects/ps-epigen/bpToolbox/snATAC/")
import splitfrags

### add space for a split later 
frag_paths = [f' {i}' for i in frag_paths]
samples = [p.split("/")[-3] for p in frag_paths]

#outdir = 'ac_splitfrags_out/'
os.makedirs(outdir, exist_ok=True)

### Specify metadata col to split by 
#celltype_col = 'celltype'

#metadata = pd.read_csv("k27ac_metadata.csv", index_col=0)
## Build sample and barcode columns to query 
## Note: run_split_frags want's explicent "sample" and "barcode" named columns 
metadata['sample'] = ["_".join(i.split("_")[:-1]) for i in metadata.index]
metadata['barcode'] = [i.split("_")[-1] for i in metadata.index]

frags = pd.DataFrame.from_dict({"samples": samples, "paths": frag_paths})
frags['comb'] = frags["samples"].astype(str) + frags["paths"].astype(str)

### Create all combinations off celltypes and frag files with itertools product() 
combinations_temp = list(iter.product(metadata[celltype_col].unique(), frags['comb']))
### Separate the sample from frag path
## len(combinations) == (celltypes * samples) - dropout
combinations = [[c[0], c[1].split(" ")[0], c[1].split(" ")[1]] for c in combinations_temp]

splitfrags.run_split_frags(combinations, 20, outdir, ct_col=celltype_col, meta=metadata)


### Merging 
os.makedirs(merge_out, exist_ok=True)
mcombos = [[f'{c}', outdir, merge_out] for c in metadata[celltype_col].unique()] 
splitfrags.run_merge_frags(mcombos, 12)



