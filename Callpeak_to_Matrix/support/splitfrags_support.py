import sys
import argparse
import numpy as np
import pandas as pd 

sys.path.append("/tscc/projects/ps-epigen/bpToolbox/snATAC/")
import splitfrags


def getargs():
    # 1. Initialize the parser
    parser = argparse.ArgumentParser(
        description="A script to filter and split fragment files",
    )

    # 2. Add Positional Arguments (Required)
    parser.add_argument("--paths_txt", help="Path to the input frags.txt")
    parser.add_argument("--groups", help="column of to split frags by in meta. e.g celltype")
    parser.add_argument("--metadata", help="path to metadata")
    parser.add_argument("--out", help="out path")
    parser.add_argument("--cores", default = 20, help="ncores for split frags. > 20 is good for large data")


        # 5. Parse the arguments
    args = parser.parse_args()
    return(args)


# gather args 
args = getargs()

metadata = pd.read_csv(args.metadata, index_col=0)
frag_paths = pd.read_csv(args.paths_txt)[1]

### add space for a split later 
frag_paths = [f' {i}' for i in frag_paths]
samples = [p.split("/")[-3] for p in frag_paths]

outdir = f'{args.out}sample_tmp/'
os.makedirs(outdir, exist_ok=True)


## Build sample and barcode columns to query 
## Note: run_split_frags want's explicent "sample" and "barcode" named columns 
### ****** This may need to be changed if barcodes are not formatted as sample_barcode 
metadata['sample'] = ["_".join(i.split("_")[:-1]) for i in metadata.index]
metadata['barcode'] = [i.split("_")[-1] for i in metadata.index]

frags = pd.DataFrame.from_dict({"samples": samples, "paths": frag_paths})
frags['comb'] = frags["samples"].astype(str) + frags["paths"].astype(str)

### Create all combinations off celltypes and frag files with itertools product() 
combinations_temp = list(iter.product(metadata[args.groups].unique(), frags['comb']))
### Separate the sample from frag path using our space
combinations = [[c[0], c[1].split(" ")[0], c[1].split(" ")[1]] for c in combinations_temp]

splitfrags.run_split_frags(combinations, args.cores, outdir, ct_col=args.groups, meta=metadata)


### Merging 
merge_out = f'{args.out}merged/' 
os.makedirs(merge_out, exist_ok=True)
mcombos = [[f'{c}', outdir, merge_out] for c in metadata[args.groups].unique()] 
splitfrags.run_merge_frags(mcombos, args.cores/2)



