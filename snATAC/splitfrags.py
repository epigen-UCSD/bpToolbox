import os 
import pandas as pd 
import numpy as np
import itertools as iter

from pathlib import Path
from functools import partial
from multiprocessing import Pool


### Metadata must have sample, celltype, barcodes
### Since we have the celltypes samples, and fragments paired 
### We can read from the fragments concurrently. Then write to independent frag files

### Functions 
def split_frags(celltype, sample, fragment_file, outpath, ct_col = 'celltype', chunksize = 10**6, metadata = None):
    
    ### Want to use barcodes for only ct / cluster 
    ctmeta = metadata[metadata[ct_col] == celltype]

    ### group barcodes by sample id's and store in dict 
    ### Could refactor to only build group_values once 
    group_values = ctmeta.groupby("sample")["barcode"].apply(list).to_dict()

    ### Extract barcodes for sample of interest 
    ## If sample exists in this dict. That means barcodes exist as values 
    if sample in group_values.keys():
        keep_barcodes = group_values[sample]

        ### Define columns
        cols = ["chrom", "start", "end", "barcode", "mapq"]
        ### Built iterable of our frag file.. basically chunked up version of file 
        reader = pd.read_csv(fragment_file, 
                            sep="\t", 
                            comment="#", 
                            names=cols,
                            chunksize=chunksize)
    
        ### Formats filename 
        output_file = f'{outpath}{celltype}_{sample}.tsv'
    
        first_chunk = True  
        for chunk in reader:
            
            # Now chunk is a DataFrame, so we can filter
            filtered = chunk[chunk["barcode"].isin(keep_barcodes)]
    
            # write
            if not filtered.empty:
                if first_chunk:
                    filtered.to_csv(output_file, sep="\t", index=False, mode="w")  # include header
                    first_chunk = False
                else:
                    filtered.to_csv(output_file, sep="\t", index=False, mode="a", header=False)  # append, no header
                    
    ### else no sample in our dict
    else: return 


def run_split_frags(combinations, p, odir, ct_col, meta):

    # wrapper func helps delclare keyword args not in combinations 
    # **************** had to change to leiden 
    # add to run_split_frags() to handle ct_col
    worker_func = partial(split_frags, ct_col=ct_col, chunksize=10**7, metadata=meta)


    # format args as 2d list for starmap
    ### append outpath to combinations
    args=[i + [odir] for i in combinations]

    with Pool(processes=p) as pool:
        pool.starmap(func=worker_func, iterable=args)


### Needs to be adapted to leiden clusters e.g 2_MM_2431_2 
def merge_frags(celltype, indir, outdir):
    tsv_files = list(Path(indir).glob("*.tsv"))
    ### extract all files containing the celltype string 
    #filtered = [p for p in tsv_files if celltype in str(p)]
    filtered = [p for p in tsv_files if f'{indir}{celltype}_' in str(p)]
    
    output_file = Path(f"{outdir}{celltype}.bed")

    #no header 
    first_file = False  # add header if first chunk
    for tsv in filtered:
        # Read in chunks
        for chunk in pd.read_csv(tsv, sep="\t", chunksize=10**7):
            chunk.to_csv(output_file, sep="\t", mode='a', index=False, header=first_file)
            first_file = False  # After the first chunk of the first file, header is written


def run_merge_frags(combinations, p):
    
    # format args as 2d list for starmap
    ### append outpath to combinations
    args=combinations

    with Pool(processes=p) as pool:
        pool.starmap(func=merge_frags, iterable=args)


def getargs():
    # Create the parser
    parser = argparse.ArgumentParser(description="Gather args to build ATAC matrices by sample/frags")

    # Add arguments
    parser.add_argument("-f", "--frag", type=str, required=True,
                        help="path to fragments")

    parser.add_argument("-r", "--regions", type=str, required=True,
                        help="path to regions")

    parser.add_argument("-o", "--outdir", type=str, required=True,
                        help="output directory")

    parser.add_argument("-m", "--meta", type=str, required=True,
                        help="metadata path")

    # Parse the arguments
    args = parser.parse_args()

    return args



if __name__ == '__main__':

    ## Get args 
    args = getargs()
    print(f'Fragments: {args.frag_paths}')
    print(f'Metadata: {args.meta}')    
    print(f'celltype_column: {args.ct_col}')
    print(f'OutDir: {args.outdir}')

    ### Init metadata as global var so don't need to add to combinations 
    ### Will be drawn on by splitfrags 
    global metadata
    #metadata = pd.read_csv("/tscc/projects/ps-epigen/users/rlan/LFNIH_Pt2/Subclustering/SCVI/multi_metadata.csv", index_col=0)
#    metadata = pd.read_csv("hsc_meta.txt", index_col=0, sep=' ')
    ## Only call peaks with multiome samples/barcodes 
    metadata = metadata[metadata['assay'] == 'multiome']
    ## Build sample and barcode columns to query 
    metadata['sample'] = ["_".join(i.split("_")[:-1]) for i in metadata.index]
    metadata['barcode'] = [i.split("_")[-1] for i in metadata.index] 

#    os.makedirs(outdir, exist_ok=True)


#    frag_paths = args.frags

    ### Format input 
#    samples = [p.split("/")[-3] for p in frag_paths]
#    samples

    frags = pd.DataFrame.from_dict({"samples": samples, "paths": frag_paths})
    frags['comb'] = frags["samples"].astype(str) + frags["paths"].astype(str)



    ## *****  
    ### Create all combinations off celltypes and frag files with itertools product() 
    combinations_temp = list(iter.product(metadata[args.cluster_col].unique(), frags['comb']))
    ### Separate the sample from frag path
    combinations = [[c[0], c[1].split(" ")[0], c[1].split(" ")[1]] for c in combinations_temp]




## Splitfrags proccessess in cluster/celltype fragment pairs. 
## 10 new celltypes mean 10 more iterations per fragment file.
## With 28 fragment files. That's 280 more times to parse the entire same fragment files with 10 celltypes 
## So even with less barcodes, adding subclusters is computationally expensive 

### Improved method might retain output structure 
### But only iterate over each frag once.
### While keeping a map of index to celltype barcodes 
### Then indexes are drawn for each celltype to write relevant files 
### use many .isins() to dict instead of one 

