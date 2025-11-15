import numpy as np
import pandas as pd 
import scanpy as sc
import episcanpy as epsc
import argparse

### Functions 
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

## add celltype param to not count all barcodes 
def build_obj_from_peaks(sample, frag, regions, outdir):
    
    tmeta = metadata[metadata['sample'] == sample]
    barcodes = tmeta['barcode']

    ### Dev version fixed it ****** 
    adata = epsc.count_matrix.peak_mtx(frag, regions, barcodes)
    
    adata.obs.index = [f'{sample}_{b}' for b in adata.obs.index]

    adata.write_h5ad(f'{outdir}{sample}.h5ad')



if __name__ == '__main__':

    args = getargs()
    print(f'Fragments: {args.frag}')    
    print(f'Regions: {args.regions}')    
    print(f'OutDir: {args.outdir}')    
    print(f'Metadata: {args.meta}')    

    global metadata
    metadata = pd.read_csv(args.meta, index_col=0)
    ### Might have to change in future to suit barcode nomenclature 
    metadata['barcode'] = [i.split("_")[-1] for i in metadata.index]
    metadata['sample'] = ["_".join(i.split("_")[:-1]) for i in metadata.index]

    sample = args.frag
    sample = sample.split('/')[-3]
    print(f'Sample: {sample}')    
    print(f'Barcodes in sample: ', len(metadata[metadata['sample'] == sample]))
    

    build_obj_from_peaks(sample, args.frag, args.regions, args.outdir)

    print('Finish')




