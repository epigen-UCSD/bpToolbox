import sys
import os
sys.path.append(os.path.abspath("/tscc/projects/ps-epigen/bpToolbox/Scenic/scripts/"))
from heatmap_utils import *


### Inputs 
proj_path = "scplus_otsu"
celltype_col = "celltypes"
condition_col = "condition"
outdir = "Reports_otsu"

os.makedirs(f"{outdir}", exist_ok = True)
scplus_mdata = mudata.read(f"{proj_path}/Snakemake/scplusmdata.h5mu")


###### Specificity Reports 
specificity_report(scplus_mdata, celltype_col, condition_col, outdir):


###### Plotting heatmaps 

### Report contents:
# heatmaps
    # gAUC / rAUC 
    # TFe / gAUC 
    # TFe / rAUC 
    # TFe / specificity 

# Autoscale for nGRN / Catagory 

### additional params for 
    # direct / extended 
    # Top N GRN

dff = make_heatmap_df(scplus_mdata, celltype_col, direct == True, onlypos == False)

heatmap_plotting(scplus_mdata, dff, motifAnno = "direct", celltype_col, outdir)
