library(Seurat)
library(ggplot2)
library(Signac)
library(argparse)
library(future)
library(tableHTML)
library(MAST)
library(sctransform)

### START OF FUNCTIONS ###
process_args <- function(){
    # create parser object
    parser <- ArgumentParser(description = "Create merged SoloTE Sobj from user selected SoloTE Mtxs")

    # Create parser arg group for our required args
    required_arg_group = parser$add_argument_group('flagged required arguments',
        'the script will fail if these args are not included')
    
    required_arg_group$add_argument("-mtx","--soloteMtx", required = TRUE,
        help="SoloTE output matrix. Options include: classtes, familytes, subfamilytes, and legacytes, locustes, but recommend legacy or locus.")

    required_arg_group$add_argument("-obj","--sobj", required = TRUE,
        help="Reference/Main seurat object (must include complete pathway); will be used to add UMAP coordinates and metacol to SoloTE sobj.")

    required_arg_group$add_argument("-col","--metaCol", required = TRUE,
        help="Metacol(s) from main seurat object to add to SoloTE sobj; either one or multiple metacols allowed.", nargs='+')

    required_arg_group$add_argument("-proj","--projName", required = TRUE,
        help="Project name used in merge function and as prefix for final merged SoloTE rds object.")

    required_arg_group$add_argument("-odir","--outputDir", required = TRUE,
    help="Output directory for solote sobj, output files, and other visualizations.")

    required_arg_group$add_argument("-umap","--umapReduc", required = TRUE,
    help="Reference UMAP reduction name to project over SoloTE UMAP.")

    args <- parser$parse_args()
    return(args)
}

# this function goes through user selected solote mtx, creates sobj for each dset/seq lib, and merges all dset sobjs into one final solote sobj
# note that this merged solote sobj does not match main sobj yet
createSoloteSobj <- function(sobj_path, solote_mtx, odir, projectName, addMetaCol, umap){
    # read in reference/main sobj
    sobj <- readRDS(sobj_path)
    DefaultAssay(sobj) <- "RNA"

    # all of the mtx directory from solote; go through each directory and make seurat objs from these
    solote_mtxs <- grep(solote_mtx, list.dirs(odir, full.names = T), value = T)
    print(solote_mtxs)

    # go through each mtx outputted from solote
    # num mtx = num of dsets/seq ids
    # make solote sobj for each dset
    # merge sobjs into one giant solote object
    sobj_list <- c()
    cellIDs <- c()
    for(mtx in solote_mtxs){
        print(mtx)
    
        # needed when merge sobjs together - added to beginning of barcode
        dsetName <- strsplit(mtx, split = "/")[[1]]
        dsetName <- strsplit(dsetName[length(dsetName)], split = paste0("_",solote_mtx,"_MATRIX"))[[1]]
        print(dsetName)
    
        # getting solote files needed to create sobj
        barcodesFile <- list.files(mtx, pattern = "barcodes", full.names = T)
        featuresFile <- list.files(mtx, pattern = "features", full.names = T)
        matrixFile <- list.files(mtx, pattern = "matrix", full.names = T)
    
        # making each ind sobj (based on seq libs/dsets)
        # note: this makes sobj with only RNA assay (counts, data) and 0 variable features
        solote_matrix <- ReadMtx(mtx = matrixFile, cells = barcodesFile, features = featuresFile, feature.column=1)
        
        dsetSobj <- CreateSeuratObject(count=solote_matrix, min.cells=3, project=dsetName)
        print(paste0("Done creating sobj for: ", dsetName))
        
        # collecting cell ids to add on barcode of merged solote obj
        cellIDs <- append(cellIDs, dsetName)
        
        sobj_list <- append(sobj_list, dsetSobj)   
    }

    # done creating individual sobjs; now merge
    mobj <- merge(x = sobj_list[[1]], y = sobj_list[2:length(sobj_list)], add.cell.ids = cellIDs, proj = projectName)
    mobj
    print("Done merging all solote sobjs")

    # join layers
    mobj <- JoinLayers(mobj)

    # take merged obj and match to main sobj (cell num and order barcodes)
    print("Reference:")
    print(sobj)

    print("pre-subset merged SoloTE sobj:")
    print(head(mobj@meta.data))

    # done merging; subset merged solote obj to match main sobj
    sub_solote <- subset(mobj, cells = colnames(sobj))
    solote_meta <- sub_solote@meta.data
    sub_solote@meta.data <- solote_meta[order(match(rownames(solote_meta), colnames(sobj))),]

    print("Subsetted SoloTE sobj:")
    print(head(sub_solote@meta.data))

    print("Main sobj:")
    print(head(sobj@meta.data))

    # looping through reference metacols and adding to solote sobj
    for(mcol in addMetaCol){
        sub_solote <- AddMetaData(sub_solote, sobj[[mcol]], paste0("Ref_", mcol))
        print(paste0("Added metacol: ", mcol))
    }

    # transferring UMAP coordinates to solote sobj
    sobjUMAP <- as.matrix(Embeddings(sobj[[umap]]))
    head(sobjUMAP)
    sub_solote@reductions[[paste0("Ref_", umap)]] <- CreateDimReducObject(embeddings = sobjUMAP,
                                                                               key = paste0("Ref_", umap, "_"),
                                                                               assay = DefaultAssay(sobj))
    print("Done adding metacol(s) and umap reduction to solote")
    print(colnames(sub_solote@meta.data))

    Idents(sub_solote) <- "orig.ident"
    sub_solote <- NormalizeData(sub_solote)

    # saving final merged solote sobj
    saveRDS(sub_solote, paste0(odir, projectName, ".rds"))
    print("RDS saved! yay")
}
### END OF FUNCTIONS ###


### MAIN ###
args <- process_args()

createSoloteSobj(sobj_path = args$sobj, solote_mtx = args$soloteMtx,
                 odir = args$outputDir, projectName = args$projName,
                 addMetaCol = args$metaCol, umap = args$umapReduc)
