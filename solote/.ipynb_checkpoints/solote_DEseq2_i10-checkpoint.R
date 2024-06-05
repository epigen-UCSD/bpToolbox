suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(DESeq2))
suppressMessages(library(enrichR))
suppressMessages(library(stringr))
suppressMessages(library(tableHTML))
suppressMessages(library(argparse))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(enrichplot))
suppressMessages(library(gridExtra))
suppressMessages(library(ReactomePA))

### Functions
# please note - in the original RNA_DEseq2_i10.R script from Ryan there are mult filter feature functions:
# filt.features, filt.features.avg, filt.features.avg_v2, and filt.features.avg_v3
# default one used is filt.features.avg_v3; some filt functions contain "lung BPD" specific conditions
# for solote i removed unique lung BPD condition statements from these functions

process_args <- function(){
    # create parser object
    parser <- ArgumentParser(description = "Create seurat object from 10x h5 and collect QC metrics")

    # Create parser arg group for our required args
    required_arg_group = parser$add_argument_group('flagged required arguments',
        'the script will fail if these args are not included')

    # These args belong to my group of required flagged arguments
    required_arg_group$add_argument("-sobj", "--seurat_object", required = TRUE,
        help="path to dataset h5 file for object creation")

    required_arg_group$add_argument("-ctcol","--celltype_column", required = TRUE,
        help="column containing id of the patient")

    required_arg_group$add_argument("-mpath","--matrix_path", required = TRUE,
        help="column containing condition group")

    required_arg_group$add_argument("-proj","--project_name", required = TRUE,
        help="name of the DE project. e.g. Lung_BPD_subD352")

    required_arg_group$add_argument("-o","--output_path", required = TRUE,
        help="path to output directory")

    required_arg_group$add_argument("-nf","--n_features", required = TRUE,
        help="n features must be observed in every donor for at least one condition group")

    # optional args
    parser$add_argument("-comp","--comparisons",
        help="path to tsv containing comparisons. Use header: Target  Reference")

    parser$add_argument("-form","--formula",
        help="a string representing the DESeq formula. Default is: ~ condition")

    
    args <- parser$parse_args()
    return(args)
}

filt.features.avg <- function(blk.matrix, meta, avg.thresh, comp_pair){
    # will store our list of features that passed theshold
    pf <- c()

    g1.donors <- rownames(meta[meta$condition == comp_pair[1],])
    g2.donors <- rownames(meta[meta$condition == comp_pair[2],])

    # Iterate over rows (features) of bulk matrix
    for (r in seq(1, nrow(blk.matrix))){
        g1.values <- as.integer(blk.matrix[r, g1.donors])
        g2.values <- as.integer(blk.matrix[r, g2.donors])

        # If there were at least n counts in all donors in one of the condition groups
        # Add that feature to the list
        if (mean(g1.values) >= avg.thresh | mean(g2.values) >= avg.thresh){
            pf <- append(pf, rownames(blk.matrix[r,]))
        }
    }
    return(pf)
}

# Avg feature filter using avg(all donors)
filt.features.avg_v2 <- function(blk.matrix, meta, avg.thresh, comp_pair){
    # will store our list of features that passed theshold
    pf <- c()
    g1.donors <- rownames(meta[meta$condition == comp_pair[1],])
    g2.donors <- rownames(meta[meta$condition == comp_pair[2],])
    
    # Iterate over rows (features) of bulk matrix
    for (r in seq(1, nrow(blk.matrix))){
                                
        values <- as.integer(blk.matrix[r, c(g1.donors, g2.donors)])
        
        if (mean(values) >= avg.thresh){
            pf <- append(pf, rownames(blk.matrix[r,]))
        }
    }

    return(pf)
}


# Avg feature filter using avg(all donors) + require at least one count from each donor group
filt.features.avg_v3 <- function(blk.matrix, meta, avg.thresh, comp_pair){
    # will store our list of features that passed theshold
    passed.features <- c()
    # vectors of donor ids pertaining to the comparison
    g1.donors <- rownames(meta[meta$condition == comp_pair[1],])
    g2.donors <- rownames(meta[meta$condition == comp_pair[2],])
    
    # Iterate over rows (features) of bulk matrix
    for (r in seq(1, nrow(blk.matrix))){
                                
        g1.values <- as.integer(blk.matrix[r, g1.donors])
        g2.values <- as.integer(blk.matrix[r, g2.donors])
        values <- c(g1.values, g2.values)
        
        g1.sum <- sum(g1.values)
        g2.sum <- sum(g2.values)
        
        if (mean(values) >= avg.thresh & g1.sum >= 1 & g2.sum >= 1){
            passed.features <- append(passed.features, rownames(blk.matrix[r,]))
        }
    }

    return(passed.features)
}



filt.features <- function(blk.matrix, meta, n, comp_pair){
    # will store our list of features that passed theshold
    pf <- c()

    g1.donors <- rownames(meta[meta$condition == comp_pair[1],])
    g2.donors <- rownames(meta[meta$condition == comp_pair[2],])
    
    # Iterate over rows (features) of bulk matrix
    for (r in seq(1, nrow(blk.matrix))){
      
        g1.values <- blk.matrix[r, g1.donors]
        g2.values <- blk.matrix[r, g2.donors]
    
        # Vectorize the query across the columns 
        g1.bool.list <- g1.values >= n
        g2.bool.list <- g2.values >= n
    
        # Sum the trues from our query
        g1.ntrue <- sum(g1.bool.list)
        g2.ntrue <- sum(g2.bool.list)
    
        # If there were at least n counts in all donors in one of the condition groups
        # Add that feature to the list
        if (g1.ntrue >= length(g1.donors) | g2.ntrue >= length(g2.donors)){
            pf <- append(pf, rownames(blk.matrix[r,]))
        }
    }
    
    return(pf)
}


pseudoBulk_DE <- function(sobj.path, matrices.path, ct.col = "seurat_clusters",
				output_dir, project_name, comps = "pairwise", formula = "~ condition", nfeatures = 0){  
  print("Initiating DESeq Pipeline")
  seuObj <- readRDS(sobj.path)

  ### Create output outdirs
  dir.create(output_dir, showWarnings = FALSE)
  dir.create(paste0(output_dir,"DESeq_out"), showWarnings = FALSE)

  if (ct.col == "seurat_clusters"){
      cell_type <- sort.int(unique(seuObj[[ct.col]])[,1])
  } else {
      cell_type <- unique(seuObj[[ct.col]])[,1]
      cell_type <- str_replace_all(cell_type, "/", ".")
  }

  ### Iterate through celltypes. Run DESeq. Iterate through comparisons. Extract results. Call Go.
  for (c in cell_type){
      print(paste("STARTING: Celltype", c))

      ### Reading files
      blk.matrix <- read.csv(paste(matrices.path, project_name, "_", c, "_bulk_matrix.csv", sep=""), sep=',',
            header=TRUE, row.names=1, stringsAsFactors=FALSE)

      metadata <- read.csv(paste(matrices.path, project_name,"_metadata.csv", sep=""), sep=',', header=TRUE,
            row.names = 1, stringsAsFactors=FALSE)

      ### subset metadata for donors in celltype blk matrix
      sub.meta <- subset(metadata, row.names(metadata) %in% colnames(blk.matrix))

      ### Extract our comparisons it iterate through. Either pairwise combo of all or user defined.
      # if comps is null, run pairwise comparisons. null argparse as input to pseudobulk_DE overwrites default "pairwise"
      if (is.null(comps)){
          # DEFAULT. Run pairwise
    	  pair_list <- combn(unique(metadata$condition), 2)
    	  num_pairs <- ncol(pair_list)
      } else {
          comp.table <- read.table(comps, header = T, sep = "\t")
          num_pairs <- nrow(comp.table)
          pair_list <- t(comp.table)
      }
	
      ### Begin iterating through comparisons to extract results for each group
      for (i in 1:num_pairs){
          cpair <- pair_list[,i]
          cond1 <- cpair[1]
          cond2 <- cpair[2]
	  
          # subset metadata for desired conditions. ### Still need to check whether enough donors for comparison
          sub.meta.comp <- sub.meta[which(sub.meta$condition==cond1 | sub.meta$condition==cond2), ]
    	  blk.mat.comp <- subset(blk.matrix, select = row.names(sub.meta.comp))
    	  # If only one condition is represented by donors existing in a celltype. Skp DESeq
    
    	  ###### Here is where features we test will be filtered by a specified N.
    	  ### These filtered matrices will be saved for convinience
    	  # if filter.features == TRUE do this. if not. pass blk.mat.comp to blk.mat.comp.filt
    	  print(paste0("identify valid features: ", cond1, "v", cond2))
    	  nfeatures <- as.integer(nfeatures)
    	  # use average feature filter. for normal filter use n. for avg use avg.thresh param
    	  pass.features <- filt.features.avg_v3(blk.matrix = blk.mat.comp, meta = sub.meta.comp, avg.thresh = nfeatures, comp_pair = cpair)
    	  
    	  blk.mat.comp.filt <- blk.mat.comp[pass.features,]
    	  print(paste0("nrow filt matrix:", nrow(blk.mat.comp.filt)))

      	  # save filtered matrix 
    	  condPair <- paste0(cond1, "-", cond2)
    	  dir.create(paste0(output_dir, "Filtered_Matrices/"), showWarnings = FALSE) 	  
          filt_mat_file <- paste0(output_dir, "Filtered_Matrices/", c, "_" , condPair, "_filt_bulk_mat.csv")
          write.csv(blk.mat.comp.filt, file = filt_mat_file, quote=FALSE)
    
      	  ### ************ else if condi = active | chronic batch add TOD as covariate 
          if (length(unique(sub.meta.comp$condition)) > 1 & length(pass.features) >= 100){
              ### If no formula is input. Only test condition. If formula is present, try to use it. If it doesn't work, only test condition
              form_used <- "~ condition"
              if(is.null(formula)){
                  # DEFAULT. Simply test each condition
                  deseq_results <- DESeqDataSetFromMatrix(countData = blk.mat.comp.filt, colData = sub.meta.comp, design = ~ condition)
              } else {
                  # try to use formula  
                  print("try to use form")
                  form_used <- formula
                  an.error.occured <- FALSE
                  tryCatch({deseq_results <- DESeqDataSetFromMatrix(countData = blk.mat.comp.filt,
                                                                    colData = sub.meta.comp, design = as.formula(formula)) },
                           error = function(e) {an.error.occured <<- TRUE})
                         
                  # if error run default
                  if (an.error.occured == TRUE){
                      # Change form used back to default since input didn't work 
                      print("error occured, use default form")	
                      form_used <- "~ condition"
                      deseq_results <- DESeqDataSetFromMatrix(countData = blk.mat.comp.filt, colData = sub.meta.comp, design = ~ condition)
                  }
               }

               deseq_results <- DESeq(deseq_results)
               ############ Declare comparison so direction of results are known ##############
               final <- results(deseq_results, contrast = c("condition", cond1, cond2))
               # keep only if both padj and log2FC is NOT NA
               final <- final[is.na(final$log2FoldChange) == FALSE, ]
               final <- final[is.na(final$padj) == FALSE, ]

               # write deseq matrix to file
               deseq_matrix <- paste(output_dir, "DESeq_out/", "DESeq2_", c, "_" , condPair, ".csv", sep="")
               write.csv(final, file=deseq_matrix, quote=FALSE)
           } else {
               # Deseq won't run if these conditions aren't met
               print("not enough replicates to run DESeq for comparison")
               print("or less than 100 features passed filter")
               # Write not enough donors to make comparison
               summarylist <- c()
               summarylist <- c(summarylist, paste(c, "COMBINED", condPair, "NA", "NA", sep = "\t"))
               #write.table(summarylist, paste(output_dir, "DEG_countSummary.txt", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
               # 240603 KD edit - i didn't write out the DEG_countSummary here since it includes DEGs and TEs
               # for solote we want to run deseq2 on entire transcriptome, but only keep the TE outputs; having this combination DEG_countSummary might be confusing
           }
      
           print(paste("DONE: Celltype ", c, " | Conditions ", cond1, " vs. ", cond2, sep=""))
      } # end going through comparisons
  } # end going through CTs  
} # end pseudobulk de func

########### MAIN ###########
## Aquire inputs
args <- process_args()

### Record input
print("****** DE Pipeline Inputs ******")
print(paste0("sobj.path: ", args$seurat_object))
print(paste0("matrices.path: ", args$matrix_path))
print(paste0("output.dir: ", args$output_path))
print(paste0("project.name: ", args$project_name))
print(paste0("celltype column: ", args$celltype_column))
print(paste0("comparison file: ", args$comparisons))
print(paste0("n features: ", args$n_features)) 
if(is.null(args$formula) == TRUE){
	print("formula: ~ condition")
} else{
	print(paste0("formula: ", args$formula))
}
if(is.null(args$comparisons) == TRUE){
	print("comparisons: pairwise")
} else{
	print("comparisons: ")
  print(read.table(args$comparisons, header = T, sep = "\t"))
}

pseudoBulk_DE(sobj.path = args$seurat_object, ct.col = args$celltype_column, matrices.path = args$matrix_path,
              output_dir = args$output_path, project_name = args$project_name,
              comps = args$comparisons, formula = args$formula, nfeatures = args$n_features)
