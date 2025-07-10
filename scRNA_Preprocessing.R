###### Clear screen ######
# cat("\014")

###### Clear plot history ######
# graphics.off()

###### Clear environment ######
# rm(list=ls())

###### Set working directory ######
# this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
# setwd(this.dir)


###### Packages ######
# install.packages("Seurat")
library(Seurat)
# install.packages("data.table")
library(data.table)
# install.packages("dplyr")
library(dplyr)
# install.packages("BiocManager")
library(BiocManager)
# BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)
# BiocManager::install("scDblFinder")
library(scDblFinder)
# BiocManager::install("decontX")
library(decontX)

###### Data read-in and preprocess ######
scRNA_Preprocess = function(matrix, Duplicate.Genes = "Keep", # Duplicate.Genes has three options: "Keep", "Merge", and "Remove".
                            Denoise = TRUE, Remove.Doublets = TRUE, # Choose whether to perform denoising and remove doublets.
                            nCells.percentage = FALSE, # Choose whether to use percentage of total cells as thereshold.
                            nCells.lowerpercent = 0, nCells.lower = 0, # Filter out genes expressed in < nCell.lowerpercent * number of cells/100 or nCell.lower cells.
                            nFeatures.MAD = TRUE, # Choose whether to use MAD methods to filter cells.
                            nFeatures.nMAD = 5, # The number of MADs below/above median to be considered as outliers.
                            nFeatures.fixed = TRUE, # Choose whether to use fixed thresholds.
                            nFeatures.lower = 0, nFeatures.upper = Inf, # Filter out cells with expressed genes < nFeatures.lower and > nFeatures.upper.
                            nCounts.MAD = TRUE, # Choose whether to use MAD methods to filter cells.
                            nCounts.nMAD = 5, # The number of MADs below/above median to be considered as outliers.
                            nCounts.fixed = TRUE, # Choose whether to use fixed thresholds.
                            nCounts.lower = 0, nCounts.upper = Inf, # Filter out cells with expressed genes < nCounts.lower and > nCounts.upper.
                            mt.MAD = TRUE, # Choose whether to use MAD methods to filter cells.
                            mt.nMAD = 3, # The number of MADs below/above median to be considered as outliers.
                            mt.fixed = TRUE, # Choose whether to use fixed thresholds or MAD methods.
                            mt.upper = 10){ # Filter out cells with percentage of reads that map to the mitochondrial genome > mt.upper.
                            
  #### Read in matrix ####
  message("Converting the matrix to sparse format...")
  matrix = as(matrix, "dgCMatrix")
  
  #### Duplicate genes ####
  if (Duplicate.Genes != "Keep"){
    gene_names = rownames(matrix)
    if (Duplicate.Genes == "Merge"){
      message("Merging duplicate genes...")
      base_gene = sub("\\.\\d+$", "", gene_names)
      gene_index_map = lapply(unique(base_gene), function(gene){
        which(base_gene == gene)
      })
      matrix = do.call(rbind, lapply(gene_index_map, function(index){
        if (length(index) > 1){
          colSums(matrix[index, , drop = FALSE])
        } else {
          matrix[index, , drop = FALSE]
        }
      }))
      rownames(matrix) = unique(base_gene)
      matrix = as(matrix, "dgCMatrix")
    } else if (Duplicate.Genes == "Remove"){
      message("Removing duplicate genes...")
      base_gene = !grepl("\\.\\d+$", gene_names)
      matrix = matrix[base_gene, ]
      matrix = as(matrix, "dgCMatrix")
    }
  }
  
  #### Remove empty cells ####
  message("Removing empty cells...")
  sce = SingleCellExperiment(list(counts = matrix))
  lib_sizes = colSums(counts(sce))
  if (sum(lib_sizes == 0) > 0){
    message(sum(lib_sizes == 0), " empty droplets detected.")
  }
  sce = sce[, lib_sizes > 0]
  
  #### Denoise ####
  if (Denoise){
    message("Denoising...")
    sce = decontX(sce) # Remove contamination from ambient RNA and technical noise.
  }
  
  #### Detect doublets ####
  if (Remove.Doublets){
    message("Detecting doublets...")
    sce = scDblFinder(sce)
  }
  
  #### Create Seurat object and QC filtering ####
  ## Filter out genes expressed in fewer than nCells.lower cells ##
  message("Filtering genes...")
  if (nCells.percentage){
    nCells.lower = ceiling(nCells.lowerpercent * dim(matrix)[2] / 100)
  }
  if (Denoise){
    srt = CreateSeuratObject(counts = round(decontXcounts(sce)), min.cells = nCells.lower)
  } else {
    srt = CreateSeuratObject(counts = round(counts(sce)), min.cells = nCells.lower)
  }
  
  
  ## Check if MAD methods is used ##
  if (nFeatures.MAD | nCounts.MAD | mt.MAD){
    is_outlier_MAD = function(vec, nMAD){
      is_outlier = vec < median(vec) - nMAD * mad(vec) | vec > median(vec) + nMAD * mad(vec)
      return(is_outlier)
    }
  }
  
  ## Filter out cells ##
  message("Filtering cells...")
  if (nFeatures.MAD){
    srt$outlier = is_outlier_MAD(srt$nFeature_RNA, nFeatures.nMAD)
    srt = subset(srt, subset = outlier == FALSE)
    srt$outlier = NULL
  }
  if (nFeatures.fixed){
    srt = subset(srt, subset = nFeature_RNA > nFeatures.lower & nFeature_RNA < nFeatures.upper)
  }
  
  if (nCounts.MAD){
    srt$outlier = is_outlier_MAD(srt$nCount_RNA, nCounts.nMAD)
    srt = subset(srt, subset = outlier == FALSE)
    srt$outlier = NULL
  }
  if (nCounts.fixed){
    srt = subset(srt, subset = nCount_RNA > nCounts.lower & nCount_RNA < nCounts.upper)
  }
  
  ## Remove doublets ##
  if (Remove.Doublets){
    srt = AddMetaData(srt, sce$scDblFinder.class, col.name = "scDblFinder.class")
    # srt = AddMetaData(srt, sce$scDblFinder.score, col.name = "scDblFinder.score")
    # srt = AddMetaData(srt, sce$scDblFinder.weighted, col.name = "scDblFinder.weighted")
    # srt = AddMetaData(srt, sce$scDblFinder.cxds_score, col.name = "scDblFinder.cxds_score")
    message("Removing doublets...")
    srt = subset(srt, subset = scDblFinder.class == "singlet")
  }
  
  ## Remove low-quality and dying cells ##
  message("Removing low-quality and dying cells...") # Low-quality and dying cells often exhibit extensive mitochondrial contamination.
  srt[["percent.mt"]] = PercentageFeatureSet(srt, pattern = "^MT-", assay = "RNA")
  if (mt.MAD){
    srt$outlier = is_outlier_MAD(srt$percent.mt, mt.nMAD)
    srt = subset(srt, subset = outlier == FALSE)
    srt$outlier = NULL
  }
  if (mt.fixed){
    srt = subset(srt, subset = percent.mt < mt.upper)
  }
  
  #### Processed matrix output ####
  # message("Saving Seurat object...")
  # matrix = srt@assays$RNA$counts
  message("Preprocess finished.")
  return(srt)
  # saveRDS(srt, paste0(output_path, "Seurat.rds"))
  # message("Seurat object available in 'Preprocessed' folder.")
  
  #### Free up memory ####
  rm(matrix, lib_sizes, sce, gene_names, base_gene, gene_index_map)
  gc()
  
}
