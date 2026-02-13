# Single-nuclei RNA-seq and matched CAG sizing - ELongATE

# Index

[Load packages](#Load_packages)

[Read Handsaker metadata](#Read_Handsaker_metadata)
    
[Read Handsaker expression matrix](#Read_Handsaker_expression_matrix)

[Filter Handsaker SPN UMAP](#Filter_Handsaker_SPN_UMAP)

[Read phase C genes](#Read_phaseC_genes)

[Learn phase model](#Learn_phase_model)
    
[Learn CAG sizing model](#Learn_CAG_sizing_model)

[Import Lee dataset](#Import_Lee_dataset)
    
[Run model on Lee MSNs dataset from Caudate](#Lee_Caudate_only)
    
[Run model on Lee MSNs dataset from Putamen](#Lee_Putamen_only)
    
[Import Paryani dataset](#Import_Paryani_dataset)    
    
[Run model on Paryani MSNs dataset from Caudate](#Paryani_Caudate_only)
    
[Run model on Paryani MSNs dataset from Accumbens](#Paryani_Accumbens_only)
    
[Plot predictions on controls](#Plot_controls)
    
[Import Xu dataset](#Import_Xu_dataset)

[Identify SPNs in Xu dataset](#Identify_SPNs_in_Xu_dataset)

[Run model on Xu dataset](#Run_model_on_Xu_dataset)

[All datasets](#All_datasets)

# Load_packages


```R
library("dplyr")
library("Seurat")
library("patchwork")
library("Matrix")
library("biomaRt")
library("ggplot2")
library("celda")
library("DoubletFinder")
library("hdf5r")
library("harmony")
library("presto")
library("neuralnet")
library("keras")
library("tensorflow")
library("glmnet")
library("caret")
library("descriptio")
library("gridExtra")
library("ggpubr")
library("weights")
library("patchwork")
library("ggplotify")
library("PRROC")
setwd("/mnt/projects/labs/CLAB/PROJECT_ELongATE/results_rebuttal/")
```

# Read_Handsaker_metadata


```R
pb_metadata_file <- "/mnt/projects/labs/CLAB/PROJECT_ELongATE/Handsaker/Analysis_bag_6_Broad_Huntingtons_Caudate_2024_PacBio_Metadata_Open/pacbio_deepdives_cells_unfiltered.txt"
pb_metadata <- read.table(pb_metadata_file, sep = "\t", header = TRUE)
sample_name <- gsub(pattern = '(.*)_\\w+', replacement = '\\1', pb_metadata$CELL_BARCODE)
CB <- gsub(pattern = '.*_(\\w+)', replacement = '\\1', pb_metadata$CELL_BARCODE)
pb_metadata <- cbind(pb_metadata, CB, sample_name)
rownames(pb_metadata) <- gsub(x = pb_metadata$CELL_BARCODE, pattern = "Caudate", replacement = "HD_Caudate")
```


```R
tmp <- lapply(split(pb_metadata, pb_metadata$DONOR), function(x) {
    if (x$DONOR[1] != "S04577") {
      x <- x[which(x$CAGLENGTH > 36), ]    
    }
    num_cells <- dim(x)[1]
    num_spns <- dim(x[which(x$CELLTYPE == "SPN"), ])[1]
    return(data.frame(num_cells, num_spns))
})

num_cells_Handsaker <- do.call(rbind, tmp)
#num_cells_Handsaker
```


```R
num_cells_Handsaker
```

# Read_Handsaker_expression_matrix


```R
h5_files <- list.files(path = "/mnt/projects/labs/CLAB/PROJECT_ELongATE/Handsaker/", pattern = "\\.h5$", full.names = TRUE)
names(h5_files) <- gsub(x = basename(h5_files), pattern = "\\.umi.*", replacement = "")
```


```R
Read_h5_file = function(h5_file, metadata) {
    #extract counts
    counts <- Read10X_h5(filename = h5_file, use.names = TRUE)
    #create a single cell experiment object
    sce <- SingleCellExperiment(list(counts = counts))
    #create a seurat object
    seuratObject <- CreateSeuratObject(counts(sce), project = names(h5_file), min.cells = 0, min.features = 0)
    seuratObject[["CB"]] <- rownames(seuratObject@meta.data)
    #print(head(rownames(metadata)))
    #print(head(seuratObject@meta.data))
    #create tmp variable and edit rownames
    tmp <- seuratObject@meta.data
    rownames(tmp) <- gsub(x = gsub(x = paste0(seuratObject@meta.data$orig.ident, "_",
                                              gsub(x = seuratObject@meta.data$CB, pattern = "_.*", replacement = "")), 
                                   pattern = "_10X_", replacement = "_"), 
                          pattern = "_merged_", replacement = "_")
    #exploit edited rownames to add metadata to Seurat object
    tmp <- cbind(tmp, metadata[rownames(tmp), ])
    #print(head(rownames(tmp)))
    #replace row names with old ones
    rownames(tmp) <- rownames(seuratObject@meta.data)
    #update original object
    seuratObject@meta.data <- tmp
    #print(head(seuratObject@meta.data))
    #print(head(rownames(seuratObject@meta.data)))
    #retain only SPN
    seuratObject <- subset(x = seuratObject, subset = CELLTYPE == "SPN")
    #run SCTransform on SPN only
    seuratObject <- SCTransform(seuratObject)
    cat(sprintf("%s; Num. features = %d; Num. SPN = %d\n", names(h5_file), dim(seuratObject)[1], dim(seuratObject)[2]))
    return(seuratObject)
}
```


```R
names(h5_files)
```


```R
# S02205_HD_Caudate_DeepDive_rxn1 <- Read_h5_file(h5_files["S02205_HD_Caudate_DeepDive_rxn1"], pb_metadata)
# S02205_HD_Caudate_DeepDive_rxn2 <- Read_h5_file(h5_files["S02205_HD_Caudate_DeepDive_rxn2"], pb_metadata)
# S02205_HD_Caudate_DeepDive_rxn3 <- Read_h5_file(h5_files["S02205_HD_Caudate_DeepDive_rxn3"], pb_metadata)
# S02205_HD_Caudate_DeepDive_rxn4 <- Read_h5_file(h5_files["S02205_HD_Caudate_DeepDive_rxn4"], pb_metadata)

# S02205 <- merge(x = S02205_HD_Caudate_DeepDive_rxn1, 
#                 y = c(S02205_HD_Caudate_DeepDive_rxn2, 
#                       S02205_HD_Caudate_DeepDive_rxn3,
#                       S02205_HD_Caudate_DeepDive_rxn4))

# saveRDS(object = S02205, file = "S02205.rds")
```


```R
# S04002_10X_rxn1 <- Read_h5_file(h5_files["S04002_10X_rxn1"], pb_metadata)
# S04002_10X_rxn2_merged <- Read_h5_file(h5_files["S04002_10X_rxn2_merged"], pb_metadata)
# S04002_10X_rxn3_merged <- Read_h5_file(h5_files["S04002_10X_rxn3_merged"], pb_metadata)
# S04002_10X_rxn4_merged <- Read_h5_file(h5_files["S04002_10X_rxn4_merged"], pb_metadata)
# S04002_rxn5 <- Read_h5_file(h5_files["S04002_rxn5"], pb_metadata)
# S04002_rxn6 <- Read_h5_file(h5_files["S04002_rxn6"], pb_metadata)
# S04002_rxn7 <- Read_h5_file(h5_files["S04002_rxn7"], pb_metadata)
# S04002_rxn8 <- Read_h5_file(h5_files["S04002_rxn8"], pb_metadata)

# S04002 <- merge(x = S04002_10X_rxn1,
#                 y = c(S04002_10X_rxn2_merged,
#                      S04002_10X_rxn3_merged,
#                      S04002_10X_rxn4_merged,
#                      S04002_rxn5,
#                      S04002_rxn6,
#                      S04002_rxn7,
#                      S04002_rxn8))

# saveRDS(object = S04002, file = "S04002.rds")
```


```R
# S04577_HD_Caudate_DeepDive_rxn1 <- Read_h5_file(h5_files["S04577_HD_Caudate_DeepDive_rxn1"], pb_metadata)
# S04577_HD_Caudate_DeepDive_rxn2 <- Read_h5_file(h5_files["S04577_HD_Caudate_DeepDive_rxn2"], pb_metadata)
# S04577_HD_Caudate_DeepDive_rxn3 <- Read_h5_file(h5_files["S04577_HD_Caudate_DeepDive_rxn3"], pb_metadata)
# S04577_HD_Caudate_DeepDive_rxn4 <- Read_h5_file(h5_files["S04577_HD_Caudate_DeepDive_rxn4"], pb_metadata)

# S04577 <- merge(x = S04577_HD_Caudate_DeepDive_rxn1,
#                 y = c(S04577_HD_Caudate_DeepDive_rxn2,
#                      S04577_HD_Caudate_DeepDive_rxn3,
#                      S04577_HD_Caudate_DeepDive_rxn4))

# saveRDS(object = S04577, file = "S04577.rds")
```


```R
# S05202_HD_Caudate_DeepDive_rxn1 <- Read_h5_file(h5_files["S05202_HD_Caudate_DeepDive_rxn1"], pb_metadata)
# S05202_HD_Caudate_DeepDive_rxn2 <- Read_h5_file(h5_files["S05202_HD_Caudate_DeepDive_rxn2"], pb_metadata)
# S05202_HD_Caudate_DeepDive_rxn3 <- Read_h5_file(h5_files["S05202_HD_Caudate_DeepDive_rxn3"], pb_metadata)

# S05202 <- merge(x = S05202_HD_Caudate_DeepDive_rxn1,
#                 y = c(S05202_HD_Caudate_DeepDive_rxn2,
#                      S05202_HD_Caudate_DeepDive_rxn3))

# saveRDS(object = S05202, file = "S05202.rds")
```


```R
# S05368_HD_Caudate_DeepDive_rxn1 <- Read_h5_file(h5_files["S05368_HD_Caudate_DeepDive_rxn1"], pb_metadata)
# S05368_HD_Caudate_DeepDive_rxn2 <- Read_h5_file(h5_files["S05368_HD_Caudate_DeepDive_rxn2"], pb_metadata)
# S05368_HD_Caudate_DeepDive_rxn3 <- Read_h5_file(h5_files["S05368_HD_Caudate_DeepDive_rxn3"], pb_metadata)

# S05368 <- merge(x = S05368_HD_Caudate_DeepDive_rxn1,
#                 y = c(S05368_HD_Caudate_DeepDive_rxn2,
#                      S05368_HD_Caudate_DeepDive_rxn3))

# saveRDS(object = S05368, file = "S05368.rds")
```


```R
# S06758_HD_Caudate_DeepDive_rxn1 <- Read_h5_file(h5_files["S06758_HD_Caudate_DeepDive_rxn1"], pb_metadata)
# S06758_HD_Caudate_DeepDive_rxn2 <- Read_h5_file(h5_files["S06758_HD_Caudate_DeepDive_rxn2"], pb_metadata)
# S06758_HD_Caudate_DeepDive_rxn3 <- Read_h5_file(h5_files["S06758_HD_Caudate_DeepDive_rxn3"], pb_metadata)
# S06758_HD_Caudate_DeepDive_rxn4 <- Read_h5_file(h5_files["S06758_HD_Caudate_DeepDive_rxn4"], pb_metadata)
# S06758_HD_Caudate_DeepDive_rxn5 <- Read_h5_file(h5_files["S06758_HD_Caudate_DeepDive_rxn5"], pb_metadata)

# S06758 <- merge(x = S06758_HD_Caudate_DeepDive_rxn1,
#                 y = c(S06758_HD_Caudate_DeepDive_rxn2,
#                      S06758_HD_Caudate_DeepDive_rxn3,
#                      S06758_HD_Caudate_DeepDive_rxn4,
#                      S06758_HD_Caudate_DeepDive_rxn5))

# saveRDS(object = S06758, file = "S06758.rds")
```


```R
# S07681_HD_Caudate_DeepDive_rxn1 <- Read_h5_file(h5_files["S07681_HD_Caudate_DeepDive_rxn1"], pb_metadata)
# S07681_HD_Caudate_DeepDive_rxn2 <- Read_h5_file(h5_files["S07681_HD_Caudate_DeepDive_rxn2"], pb_metadata)
# S07681_HD_Caudate_DeepDive_rxn3 <- Read_h5_file(h5_files["S07681_HD_Caudate_DeepDive_rxn3"], pb_metadata)
# S07681_HD_Caudate_DeepDive_rxn4 <- Read_h5_file(h5_files["S07681_HD_Caudate_DeepDive_rxn4"], pb_metadata)
# S07681_HD_Caudate_DeepDive_rxn5 <- Read_h5_file(h5_files["S07681_HD_Caudate_DeepDive_rxn5"], pb_metadata)
# S07681_HD_Caudate_DeepDive_rxn6 <- Read_h5_file(h5_files["S07681_HD_Caudate_DeepDive_rxn6"], pb_metadata)

# S07681 <- merge(x = S07681_HD_Caudate_DeepDive_rxn1,
#                 y = c(S07681_HD_Caudate_DeepDive_rxn2,
#                      S07681_HD_Caudate_DeepDive_rxn3,
#                      S07681_HD_Caudate_DeepDive_rxn4,
#                      S07681_HD_Caudate_DeepDive_rxn5,
#                      S07681_HD_Caudate_DeepDive_rxn6))

# saveRDS(object = S07681, file = "S07681.rds")
```


```R
# S09619_HD_Caudate_DeepDive_rxn1 <- Read_h5_file(h5_files["S09619_HD_Caudate_DeepDive_rxn1"], pb_metadata)
# S09619_HD_Caudate_DeepDive_rxn2 <- Read_h5_file(h5_files["S09619_HD_Caudate_DeepDive_rxn2"], pb_metadata)
# S09619_HD_Caudate_DeepDive_rxn3 <- Read_h5_file(h5_files["S09619_HD_Caudate_DeepDive_rxn3"], pb_metadata)
# S09619_HD_Caudate_DeepDive_rxn4 <- Read_h5_file(h5_files["S09619_HD_Caudate_DeepDive_rxn4"], pb_metadata)

# S09619 <- merge(x = S09619_HD_Caudate_DeepDive_rxn1,
#                 y = c(S09619_HD_Caudate_DeepDive_rxn2,
#                      S09619_HD_Caudate_DeepDive_rxn3,
#                      S09619_HD_Caudate_DeepDive_rxn4))

# saveRDS(object = S09619, file = "S09619.rds")
```


```R
# S12365_HD_Caudate_DeepDive_rxn1 <- Read_h5_file(h5_files["S12365_HD_Caudate_DeepDive_rxn1"], pb_metadata)
# S12365_HD_Caudate_DeepDive_rxn2 <- Read_h5_file(h5_files["S12365_HD_Caudate_DeepDive_rxn2"], pb_metadata)
# S12365_HD_Caudate_DeepDive_rxn3 <- Read_h5_file(h5_files["S12365_HD_Caudate_DeepDive_rxn3"], pb_metadata)
# S12365_HD_Caudate_DeepDive_rxn4 <- Read_h5_file(h5_files["S12365_HD_Caudate_DeepDive_rxn4"], pb_metadata)

# S12365 <- merge(x = S12365_HD_Caudate_DeepDive_rxn1,
#                 y = c(S12365_HD_Caudate_DeepDive_rxn2,
#                      S12365_HD_Caudate_DeepDive_rxn3,
#                      S12365_HD_Caudate_DeepDive_rxn4))

# saveRDS(object = S12365, file = "S12365.rds")
```


```R
# SPN_all <- merge(x = S02205_HD_Caudate_DeepDive_rxn1, 
#                 y = c(S02205_HD_Caudate_DeepDive_rxn2, 
#                       S02205_HD_Caudate_DeepDive_rxn3,
#                       S02205_HD_Caudate_DeepDive_rxn4,
#                       S04002_10X_rxn1,
#                       S04002_10X_rxn2_merged,
#                       S04002_10X_rxn3_merged,
#                       S04002_10X_rxn4_merged,
#                       S04002_rxn5,
#                       S04002_rxn6,
#                       S04002_rxn7,
#                       S04002_rxn8,
#                       S04577_HD_Caudate_DeepDive_rxn1,
#                       S04577_HD_Caudate_DeepDive_rxn2,
#                       S04577_HD_Caudate_DeepDive_rxn3,
#                       S04577_HD_Caudate_DeepDive_rxn4,
#                       S05202_HD_Caudate_DeepDive_rxn1,
#                       S05202_HD_Caudate_DeepDive_rxn2,
#                       S05202_HD_Caudate_DeepDive_rxn3,
#                       S05368_HD_Caudate_DeepDive_rxn1,
#                       S05368_HD_Caudate_DeepDive_rxn2,
#                       S05368_HD_Caudate_DeepDive_rxn3,
#                       S06758_HD_Caudate_DeepDive_rxn1,
#                       S06758_HD_Caudate_DeepDive_rxn2,
#                       S06758_HD_Caudate_DeepDive_rxn3,
#                       S06758_HD_Caudate_DeepDive_rxn4,
#                       S06758_HD_Caudate_DeepDive_rxn5,
#                       S07681_HD_Caudate_DeepDive_rxn1,
#                       S07681_HD_Caudate_DeepDive_rxn2,
#                       S07681_HD_Caudate_DeepDive_rxn3,
#                       S07681_HD_Caudate_DeepDive_rxn4,
#                       S07681_HD_Caudate_DeepDive_rxn5,
#                       S07681_HD_Caudate_DeepDive_rxn6,
#                       S09619_HD_Caudate_DeepDive_rxn1,
#                       S09619_HD_Caudate_DeepDive_rxn2,
#                       S09619_HD_Caudate_DeepDive_rxn3,
#                       S09619_HD_Caudate_DeepDive_rxn4,
#                       S12365_HD_Caudate_DeepDive_rxn1,
#                       S12365_HD_Caudate_DeepDive_rxn2,
#                       S12365_HD_Caudate_DeepDive_rxn3,
#                       S12365_HD_Caudate_DeepDive_rxn4))

# #Add metadata
# SPN_all[["percent.mt"]] <- PercentageFeatureSet(SPN_all, pattern = "^MT-")
# SPN_all[["SAMPLE"]] <- factor(unlist(lapply(strsplit(SPN_all@meta.data$sample_name, "_"), '[[', 1)))
# SPN_all[["EXP_CAGLENGTH"]] <- NA
# SPN_all[["PHASE"]] <- NA
# SPN_all@meta.data$EXP_CAGLENGTH[which(SPN_all@meta.data$CAGLENGTH > 37)] <- SPN_all@meta.data$CAGLENGTH[which(SPN_all@meta.data$CAGLENGTH > 37)]
# SPN_all@meta.data$PHASE[which(SPN_all@meta.data$EXP_CAGLENGTH < 150)] <- "A-B"
# SPN_all@meta.data$PHASE[which(SPN_all@meta.data$EXP_CAGLENGTH >= 150)] <- "C-D-E"

# head(SPN_all@meta.data)

# # SPN_all <- merge(x = S02205,
# #                 y = c(S04002,
# #                      S04577,
# #                      S05202,
# #                      S05368,
# #                      S06758,
# #                      S07681,
# #                      S09619,
# #                      S12365),
# #                 add.cell.ids = c("S02205", "S04002", "S04577", "S05202", "S05368", "S06758", "S07681", "S09619", "S12365"),
# #                  project = "snRNAseqCAG")

# saveRDS(object = SPN_all, file = "/mnt/projects/labs/CLAB/PROJECT_ELongATE/SPN_all_raw.rds")
```


```R
#load dataset with SPNs from all samples
SPN_all <- readRDS("/mnt/projects/labs/CLAB/PROJECT_ELongATE/SPN_all_raw.rds")
```


```R
CONDITION <- factor(c("HD", "HD", "CTRL", "HD", "HD_PREM", "HD", "HD", "HD_PREM", "HD"))
names(CONDITION) <- c("S02205", "S04002", "S04577", "S05202", "S05368", "S06758", "S07681" , "S09619", "S12365")
```


```R
SPN_all@meta.data$CONDITION <- CONDITION[SPN_all@meta.data$DONOR]
```


```R
SPN_all@meta.data$PHASE[which(SPN_all@meta.data$CONDITION == "CTRL")] <- "A-B"
```


```R
table(SPN_all@meta.data$DONOR)
```


```R
str(SPN_all@meta.data)
```


```R
sort(table(SPN_all@meta.data$orig.ident))
#sort(table(SPN_all@meta.data$orig.ident[grep(x = SPN_all@meta.data$orig.ident, pattern = "S05202")]))
```


```R
#find total number of SPNs and ratio of SPNs with Wt-allele only
tot <- table(unlist(lapply(strsplit(SPN_all@meta.data$sample_name, "_"), '[[', 1)))
num_WT <- table(unlist(lapply(strsplit(SPN_all@meta.data[which(SPN_all@meta.data$CAGLENGTH < 36), "sample_name"], "_"), '[[', 1)))
num_HD <- table(unlist(lapply(strsplit(SPN_all@meta.data[which(SPN_all@meta.data$CAGLENGTH > 36), "sample_name"], "_"), '[[', 1)))
num_HD["S04577"] <- 0
num_HD_exp <- table(unlist(lapply(strsplit(SPN_all@meta.data[which(SPN_all@meta.data$CAGLENGTH > 150), "sample_name"], "_"), '[[', 1)))
num_HD_exp["S04577"] <- 0

#tot
#num_WT/tot[names(num_WT)]
#num_HD/tot[names(num_HD)]
#num_HD_exp/tot[names(num_HD_exp)]
#num_HD_exp/num_HD[names(num_HD_exp)]

df <- data.frame(DONOR = names(tot),
                 CONDITION = CONDITION[names(tot)],
                 NUM_SPN = as.numeric(tot),
                 NUM_HD_SPN = as.numeric(num_HD[names(tot)]),
                 FRACT_HD_SPN = sprintf("%.2f%%", 100*as.numeric(num_HD[names(tot)])/as.numeric(tot)),
                 NUM_HD_SPN_CDE = as.numeric(num_HD_exp[names(tot)]),
                 FRACT_HD_SPN_CDE = sprintf("%.2f%%", 100*as.numeric(num_HD_exp[names(tot)])/as.numeric(num_HD[names(tot)])))
df["S04577", "FRACT_HD_SPN_CDE"] <- "0.00%"
#df
df$CONDITION <- factor(df$CONDITION, levels = c("CTRL", "HD_PREM", "HD"))
#df
df_sort <- df[order(df$CONDITION), ]
df_sort
```


```R
df_sort$FRACT_HD_SPN_CDE
mean(as.numeric(gsub(x = df_sort$FRACT_HD_SPN_CDE[-c(1)], pattern = "%", replacement = "")))
mean(as.numeric(gsub(x = df_sort$FRACT_HD_SPN_CDE[-c(1, 2, 3)], pattern = "%", replacement = "")))
```


```R
# head(S02205@meta.data)
# head(S04002@meta.data)
# head(S04577@meta.data)
# head(S05202@meta.data)
# head(S05368@meta.data)
# head(S06758@meta.data)
# head(S07681@meta.data)
# head(S09619@meta.data)
# head(S12365@meta.data)
head(SPN_all@meta.data)
SPN_all
```

# Filter_Handsaker_SPN_UMAP


```R
# #Run SCTransform
# SPN_all <- SCTransform(SPN_all, variable.features.n = 3000,  ncells = 5000)
# #Run PCA
# SPN_all <- RunPCA(SPN_all, features = VariableFeatures(object = SPN_all), npcs = 50, reduction.name = "pca")
# ElbowPlot(SPN_all, reduction = "pca", ndims = 50)
# num_dims <- 40
# #Run UMAP
# SPN_all <- RunUMAP(SPN_all, dims = 1:num_dims, reduction = "pca", reduction.name = "umap", reduction.key = "umap")
```


```R
# #do clustering
# SPN_all <- FindNeighbors(SPN_all, reduction = "pca", dims = 1:num_dims)
# SPN_all <- FindClusters(SPN_all, resolution = 0.01, method = 4)
# DimPlot(SPN_all, reduction = "umap", group.by = "seurat_clusters", pt.size = 1)
# ggsave("UMAP_Handsaker_clusters.pdf", width = 8, height = 8)
```


```R
# table(SPN_all@meta.data$seurat_clusters)
```


```R
# #find gene markers for the two clusters
# SPN_all <- PrepSCTFindMarkers(SPN_all, assay = "SCT", verbose = TRUE)

# # find markers for every cluster compared to all remaining cells, report only the positive ones
# SPN_all.markers <- FindAllMarkers(SPN_all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")

# #find the top marker for each cluster and plot their expression values
# SPN_all_top_marker <- SPN_all.markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 1, order_by = avg_log2FC)

####################################
# num_genes_plot <- 10
# if (length(SPN_all_top_marker) < num_genes_plot) {
#   chunks_genes <- list(1:length(SPN_all_top_marker))
# } else {
#   chunks_genes <- split(1:length(SPN_all_top_marker), ceiling(seq(from = 1, to = length(SPN_all_top_marker))/num_genes_plot))
# }

# #plot normalized counts
# for (i in 1:length(chunks_genes)) {
#     p <- VlnPlot(SPN_all, features = SPN_all_top_marker$gene[chunks_genes[[i]]]) 
#     #+ ggtitle("Normalized counts for top markers")
#     plot(p)
#     ggsave(paste0("Handsaker_Normalized_MarkerGeneCounts_chunk", i, ".pdf"), device = pdf, height = 8, width = 8)
# }
```


```R
# saveRDS(object = SPN_all, file = "/mnt/projects/labs/CLAB/PROJECT_ELongATE/SPN_all.rds")
```


```R
SPN_all <- readRDS("/mnt/projects/labs/CLAB/PROJECT_ELongATE/SPN_all.rds")
```


```R
table(SPN_all@meta.data$seurat_clusters)
```


```R
#extract top 10 markers for each cluster and plot heatmap
# find markers for every cluster compared to all remaining cells, report only the positive ones
SPN_all.markers <- FindAllMarkers(SPN_all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")

#find the top marker for each cluster and plot their expression values
SPN_all_top_marker <- SPN_all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC)
SPN_all_top_10_markers <- SPN_all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

p <- DoHeatmap(SPN_all, features = SPN_all_top_10_markers$gene) + NoLegend() + ggtitle("Heatmap with top 10 marker genes per cluster")
plot(p)
ggsave("Handsaker_Heatmap_top10MarkerGeneCounts.pdf", device = pdf, height = 8, width = 16)
```


```R
SPN_all_top_10_markers
```


```R
#do plots
options(repr.plot.width=20, repr.plot.height=20) 
DimPlot(SPN_all, reduction = "umap", group.by = "seurat_clusters", pt.size = 1)
ggsave("UMAP_Handsaker_clusters.pdf", width = 8, height = 8)
FeaturePlot(SPN_all, reduction = "umap", features = "CAGLENGTH", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Handsaker_CAGLENGTH.pdf", width = 8, height = 8)
FeaturePlot(SPN_all, reduction = "umap", features = "EXP_CAGLENGTH", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Handsaker_EXP_CAGLENGTH.pdf", width = 8, height = 8)
DimPlot(SPN_all, reduction = "umap", group.by = "PHASE", pt.size = 1)
ggsave("UMAP_Handsaker_PHASE.pdf", width = 8, height = 8)
FeaturePlot(SPN_all, reduction = "umap", features = "nCount_RNA", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Handsaker_NCOUNTS.pdf", width = 8, height = 8)
FeaturePlot(SPN_all, reduction = "umap", features = "nFeature_RNA", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Handsaker_NFEATURES.pdf", width = 8, height = 8)
FeaturePlot(SPN_all, reduction = "umap", features = "percent.mt", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Handsaker_PERCMT.pdf")
DimPlot(SPN_all, group.by = "SAMPLE", reduction = "umap", pt.size = 1)
ggsave("UMAP_Handsaker_splitBySampleMerged.pdf", width = 8, height = 8)
DimPlot(SPN_all, split.by = "SAMPLE", reduction = "umap", pt.size = 1, ncol = 3, group.by = "SAMPLE") + NoLegend()
ggsave("UMAP_Handsaker_splitBySample.pdf", width = 20, height = 20)
```


```R
#read phase C-D file with all genes
phaseCD_all_file <- "/mnt/projects/labs/CLAB/PROJECT_ELongATE/PhaseCD_all.txt"

phaseCD_genes_all <- read.table(phaseCD_all_file, sep = "\t", header = TRUE)
rownames(phaseCD_genes_all) <- phaseCD_genes_all$Gene

phaseCD_genes_all[intersect(rownames(phaseCD_genes_all), SPN_all_top_10_markers[which(SPN_all_top_10_markers$cluster == 1), ]$gene), ]
```


```R
#create file with only phase C-, C+ and D genes
phaseCD_file <- "/mnt/projects/labs/CLAB/PROJECT_ELongATE/PhaseCD.txt"
system(command = paste0("cat ", phaseCD_all_file, " | awk 'BEGIN {FS=\"\t\"} { if ($8 != NA) print }' > ", phaseCD_file))
```


```R
#avoid discarding cluster 1; the following code reports the code for filtering out cluster 1
# #filter only SPNs from cluster 0
# SPN_all_unfiltered <- SPN_all
# SPN_filt <- subset(x = SPN_all, subset = seurat_clusters == "0")
# saveRDS(object = SPN_filt, file = "/mnt/projects/labs/CLAB/PROJECT_ELongATE/SPN_filtered.rds")
# SPN_filt <- readRDS("/mnt/projects/labs/CLAB/PROJECT_ELongATE/SPN_filtered.rds")
```


```R
SPN_filt <- SPN_all
```


```R
num_dims <- 40
SPN_filt <- RunPCA(SPN_filt, features = VariableFeatures(object = SPN_filt), npcs = 50, reduction.name = "pca")
SPN_filt <- RunUMAP(SPN_filt, dims = 1:num_dims, reduction = "pca", reduction.name = "umap", reduction.key = "umap")
```


```R
#do plots for CAG size of HD samples
p <- ggplot(SPN_filt@meta.data, aes(x = CAGLENGTH)) +
geom_histogram(binwidth = 5, alpha = 0.5, position = "identity", fill = "#F8766D") +
geom_vline(aes(xintercept = 150), colour="black", linetype = 2) +
theme_classic()
#   theme(legend.title = element_text(size = 10),
#         legend.text = element_text(size = 10),
#         axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
#         axis.text.y = element_text(size = 10),
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         panel.background = element_blank())
plot(p)
ggsave("Handsaker_CAGLENGTH_hist.pdf", width = 8, height = 8)

p <- ggplot(SPN_filt@meta.data, aes(x = EXP_CAGLENGTH)) +
geom_histogram(binwidth = 5, alpha = 0.5, position = "identity", fill = "#F8766D") +
geom_vline(aes(xintercept = 150), colour="black", linetype = 2) +
theme_classic()
#   theme(legend.title = element_text(size = 10),
#         legend.text = element_text(size = 10),
#         axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
#         axis.text.y = element_text(size = 10),
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         panel.background = element_blank())
plot(p)
ggsave("Handsaker_CAGLENGTH_EXP_hist.pdf", width = 8, height = 8)
p <- ggplot(SPN_filt@meta.data, aes(x = EXP_CAGLENGTH, fill = SAMPLE)) +
geom_histogram(binwidth = 5, alpha = 0.5, position = "identity") + 
geom_vline(aes(xintercept = 150), colour="black", linetype = 2) +
theme_classic()
#   theme(legend.title = element_text(size = 10),
#         legend.text = element_text(size = 10),
#         axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
#         axis.text.y = element_text(size = 10),
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         panel.background = element_blank())
plot(p)
ggsave("Handsaker_CAGLENGTH_hist_splitBySample.pdf", width = 8, height = 8)

p <- ggplot(SPN_filt@meta.data, aes(x = EXP_CAGLENGTH, fill = SAMPLE)) +
geom_histogram(binwidth = 5, alpha = 0.5, position = "identity") + 
geom_vline(aes(xintercept = 150), colour="black", linetype = 2) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank())
plot(p)
ggsave("Handsaker_CAGLENGTH_EXP_hist_splitBySample.pdf", width = 8, height = 8)

p <- ggplot(SPN_filt@meta.data, aes(x = EXP_CAGLENGTH)) +
geom_histogram(binwidth = 5, alpha = 0.5, position = "identity", fill = "#F8766D") +
xlim(c(150, 1000)) +
theme_classic()
#   theme(legend.title = element_text(size = 10),
#         legend.text = element_text(size = 10),
#         axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
#         axis.text.y = element_text(size = 10),
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         panel.background = element_blank())
plot(p)
ggsave("Handsaker_CAGLENGTH_hist_phaseC.pdf", width = 8, height = 8)

p <- ggplot(SPN_filt@meta.data, aes(x = EXP_CAGLENGTH, fill = SAMPLE)) +
geom_histogram(binwidth = 5, alpha = 0.5, position = "identity") + 
xlim(c(150, 1000)) +
theme_classic()
#   theme(legend.title = element_text(size = 10),
#         legend.text = element_text(size = 10),
#         axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
#         axis.text.y = element_text(size = 10),
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         panel.background = element_blank())
plot(p)
ggsave("Handsaker_CAGLENGTH_hist_phaseC_splitBySample.pdf", width = 8, height = 8)
```


```R
Handsaker_donor_metadata <- read.table("/mnt/projects/labs/CLAB/PROJECT_ELongATE/Handsaker/donor_metadata.txt", header = TRUE)
Handsaker_donor_metadata <- Handsaker_donor_metadata[which(Handsaker_donor_metadata$SID %in% unique(SPN_filt$DONOR)), ]
rownames(Handsaker_donor_metadata) <- Handsaker_donor_metadata$SID
Handsaker_donor_metadata
```


```R
tmp_Handsaker <- do.call(rbind,lapply(split(SPN_filt@meta.data$PHASE, SPN_filt@meta.data$SAMPLE), function(x) {
    tab <- table(x)
    if (length(tab) == 1) {
        empty_class <- setdiff(c("A-B", "C-D-E"), names(tab))
        tab <- c(tab, 0)
        names(tab)[2] <- empty_class
    }
    return(tab)
}))

#tmp_Handsaker
Handsaker_fract_CDE <- data.frame(DATASET = "Handsaker_Caudate", SAMPLE = rownames(tmp_Handsaker), FRACT_CDE=100*tmp_Handsaker[, 2]/(tmp_Handsaker[, 1] + tmp_Handsaker[, 2]), NUM_SPN=tmp_Handsaker[, 1] + tmp_Handsaker[, 2])
Handsaker_fract_CDE$GRADE <- Handsaker_donor_metadata[rownames(Handsaker_fract_CDE), "VS_Grade"]
Handsaker_fract_CDE$FRACT_SPN <- num_cells_Handsaker[rownames(Handsaker_fract_CDE), "num_spns"]/num_cells_Handsaker[rownames(Handsaker_fract_CDE), "num_cells"]
Handsaker_fract_CDE$GRADE <- factor(Handsaker_fract_CDE$GRADE, levels = c("n/a", "HD-0", "HD-1", "HD-2", "HD-3"))
Handsaker_fract_CDE
```


```R
Handsaker_fract_CDE_HD <- Handsaker_fract_CDE[which(Handsaker_fract_CDE$SAMPLE %in% rownames(Handsaker_donor_metadata)[which(Handsaker_donor_metadata$Status == "Case")]), ]
#Handsaker_fract_CDE_CTRL <- Handsaker_fract_CDE[which(Handsaker_fract_CDE$SAMPLE %in% rownames(Handsaker_donor_metadata)[which(Handsaker_donor_metadata$Status == "Control")]), ]

ggplot(Handsaker_fract_CDE_HD, aes(x = GRADE, y = FRACT_CDE, size = FRACT_SPN)) +
geom_point(aes(size = FRACT_SPN), color = "#F8766D", alpha = 0.8, position = position_jitter(width = 0, height = 0)) +
#scale_color_manual(values =c("#F8766D")) +
theme_classic() +
xlab("HD Grade") +
ylim(c(0, 25)) +
ylab("Fraction cells in C-D-E phase (%)")
ggsave("Handsaker_HD_fraction_CDE.pdf", width = 6, height = 6)
```

# Read_phaseC_genes


```R
phaseC_file <- "/mnt/projects/labs/CLAB/PROJECT_ELongATE/PhaseCD.txt"

phaseC_genes_all <- read.table(phaseC_file, sep = "\t", header = TRUE)
phaseC_genes <- phaseC_genes_all$Gene
```


```R
length(phaseC_genes)
```


```R
phaseC_plus_genes <- phaseC_genes_all[which(phaseC_genes_all$Phase.C.effect > 0), "Gene"]
phaseC_minus_genes <- phaseC_genes_all[which(phaseC_genes_all$Phase.C.effect < 0), "Gene"]
```


```R
length(phaseC_plus_genes)
length(phaseC_minus_genes)
```


```R
#plot UMAP using only PhaseC genes as a test
SPN_filt <- RunPCA(SPN_filt, features = phaseC_genes, npcs = 50, reduction.name = "pca_phaseC")
SPN_filt <- RunUMAP(SPN_filt, dims = 1:num_dims, reduction = "pca_phaseC", reduction.name = "umap_phaseC", reduction.key = "umap_phaseC")
FeaturePlot(SPN_filt, reduction = "umap_phaseC", features = "CAGLENGTH", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Handsaker_phaseC_genesOnly_CAGLENGTH.pdf")
DimPlot(SPN_filt, reduction = "umap_phaseC", group.by = "PHASE", pt.size = 1)
ggsave("UMAP_Handsaker_phaseC_genesOnly_PHASE.pdf")
```

# Learn_phase_model


```R
#extract phase C genes that are expressed in the dataset
phaseC_genes <- intersect(phaseC_genes, rownames(SPN_filt@assays$SCT$data))
phaseC_plus_genes <- intersect(phaseC_plus_genes, phaseC_genes)
phaseC_minus_genes <- intersect(phaseC_minus_genes, phaseC_genes)

#extract SPNs with a phase
PHASE_SPN <- SPN_filt@meta.data[, "PHASE"]

#extract counts for Phase C genes
phaseC_counts_SPN_phased <- SPN_filt@assays$SCT$data[phaseC_genes, which(!is.na(PHASE_SPN))]
phaseC_counts_SPN <- SPN_filt@assays$SCT$data[phaseC_genes, ]

#create dataframe with PHASE and counts of phase C genes
training_test_data_phase <- data.frame(PHASE = PHASE_SPN, DONOR = SPN_filt@meta.data$DONOR, t(phaseC_counts_SPN))

#Phase C genes expression for all SPNs
phase_all <- training_test_data_phase[, -c(1, 2)]
```


```R
learn_phase_model = function(data, class, lambda, modelTrained = NULL, thr = NULL, titlePR = NULL, thr_stats = NULL) {
  if (length(modelTrained) == 0) {
    model <- cv.glmnet(x = data,
                       y = class,
                       family = "binomial",
                       type.measure = "class")
  } else {
    model <- modelTrained
  }
  #run phase model
  y_pred <- predict(model, 
                    newx = data,
                    s = lambda,
                    type = "response")

  #plot PR curve  
  pr <- pr.curve(scores.class0 = y_pred[class == "C-D-E"],
                 scores.class1 = y_pred[class == "A-B"],
                 curve = TRUE)
  plot(pr, main = titlePR)
  perf <- assess_metrics(predictions = y_pred, best_t = thr, y_true = class, thresholds = thr_stats)
  return(list(model, pr, y_pred, perf))
}
```


```R
assess_metrics <- function(predictions, best_t, y_true, thresholds = NULL) {
  #sort probabilities to assess precision, recall, F1 and accuracy for each threshold
  if (length(thresholds) == 0) {
    thresholds <- sort(unique(predictions), decreasing = TRUE)
  }
  tmp <- lapply(thresholds, function(t, pred, class) {
    #find predicted class using t as a threshold
    y_hat <- as.integer(pred[, 1] >= t)
    #counts TP, FP, FN and TN
    TP <- sum(y_hat == 1 & class == "C-D-E")
    FP <- sum(y_hat == 1 & class == "A-B")
    FN <- sum(y_hat == 0 & class == "C-D-E")
    TN <- sum(y_hat == 0 & class == "A-B")
    #assess metrics
    precision <- ifelse(TP + FP == 0, NA, TP / (TP + FP))
    recall    <- ifelse(TP + FN == 0, NA, TP / (TP + FN))
    F1        <- ifelse(is.na(precision + recall), NA,
                      2 * precision * recall / (precision + recall))
    acc       <- ifelse(TP + FP + TN + FN == 0, NA, (TP + TN) / (TP + FP + TN + FN))
    #save metrics to df
    metrics <- data.frame(thr = t, P = precision, R = recall, F1 = F1, acc = acc)
    return(metrics)
  }, pred = predictions, class = y_true)
  #convert list of df to df   
  metrics_all <- do.call("rbind", tmp)
  #choose which metric should be optimized
  ##optimize F1 score
  #best_thr <- thresholds[which.max(metrics[, "F1"])]
  #optimize recall with minimum target precision
  target_precision <- 0.99
  #if best_t was not already set on the training set, assess it
  if (length(best_t) == 0) {
    valid <- which(metrics_all[, "P"] >= target_precision)
    best_thr <- thresholds[valid[which.max(metrics_all[valid, "R"])]]      
  } else {
    valid <- which.min(abs(metrics_all[, "thr"] - best_t))
    best_thr <- thresholds[valid]
  }
  ##optimize precision with minimum target recall
  #target_recall <- 0.95
  #valid <- which(metrics[, "R"] >= target_recall)
  #best_thr <- thresholds[valid[which.max(metrics[valid, "P"])]]
  metrics_final <- metrics_all[which(metrics_all$thr == best_thr)[1], ]
  return(metrics_final)
}
```


```R
#set this flag to TRUE to perform training on n-1 samples and test on the remaining sample; otherwise, training and test are randomly sampled
leave_one_out_flag <- TRUE
#initialize variables
models_phase <- list()
PR_curves_phase_training_noNA <- list()
preds_phase_training_noNA <- list()
perf_phase_training_noNA <- list()
PR_curves_phase_test_noNA <- list()
preds_phase_test_noNA <- list()
perf_phase_test_noNA <- list()
res_phase_training_noNA <- list()
res_phase_test_noNA <- list()
res_phase_test <- list()

#cycle across donors
for (i in 1:length(unique(training_test_data_phase$DONOR))) { 
    set.seed(i)
    #training on n-1 samples and test on the remaining sample
    if (leave_one_out_flag) {
      ind_training_phase <- which(training_test_data_phase$DONOR != unique(training_test_data_phase$DONOR)[i])
      ind_test_phase <- which(training_test_data_phase$DONOR == unique(training_test_data_phase$DONOR)[i])
      it <- paste0("Iteration with held-out donor: ", unique(training_test_data_phase$DONOR)[i])
    #training on random 70% of cells and test on the remaining 30%
    } else {
      ind_training_phase <- sample(x = 1:round(length(PHASE_SPN)), size = round(length(PHASE_SPN)*0.7), replace = FALSE)
      ind_test_phase <- setdiff(1:round(length(PHASE_SPN)), ind_training_phase)
      it <- paste0("Iteration #", i)
    }
  #select data for training
  training_data_phase <- training_test_data_phase[ind_training_phase, ]
  test_data_phase <- training_test_data_phase[ind_test_phase, ]

  #select data for training excluding those with NA phase
  training_data_phase_noNA <- training_data_phase[which(!is.na(training_data_phase$PHASE)), ]
  test_data_phase_noNA <- test_data_phase[which(!is.na(test_data_phase$PHASE)), ]

  phase_training <- training_data_phase[, -c(1, 2)]
  phase_test <- test_data_phase[, -c(1, 2)]

  phase_training_noNA <- training_data_phase_noNA[, -c(1, 2)]
  phase_test_noNA <- test_data_phase_noNA[, -c(1, 2)]

  tmp_training <- learn_phase_model(data = as.matrix(phase_training_noNA),
                                    class = factor(training_data_phase_noNA[, 1], levels = c("A-B", "C-D-E")),
                                    lambda = "lambda.min",
                                    titlePR = paste0("Training set: ", it))
    
  models_phase[[i]] <- tmp_training[[1]]
  PR_curves_phase_training_noNA[[i]] <- tmp_training[[2]]
  preds_phase_training_noNA[[i]] <- tmp_training[[3]]
  perf_phase_training_noNA[[i]] <- tmp_training[[4]]
    
  predicted_phase_training_noNA <- rep("A-B", dim(preds_phase_training_noNA[[i]])[1])
  predicted_phase_training_noNA[which(preds_phase_training_noNA[[i]] > perf_phase_training_noNA[[i]]$thr)] <- "C-D-E"
  
  res_phase_training_noNA[[i]] <- data.frame(pred = preds_phase_training_noNA[[i]], predicted_phase_training_noNA, PHASE = training_data_phase_noNA[, 1])

  #acc_training <- (conf_mat_training[1, 1] + conf_mat_training[2, 2])/(sum(conf_mat_training))

  cat(sprintf("\n%s\n", it))
  cat(sprintf("\nOptimal threshold: %.3f\n", perf_phase_training_noNA[[i]]$thr))
  cat(sprintf("Phase model performances on training set\nAUC: %.3f\nThreshold: %.3f - Precision: %.3f - Recall: %.3f - F1 score: %.3f - Accuracy: %.3f\n", PR_curves_training_noNA[[i]]$auc.integral, perf_training_noNA[[i]]$thr, perf_training_noNA[[i]]$P, perf_training_noNA[[i]]$R, perf_training_noNA[[i]]$F1, perf_training_noNA[[i]]$acc))    
  cat(sprintf("Confusion matrix on training set\n"))
  conf_mat_training <- table(res_phase_training_noNA[[i]][, c(2, 3)])
  print(conf_mat_training)

  tmp_test <- learn_phase_model(data = as.matrix(phase_test_noNA),
                                class = factor(test_data_phase_noNA[, 1], levels = c("A-B", "C-D-E")),
                                lambda = "lambda.min",
                                modelTrained = models_phase[[i]],
                                thr = perf_training_noNA[[i]]$thr,
                                titlePR = paste0("Test set: ", it),
                                thr_stats = preds_phase_training_noNA[[i]]) 

  #models_phase[[i]] <- tmp_test[[1]]
  PR_curves_phase_test_noNA[[i]] <- tmp_test[[2]]
  preds_phase_test_noNA[[i]] <- tmp_test[[3]]
  perf_phase_test_noNA[[i]] <- tmp_test[[4]]
    
  predicted_phase_test_noNA <- rep("A-B", dim(preds_phase_test_noNA[[i]])[1])
  predicted_phase_test_noNA[which(preds_phase_test_noNA[[i]] > perf_phase_test_noNA[[i]]$thr)] <- "C-D-E"
    
  res_phase_test_noNA[[i]] <- data.frame(pred = preds_phase_test_noNA[[i]], predicted_phase_test_noNA, PHASE = test_data_phase_noNA[, 1])
  
  #run phase model on all test data (also those with NAs)
  #assign predicted PHASE to all SPNs
  preds_phase_test <- predict(object = models_phase[[i]], 
                        newx = as.matrix(test_data_phase[, -c(1, 2)]),
                        s = "lambda.min",
                        type = "response")
    
   predicted_phase_test <- rep("A-B", dim(preds_phase_test)[1])
   predicted_phase_test[which(preds_phase_test > perf_phase_test_noNA[[i]]$thr)] <- "C-D-E"
  
  res_phase_test[[i]] <- data.frame(pred = preds_phase_test, predicted_phase_test, PHASE = test_data_phase[, 1])
 
  #acc_test <- (conf_mat_test[1, 1] + conf_mat_test[2, 2])/(sum(conf_mat_test))
    
  cat(sprintf("Phase model performances on test set\nAUC: %.3f\nThreshold: %.3f - Precision: %.3f - Recall: %.3f - F1 score: %.3f - Accuracy: %.3f\n", PR_curves_test_noNA[[i]]$auc.integral, perf_test_noNA[[i]]$thr, perf_test_noNA[[i]]$P, perf_test_noNA[[i]]$R, perf_test_noNA[[i]]$F1, perf_test_noNA[[i]]$acc))
  cat(sprintf("Confusion matrix on test set\n"))
  conf_mat_test_noNA <- table(res_phase_test_noNA[[i]][, c(2, 3)])
  print(conf_mat_test_noNA)

  conf_mat_test <- table(res_phase_test[[i]][, c(2, 3)])
  print(conf_mat_test)
    
}
```


```R
saveRDS(object = models_phase, file = "Phase_models_leaveOneOut.Rds")
saveRDS(object = PR_curves_phase_training_noNA, file = "PR_curves_training.Rds")
#saveRDS(object = preds_phase_training_noNA, file = "Predictions_training.Rds")
saveRDS(object = perf_phase_training_noNA, file = "Performances_training.Rds")
saveRDS(object = PR_curves_phase_test_noNA, file = "PR_curves_training.Rds")
#saveRDS(object = preds_phase_test_noNA, file = "Predictions_test.Rds")
saveRDS(object = perf_phase_test_noNA, file = "Performances_test.Rds")
saveRDS(object = res_phase_training_noNA, file = "Results_training_noNA.Rds")
saveRDS(object = res_phase_test_noNA, file = "Results_test_noNA.Rds")
saveRDS(object = res_phase_test, file = "Results_test.Rds")
```


```R
models_phase <- readRDS("Phase_models_leaveOneOut.Rds")
PR_curves_phase_training_noNA <- readRDS("PR_curves_training.Rds")
perf_phase_training_noNA <- readRDS("Performances_training.Rds")
PR_curves_phase_test_noNA <- readRDS("PR_curves_training.Rds")
perf_phase_test_noNA <- readRDS("Performances_test.Rds")
res_phase_training_noNA <- readRDS("Results_training_noNA.Rds")
res_phase_test_noNA <- readRDS("Results_test_noNA.Rds")
res_phase_test <- readRDS("Results_test.Rds")
```


```R
names(res_phase_test) <- unique(training_test_data_phase$DONOR)
```


```R
names(res_phase_test_noNA) <- unique(training_test_data_phase$DONOR)
```


```R
lapply(res_phase_test, function(x) table(x$predicted_phase_test))
```


```R
#add predicted probability and predicted phase in leave-one-out test 
for (i in 1:length(unique(SPN_filt@meta.data$DONOR))) {
    ind_curr <- which(SPN_filt@meta.data$DONOR == unique(SPN_filt@meta.data$DONOR)[i])
    SPN_filt@meta.data[ind_curr, "PREDICTED_PHASE"] <- res_phase_test[[unique(SPN_filt@meta.data$DONOR)[i]]][rownames(SPN_filt@meta.data[ind_curr, ]), "predicted_phase_test"]
    SPN_filt@meta.data[ind_curr, "PROB_PHASE_CDE"] <- res_phase_test[[unique(SPN_filt@meta.data$DONOR)[i]]][rownames(SPN_filt@meta.data[ind_curr, ]), "lambda.min"]
}
```


```R
DimPlot(SPN_filt, reduction = "umap", group.by = "PREDICTED_PHASE", pt.size = 1)
ggsave("UMAP_Handsaker_PREDICTED_PHASE.pdf", width = 8, height = 8)
FeaturePlot(SPN_filt, reduction = "umap", features = "PROB_PHASE_CDE", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Handsaker_PROB_PHASE_CDE.pdf", width = 8, height = 8)
```


```R
#full dataset 
training_test_data_phase_noNA <- training_test_data_phase[which(!is.na(training_test_data_phase$PHASE)), ]
phase_training_test_data_noNA <- training_test_data_phase_noNA[, -c(1, 2)]

#learn model using full dataset
tmp_phase_training_test <- learn_phase_model(data = as.matrix(phase_training_test_data_noNA),
                                class = factor(training_test_data_phase_noNA[, 1], levels = c("A-B", "C-D-E")),
                                lambda = "lambda.min",
                                titlePR = paste0("Full training dataset")) 

#extract model and optimal threshold
model_phase <- tmp_phase_training_test[[1]]
PR_curves_phase_training_test_noNA <- tmp_phase_training_test[[2]]
preds_phase_training_test_noNA <- tmp_phase_training_test[[3]]
perf_phase_training_test_noNA <- tmp_phase_training_test[[4]]

#assign predicted PHASE to all SPNs
predicted_phase_training_test_noNA <- rep("A-B", dim(preds_phase_training_test_noNA)[1])
predicted_phase_training_test_noNA[which(preds_phase_training_test_noNA > perf_phase_training_test_noNA$thr)] <- "C-D-E"
  
res_phase_training_test <- data.frame(pred = preds_phase_training_test_noNA, predicted_phase_training_test_noNA, PHASE = training_test_data_phase_noNA[, 1])
```


```R
saveRDS(object = model_phase, file = "Phase_model_fullDataset.Rds")
saveRDS(object = PR_curves_phase_training_test_noNA, file = "PR_curves_training_test.Rds")
#saveRDS(object = preds_phase_training_test_noNA, file = "Predictions_training_test.Rds")
saveRDS(object = perf_phase_training_test_noNA, file = "Performances_training_test.Rds")
saveRDS(object = res_phase_training_test, file = "Results_training_test.Rds")
```


```R
length(predicted_phase_test_noNA)
dim(phase_test)
dim(test_data_phase)
```


```R
names(PR_curves_test_noNA) <- unique(training_test_data_phase$DONOR)
names(PR_curves_test_noNA)
```


```R
model_phase <- readRDS("Phase_model_fullDataset.Rds")
perf_phase_training_test_noNA <-  readRDS("Performances_training_test.Rds")
perf_phase_training_test_noNA
```


```R
class_thr <- perf_training_test_noNA$thr
lambda_val <- "lambda.min"

#run phase model on all SPNs
predicted_phase_all_wprob <- predict(model_phase, newx = as.matrix(phase_all),
                                     s = lambda_val,
                                     type = "response")
predicted_phase_all <- rep("A-B", dim(predicted_phase_all_wprob)[1])  
predicted_phase_all[which(predicted_phase_all_wprob > class_thr)] <- "C-D-E"
```


```R
head(phase_all)
```


```R
table(predicted_phase_training_test_noNA)
table(predicted_phase_all)
```


```R
#extract genes selected by the phase model
tmp_coef_phase <- coef(model_phase, s = lambda_val)
selected_features_phase <- names(tmp_coef_phase[tmp_coef_phase@i + 1, ])[-1]
length(selected_features_phase)
```


```R
#head(res_training)
length(which(res_training_test$predicted_phase == "A-B"))
length(which(res_training_test$predicted_phase == "C-D-E"))
length(which(res_training_test$PHASE == "A-B"))
length(which(res_training_test$PHASE == "C-D-E"))
```


```R
#add predicted phase to HD samples and plot UMAP
#SPN_filt@meta.data[, "PREDICTED_PHASE"] <- predicted_phase_all
```


```R
SPN_filt@meta.data$CONDITION_BIN <- SPN_filt@meta.data$CONDITION
SPN_filt@meta.data$CONDITION_BIN[which(SPN_filt@meta.data$CONDITION_BIN == "HD_PREM")] <- "HD"
```


```R
table(SPN_filt@meta.data$PHASE, SPN_filt@meta.data$PREDICTED_PHASE)
```


```R
(8580+214)/(8580+214+19+560)
```


```R
FeaturePlot(SPN_filt, reduction = "umap", features = "PROB_PHASE_CDE", cols = c("blue", "red"), split.by = "CONDITION_BIN", pt.size = 1)
ggsave("UMAP_Handsaker_PROB_PHASE_CDE.pdf", width = 8, height = 8)
```


```R
#add metadata
#assess fraction of SPNs in CDE Phase
tmp <- lapply(split(SPN_filt$PHASE, SPN_filt$DONOR), function(x) {
    if (length(unique(x)) > 1) {
        val <- table(x)[2]/(table(x)[1] + table(x)[2])
    } else {
        val <- 0
    }
    names(val) <- "FRACT_CDE"; return(val)})
tmp2 <- do.call(rbind, tmp)

#assess fraction of SPNs predicted to be in CDE Phase in leave-one-out predictions
names(res_test) <- unique(training_test_data_phase$DONOR)
tmp3 <- lapply(res_test, function(x) {
     if (length(unique(x$predicted_phase_test_noNA)) > 1) {
        val <- table(x$predicted_phase_test_noNA)[2]/(table(x$predicted_phase_test_noNA)[1] + table(x$predicted_phase_test_noNA)[2])
    } else {
        val <- 0
    }
    names(val) <- "PRED_FRACT_CDE"; return(val)})
tmp4 <- do.call("rbind", tmp3)

tmp5 <- lapply(split(SPN_filt$PHASE, SPN_filt$DONOR), function(x) {
    val = length(x)
    names(val) <- "NUM_SPN"; return(val)})
tmp6 <- do.call(rbind, tmp5)

GERM_CAGEXP <- as.numeric(gsub(x = Handsaker_donor_metadata$CAG, pattern = ".*/", replacement = ""))
names(GERM_CAGEXP) <- rownames(Handsaker_donor_metadata)

Handsaker_donor_metadata_full <- cbind(FRACT_CDE = tmp2, PRED_FRACT_CDE = tmp4[rownames(tmp2), ], NUM_SPN = tmp6[rownames(tmp2), ], GERM_CAGEXP = GERM_CAGEXP[rownames(tmp2)], Handsaker_donor_metadata[rownames(tmp2), ])

Handsaker_donor_metadata_full$CAP <- as.numeric(Handsaker_donor_metadata_full$CAP)
Handsaker_donor_metadata_full$PRED_FRACT_CDE[which(is.na(Handsaker_donor_metadata_full$PRED_FRACT_CDE))] <- 0
Handsaker_donor_metadata_full$FRACT_CDE[which(is.na(Handsaker_donor_metadata_full$FRACT_CDE))] <- 0
Handsaker_donor_metadata_full
```


```R
str(Handsaker_donor_metadata_full)
```


```R
wtd.cor(x = Handsaker_donor_metadata_HD$GERM_CAGEXP, y = Handsaker_donor_metadata_HD$FRACT_CDE, weight = Handsaker_donor_metadata_HD$NUM_SPN)
wtd.cor(x = Handsaker_donor_metadata_HD$CAP, y = Handsaker_donor_metadata_HD$FRACT_CDE, weight = Handsaker_donor_metadata_HD$NUM_SPN)
wtd.cor(x = Handsaker_donor_metadata_HD$FRACT_CDE, y = Handsaker_donor_metadata_HD$PRED_FRACT_CDE, weight = Handsaker_donor_metadata_HD$NUM_SPN)
```


```R
Calculate_confInt = function(r, w) {
  #n is the number of samples
  n <- length(w)
    
  #Fisher transformation
  z <- 0.5 * log((1 + r) / (1 - r))
    
  #assess effective n
  n_eff <- (sum(w)^2) / sum(w^2)

  #Standard error
  #SE_z <- 1 / sqrt(n_eff - 3)
  SE_z <- 1 / sqrt(n - 3)

  #Z-scale CI
  z_lower <- z - 1.96 * SE_z
  z_upper <- z + 1.96 * SE_z

  #r-scale CI
  r_lower <- (exp(2 * z_lower) - 1) / (exp(2 * z_lower) + 1)
  r_upper <- (exp(2 * z_upper) - 1) / (exp(2 * z_upper) + 1)

  c(r_lower, r, r_upper) 
}
```


```R
Calculate_confInt(r = wtd.cor(x = Handsaker_donor_metadata_HD$GERM_CAGEXP, y = Handsaker_donor_metadata_HD$FRACT_CDE, weight = Handsaker_donor_metadata_HD$NUM_SPN)[, "correlation"], w = Handsaker_donor_metadata_HD$NUM_SPN)
Calculate_confInt(r = wtd.cor(x = Handsaker_donor_metadata_HD$CAP, y = Handsaker_donor_metadata_HD$FRACT_CDE, weight = Handsaker_donor_metadata_HD$NUM_SPN)[, "correlation"], w = Handsaker_donor_metadata_HD$NUM_SPN)
Calculate_confInt(r = wtd.cor(x = Handsaker_donor_metadata_HD$FRACT_CDE, y = Handsaker_donor_metadata_HD$PRED_FRACT_CDE, weight = Handsaker_donor_metadata_HD$NUM_SPN)[, "correlation"], w = Handsaker_donor_metadata_HD$NUM_SPN)
```


```R
wtd.cor(x = Handsaker_donor_metadata_full$GERM_CAGEXP, y = Handsaker_donor_metadata_full$FRACT_CDE, weight = Handsaker_donor_metadata_full$NUM_SPN)
wtd.cor(x = Handsaker_donor_metadata_full$CAP, y = Handsaker_donor_metadata_full$FRACT_CDE, weight = Handsaker_donor_metadata_full$NUM_SPN)
wtd.cor(x = Handsaker_donor_metadata_full$FRACT_CDE, y = Handsaker_donor_metadata_full$PRED_FRACT_CDE, weight = Handsaker_donor_metadata_full$NUM_SPN)
```


```R
Calculate_confInt(r = wtd.cor(x = Handsaker_donor_metadata_full$GERM_CAGEXP, y = Handsaker_donor_metadata_full$FRACT_CDE, weight = Handsaker_donor_metadata_full$NUM_SPN)[, "correlation"], w = Handsaker_donor_metadata_full$NUM_SPN)
Calculate_confInt(r = wtd.cor(x = Handsaker_donor_metadata_full$CAP, y = Handsaker_donor_metadata_full$FRACT_CDE, weight = Handsaker_donor_metadata_full$NUM_SPN)[, "correlation"], w = Handsaker_donor_metadata_full$NUM_SPN)
Calculate_confInt(r = wtd.cor(x = Handsaker_donor_metadata_full$FRACT_CDE, y = Handsaker_donor_metadata_full$PRED_FRACT_CDE, weight = Handsaker_donor_metadata_full$NUM_SPN)[, "correlation"], w = Handsaker_donor_metadata_full$NUM_SPN)
```


```R
#plot fraction of SPN in CDE phase in the test set vs CAP score
cor_CAP_fractCDE <- weighted.cor(x = Handsaker_donor_metadata_full$FRACT_CDE, y = Handsaker_donor_metadata_full$CAP, weights = Handsaker_donor_metadata_full$NUM_SPN, method = "pearson", na.rm=TRUE)
ggplot(Handsaker_donor_metadata_full, aes(x=CAP, y=100*FRACT_CDE, size = NUM_SPN)) + 
  geom_point(color = "blue") +
  geom_smooth(method=lm, se=TRUE, mapping = aes(weight = NUM_SPN)) +
  theme_classic() +
  xlab("CAP score") +
  ylab("Fraction cells in C-D-E phase (%)")
sprintf("Pearson corr = %.2f", cor_CAP_fractCDE)
ggsave("Handsaker_fraction_SPN_CDE_phase_CAP.pdf", width = 8, height = 8)
```


```R
#plot fraction of SPN in CDE phase in the test set vs num. CAG in the germline
cor_numCAGGerm_fractCDE <- weighted.cor(x = Handsaker_donor_metadata_full$GERM_CAGEXP, y = Handsaker_donor_metadata_full$FRACT_CDE, weights = Handsaker_donor_metadata_full$NUM_SPN, method = "pearson", na.rm=TRUE)
cor_numCAGGerm_fractCDE
#ggplot(Handsaker_donor_metadata_full, aes(x=GERM_CAGEXP, y=100*FRACT_CDE, size = NUM_SPN)) + 
ggplot(Handsaker_donor_metadata_HD, aes(x=GERM_CAGEXP, y=100*FRACT_CDE, size = NUM_SPN)) + 
  geom_point(color = "blue") +
  geom_smooth(method=lm, se=TRUE, mapping = aes(weight = NUM_SPN)) +
  theme_classic() +
   xlab("Num. CAG germline") +
  ylab("Fraction cells in C-D-E phase (%)")

print(sprintf("Pearson corr = %.2f", cor_numCAGGerm_fractCDE))
ggsave("Handsaker_numCAGGermline_vs_fraction_SPN_CDE.pdf", width = 8, height = 8)
```


```R
#plot fraction of SPN in CDE phase in the test set
cor_predFractCDE_fractCDE <- weighted.cor(x = Handsaker_donor_metadata_full$FRACT_CDE, y = Handsaker_donor_metadata_full$PRED_FRACT_CDE, weights = Handsaker_donor_metadata_full$NUM_SPN, method = "pearson", na.rm=TRUE)
#ggplot(Handsaker_donor_metadata_full, aes(x=100*FRACT_CDE, y=100*PRED_FRACT_CDE, size = NUM_SPN)) + 
ggplot(Handsaker_donor_metadata_HD, aes(x=100*FRACT_CDE, y=100*PRED_FRACT_CDE, size = NUM_SPN)) + 
  geom_point(color = "blue") +
  geom_smooth(method=lm, se=TRUE, mapping = aes(weight = NUM_SPN)) +
  theme_classic() +
  xlab("Fraction cells in C-D-E phase (%)") +
  ylab("Predicted fraction cells in C-D-E phase (%)")

sprintf("Pearson corr = %.2f", cor_predFractCDE_fractCDE)
ggsave("Handsaker_fraction_SPN_CDE_phase_test_set.pdf", width = 8, height = 8)
```

# Define training and test set for CAG sizing model


```R
CAGLENGTH_SPN_sized <- SPN_filt@meta.data[, "CAGLENGTH"]
training_test_data <- data.frame(CAGLENGTH_SPN_sized, DONOR = SPN_filt@meta.data[, "DONOR"], t(phaseC_counts_SPN))
#extract filtered training set with CAG > CAG_thr
expr_C_training_test <- training_test_data[, -c(1, 2)]
CAG_THR <- 150
#select all data
training_test_data_filtered <- training_test_data[which(training_test_data$CAGLENGTH_SPN_sized > CAG_THR), ]
expr_C_training_test_filtered <- training_test_data_filtered[, -c(1, 2)]
```


```R
length(CAGLENGTH_SPN_sized)
dim(training_test_data)
```


```R
summary(CAGLENGTH_SPN_sized)
length(sort(CAGLENGTH_SPN_sized))
```

# Learn_CAG_sizing_model


```R
#plot correlation between CAG sizing measured values and predictions
Evaluate_correlation <- function(Data, x, y, quant_thr = 1, abl = TRUE, notes = "") {
  #correlation <- cor(Data[, x], Data[, y], method = "pearson", use = "complete.obs")
  correlation <- cor(Data[, x], Data[, y], method = "spearman", use = "complete.obs")
  
  pdf(paste0("AccuracyPred_", notes, ".pdf"))
  if (abl) {
    smoothScatter(Data[, x], Data[, y], xlab = x, ylab = y, cex.main = 0.6, cex.axis = 0.8,
                xlim = c(0, quantile(c(Data[, x], Data[, y]), quant_thr, na.rm = TRUE)), 
                ylim = c(0, quantile(c(Data[, x], Data[, y]), quant_thr, na.rm = TRUE)))
    abline(0, 1, col = "green")
  } else {
      smoothScatter(Data[, x], Data[, y], xlab = x, ylab = y, cex.main = 0.6, cex.axis = 0.8,
                xlim = c(0, quantile(c(Data[, x], Data[, y]), quant_thr, na.rm = TRUE)), 
                ylim = c(0, quantile(c(Data[, x], Data[, y]), quant_thr, na.rm = TRUE)))
  }
  dev.off()
  
  #smoothScatter(Data[, x], Data[, y], xlab = x, main = paste0(y, " VS ", x, "\n r Pearson. = ",
  smoothScatter(Data[, x], Data[, y], xlab = x, main = paste0(y, " VS ", x, "\n r Spearman = ",
                sprintf("%.2f", correlation), "\n ", notes, " - n = ", dim(Data)[1]), ylab = y, cex.main = 0.6, cex.axis = 0.8,
                xlim = c(0, quantile(c(Data[, x], Data[, y]), quant_thr, na.rm = TRUE)), ylim = c(0, quantile(c(Data[, x], Data[, y]), quant_thr, na.rm = TRUE)))
 
  if (abl) {
    abline(0, 1, col = "green")
  }
    
  #print(paste0(y, " VS ", x, "\n r Pearson. = ", 
  print(paste0(y, " VS ", x, "\n r Spearman = ", 
              sprintf("%.2f", correlation), "\n ", notes, " - n = ", dim(Data)[1]))
  
  ind_noNA <- which(apply(Data, 1, function(x) !any(is.na(x))))
  Data_filt <- Data[ind_noNA, ]
  df <- data.frame(CAG = c(Data_filt[, x], Data_filt[, y]), method = c(rep("Measured", dim(Data_filt)[1]), rep("Predicted", dim(Data_filt)[1])))
  
  p <- ggplot(df, aes(x= CAG, fill = method, color = method)) +
  geom_density(alpha = 0.5, adjust = 1) +
  #geom_histogram(binwidth = 5, alpha = 0.5, position = "identity") +
  scale_color_manual(values = c("purple", "green")) +
  scale_fill_manual(values = c("purple", "green")) +
  xlab("Num. CAG") + ylab("Num. SPNs") +
  theme_classic()
#   theme(legend.title = element_text(size = 10),
#         legend.text = element_text(size = 10),
#         axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
#         axis.text.y = element_text(size = 10),
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         panel.background = element_blank())

print(gsub(x = notes, pattern = "_", replacement = " "))
  #  +  geom_vline(aes(xintercept = 150), colour="black", linetype = 2)
  
  plot(p)
  
  ggsave(paste0("AccuracyPred_hist_", notes, ".pdf"), width = 8, height = 8)
  #return(p)

  #residuals
  residuals <- Data[, x] - Data[, y]
  
  p1 <- plot(Data[, y], residuals,
     xlab = "Predicted",
     ylab = "Residuals",
     main = "Residuals vs Predictions")
  abline(h = 0, col = "red")

  print(p1)
  pdf(paste0("Residuals_vs_Predictions_CAGsizingModel_", notes, ".pdf"), width = 8, height = 8)
  print(p1)
  dev.off()
         
  hist(residuals, breaks = 30, main = paste0("Residuals distribution ", notes, " Mean residuals = ", sprintf("%.3f", mean(residuals, na.rm = TRUE)))) 
  pdf(paste0("Residuals_hist_CAGsizingModel_", notes, ".pdf"), width = 8, height = 8)
  hist(residuals, breaks = 30, main = paste0("Residuals distribution ", notes, " Mean residuals = ", sprintf("%.3f", mean(residuals, na.rm = TRUE))))
  dev.off()
                          
  RMSE <- sqrt(mean(residuals^2))
  MAE  <- mean(abs(residuals))
  R2   <- 1 - sum(residuals^2) / sum((Data[, y] - mean(Data[, y]))^2)
  corr_test <- cor.test(Data[, x], Data[, y], method = "spearman", use = "complete.obs")
  return(list(RMSE, MAE, R2, corr_test))
}
```


```R
#set this flag to TRUE to perform training on n-1 samples and test on the remaining sample; otherwise, training and test are randomly sampled
leave_one_out_flag <- TRUE
#initialize variables
models_CAGsizing <- list()
res_CAGsizing_training_filtered <- list()
res_CAGsizing_training_expPred <- list()
perf_CAGsizing_training_filtered <- list()
perf_CAGsizing_training_expPred <- list()
res_CAGsizing_test <- list()
res_CAGsizing_test_expPred <- list()
res_CAGsizing_test_expMeas <- list()
perf_CAGsizing_test_expPred <- list()
perf_CAGsizing_test_expMeas <- list()

#cycle across donors
for (i in 1:length(unique(training_test_data$DONOR))) {
    set.seed(i)
    #training on n-1 samples and test on the remaining sample
    if (leave_one_out_flag) {
      ind_training <- which(training_test_data$DONOR != unique(training_test_data$DONOR)[i])
      ind_test <- which(training_test_data$DONOR == unique(training_test_data$DONOR)[i])
      it <- paste0("Held_out_donor_", unique(training_test_data$DONOR)[i])
    #training on random 70% of cells and test on the remaining 30%
    } else {
      ind_training <- sample(x = 1:round(length(SPN_filt@assays$SCT$data)), size = round(length(SPN_filt@assays$SCT$data)*0.7), replace = FALSE)
      ind_test <- setdiff(1:round(length(SPN_filt@assays$SCT$data)), ind_training)
      it <- paste0("Iteration_", i)
    }
  #select training data  
  training_data <-  training_test_data[ind_training, ]
  expr_C_training <- training_data[, -c(1, 2)] 
  training_data_filtered <- training_data[which(training_data$CAGLENGTH_SPN_sized >= CAG_THR), ]
  expr_C_training_filtered <- training_data_filtered[, -c(1, 2)] 
  #select test data
  test_data <- training_test_data[ind_test, ]
  expr_C_test <- test_data[, -c(1, 2)] 
  test_data_filtered <- test_data[which(test_data$CAGLENGTH_SPN_sized >= CAG_THR), ]
  expr_C_test_filtered <- test_data_filtered[, -c(1, 2)] 
    
  #learn CAG sizing model
  fit_glm <- cv.glmnet(x = as.matrix(expr_C_training_filtered),
                       y = as.vector(training_data_filtered[, 1]),
                       family = "gaussian")
  
  #plot(fit_glm)
  lambda_val = "lambda.min"
  
  models_CAGsizing[[i]] <- fit_glm
  names(models_CAGsizing)[i] <- it
  tmp_coef <- coef(fit_glm, s = lambda_val)
  selected_features <- names(tmp_coef[tmp_coef@i + 1, ])[-1]
    
  #find number of phase C genes expressed in each SPN
  num_expr_genes_training <- apply(as.matrix(expr_C_training[, selected_features]), 1, function(x) length(which(x > 0)))
  num_expr_genes_training_filtered <- apply(as.matrix(expr_C_training_filtered[, selected_features]), 1, function(x) length(which(x > 0)))
  num_expr_genes_test <- apply(as.matrix(expr_C_test[, selected_features]), 1, function(x) length(which(x > 0)))

  #predict CAG size on filtered training set
  pred_training_filtered <- predict(fit_glm,
                                    newx = as.matrix(expr_C_training_filtered),
                                    s = lambda_val,
                                    type = "response")
  #predict CAG size on training set
  pred_training <- predict(fit_glm, newx = as.matrix(expr_C_training), s = lambda_val)
  #predict CAG size on test set
  pred_test <- predict(fit_glm,
                       newx = as.matrix(expr_C_test),
                       s = lambda_val,
                       type = "response")
  
  #create dataframes with model predictions
  #model predictions for SPNs in the training set with > 150 measured CAG repeats
  results_training_filtered <- data.frame(CAG_measured = training_data_filtered[, 1],
                                          CAG_predicted = as.numeric(pred_training_filtered[, 1]),
                                          phase = SPN_filt@meta.data[intersect(ind_training, which(SPN_filt@meta.data[, "PHASE"] == "C-D-E")), "PHASE"],
                                          pred_phase = SPN_filt@meta.data[intersect(ind_training, which(SPN_filt@meta.data[, "PHASE"] == "C-D-E")), "PREDICTED_PHASE"],
                                          num_expr_genes = num_expr_genes_training_filtered)
  
  res_training_filtered[[i]] <- results_training_filtered
  
  #model predictions for SPNs in the training set with > 150 predicted CAG repeats 
  results_training <- data.frame(CAG_measured = training_data[, 1],
                                 CAG_predicted = as.numeric(pred_training),
                                 phase = SPN_filt@meta.data[ind_training, "PHASE"],
                                 pred_phase = SPN_filt@meta.data[ind_training, "PREDICTED_PHASE"],
                                 num_expr_genes = num_expr_genes_training)
  results_training_expPred <- results_training[which(results_training$pred_phase == "C-D-E"), ]
  results_training_expPred[which(results_training_expPred$CAG_measured < 36), "CAG_measured"] <- NA   
  res_CAGsizing_training_expPred[[i]] <- results_training_expPred
  
  #model predictions for SPNs in the test set with > 150 predicted CAG repeats
  results_test <- data.frame(CAG_measured = test_data[, 1],
                             CAG_predicted = as.numeric(pred_test),
                             phase = SPN_filt@meta.data[ind_test, "PHASE"],
                             pred_phase = SPN_filt@meta.data[ind_test, "PREDICTED_PHASE"],
                             num_expr_genes = num_expr_genes_test)
  
  res_CAGsizing_test[[i]] <- results_test
  results_test_expPred <- results_test[which(results_test$pred_phase == "C-D-E"), ]
  results_test_expPred[which(results_test_expPred$CAG_measured < 36), "CAG_measured"] <- NA                      
  res_CAGsizing_test_expPred[[i]] <- results_test_expPred
                               
  #model predictions for SPNs in the test set with > 150 measured CAG repeats                             
  results_test_expMeas <- results_test[which(results_test$phase == "C-D-E"), ]
  results_test_expMeas[which(results_test_expMeas$CAG_measured < 36), "CAG_measured"] <- NA
  res_CAGsizing_test_expMeas[[i]] <- results_test_expMeas
  
  #assess correlation and model performances on training set with SPNs measured to be in C-D-E phase
  perf_CAGsizing_training_filtered[[i]] <- Evaluate_correlation(results_training_filtered, "CAG_measured", "CAG_predicted", notes = paste0("Training_set_Filtered_SPN_", it), quant_thr = 0.999)
  #assess correlation and model performances on training set with SPNs predicted to be in C-D-E phase
  perf_CAGsizing_training_expPred[[i]] <- Evaluate_correlation(results_training_expPred, "CAG_measured", "CAG_predicted", notes = paste0("Training_set_SPN_expPred_", it), quant_thr = 0.999)
  
  #assess correlation and model performances on training set with SPNs measured to be in C-D-E phase
  if (length(which(results_test_expMeas$phase == "C-D-E")) > 1) { 
    perf_CAGsizing_test_expMeas[[i]] <- Evaluate_correlation(results_test_expMeas, "CAG_measured", "CAG_predicted", notes = paste0("Test_set_SPN_expMeas_", it), quant_thr = 0.999)
  } else {
    cat(sprintf("Not assessing performances on test set for SPNs measured in C-D-E phase\n"))
    perf_CAGsizing_test_expMeas[[i]] <- NA
  }
  #assess correlation and model performances on test set with SPNs predicted to be in C-D-E phase
  if (length(which(results_test_expPred$pred_phase == "C-D-E")) > 1 && length(which(results_test_expPred$phase == "C-D-E")) > 1) {
    perf_CAGsizing_test_expPred[[i]] <- Evaluate_correlation(results_test_expPred, "CAG_measured", "CAG_predicted", notes = paste0("Test_set_SPN_expPred_", it), quant_thr = 0.999)
  } else {
    cat(sprintf("Not assessing performances on test set for SPNs predicted in C-D-E phase\n"))
    perf_CAGsizing_test_expPred[[i]] <- NA
  }
}
```


```R
names(res_CAGsizing_test_expMeas) <- unique(training_test_data$DONOR)
```


```R
#save variables to Rds files
saveRDS(object = models_CAGsizing, file = "CAGsizing_models_leaveOneOut.Rds")
saveRDS(object = res_CAGsizing_training_filtered, file = "CAGsizing_results_training_filtered.Rds")
saveRDS(object = res_CAGsizing_training_expPred, file = "CAGsizing_results_training_expPred.Rds")
saveRDS(object = res_CAGsizing_test, file = "CAGsizing_results_test.Rds")
saveRDS(object = res_CAGsizing_test_expPred, file = "CAGsizing_results_test_expPred.Rds")
saveRDS(object = res_CAGsizing_test_expMeas, file = "CAGsizing_results_test_expMeas.Rds")
names(perf_CAGsizing_training_filtered) <- names(models_CAGsizing)
saveRDS(object = perf_CAGsizing_training_filtered, file = "CAGsizing_performances_training_filtered.Rds")
names(perf_CAGsizing_training_expPred) <- names(models_CAGsizing)
saveRDS(object = perf_CAGsizing_training_expPred, file = "CAGsizing_performances_training_expPred.Rds")
names(perf_CAGsizing_test_expPred) <- names(models_CAGsizing)
saveRDS(object = perf_CAGsizing_test_expPred, file = "CAGsizing_performances_test_expPred.Rds")
names(perf_CAGsizing_test_expMeas) <- names(models_CAGsizing)
saveRDS(object = perf_CAGsizing_test_expMeas, file = "CAGsizing_performances_test_expMeas.Rds")
```


```R
names(models_CAGsizing)
```


```R
names(res_CAGsizing_test) <- unique(SPN_filt@meta.data$DONOR)
```


```R
table(training_test_data_filtered$DONOR)
```


```R
#add predicted probability and predicted phase in leave-one-out test 
for (i in 1:length(unique(SPN_filt@meta.data$DONOR))) {
    ind_curr <- which(SPN_filt@meta.data$DONOR == unique(SPN_filt@meta.data$DONOR)[i])
    SPN_filt@meta.data[ind_curr, "PRED_CAGLENGTH"] <- res_CAGsizing_test[[unique(SPN_filt@meta.data$DONOR)[i]]][rownames(SPN_filt@meta.data[ind_curr, ]), "CAG_predicted"]
}
```


```R
#select all data
training_test_data_filtered <- training_test_data[which(training_test_data$CAGLENGTH_SPN_sized >= CAG_THR), ]
expr_C_training_test_filtered <- training_test_data_filtered[, -c(1, 2)]

#learn CAG sizing model
fit_glm <- cv.glmnet(as.matrix(expr_C_training_test_filtered),
                     as.vector(training_test_data_filtered[, 1]))
#extract features
tmp_coef <- coef(fit_glm, s = lambda_val)
selected_features <- names(tmp_coef[tmp_coef@i + 1, ])[-1]
```


```R
#run CAG sizing model on all SPNs
pred_all <- predict(fit_glm,
                    newx = as.matrix(expr_C_training_test),
                    s = lambda_val,
                    type = "response")
```


```R
head(pred_all)
```


```R
head(predicted_phase_all_wprob)
```


```R
#saveRDS(fit_glm, "CAG_sizing_model.rds")
```


```R
fit_glm <- readRDS("CAG_sizing_model.rds")
```


```R
#dim(training_data[, -1])
dim(expr_C_training)
dim(expr_C_training_filtered)
```


```R
selected_features
length(selected_features)
```


```R
#training_data_filtered_sorted
coef(fit_glm, s = lambda_val)[, 1][selected_features]
```


```R
#sort SPNs of filtered training set by CAG length
training_test_data_filtered_sorted <- training_test_data_filtered[order(training_test_data_filtered$CAGLENGTH_SPN_sized), ]
head(training_test_data_filtered_sorted)
#plot(training_data_filtered_sorted[, 1], apply(training_data_filtered_sorted[, intersect(colnames(training_data_filtered_sorted), phaseC_plus_genes)], 1, mean))
#plot(training_data_filtered_sorted[, 1], apply(training_data_filtered_sorted[, intersect(colnames(training_data_filtered_sorted), phaseC_minus_genes)], 1, mean))

#plot average expression of all phase C+ and phase C- genes
plot(training_test_data_filtered_sorted[, 1], apply(training_test_data_filtered_sorted[, intersect(intersect(colnames(training_test_data_filtered_sorted), phaseC_plus_genes), selected_features)], 1, mean))
plot(training_test_data_filtered_sorted[, 1], apply(training_test_data_filtered_sorted[, intersect(intersect(colnames(training_test_data_filtered_sorted), phaseC_minus_genes), selected_features)], 1, mean))
```


```R
#combine the predictions of the two models
df <- data.frame(CAGLENGTH_SPN_sized=training_test_data[, "CAGLENGTH_SPN_sized"], score = pred_all[, 1], PROB_PHASE_CDE = predicted_phase_all_wprob[, 1], PRED_PHASE = predicted_phase_all)
df <- df[which(df$CAGLENGTH_SPN_sized > 36), ]
df_sorted <- df[order(df$CAGLENGTH_SPN_sized), ]
df_sorted$index <- seq_along(df_sorted$score)
```


```R
#extract features selected by CAG sizing model and associated coefficients
selected_features_plus <-  intersect(intersect(colnames(training_data_filtered_sorted), phaseC_plus_genes), selected_features)
selected_features_minus <-  intersect(intersect(colnames(training_data_filtered_sorted), phaseC_minus_genes), selected_features)
coef_selected_features_plus <- coef(fit_glm, s = lambda_val)[, 1][selected_features_plus]
coef_selected_features_minus <- coef(fit_glm, s = lambda_val)[, 1][selected_features_minus]
df_plus_training_filtered <- data.frame(CAGLENGTH_SPN_sized=training_data_filtered_sorted[, "CAGLENGTH_SPN_sized"], score = apply(mapply(function(x, y) x*y, training_data_filtered_sorted[, selected_features_plus], coef_selected_features_plus), 1, function(x) sum(x)) + coef(fit_glm, s = lambda_val)[1, 1]/2, signature = "plus")
df_minus_training_filtered <- data.frame(CAGLENGTH_SPN_sized=training_data_filtered_sorted[, "CAGLENGTH_SPN_sized"], score = apply(mapply(function(x, y) x*y, training_data_filtered_sorted[, selected_features_minus], coef_selected_features_minus), 1, function(x) sum(x) + coef(fit_glm, s = lambda_val)[1, 1]/2), signature = "minus")

#extract features selected by CAG sizing model and obtain score for SPNs
training_test_data_w_predPhase <- cbind(training_test_data, PROB_PHASE_CDE=predicted_phase_all)
training_test_data_w_predPhase <- training_test_data_w_predPhase[which(training_test_data_w_predPhase$CAGLENGTH_SPN_sized > 36), ]
df_plus <- data.frame(CAGLENGTH_SPN_sized=training_test_data_w_predPhase[, "CAGLENGTH_SPN_sized"], score = apply(mapply(function(x, y) x*y, training_test_data_w_predPhase[, selected_features_plus], coef_selected_features_plus), 1, function(x) sum(x)) + coef(fit_glm, s = lambda_val)[1, 1]/2, signature = "plus")
df_minus <- data.frame(CAGLENGTH_SPN_sized=training_test_data_w_predPhase[, "CAGLENGTH_SPN_sized"], score = apply(mapply(function(x, y) x*y, training_test_data_w_predPhase[, selected_features_minus], coef_selected_features_minus), 1, function(x) sum(x)) + coef(fit_glm, s = lambda_val)[1, 1]/2, signature = "minus")
#extract all Phase C+ and Phase C- genes and obtain C+ and C- scores for SPNs from test set
df_plus_allPhaseCgenes <- data.frame(CAGLENGTH_SPN_sized=training_test_data_w_predPhase[, "CAGLENGTH_SPN_sized"], score = apply(training_test_data_w_predPhase[, gsub(x = phaseC_plus_genes, pattern = "-", replacement = ".")], 1, function(x) mean(x)), signature = "plus")
df_minus_allPhaseCgenes <- data.frame(CAGLENGTH_SPN_sized=training_test_data_w_predPhase[, "CAGLENGTH_SPN_sized"], score = apply(training_test_data_w_predPhase[, gsub(x = phaseC_minus_genes, pattern = "-", replacement = ".")], 1, function(x) mean(x)), signature = "minus")                                                                                                                  
```


```R
df_allPhaseCgenes <- data.frame(CAGLENGTH_SPN_sized=training_test_data_w_predPhase[, "CAGLENGTH_SPN_sized"], avg_expr_PhaseC = c(df_minus_allPhaseCgenes$score, df_plus_allPhaseCgenes$score))
df_allPhaseCgenes_sorted <- df_allPhaseCgenes[order(df_allPhaseCgenes$CAGLENGTH_SPN_sized), ]
df_allPhaseCgenes_sorted$PHASE <- "A-B"
df_allPhaseCgenes_sorted$PHASE[which(df_allPhaseCgenes_sorted$CAGLENGTH_SPN_sized > 150)] <- "C-D-E"
df_allPhaseCgenes_sorted$index <- seq_along(df_allPhaseCgenes_sorted$avg_expr_PhaseC)
```


```R
ind_150_all <- df_allPhaseCgenes_sorted[which(df_allPhaseCgenes_sorted$CAGLENGTH_SPN_sized > 150), "index"][1]
p1 <- ggplot(df_allPhaseCgenes_sorted, aes(x = index, y = avg_expr_PhaseC, color = PHASE)) +
geom_point(aes(x = index, y = avg_expr_PhaseC, color = PHASE), size = 1) +
theme_classic() +
# theme(legend.title = element_text(size = 10),
#         legend.text = element_text(size = 10),
#         #axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size = 10),
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         panel.background = element_blank()) +
geom_vline(aes(xintercept = ind_150_all), colour="black", linetype = 2) +
xlab("SPNs sorted by CAG length") +
ylab("Average expression Phase C genes")

plot(p1)
ggsave("Handsaker_HD_PhaseC_expression.pdf", height = 8, width = 8)
```


```R
ind_150 <- df_sorted[which(df_sorted$CAGLENGTH_SPN_sized > 150), "index"][1]
p2 <- ggplot(df_sorted, aes(x = index)) +
geom_point(aes(y = score, color = PRED_PHASE, alpha = PROB_PHASE_CDE), size = 1) +
geom_point(aes(y = CAGLENGTH_SPN_sized), size = 0.1) +
theme_classic() +
# theme(legend.title = element_text(size = 10),
#         legend.text = element_text(size = 10),
#         #axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size = 10),
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         panel.background = element_blank()) +
geom_vline(aes(xintercept = ind_150), colour="black", linetype = 2) +
#geom_hline(aes(yintercept = 150), colour="black", linetype = 2) +
xlab("SPNs sorted by CAG length") +
ylab("Number of CAG repeats")
plot(p2)

ggsave("Handsaker_Phase_and_CAG_sizing_model.pdf", height = 8, width = 8)
```


```R
ind_150 <- df_sorted[which(df_sorted$CAGLENGTH_SPN_sized > 150), "index"][1]
p2_alt <- ggplot(df_sorted, aes(x = index)) +
geom_point(aes(y = PROB_PHASE_CDE, color = PRED_PHASE, alpha = score), size = 1) +
geom_point(aes(y = ifelse(CAGLENGTH_SPN_sized > 150, 1, 0))) +
theme_classic() +
# theme(legend.title = element_text(size = 10),
#         legend.text = element_text(size = 10),
#         #axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size = 10),
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         panel.background = element_blank()) +
xlab("SPNs sorted by CAG length") +
ylab("Probability of Phase C-D-E")
plot(p2_alt)

ggsave("Handsaker_Phase_and_CAG_sizing_model_alternative.pdf", height = 8, width = 8)
```


```R
df_sorted_filtered <- df_sorted[which(df_sorted$CAGLENGTH_SPN_sized > 150), ]
ind_150_filt <- df_sorted_filtered[which(df_sorted_filtered$CAGLENGTH_SPN_sized > 150), "index"][1]
p3 <- ggplot(df_sorted_filtered, aes(x = index)) +
geom_point(aes(y = score, color = PRED_PHASE, alpha = PROB_PHASE_CDE), size = 1) +
geom_point(aes(y = CAGLENGTH_SPN_sized), size = 0.1) +
theme_classic() +
# theme(legend.title = element_text(size = 10),
#         legend.text = element_text(size = 10),
#         #axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size = 10),
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         panel.background = element_blank()) +
xlim(c(ind_150_filt, max(df_sorted_filtered$index))) +
xlab("SPNs sorted by CAG length") +
ylab("Number of CAG repeats")
plot(p3)

ggsave("Handsaker_Phase_and_CAG_sizing_model_expOnly.pdf", height = 8, width = 8)
```


```R
#find number of phase C genes expressed in each SPN
num_expr_genes_training_filtered <- apply(as.matrix(expr_C_training_filtered[, selected_features]), 1, function(x) length(which(x > 0)))
num_expr_genes_training <- apply(as.matrix(expr_C_training[, selected_features]), 1, function(x) length(which(x > 0)))
num_expr_genes_test <- apply(as.matrix(expr_C_test[, selected_features]), 1, function(x) length(which(x > 0)))

summary(num_expr_genes_training_filtered)
summary(num_expr_genes_training)
summary(num_expr_genes_test)
```


```R
#check whether the number of expressed genes has an effect on prediction accuracy
ggplot(results_training_filtered, aes(x = CAG_measured, y = CAG_predicted, color = num_expr_genes)) +
geom_point() +
xlim(c(0, 800)) + ylim(c(0, 800)) +
geom_abline()
```


```R
SPN_HD <- subset(x = SPN_filt, subset = CONDITION != "CTRL")
SPN_CTRL <- subset(x = SPN_filt, subset = CONDITION == "CTRL")
```

# Import_Lee_dataset


```R
# #import Lee counts
# Lee_counts <- Read10X("/mnt/projects/labs/CLAB/PROJECT_ELongATE/Lee/data/renamed_data")
# #create single-cell experiment
# Lee_sce <- SingleCellExperiment(list(counts = Lee_counts))
# #create Seurat object
# Lee <- CreateSeuratObject(counts(Lee_sce), project = "Lee", min.cells = 0, min.features = 0)
# #read metadata
# Lee_metadata <- read.table("/mnt/projects/labs/CLAB/PROJECT_ELongATE/Lee/data/renamed_data/GSE152058_human_snRNA_processed_coldata.tsv", header = TRUE, sep = "\t")
# #Lee_metadata
# #table(Lee_metadata$CellType)
# rownames(Lee_metadata) <- Lee_metadata$Barcode
# Lee@meta.data <- cbind(Lee@meta.data, Lee_metadata[rownames(Lee@meta.data), ])
# #head(Lee@meta.data)
# Lee[["percent.mt"]] <- PercentageFeatureSet(Lee, pattern = "^MT-")
```


```R
#read metadata
Lee_metadata <- read.table("/mnt/projects/labs/CLAB/PROJECT_ELongATE/Lee/data/renamed_data/GSE152058_human_snRNA_processed_coldata.tsv", header = TRUE, sep = "\t")
#Lee_metadata
table(Lee_metadata$CellType)
```


```R
#find the number of cells for each MSN subtype
tmp <- lapply(split(Lee_metadata, Lee_metadata$NBB_ID), function(x) {
    x_caudate <- x[which(x$Region == "Caudate"), ]
    x_putamen <- x[which(x$Region == "Putamen"), ]
    num_cells_caudate <- data.frame(num_cells_caudate = dim(x_caudate)[1])
    num_D1_spns_caudate <- data.frame(num_D1_SPN_caudate = dim(x[grep(x = x_caudate$CellType, pattern = "D1_MSN"), ])[1])
    num_D2_spns_caudate <- data.frame(num_D2_SPN_caudate = dim(x[grep(x = x_caudate$CellType, pattern = "D2_MSN"), ])[1])
    num_cells_putamen <- data.frame(num_cells_putamen = dim(x_putamen)[1])
    num_D1_spns_putamen <- data.frame(num_D1_SPN_putamen = dim(x[grep(x = x_putamen$CellType, pattern = "D1_MSN"), ])[1])
    num_D2_spns_putamen <- data.frame(num_D2_SPN_putamen = dim(x[grep(x = x_putamen$CellType, pattern = "D1_MSN"), ])[1])
    return(data.frame(num_cells_caudate, num_D1_spns_caudate, num_D2_spns_caudate, num_cells_putamen, num_D1_spns_putamen, num_D2_spns_putamen))
})

num_cells_Lee <- do.call(rbind, tmp)
num_cells_Lee
```


```R
# #subset MSNs only
# Lee_MSN <- subset(x = Lee, subset = CellType == "D1_MSN" | CellType == "D2_MSN")
```

# Lee_Caudate_only


```R
# #subset caudate only
# Lee_MSN_Caudate <- subset(x = Lee_MSN, subset = Region == "Caudate")
```


```R
# #Run SCTransform
# Lee_MSN_Caudate <- SCTransform(Lee_MSN_Caudate, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
# #Lee_MSN_Caudate <- SCTransform(Lee_MSN_Caudate)
# cat(sprintf("Num. features = %d; Num. SPN = %d\n", dim(Lee_MSN_Caudate)[1], dim(Lee_MSN_Caudate)[2]))
# #Run PCA and UMAP
# Lee_MSN_Caudate <- RunPCA(Lee_MSN_Caudate, features = VariableFeatures(object = Lee_MSN_Caudate), npcs = 50, reduction.name = "pca")
# ElbowPlot(Lee_MSN_Caudate, reduction = "pca", ndims = 50)
# num_dims <- 40
# Lee_MSN_Caudate <- RunUMAP(Lee_MSN_Caudate, dims = 1:num_dims, reduction = "pca", reduction.name = "umap", reduction.key = "umap")
# DimPlot(Lee_MSN_Caudate, group.by = "Batch", reduction = "umap", pt.size = 1)
# Lee_MSN_Caudate@meta.data$Batch <- factor(Lee_MSN_Caudate@meta.data$Batch)
# Lee_MSN_Caudate <- RunHarmony(object = Lee_MSN_Caudate, reduction.use = "pca", group.by.vars = "Batch", 
#                            reduction.save = "harmonyPca", assay = "SCT", verbose = FALSE, normalization.method = "SCT")
# Lee_MSN_Caudate <- RunUMAP(Lee_MSN_Caudate, dims = 1:num_dims, reduction = "harmonyPca", reduction.name = "harmony")
# DimPlot(Lee_MSN_Caudate, group.by = "Batch", reduction = "harmony", pt.size = 1)
```


```R
# saveRDS(object = Lee_MSN_Caudate, file = "/mnt/projects/labs/CLAB/PROJECT_ELongATE/Lee/data/renamed_data/Lee_MSN_Caudate.rds")
```


```R
Lee_MSN_Caudate <- readRDS("/mnt/projects/labs/CLAB/PROJECT_ELongATE/Lee/data/renamed_data/Lee_MSN_Caudate.rds")
```


```R
#plots
options(repr.plot.width=20, repr.plot.height=20) 
FeaturePlot(Lee_MSN_Caudate, reduction = "harmony", features = "nCount_RNA", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Lee_MSN_Caudate_NCOUNTS.pdf", height = 8, width = 8)
FeaturePlot(Lee_MSN_Caudate, reduction = "harmony", features = "nFeature_RNA", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Lee_MSN_Caudate_NFEATURES.pdf", height = 8, width = 8)
FeaturePlot(Lee_MSN_Caudate, reduction = "harmony", features = "percent.mt", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Lee_MSN_Caudate_PERCMT.pdf", height = 8, width = 8)
DimPlot(Lee_MSN_Caudate, group.by = "CellType", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Lee_MSN_Caudate_CELLTYPE.pdf", height = 8, width = 8)
DimPlot(Lee_MSN_Caudate, group.by = "Batch", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Lee_MSN_Caudate_BATCH.pdf", width = 8, height = 8)
DimPlot(Lee_MSN_Caudate, group.by = "NBB_ID", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Lee_MSN_Caudate_NBBID.pdf", width = 8, height = 8)
DimPlot(Lee_MSN_Caudate, group.by = "Region", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Lee_MSN_Caudate_REGION.pdf", width = 8, height = 8)
```


```R
#extract phase C genes in the model that are also expressed in the dataset
phaseC_genes_Lee_MSN_Caudate <- intersect(colnames(phase_all), rownames(Lee_MSN_Caudate))
phaseC_genes_notExpLee_MSN_Caudate_names <- setdiff(colnames(phase_all), rownames(Lee_MSN_Caudate))
#phaseC_genes_notExpLee_MSN_Caudate_names
phaseC_genes_notExpLee_MSN_Caudate <- matrix(data = 0, nrow = length(phaseC_genes_notExpLee_MSN_Caudate_names), ncol = dim(Lee_MSN_Caudate@assays$SCT$data)[2])
rownames(phaseC_genes_notExpLee_MSN_Caudate) <- phaseC_genes_notExpLee_MSN_Caudate_names
colnames(phaseC_genes_notExpLee_MSN_Caudate) <- colnames(Lee_MSN_Caudate@assays$SCT$data)
phase_counts_Lee_MSN_Caudate <- Lee_MSN_Caudate@assays$SCT$data[phaseC_genes_Lee_MSN_Caudate, ]
phase_counts_Lee_MSN_Caudate <- rbind(phase_counts_Lee_MSN_Caudate, phaseC_genes_notExpLee_MSN_Caudate)
```


```R
phase_counts_Lee_MSN_Caudate <- phase_counts_Lee_MSN_Caudate[colnames(phase_all), ]
phase_test_Lee_MSN_Caudate <- data.frame(t(phase_counts_Lee_MSN_Caudate))
#head(phase_test_Lee_MSN_Caudate)
```


```R
class_thr <- perf_training_test_noNA$thr
predicted_phase_test_Lee_MSN_Caudate_wprob <- predict(model_phase, newx = as.matrix(phase_test_Lee_MSN_Caudate), s = lambda_val, type = "response") #glmnet
predicted_phase_test_Lee_MSN_Caudate <- rep("A-B", dim(predicted_phase_test_Lee_MSN_Caudate_wprob)[1])  
predicted_phase_test_Lee_MSN_Caudate[which(predicted_phase_test_Lee_MSN_Caudate_wprob > class_thr)] <- "C-D-E"
Lee_MSN_Caudate[["PREDICTED_PHASE"]] <- NA
Lee_MSN_Caudate@meta.data[rownames(predicted_phase_test_Lee_MSN_Caudate_wprob), "PREDICTED_PHASE"] <- predicted_phase_test_Lee_MSN_Caudate
Lee_MSN_Caudate[["PROB_PHASE_CDE"]] <- NA
Lee_MSN_Caudate@meta.data[rownames(predicted_phase_test_Lee_MSN_Caudate_wprob), "PROB_PHASE_CDE"] <- predicted_phase_test_Lee_MSN_Caudate_wprob
```


```R
lapply(split(Lee_MSN_Caudate@meta.data$PREDICTED_PHASE, Lee_MSN_Caudate@meta.data$Condition), function(x) {
    table(x)
})
```


```R
#Lee_MSN_Caudate@meta.data$Grade <- gsub(x = Lee_MSN_Caudate@meta.data$Grade, pattern = "HD", replacement = "")
table(Lee_MSN_Caudate@meta.data$Grade)
tmp <- gsub(x = gsub(x = paste0("HD", Lee_MSN_Caudate@meta.data$Grade), pattern = "HDControl", replacement = "Control"), pattern = "HDHD", replacement = "HD")
table(tmp)
```


```R
#do plots
DimPlot(Lee_MSN_Caudate, group.by = "PREDICTED_PHASE", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Lee_MSN_Caudate_PREDPHASE.pdf", width = 8, height = 8)
DimPlot(Lee_MSN_Caudate, group.by = "PREDICTED_PHASE", reduction = "harmony", split.by = "Condition", pt.size = 1)
ggsave("UMAP_Lee_MSN_Caudate_PREDPHASE_splitByCondition.pdf", width = 8, height = 8)
FeaturePlot(Lee_MSN_Caudate, features = "PROB_PHASE_CDE", cols = c("blue", "red"), pt.size = 1, reduction = "harmony")
ggsave("UMAP_Lee_MSN_Caudate_PROBPREDPHASE.pdf", width = 8, height = 8)
FeaturePlot(Lee_MSN_Caudate, features = "PROB_PHASE_CDE", cols = c("blue", "red"), pt.size = 1, split.by = "Condition", reduction = "harmony")
ggsave("UMAP_Lee_MSN_Caudate_PROBPREDPHASE_splitByCondition.pdf", width = 8, height = 8)
tmp <- gsub(x = gsub(x = paste0("HD", Lee_MSN_Caudate@meta.data$Grade), pattern = "HDControl", replacement = "Control"), pattern = "HDHD", replacement = "HD")
Lee_MSN_Caudate@meta.data$Grade <- factor(tmp, levels = c("Control", "HD2", "HD3", "HD4"))
DimPlot(Lee_MSN_Caudate, group.by = "Grade", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Lee_MSN_Caudate_GRADE.pdf", width = 8, height = 8)
DimPlot(Lee_MSN_Caudate, group.by = "Grade", reduction = "harmony", split.by = "Condition", pt.size = 1)
ggsave("UMAP_Lee_MSN_Caudate_GRADE_splitByCondition.pdf", width = 8, height = 8)
```


```R
#create metadata file for plotting
Lee_Caudate_Grade_CDE_tmp <- lapply(split(Lee_MSN_Caudate@meta.data, Lee_MSN_Caudate@meta.data$NBB_ID), function(x) {
    GRADE <- factor(unique(x$Grade), levels = c("Control", "HD1", "HD2", "HD3", "HD4"))
    SN <- unique(x$NBB_ID)
    #NUM_MSN <- length(x$PREDICTED_PHASE)
    tmp <- lapply(split(x, x$CellType), function(y) {
        NUM_MSN_CT <- length(y$PREDICTED_PHASE)
        NUM_MSN_CT_SQ <- NUM_MSN_CT**2
        CELL_TYPE <- y$CellType
        if (length(unique(y$PREDICTED_PHASE)) > 1) {
          val <- 100*table(y$PREDICTED_PHASE)[2]/(table(y$PREDICTED_PHASE)[1] + table(y$PREDICTED_PHASE)[2])
        } else {
          val <- 0
        }
        names(val) <- "Fract_CDE";
        Fract_CDE <- rep(val, NUM_MSN_CT_SQ);
    return(data.frame(NUM_MSN_CT_SQ, CELL_TYPE, Fract_CDE))})
    df <- do.call(rbind, tmp)
    #print(df)
    #NUM_MSN <- length(x$PREDICTED_PHASE)
    #Fract_CDE <- 100*length(which(x$PREDICTED_PHASE == "C-D-E"))/length(x$PREDICTED_PHASE)
    #CELL_TYPE <- x$sub_type_4
    #return(df)
    #return(data.frame(SAMPLE = rep(SN, df$NUM_MSN_CT_SQ), GRADE = rep(GRADE, df$NUM_MSN_CT_SQ), CELL_TYPE = df$CELL_TYPE, Fract_CDE = df$Fract_CDE, NUM_MSN_SQ=df$NUM_MSN_CT_SQ))
    return(data.frame(SAMPLE = SN, GRADE = GRADE, CELL_TYPE = df$CELL_TYPE, Fract_CDE = df$Fract_CDE, NUM_MSN_SQ=df$NUM_MSN_CT_SQ))
})
Lee_Caudate_Grade_CDE <- do.call(rbind, Lee_Caudate_Grade_CDE_tmp)
Lee_Caudate_Grade_CDE <- Lee_Caudate_Grade_CDE[order(Lee_Caudate_Grade_CDE$GRADE), ]
Lee_Caudate_Grade_CDE$SAMPLE <- factor(Lee_Caudate_Grade_CDE$SAMPLE, levels = unique(Lee_Caudate_Grade_CDE$SAMPLE))
#head(Lee_Caudate_Grade_CDE)
uni_Lee_Caudate_Grade_CDE <- unique(Lee_Caudate_Grade_CDE)
uni_Lee_Caudate_Grade_CDE$NUM_SPN <- sqrt(uni_Lee_Caudate_Grade_CDE$NUM_MSN_SQ)
#uni_Lee_Caudate_Grade_CDE
```


```R
#num_cells_Lee
Fract_SPN_Caudate <- c(num_cells_Lee["3345", "num_D1_SPN_caudate"]/num_cells_Lee["3345", "num_cells_caudate"],
                       num_cells_Lee["3345", "num_D2_SPN_caudate"]/num_cells_Lee["3345", "num_cells_caudate"],
                       num_cells_Lee["4294", "num_D1_SPN_caudate"]/num_cells_Lee["4294", "num_cells_caudate"],
                       num_cells_Lee["4294", "num_D2_SPN_caudate"]/num_cells_Lee["4294", "num_cells_caudate"],
                       num_cells_Lee["4308", "num_D1_SPN_caudate"]/num_cells_Lee["4308", "num_cells_caudate"],
                       num_cells_Lee["4308", "num_D2_SPN_caudate"]/num_cells_Lee["4308", "num_cells_caudate"],
                       num_cells_Lee["4494", "num_D1_SPN_caudate"]/num_cells_Lee["4494", "num_cells_caudate"],
                       num_cells_Lee["4494", "num_D2_SPN_caudate"]/num_cells_Lee["4494", "num_cells_caudate"],
                       num_cells_Lee["4621", "num_D1_SPN_caudate"]/num_cells_Lee["4621", "num_cells_caudate"],
                       num_cells_Lee["4621", "num_D2_SPN_caudate"]/num_cells_Lee["4621", "num_cells_caudate"],
                       num_cells_Lee["A39R", "num_D1_SPN_caudate"]/num_cells_Lee["A39R", "num_cells_caudate"],
                       num_cells_Lee["A39R", "num_D2_SPN_caudate"]/num_cells_Lee["A39R", "num_cells_caudate"],
                       num_cells_Lee["A47L", "num_D1_SPN_caudate"]/num_cells_Lee["A47L", "num_cells_caudate"],
                       num_cells_Lee["A47L", "num_D2_SPN_caudate"]/num_cells_Lee["A47L", "num_cells_caudate"],
                       num_cells_Lee["3730", "num_D1_SPN_caudate"]/num_cells_Lee["3730", "num_cells_caudate"],
                       num_cells_Lee["3730", "num_D2_SPN_caudate"]/num_cells_Lee["3730", "num_cells_caudate"],
                       num_cells_Lee["3881", "num_D1_SPN_caudate"]/num_cells_Lee["3881", "num_cells_caudate"],
                       num_cells_Lee["3881", "num_D2_SPN_caudate"]/num_cells_Lee["3881", "num_cells_caudate"],
                       num_cells_Lee["2030", "num_D1_SPN_caudate"]/num_cells_Lee["2030", "num_cells_caudate"],
                       num_cells_Lee["2030", "num_D2_SPN_caudate"]/num_cells_Lee["2030", "num_cells_caudate"],
                       num_cells_Lee["2952", "num_D1_SPN_caudate"]/num_cells_Lee["2952", "num_cells_caudate"],
                       num_cells_Lee["2952", "num_D2_SPN_caudate"]/num_cells_Lee["2952", "num_cells_caudate"],
                       num_cells_Lee["4254", "num_D1_SPN_caudate"]/num_cells_Lee["4254", "num_cells_caudate"],
                       num_cells_Lee["4254", "num_D2_SPN_caudate"]/num_cells_Lee["4254", "num_cells_caudate"],
                       num_cells_Lee["2665", "num_D1_SPN_caudate"]/num_cells_Lee["2665", "num_cells_caudate"],
                       num_cells_Lee["2665", "num_D2_SPN_caudate"]/num_cells_Lee["2665", "num_cells_caudate"])
uni_Lee_Caudate_Grade_CDE$Fract_SPN <- Fract_SPN_Caudate
uni_Lee_Caudate_Grade_CDE
```


```R
unique(uni_Lee_Caudate_Grade_CDE[, c("SAMPLE", "GRADE")])
```


```R
#plot fraction of SPN in C-D-E phase by HD grade
uni_Lee_Caudate_Grade_CDE_HD <- uni_Lee_Caudate_Grade_CDE[which(uni_Lee_Caudate_Grade_CDE$GRADE != "Control"), ]
#uni_Lee_Caudate_Grade_CDE_CTRL <- uni_Lee_Caudate_Grade_CDE[which(uni_Lee_Caudate_Grade_CDE$GRADE == "Control"), ]

ggplot(uni_Lee_Caudate_Grade_CDE_HD, aes(x = GRADE, y = Fract_CDE, color = CELL_TYPE)) +
geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.1, height = 0)) +
scale_color_manual(values = c("#F8766D", "#619CFF")) +
theme_classic() +
xlab("HD GRADE") +
ylim(c(0, 100)) +
ylab("Fraction cells in C-D-E phase (%)") 
ggsave("Lee_Caudate_HD_GRADE_vs_fraction_SPN_CDE.pdf", width = 6, height = 6)
```


```R
# ggplot(uni_Lee_Caudate_Grade_CDE[which(uni_Lee_Caudate_Grade_CDE$NUM_SPN > 50), ], aes(x = GRADE, y = Fract_CDE, color = CELL_TYPE)) +
# geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.1, height = 0)) +
# scale_color_manual(values = c("#F8766D", "#619CFF")) +
# theme_classic() +
# xlab("HD GRADE") +
# ylim(c(0, 100)) +
# ylab("Fraction cells in C-D-E phase (%)")
# ggsave("Lee_Caudate_GRADE_vs_fraction_SPN_CDE_gt50.pdf", width = 6, height = 6)
```


```R
#subset HD
Lee_HD_MSN_Caudate <- subset(x = Lee_MSN_Caudate, subset = Condition == "HD")
```

# Lee_Putamen_only


```R
# #subset caudate only
# Lee_MSN_Putamen <- subset(x = Lee_MSN, subset = Region == "Putamen")
```


```R
# #Run SCTransform
# Lee_MSN_Putamen <- SCTransform(Lee_MSN_Putamen, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
# #Lee_MSN_Putamen <- SCTransform(Lee_MSN_Putamen)
# cat(sprintf("Num. features = %d; Num. SPN = %d\n", dim(Lee_MSN_Putamen)[1], dim(Lee_MSN_Putamen)[2]))
# #Run PCA and UMAP for HD donors
# Lee_MSN_Putamen <- RunPCA(Lee_MSN_Putamen, features = VariableFeatures(object = Lee_MSN_Putamen), npcs = 50, reduction.name = "pca")
# ElbowPlot(Lee_MSN_Putamen, reduction = "pca", ndims = 50)
# num_dims <- 40
# Lee_MSN_Putamen <- RunUMAP(Lee_MSN_Putamen, dims = 1:num_dims, reduction = "pca", reduction.name = "umap", reduction.key = "umap")
# DimPlot(Lee_MSN_Putamen, group.by = "Batch", reduction = "umap", pt.size = 1)
# Lee_MSN_Putamen@meta.data$Batch <- factor(Lee_MSN_Putamen@meta.data$Batch)
# Lee_MSN_Putamen <- RunHarmony(object = Lee_MSN_Putamen, reduction.use = "pca", group.by.vars = "Batch", 
#                            reduction.save = "harmonyPca", assay = "SCT", verbose = FALSE, normalization.method = "SCT")
# Lee_MSN_Putamen <- RunUMAP(Lee_MSN_Putamen, dims = 1:num_dims, reduction = "harmonyPca", reduction.name = "harmony")
# DimPlot(Lee_MSN_Putamen, group.by = "Batch", reduction = "harmony", pt.size = 1)
```


```R
# saveRDS(object = Lee_MSN_Putamen, file = "/mnt/projects/labs/CLAB/PROJECT_ELongATE/Lee/data/renamed_data/Lee_MSN_Putamen.rds")
```


```R
Lee_MSN_Putamen <- readRDS("/mnt/projects/labs/CLAB/PROJECT_ELongATE/Lee/data/renamed_data/Lee_MSN_Putamen.rds")
```


```R
#do plots
options(repr.plot.width=20, repr.plot.height=20) 
FeaturePlot(Lee_MSN_Putamen, reduction = "harmony", features = "nCount_RNA", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Lee_MSN_Putamen_NCOUNTS.pdf", height = 8, width = 8)
FeaturePlot(Lee_MSN_Putamen, reduction = "harmony", features = "nFeature_RNA", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Lee_MSN_Putamen_NFEATURES.pdf", height = 8, width = 8)
FeaturePlot(Lee_MSN_Putamen, reduction = "harmony", features = "percent.mt", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Lee_MSN_Putamen_PERCMT.pdf", height = 8, width = 8)
DimPlot(Lee_MSN_Putamen, group.by = "CellType", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Lee_MSN_Putamen_CELLTYPE.pdf", height = 8, width = 8)
DimPlot(Lee_MSN_Putamen, group.by = "Batch", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Lee_MSN_Putamen_BATCH.pdf", width = 8, height = 8)
DimPlot(Lee_MSN_Putamen, group.by = "NBB_ID", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Lee_MSN_Putamen_NBBID.pdf", width = 8, height = 8)
DimPlot(Lee_MSN_Putamen, group.by = "Region", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Lee_MSN_Putamen_REGION.pdf", width = 8, height = 8)
```


```R
#extract phase C genes in the model that are also expressed in the dataset
phaseC_genes_Lee_MSN_Putamen <- intersect(colnames(phase_all), rownames(Lee_MSN_Putamen))
phaseC_genes_notExpLee_MSN_Putamen_names <- setdiff(colnames(phase_all), rownames(Lee_MSN_Putamen))
#phaseC_genes_notExpLee_MSN_Putamen_names
phaseC_genes_notExpLee_MSN_Putamen <- matrix(data = 0, nrow = length(phaseC_genes_notExpLee_MSN_Putamen_names), ncol = dim(Lee_MSN_Putamen@assays$SCT$data)[2])
rownames(phaseC_genes_notExpLee_MSN_Putamen) <- phaseC_genes_notExpLee_MSN_Putamen_names
colnames(phaseC_genes_notExpLee_MSN_Putamen) <- colnames(Lee_MSN_Putamen@assays$SCT$data)
phase_counts_Lee_MSN_Putamen <- Lee_MSN_Putamen@assays$SCT$data[phaseC_genes_Lee_MSN_Putamen, ]
phase_counts_Lee_MSN_Putamen <- rbind(phase_counts_Lee_MSN_Putamen, phaseC_genes_notExpLee_MSN_Putamen)
```


```R
phase_counts_Lee_MSN_Putamen <- phase_counts_Lee_MSN_Putamen[colnames(phase_all), ]
phase_test_Lee_MSN_Putamen <- data.frame(t(phase_counts_Lee_MSN_Putamen))
#head(phase_test_Lee_MSN_Putamen)
```


```R
class_thr <- perf_training_test_noNA$thr
#class_thr <- 0.1
predicted_phase_test_Lee_MSN_Putamen_wprob <- predict(model_phase,
                                                      newx = as.matrix(phase_test_Lee_MSN_Putamen),
                                                      s = lambda_val,
                                                      type = "response")
predicted_phase_test_Lee_MSN_Putamen <- rep("A-B", dim(predicted_phase_test_Lee_MSN_Putamen_wprob)[1])  
predicted_phase_test_Lee_MSN_Putamen[which(predicted_phase_test_Lee_MSN_Putamen_wprob > class_thr)] <- "C-D-E"
Lee_MSN_Putamen[["PREDICTED_PHASE"]] <- NA
Lee_MSN_Putamen@meta.data[rownames(predicted_phase_test_Lee_MSN_Putamen_wprob), "PREDICTED_PHASE"] <- predicted_phase_test_Lee_MSN_Putamen
Lee_MSN_Putamen[["PROB_PHASE_CDE"]] <- NA
Lee_MSN_Putamen@meta.data[rownames(predicted_phase_test_Lee_MSN_Putamen_wprob), "PROB_PHASE_CDE"] <- predicted_phase_test_Lee_MSN_Putamen_wprob
```


```R
lapply(split(Lee_MSN_Putamen@meta.data$PREDICTED_PHASE, Lee_MSN_Putamen@meta.data$Condition), function(x) {
    table(x)
})
```


```R
# lapply(split(Lee_MSN_Putamen@meta.data$PROB_PHASE_CDE, Lee_MSN_Putamen@meta.data$Condition), function(x) {
#     quantile(x, probs = seq(0.8, 1, length.out = 101))
# })
```


```R
#summary(apply(phase_test_Lee_MSN_Putamen, 1, function(x) length(which(x > 0))))
table(predicted_phase_test_Lee_MSN_Putamen)
```


```R
#Lee_MSN_Putamen@meta.data$Grade <- gsub(x = Lee_MSN_Putamen@meta.data$Grade, pattern = "HD", replacement = "")
table(Lee_MSN_Putamen@meta.data$Grade)
tmp <- gsub(x = gsub(x = paste0("HD", Lee_MSN_Putamen@meta.data$Grade), pattern = "HDControl", replacement = "Control"), pattern = "HDHD", replacement = "HD")
table(tmp)
```


```R
#do plots
DimPlot(Lee_MSN_Putamen, group.by = "PREDICTED_PHASE", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Lee_MSN_Putamen_PREDPHASE.pdf", width = 8, height = 8)
DimPlot(Lee_MSN_Putamen, group.by = "PREDICTED_PHASE", reduction = "harmony", split.by = "Condition", pt.size = 1)
ggsave("UMAP_Lee_MSN_Putamen_PREDPHASE_splitByCondition.pdf", width = 8, height = 8)
FeaturePlot(Lee_MSN_Putamen, features = "PROB_PHASE_CDE", cols = c("blue", "red"), pt.size = 1, reduction = "harmony")
ggsave("UMAP_Lee_MSN_Putamen_PROBPREDPHASE.pdf", width = 8, height = 8)
FeaturePlot(Lee_MSN_Putamen, features = "PROB_PHASE_CDE", cols = c("blue", "red"), pt.size = 1, split.by = "Condition", reduction = "harmony")
ggsave("UMAP_Lee_MSN_Putamen_PROBPREDPHASE_splitByCondition.pdf", width = 8, height = 8)
tmp <- gsub(x = gsub(x = paste0("HD", Lee_MSN_Putamen@meta.data$Grade), pattern = "HDControl", replacement = "Control"), pattern = "HDHD", replacement = "HD")
Lee_MSN_Putamen@meta.data$Grade <- factor(tmp, levels = c("Control", "HD2", "HD3", "HD4"))
DimPlot(Lee_MSN_Putamen, group.by = "Grade", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Lee_MSN_Putamen_GRADE.pdf", width = 8, height = 8)
DimPlot(Lee_MSN_Putamen, group.by = "Grade", reduction = "harmony", split.by = "Condition", pt.size = 1)
ggsave("UMAP_Lee_MSN_Putamen_GRADE_splitByCondition.pdf", width = 8, height = 8)
```


```R
#create metadata file for plotting
Lee_Putamen_Grade_CDE_tmp <- lapply(split(Lee_MSN_Putamen@meta.data, Lee_MSN_Putamen@meta.data$NBB_ID), function(x) {
    GRADE <- factor(unique(x$Grade), levels = c("Control", "HD1", "HD2", "HD3", "HD4"))
    SN <- unique(x$NBB_ID)
    #NUM_MSN <- length(x$PREDICTED_PHASE)
    tmp <- lapply(split(x, x$CellType), function(y) {
        NUM_MSN_CT <- length(y$PREDICTED_PHASE)
        NUM_MSN_CT_SQ <- NUM_MSN_CT**2
        CELL_TYPE <- y$CellType
        if (length(unique(y$PREDICTED_PHASE)) > 1) {
          val <- 100*table(y$PREDICTED_PHASE)[2]/(table(y$PREDICTED_PHASE)[1] + table(y$PREDICTED_PHASE)[2])
        } else {
          val <- 0
        }
        names(val) <- "Fract_CDE";
        Fract_CDE <- rep(val, NUM_MSN_CT_SQ);
    return(data.frame(NUM_MSN_CT_SQ, CELL_TYPE, Fract_CDE))})
    df <- do.call(rbind, tmp)
    #print(df)
    #NUM_MSN <- length(x$PREDICTED_PHASE)
    #Fract_CDE <- 100*length(which(x$PREDICTED_PHASE == "C-D-E"))/length(x$PREDICTED_PHASE)
    #CELL_TYPE <- x$sub_type_4
    #return(df)
    #return(data.frame(SAMPLE = rep(SN, df$NUM_MSN_CT_SQ), GRADE = rep(GRADE, df$NUM_MSN_CT_SQ), CELL_TYPE = df$CELL_TYPE, Fract_CDE = df$Fract_CDE, NUM_MSN_SQ=df$NUM_MSN_CT_SQ))
    return(data.frame(SAMPLE = SN, GRADE = GRADE, CELL_TYPE = df$CELL_TYPE, Fract_CDE = df$Fract_CDE, NUM_MSN_SQ=df$NUM_MSN_CT_SQ))
})
Lee_Putamen_Grade_CDE <- do.call(rbind, Lee_Putamen_Grade_CDE_tmp)
Lee_Putamen_Grade_CDE <- Lee_Putamen_Grade_CDE[order(Lee_Putamen_Grade_CDE$GRADE), ]
Lee_Putamen_Grade_CDE$SAMPLE <- factor(Lee_Putamen_Grade_CDE$SAMPLE, levels = unique(Lee_Putamen_Grade_CDE$SAMPLE))
#head(Lee_Putamen_Grade_CDE)
uni_Lee_Putamen_Grade_CDE <- unique(Lee_Putamen_Grade_CDE)
uni_Lee_Putamen_Grade_CDE$NUM_SPN <- sqrt(uni_Lee_Putamen_Grade_CDE$NUM_MSN_SQ)
#uni_Lee_Putamen_Grade_CDE
```


```R
Fract_SPN_Putamen <- c(num_cells_Lee["4294", "num_D1_SPN_putamen"]/num_cells_Lee["4294", "num_cells_putamen"],
                       num_cells_Lee["4294", "num_D2_SPN_putamen"]/num_cells_Lee["4294", "num_cells_putamen"],
                       num_cells_Lee["4308", "num_D1_SPN_putamen"]/num_cells_Lee["4308", "num_cells_putamen"],
                       num_cells_Lee["4308", "num_D2_SPN_putamen"]/num_cells_Lee["4308", "num_cells_putamen"],
                       num_cells_Lee["4494", "num_D1_SPN_putamen"]/num_cells_Lee["4494", "num_cells_putamen"],
                       num_cells_Lee["4494", "num_D2_SPN_putamen"]/num_cells_Lee["4494", "num_cells_putamen"],
                       num_cells_Lee["A39R", "num_D1_SPN_putamen"]/num_cells_Lee["A39R", "num_cells_putamen"],
                       num_cells_Lee["A39R", "num_D2_SPN_putamen"]/num_cells_Lee["A39R", "num_cells_putamen"],
                       num_cells_Lee["A47L", "num_D1_SPN_putamen"]/num_cells_Lee["A47L", "num_cells_putamen"],
                       num_cells_Lee["A47L", "num_D2_SPN_putamen"]/num_cells_Lee["A47L", "num_cells_putamen"],
                       num_cells_Lee["3881", "num_D1_SPN_putamen"]/num_cells_Lee["3881", "num_cells_putamen"],
                       num_cells_Lee["3881", "num_D2_SPN_putamen"]/num_cells_Lee["3881", "num_cells_putamen"],
                       num_cells_Lee["4225", "num_D1_SPN_putamen"]/num_cells_Lee["4225", "num_cells_putamen"],
                       num_cells_Lee["4225", "num_D2_SPN_putamen"]/num_cells_Lee["4225", "num_cells_putamen"],
                       num_cells_Lee["2952", "num_D1_SPN_putamen"]/num_cells_Lee["2952", "num_cells_putamen"],
                       num_cells_Lee["2952", "num_D2_SPN_putamen"]/num_cells_Lee["2952", "num_cells_putamen"],
                       num_cells_Lee["4254", "num_D1_SPN_putamen"]/num_cells_Lee["4254", "num_cells_putamen"],
                       num_cells_Lee["4254", "num_D2_SPN_putamen"]/num_cells_Lee["4254", "num_cells_putamen"],
                       num_cells_Lee["2665", "num_D1_SPN_putamen"]/num_cells_Lee["2665", "num_cells_putamen"],
                       #num_cells_Lee["2665", "num_D2_SPN_putamen"]/num_cells_Lee["2665", "num_cells_putamen"],
                       num_cells_Lee["2903", "num_D1_SPN_putamen"]/num_cells_Lee["2903", "num_cells_putamen"],
                       num_cells_Lee["2903", "num_D2_SPN_putamen"]/num_cells_Lee["2665", "num_cells_putamen"])
uni_Lee_Putamen_Grade_CDE$Fract_SPN <- Fract_SPN_Putamen
uni_Lee_Putamen_Grade_CDE
```


```R
unique(uni_Lee_Putamen_Grade_CDE[, c("SAMPLE", "GRADE")])
```


```R
#plot fraction of SPN in C-D-E phase by HD grade
uni_Lee_Putamen_Grade_CDE_HD <- uni_Lee_Putamen_Grade_CDE[which(uni_Lee_Putamen_Grade_CDE$GRADE != "Control"), ]
#uni_Lee_Putamen_Grade_CDE_CTRL <- uni_Lee_Putamen_Grade_CDE[which(uni_Lee_Putamen_Grade_CDE$GRADE == "Control"), ]

ggplot(uni_Lee_Putamen_Grade_CDE_HD, aes(x = GRADE, y = Fract_CDE, color = CELL_TYPE)) +
geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.1, height = 0)) +
scale_color_manual(values = c("#F8766D", "#619CFF")) +
theme_classic() +
xlab("HD GRADE") +
ylim(c(0, 100)) +
ylab("Fraction cells in C-D-E phase (%)") 
ggsave("Lee_Putamen_HD_GRADE_vs_fraction_SPN_CDE.pdf", width = 6, height = 6)
```


```R
# ggplot(uni_Lee_Putamen_Grade_CDE[which(uni_Lee_Putamen_Grade_CDE$NUM_SPN > 50), ], aes(x = GRADE, y = Fract_CDE, color = CELL_TYPE)) +
# geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.1, height = 0)) +
# scale_color_manual(values = c("#F8766D", "#619CFF")) +
# theme_classic() +
# xlab("HD GRADE") +
# ylim(c(0, 100)) +
# ylab("Fraction cells in C-D-E phase (%)")
# ggsave("Lee_Putamen_GRADE_vs_fraction_SPN_CDE_gt50.pdf", width = 6, height = 6)
```


```R
str(Lee_Caudate_Grade_CDE)
head(Lee_Caudate_Grade_CDE)
str(Lee_Putamen_Grade_CDE)
head(Lee_Putamen_Grade_CDE)
```


```R
colnames(Lee_MSN_Caudate@meta.data)
colnames(Lee_MSN_Putamen@meta.data)
```


```R
#subset HD
Lee_HD_MSN_Putamen <- subset(x = Lee_MSN_Putamen, subset = Condition == "HD")
```

# Import_Paryani_dataset


```R
Paryani <- readRDS("/mnt/projects/labs/CLAB/PROJECT_ELongATE/Paryani/neuron_hd_obj_acc_caud_paper.rds")
```


```R
head(Paryani@meta.data)
table(Paryani@meta.data$broad_type)
table(Paryani@meta.data$sub_type_4)
table(Paryani@meta.data$Condition)
table(Paryani@meta.data$CAG)
table(Paryani@meta.data$Batch)
```


```R
#add metadata
tmp <- lapply(split(Paryani@meta.data, Paryani@meta.data$Donor), function(x) {
    x_caudate <- x[which(x$Region == "Caudate"), ]
    x_accumbens <- x[which(x$Region == "Accumbens"), ]
    num_cells_caudate <- data.frame(num_cells_caudate = dim(x_caudate)[1])
    num_dspn1_caudate <- data.frame(num_dSPN1_caudate = dim(x[grep(x = x_caudate$sub_type_4, pattern = "dSPN_1"), ])[1])
    num_dspn2_caudate <- data.frame(num_dSPN2_caudate = dim(x[grep(x = x_caudate$sub_type_4, pattern = "dSPN_2"), ])[1])
    num_ispn1_caudate <- data.frame(num_iSPN1_caudate = dim(x[grep(x = x_caudate$sub_type_4, pattern = "iSPN_1"), ])[1])
    num_ispn2_caudate <- data.frame(num_iSPN2_caudate = dim(x[grep(x = x_caudate$sub_type_4, pattern = "iSPN_2"), ])[1])
    num_cells_accumbens <- data.frame(num_cells_accumbens = dim(x_accumbens)[1])
    num_dspn1_accumbens <- data.frame(num_dSPN1_accumbens = dim(x[grep(x = x_accumbens$sub_type_4, pattern = "dSPN_1"), ])[1])
    num_dspn2_accumbens <- data.frame(num_dSPN2_accumbens = dim(x[grep(x = x_accumbens$sub_type_4, pattern = "dSPN_2"), ])[1])
    num_ispn1_accumbens <- data.frame(num_iSPN1_accumbens = dim(x[grep(x = x_accumbens$sub_type_4, pattern = "iSPN_1"), ])[1])
    num_ispn2_accumbens <- data.frame(num_iSPN2_accumbens = dim(x[grep(x = x_accumbens$sub_type_4, pattern = "iSPN_2"), ])[1])
    return(data.frame(num_cells_caudate, num_dspn1_caudate, num_dspn2_caudate, num_ispn1_caudate, num_ispn2_caudate, num_cells_accumbens, num_dspn1_accumbens, num_dspn2_accumbens, num_ispn1_accumbens, num_ispn2_accumbens))
})

num_cells_Paryani <- do.call(rbind, tmp)
#num_cells_Paryani
```


```R
#subset SPNs
Paryani_MSN <- subset(x = Paryani, subset = sub_type_4 == "dSPN_1" | sub_type_4 == "dSPN_2" | sub_type_4 == "iSPN_1" | sub_type_4 == "iSPN_2")
```

# Paryani_Caudate_only


```R
# Paryani_MSN_Caudate <- subset(x = Paryani_MSN, subset = Region == "Caudate")
```


```R
# #Run SCTransform
# Paryani_MSN_Caudate <- SCTransform(Paryani_MSN_Caudate, vars.to.regress = c("nCount_RNA", "nFeature_RNA"))
# #Run PCA and UMAP
# Paryani_MSN_Caudate <- RunPCA(Paryani_MSN_Caudate, features = VariableFeatures(object = Paryani_MSN_Caudate), npcs = 50, reduction.name = "pca")
# ElbowPlot(Paryani_MSN, reduction = "pca", ndims = 50)
# num_dims <- 40
# Paryani_MSN_Caudate <- RunUMAP(Paryani_MSN_Caudate, dims = 1:num_dims, reduction = "pca", reduction.name = "umap", reduction.key = "umap")
# DimPlot(Paryani_MSN_Caudate, group.by = "Batch", reduction = "umap", pt.size = 3)
# Paryani_MSN_Caudate <- RunHarmony(object = Paryani_MSN_Caudate, reduction.use = "pca", group.by.vars = "Batch", 
#                            reduction.save = "harmonyPca", assay = "SCT", verbose = FALSE, normalization.method = "SCT")
# Paryani_MSN_Caudate <- RunUMAP(Paryani_MSN_Caudate, dims = 1:num_dims, reduction = "harmonyPca", reduction.name = "harmony", reduction.key = "harmony")
```


```R
#save object
# saveRDS(object = Paryani_MSN_Caudate, file = "Paryani_MSN_Caudate.rds")
```


```R
#load dataset
Paryani_MSN_Caudate <- readRDS("Paryani_MSN_Caudate.rds")
```


```R
#do plots
FeaturePlot(Paryani_MSN_Caudate, reduction = "harmony", features = "nCount_RNA", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Paryani_MSN_Caudate_NCOUNTS.pdf", width = 8, height = 8)
FeaturePlot(Paryani_MSN_Caudate, reduction = "harmony", features = "nFeature_RNA", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Paryani_MSN_Caudate_NFEATURES.pdf", width = 8, height = 8)
Paryani_MSN_Caudate@meta.data$GERM_EXPCAGLENGTH <- as.numeric(gsub(x = Paryani_MSN_Caudate@meta.data$CAG, pattern = "-.*", replacement = ""))
FeaturePlot(Paryani_MSN_Caudate, reduction = "harmony", features = "GERM_EXPCAGLENGTH", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Paryani_MSN_Caudate_GERM_EXPCAGLENGTH.pdf", width = 8, height = 8)
DimPlot(Paryani_MSN_Caudate, group.by = "sub_type_4", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Caudate_CELLTYPE.pdf", width = 8, height = 8)
DimPlot(Paryani_MSN_Caudate, group.by = "CAG", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Caudate_CAG.pdf", width = 8, height = 8)
DimPlot(Paryani_MSN_Caudate, group.by = "Batch", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Caudate_BATCH.pdf", width = 8, height = 8)
DimPlot(Paryani_MSN_Caudate, group.by = "Region", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Caudate_REGION.pdf", width = 8, height = 8)
DimPlot(Paryani_MSN_Caudate, group.by = "Donor", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Caudate_DONOR.pdf", width = 8, height = 8)
Paryani_MSN_Caudate@meta.data$Grade <- factor(Paryani_MSN_Caudate@meta.data$Grade, levels = c("HDJ", "HD4", "HD3", "HD2", "HD1", "Control"))
DimPlot(Paryani_MSN_Caudate, group.by = "Grade", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Caudate_GRADE.pdf", width = 8, height = 8)
```


```R
#extract phase C genes in the model that are also expressed in the dataset
phaseC_genes_Paryani_MSN_Caudate <- intersect(colnames(phase_all), rownames(Paryani_MSN_Caudate))
phaseC_genes_notExpParyani_MSN_Caudate_names <- setdiff(colnames(phase_all), rownames(Paryani_MSN_Caudate))
phaseC_genes_notExpParyani_MSN_Caudate <- matrix(data = 0, nrow = length(phaseC_genes_notExpParyani_MSN_Caudate_names), ncol = dim(Paryani_MSN_Caudate@assays$SCT$data)[2])
rownames(phaseC_genes_notExpParyani_MSN_Caudate) <- phaseC_genes_notExpParyani_MSN_Caudate_names
colnames(phaseC_genes_notExpParyani_MSN_Caudate) <- colnames(Paryani_MSN_Caudate@assays$SCT$data)
phase_counts_Paryani_MSN_Caudate <- Paryani_MSN_Caudate@assays$SCT$data[phaseC_genes_Paryani_MSN_Caudate, ]
phase_counts_Paryani_MSN_Caudate <- rbind(phase_counts_Paryani_MSN_Caudate, phaseC_genes_notExpParyani_MSN_Caudate)
```


```R
phase_counts_Paryani_MSN_Caudate <- phase_counts_Paryani_MSN_Caudate[colnames(phase_test), ]
phase_test_Paryani_MSN_Caudate <- data.frame(t(phase_counts_Paryani_MSN_Caudate))
#head(phase_test_Paryani_MSN_Caudate)
```


```R
#predict phase
class_thr <- perf_training_test_noNA$thr
predicted_phase_test_Paryani_MSN_Caudate_wprob <- predict(model_phase, newx = as.matrix(phase_test_Paryani_MSN_Caudate), s = lambda_val, type = "response") #glmnet
predicted_phase_test_Paryani_MSN_Caudate <- rep("A-B", dim(predicted_phase_test_Paryani_MSN_Caudate_wprob)[1])  
predicted_phase_test_Paryani_MSN_Caudate[which(predicted_phase_test_Paryani_MSN_Caudate_wprob > class_thr)] <- "C-D-E" #glmnet
Paryani_MSN_Caudate[["PREDICTED_PHASE"]] <- NA
Paryani_MSN_Caudate@meta.data[rownames(predicted_phase_test_Paryani_MSN_Caudate_wprob), "PREDICTED_PHASE"] <- predicted_phase_test_Paryani_MSN_Caudate
Paryani_MSN_Caudate[["PROB_PHASE_CDE"]] <- NA
Paryani_MSN_Caudate@meta.data[rownames(predicted_phase_test_Paryani_MSN_Caudate_wprob), "PROB_PHASE_CDE"] <- predicted_phase_test_Paryani_MSN_Caudate_wprob
```


```R
table(predicted_phase_test_Paryani_MSN_Caudate)
summary(predicted_phase_test_Paryani_MSN_Caudate_wprob)
colnames(Paryani_MSN_Caudate@meta.data)
```


```R
table(Paryani_MSN_Caudate@meta.data$Grade)
```


```R
#do plots
DimPlot(Paryani_MSN_Caudate, group.by = "PREDICTED_PHASE", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Caudate_PREDPHASE.pdf", width = 8, height = 8)
DimPlot(Paryani_MSN_Caudate, group.by = "PREDICTED_PHASE", reduction = "harmony", split.by = "Condition", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Caudate_PREDPHASE_splitByCondition.pdf", width = 8, height = 8)
FeaturePlot(Paryani_MSN_Caudate, features = "PROB_PHASE_CDE", cols = c("blue", "red"), pt.size = 1, reduction = "harmony")
ggsave("UMAP_Paryani_MSN_Caudate_PROBPREDPHASE.pdf", width = 8, height = 8)
FeaturePlot(Paryani_MSN_Caudate, features = "PROB_PHASE_CDE", cols = c("blue", "red"), pt.size = 1, split.by = "Condition", reduction = "harmony")
ggsave("UMAP_Paryani_MSN_Caudate_PROBPREDPHASE_splitByCondition.pdf", width = 8, height = 8)
DimPlot(Paryani_MSN_Caudate, group.by = "Donor", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Caudate_DONOR.pdf", width = 8, height = 8)
Paryani_MSN_Caudate@meta.data$Grade <- factor(Paryani_MSN_Caudate@meta.data$Grade, levels = c("Control", "HD1", "HD2", "HD3", "HD4", "HDJ"))
DimPlot(Paryani_MSN_Caudate, group.by = "Grade", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Caudate_GRADE.pdf", width = 8, height = 8)
DimPlot(Paryani_MSN_Caudate, group.by = "Grade", reduction = "harmony", split.by = "Condition", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Caudate_GRADE_splitByCondition.pdf", width = 8, height = 8)
```


```R
#plot the fraction of SPN in C-D-E phase vs the number of CAG in the germline
fract_CDE_tmp_Caudate <- lapply(split(Paryani_MSN_Caudate@meta.data, Paryani_MSN_Caudate@meta.data$CAG), function(x) {
    num_SPN <- dim(x)[1]
    num_CDE <- length(which(x$PREDICTED_PHASE == "C-D-E"))
    fraction_CDE <- num_CDE/num_SPN
    num_CAG <- gsub(x = x$CAG[1], pattern = "-.*", replacement = "")
    #cat(sprintf("Fraction MSNs in C-D-E phase: %.2f\n", fraction_CDE*100))
    return(rbind(fraction_CDE, num_SPN, num_CAG))
})
fract_CDE_Caudate <- t(do.call(cbind, fract_CDE_tmp_Caudate))
fract_CDE_Paryani_Caudate <- data.frame(num_CAG_germ = as.numeric(fract_CDE_Caudate[, 3]), fract_CDE = as.numeric(fract_CDE_Caudate[, 1]), num_SPN = as.numeric(fract_CDE_Caudate[, 2]))
#cor_numCAGGerm_fractCDE_Caudate <- cor(x = fract_CDE_Caudate$num_CAG_germ, y = fract_CDE_Caudate$fract_CDE, method = "pearson", use="complete.obs")

cor_numCAGGerm_fractCDE_Paryani_Caudate <- weighted.cor(x = fract_CDE_Paryani_Caudate$num_CAG_germ, y = fract_CDE_Paryani_Caudate$fract_CDE, weights = fract_CDE_Paryani_Caudate$num_SPN, method = "pearson", na.rm = TRUE)

print(cor_numCAGGerm_fractCDE_Paryani_Caudate)
#plot(x = fract_CDE_Caudate$num_CAG_germ, y = fract_CDE_Caudate$fract_CDE)

ggplot(fract_CDE_Paryani_Caudate, aes(x=num_CAG_germ, y=100*fract_CDE)) + 
  geom_point(col = "blue", aes(size = num_SPN))+
  geom_smooth(method=lm, se = TRUE, mapping = aes(weight = num_SPN)) +
  theme_classic() +
  xlab("Num. CAG germline") +
  ylab("Fraction cells in C-D-E phase (%)")
  #+ ggtitle(sprintf("Sp. corr = %.2f", cor_numCAGGerm_fractCDE_Caudate))
  #+ xlim(c(36, 70)) + ylim(c(0, 1))
ggsave("Paryani_NumCAGGermline_FractCDE_Caudate.pdf", width = 8, height = 8)
```


```R
wtd.cor(x = fract_CDE_Paryani_Caudate$num_CAG_germ, y = fract_CDE_Paryani_Caudate$fract_CDE, weight = fract_CDE_Paryani_Caudate$num_SPN)
```


```R
Calculate_confInt(r = wtd.cor(x = fract_CDE_Paryani_Caudate$num_CAG_germ, y = fract_CDE_Paryani_Caudate$fract_CDE, weight = fract_CDE_Paryani_Caudate$num_SPN)[, "correlation"], w = fract_CDE_Paryani_Caudate$num_SPN)
```


```R
#create metadata for plotting
Paryani_Caudate_Grade_CDE_tmp <- lapply(split(Paryani_MSN_Caudate@meta.data, Paryani_MSN_Caudate@meta.data$Donor), function(x) {
    GRADE <- factor(unique(x$Grade), levels = c("Control", "HD1", "HD2", "HD3", "HD4", "HDJ"))
    SN <- unique(x$Donor)
    #NUM_MSN <- length(x$PREDICTED_PHASE)
    tmp <- lapply(split(x, x$sub_type_4), function(y) {
        NUM_MSN_CT <- length(y$PREDICTED_PHASE)
        NUM_MSN_CT_SQ <- NUM_MSN_CT**2
        CELL_TYPE <- y$sub_type_4
        if (length(unique(y$PREDICTED_PHASE)) > 1) {
          val <- 100*table(y$PREDICTED_PHASE)[2]/(table(y$PREDICTED_PHASE)[1] + table(y$PREDICTED_PHASE)[2])
        } else {
          val <- 0
        }
        names(val) <- "Fract_CDE";
        Fract_CDE <- rep(val, NUM_MSN_CT_SQ);
    return(data.frame(NUM_MSN_CT_SQ, CELL_TYPE, Fract_CDE))})
    df <- do.call(rbind, tmp)
    #print(df)
    #NUM_MSN <- length(x$PREDICTED_PHASE)
    #Fract_CDE <- 100*length(which(x$PREDICTED_PHASE == "C-D-E"))/length(x$PREDICTED_PHASE)
    #CELL_TYPE <- x$sub_type_4
    #return(df)
    #return(data.frame(SAMPLE = rep(SN, df$NUM_MSN_CT_SQ), GRADE = rep(GRADE, df$NUM_MSN_CT_SQ), CELL_TYPE = df$CELL_TYPE, Fract_CDE = df$Fract_CDE, NUM_MSN_SQ=df$NUM_MSN_CT_SQ))
    return(data.frame(SAMPLE = SN, GRADE = GRADE, CELL_TYPE = df$CELL_TYPE, Fract_CDE = df$Fract_CDE, NUM_MSN_SQ=df$NUM_MSN_CT_SQ))
})
Paryani_Caudate_Grade_CDE <- do.call(rbind, Paryani_Caudate_Grade_CDE_tmp)
Paryani_Caudate_Grade_CDE <- Paryani_Caudate_Grade_CDE[order(Paryani_Caudate_Grade_CDE$GRADE), ]
Paryani_Caudate_Grade_CDE$SAMPLE <- factor(Paryani_Caudate_Grade_CDE$SAMPLE, levels = unique(Paryani_Caudate_Grade_CDE$SAMPLE))
head(Paryani_Caudate_Grade_CDE)
```


```R
uni_Paryani_Caudate_Grade_CDE <- unique(Paryani_Caudate_Grade_CDE)
uni_Paryani_Caudate_Grade_CDE$NUM_SPN <- sqrt(uni_Paryani_Caudate_Grade_CDE$NUM_MSN_SQ)
#uni_Paryani_Caudate_Grade_CDE
```


```R
lapply(split(Paryani_MSN_Caudate@meta.data[rownames(predicted_phase_test_Paryani_MSN_Caudate_wprob), "PREDICTED_PHASE"], Paryani_MSN_Caudate@meta.data[rownames(predicted_phase_test_Paryani_MSN_Caudate_wprob), "Donor"]), function(x) table(x))
lapply(split(Paryani_MSN_Caudate@meta.data[rownames(predicted_phase_test_Paryani_MSN_Caudate_wprob), "PREDICTED_PHASE"], Paryani_MSN_Caudate@meta.data[rownames(predicted_phase_test_Paryani_MSN_Caudate_wprob), "Grade"]), function(x) table(x))
```


```R
#add metadata
Fract_SPN_Caudate <- c(num_cells_Paryani["T-4812", "num_dSPN1_caudate"]/num_cells_Paryani["T-4812", "num_cells_caudate"],
                       num_cells_Paryani["T-4812", "num_dSPN2_caudate"]/num_cells_Paryani["T-4812", "num_cells_caudate"],
                       num_cells_Paryani["T-4812", "num_iSPN1_caudate"]/num_cells_Paryani["T-4812", "num_cells_caudate"],
                       num_cells_Paryani["T-4812", "num_iSPN2_caudate"]/num_cells_Paryani["T-4812", "num_cells_caudate"],
                       
                       num_cells_Paryani["T-4915", "num_dSPN1_caudate"]/num_cells_Paryani["T-4915", "num_cells_caudate"],
                       num_cells_Paryani["T-4915", "num_dSPN2_caudate"]/num_cells_Paryani["T-4915", "num_cells_caudate"],
                       num_cells_Paryani["T-4915", "num_iSPN1_caudate"]/num_cells_Paryani["T-4915", "num_cells_caudate"],
                       #num_cells_Paryani["T-4915", "num_iSPN2_caudate"]/num_cells_Paryani["T-4915", "num_cells_caudate"],
                       
                       num_cells_Paryani["T-5227", "num_dSPN1_caudate"]/num_cells_Paryani["T-5227", "num_cells_caudate"],
                       num_cells_Paryani["T-5227", "num_dSPN2_caudate"]/num_cells_Paryani["T-5227", "num_cells_caudate"],
                       num_cells_Paryani["T-5227", "num_iSPN1_caudate"]/num_cells_Paryani["T-5227", "num_cells_caudate"],
                       num_cells_Paryani["T-5227", "num_iSPN2_caudate"]/num_cells_Paryani["T-5227", "num_cells_caudate"],
                       
                       num_cells_Paryani["T-5596", "num_dSPN1_caudate"]/num_cells_Paryani["T-5596", "num_cells_caudate"],                       
                       num_cells_Paryani["T-5596", "num_dSPN2_caudate"]/num_cells_Paryani["T-5596", "num_cells_caudate"],
                       num_cells_Paryani["T-5596", "num_iSPN1_caudate"]/num_cells_Paryani["T-5596", "num_cells_caudate"],
                       #num_cells_Paryani["T-5596", "num_iSPN2_caudate"]/num_cells_Paryani["T-5596", "num_cells_caudate"],
                       
                       num_cells_Paryani["T-5700", "num_dSPN1_caudate"]/num_cells_Paryani["T-5700", "num_cells_caudate"],
                       num_cells_Paryani["T-5700", "num_dSPN2_caudate"]/num_cells_Paryani["T-5700", "num_cells_caudate"],
                       num_cells_Paryani["T-5700", "num_iSPN1_caudate"]/num_cells_Paryani["T-5700", "num_cells_caudate"],
                       #num_cells_Paryani["T-5700", "num_iSPN2_caudate"]/num_cells_Paryani["T-5700", "num_cells_caudate"],
                       
                       num_cells_Paryani["T-4285", "num_dSPN1_caudate"]/num_cells_Paryani["T-4285", "num_cells_caudate"],
                       num_cells_Paryani["T-4285", "num_dSPN2_caudate"]/num_cells_Paryani["T-4285", "num_cells_caudate"],
                       num_cells_Paryani["T-4285", "num_iSPN1_caudate"]/num_cells_Paryani["T-4285", "num_cells_caudate"],
                       num_cells_Paryani["T-4285", "num_iSPN2_caudate"]/num_cells_Paryani["T-4285", "num_cells_caudate"],
                       
                       num_cells_Paryani["T-5266", "num_dSPN1_caudate"]/num_cells_Paryani["T-5266", "num_cells_caudate"],
                       #num_cells_Paryani["T-5266", "num_dSPN2_caudate"]/num_cells_Paryani["T-5266", "num_cells_caudate"],
                       num_cells_Paryani["T-5266", "num_iSPN1_caudate"]/num_cells_Paryani["T-5266", "num_cells_caudate"],
                       num_cells_Paryani["T-5266", "num_iSPN2_caudate"]/num_cells_Paryani["T-5266", "num_cells_caudate"],
                      
                       num_cells_Paryani["T-5798", "num_dSPN1_caudate"]/num_cells_Paryani["T-5798", "num_cells_caudate"],
                       num_cells_Paryani["T-5798", "num_dSPN2_caudate"]/num_cells_Paryani["T-5798", "num_cells_caudate"],
                       num_cells_Paryani["T-5798", "num_iSPN1_caudate"]/num_cells_Paryani["T-5798", "num_cells_caudate"],
                       num_cells_Paryani["T-5798", "num_iSPN2_caudate"]/num_cells_Paryani["T-5798", "num_cells_caudate"],
                       
                       num_cells_Paryani["T-5263", "num_dSPN1_caudate"]/num_cells_Paryani["T-5263", "num_cells_caudate"],
                       num_cells_Paryani["T-5263", "num_dSPN2_caudate"]/num_cells_Paryani["T-5263", "num_cells_caudate"],
                       num_cells_Paryani["T-5263", "num_iSPN1_caudate"]/num_cells_Paryani["T-5263", "num_cells_caudate"],
                       num_cells_Paryani["T-5263", "num_iSPN2_caudate"]/num_cells_Paryani["T-5263", "num_cells_caudate"],
                       
                       num_cells_Paryani["T-5286", "num_dSPN1_caudate"]/num_cells_Paryani["T-5286", "num_cells_caudate"],
                       num_cells_Paryani["T-5286", "num_dSPN2_caudate"]/num_cells_Paryani["T-5286", "num_cells_caudate"],
                       num_cells_Paryani["T-5286", "num_iSPN1_caudate"]/num_cells_Paryani["T-5286", "num_cells_caudate"],
                       num_cells_Paryani["T-5286", "num_iSPN2_caudate"]/num_cells_Paryani["T-5286", "num_cells_caudate"],
                      
                       num_cells_Paryani["T-5693", "num_dSPN1_caudate"]/num_cells_Paryani["T-5693", "num_cells_caudate"],
                       num_cells_Paryani["T-5693", "num_dSPN2_caudate"]/num_cells_Paryani["T-5693", "num_cells_caudate"],
                       num_cells_Paryani["T-5693", "num_iSPN1_caudate"]/num_cells_Paryani["T-5693", "num_cells_caudate"],
                       num_cells_Paryani["T-5693", "num_iSPN2_caudate"]/num_cells_Paryani["T-5693", "num_cells_caudate"],
                      
                       num_cells_Paryani["T-4161", "num_dSPN1_caudate"]/num_cells_Paryani["T-4161", "num_cells_caudate"],
                       num_cells_Paryani["T-4161", "num_dSPN2_caudate"]/num_cells_Paryani["T-4161", "num_cells_caudate"],
                       num_cells_Paryani["T-4161", "num_iSPN1_caudate"]/num_cells_Paryani["T-4161", "num_cells_caudate"],
                       num_cells_Paryani["T-4161", "num_iSPN2_caudate"]/num_cells_Paryani["T-4161", "num_cells_caudate"],
                      
                       num_cells_Paryani["T-5493", "num_dSPN1_caudate"]/num_cells_Paryani["T-5493", "num_cells_caudate"],
                       num_cells_Paryani["T-5493", "num_dSPN2_caudate"]/num_cells_Paryani["T-5493", "num_cells_caudate"],
                       num_cells_Paryani["T-5493", "num_iSPN1_caudate"]/num_cells_Paryani["T-5493", "num_cells_caudate"],
                       num_cells_Paryani["T-5493", "num_iSPN2_caudate"]/num_cells_Paryani["T-5493", "num_cells_caudate"],
                      
                       num_cells_Paryani["T-5575", "num_dSPN1_caudate"]/num_cells_Paryani["T-5575", "num_cells_caudate"],
                       num_cells_Paryani["T-5575", "num_dSPN2_caudate"]/num_cells_Paryani["T-5575", "num_cells_caudate"],
                       num_cells_Paryani["T-5575", "num_iSPN1_caudate"]/num_cells_Paryani["T-5575", "num_cells_caudate"],
                       num_cells_Paryani["T-5575", "num_iSPN2_caudate"]/num_cells_Paryani["T-5575", "num_cells_caudate"],
                      
                       num_cells_Paryani["T-4273", "num_dSPN1_caudate"]/num_cells_Paryani["T-4273", "num_cells_caudate"],
                       #num_cells_Paryani["T-4273", "num_dSPN2_caudate"]/num_cells_Paryani["T-4273", "num_cells_caudate"],
                       num_cells_Paryani["T-4273", "num_iSPN1_caudate"]/num_cells_Paryani["T-4273", "num_cells_caudate"],
                       num_cells_Paryani["T-4273", "num_iSPN2_caudate"]/num_cells_Paryani["T-4273", "num_cells_caudate"],
                      
                       num_cells_Paryani["T-4638", "num_dSPN1_caudate"]/num_cells_Paryani["T-4638", "num_cells_caudate"],
                       #num_cells_Paryani["T-4638", "num_dSPN2_caudate"]/num_cells_Paryani["T-4638", "num_cells_caudate"],
                       num_cells_Paryani["T-4638", "num_iSPN1_caudate"]/num_cells_Paryani["T-4638", "num_cells_caudate"],
                       num_cells_Paryani["T-4638", "num_iSPN2_caudate"]/num_cells_Paryani["T-4638", "num_cells_caudate"],
                      
                       num_cells_Paryani["T-5534", "num_dSPN1_caudate"]/num_cells_Paryani["T-5534", "num_cells_caudate"],
                       num_cells_Paryani["T-5534", "num_dSPN2_caudate"]/num_cells_Paryani["T-5534", "num_cells_caudate"],
                       num_cells_Paryani["T-5534", "num_iSPN1_caudate"]/num_cells_Paryani["T-5534", "num_cells_caudate"],
                       num_cells_Paryani["T-5534", "num_iSPN2_caudate"]/num_cells_Paryani["T-5534", "num_cells_caudate"],
                      
                       num_cells_Paryani["T-5714", "num_dSPN1_caudate"]/num_cells_Paryani["T-5714", "num_cells_caudate"],
                       num_cells_Paryani["T-5714", "num_dSPN2_caudate"]/num_cells_Paryani["T-5714", "num_cells_caudate"],
                       num_cells_Paryani["T-5714", "num_iSPN1_caudate"]/num_cells_Paryani["T-5714", "num_cells_caudate"],
                       num_cells_Paryani["T-5714", "num_iSPN2_caudate"]/num_cells_Paryani["T-5714", "num_cells_caudate"])
#Fract_SPN_Caudate
uni_Paryani_Caudate_Grade_CDE$Fract_SPN <- Fract_SPN_Caudate
#uni_Paryani_Caudate_Grade_CDE
```


```R
# ggplot(Paryani_Caudate_Grade_CDE, aes(x = SAMPLE, y = Fract_CDE, color = CELL_TYPE)) +
# geom_boxplot(varwidth = TRUE, lwd = 5, position = position_dodge2(preserve = "single")) +
# scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
# theme_classic() +
# xlab("SAMPLE") +
# ylim(c(0, 100)) +
# ylab("Fraction cells in C-D-E phase (%)") 
# ggsave("Paryani_Caudate_CELLTYPE_vs_fraction_SPN_CDE.pdf", width = 8, height = 8)

uni_Paryani_Caudate_Grade_CDE_HD <- uni_Paryani_Caudate_Grade_CDE[which(uni_Paryani_Caudate_Grade_CDE$GRADE != "Control"), ]
#uni_Paryani_Caudate_Grade_CDE_CTRL <- uni_Paryani_Caudate_Grade_CDE[which(uni_Paryani_Caudate_Grade_CDE$GRADE == "Control"), ]

#plot fraction of SPNs in C-D-E phase by HD grade
ggplot(uni_Paryani_Caudate_Grade_CDE_HD, aes(x = GRADE, y = Fract_CDE, color = CELL_TYPE)) +
geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.1, height = 0)) +
scale_color_manual(values =c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
theme_classic() +
xlab("HD GRADE") +
#ylim(c(0, 100)) +
ylab("Fraction cells in C-D-E phase (%)") 
ggsave("Paryani_Caudate_HD_GRADE_vs_fraction_SPN_CDE.pdf", width = 6, height = 6)
```


```R
unique(uni_Paryani_Caudate_Grade_CDE[, c("SAMPLE", "GRADE")])
```


```R
# ggplot(uni_Paryani_Caudate_Grade_CDE[which(uni_Paryani_Caudate_Grade_CDE$NUM_SPN > 30), ], aes(x = GRADE, y = Fract_CDE, color = CELL_TYPE)) +
# geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.1, height = 0)) +
# scale_color_manual(values =c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
# theme_classic() +
# xlab("HD GRADE") +
# #ylim(c(0, 100)) +
# ylab("Fraction cells in C-D-E phase (%)") 
# ggsave("Paryani_Caudate_GRADE_vs_fraction_SPN_CDE_gt30.pdf", width = 6, height = 6)
```


```R
#subset HD
Paryani_HD_MSN_Caudate <- subset(x = Paryani_MSN_Caudate, subset = Condition == "HD")
```

# Paryani_Accumbens_only


```R
# Paryani_MSN_Accumbens <- subset(x = Paryani_MSN, subset = Region == "Accumbens")
```


```R
# #Run SCTransform for HD donors
# Paryani_MSN_Accumbens <- SCTransform(Paryani_MSN_Accumbens, vars.to.regress = c("nCount_RNA", "nFeature_RNA"))
# #Run PCA and UMAP for HD donors
# Paryani_MSN_Accumbens <- RunPCA(Paryani_MSN_Accumbens, features = VariableFeatures(object = Paryani_MSN_Accumbens), npcs = 50, reduction.name = "pca")
# ElbowPlot(Paryani_MSN, reduction = "pca", ndims = 50)
# num_dims <- 40
# Paryani_MSN_Accumbens <- RunUMAP(Paryani_MSN_Accumbens, dims = 1:num_dims, reduction = "pca", reduction.name = "umap", reduction.key = "umap")
# DimPlot(Paryani_MSN_Accumbens, group.by = "Batch", reduction = "umap", pt.size = 1)
# Paryani_MSN_Accumbens <- RunHarmony(object = Paryani_MSN_Accumbens, reduction.use = "pca", group.by.vars = "Batch", 
#                            reduction.save = "harmonyPca", assay = "SCT", verbose = FALSE, normalization.method = "SCT")
# Paryani_MSN_Accumbens <- RunUMAP(Paryani_MSN_Accumbens, dims = 1:num_dims, reduction = "harmonyPca", reduction.name = "harmony", reduction.key = "harmony")
```


```R
# saveRDS(object = Paryani_MSN_Accumbens, file = "Paryani_MSN_Accumbens.rds")
```


```R
Paryani_MSN_Accumbens <- readRDS("Paryani_MSN_Accumbens.rds")
```


```R
#do plots
FeaturePlot(Paryani_MSN_Accumbens, reduction = "harmony", features = "nCount_RNA", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Paryani_MSN_Accumbens_NCOUNTS.pdf", width = 8, height = 8)
FeaturePlot(Paryani_MSN_Accumbens, reduction = "harmony", features = "nFeature_RNA", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Paryani_MSN_Accumbens_NFEATURES.pdf", width = 8, height = 8)
Paryani_MSN_Accumbens@meta.data$GERM_EXPCAGLENGTH <- as.numeric(gsub(x = Paryani_MSN_Accumbens@meta.data$CAG, pattern = "-.*", replacement = ""))
FeaturePlot(Paryani_MSN_Accumbens, reduction = "harmony", features = "GERM_EXPCAGLENGTH", cols = c("blue", "red"), pt.size = 1)
ggsave("UMAP_Paryani_MSN_Accumbens_GERM_EXPCAGLENGTH.pdf", width = 8, height = 8)
DimPlot(Paryani_MSN_Accumbens, group.by = "sub_type_4", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Accumbens_CELLTYPE.pdf", width = 8, height = 8)
DimPlot(Paryani_MSN_Accumbens, group.by = "CAG", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Accumbens_CAG.pdf", width = 8, height = 8)
DimPlot(Paryani_MSN_Accumbens, group.by = "Batch", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Accumbens_BATCH.pdf", width = 8, height = 8)
DimPlot(Paryani_MSN_Accumbens, group.by = "Region", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Accumbens_REGION.pdf", width = 8, height = 8)
DimPlot(Paryani_MSN_Accumbens, group.by = "Donor", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Accumbens_DONOR.pdf", width = 8, height = 8)
Paryani_MSN_Accumbens@meta.data$Grade <- factor(Paryani_MSN_Accumbens@meta.data$Grade, levels = c("Control", "HDJ", "HD4", "HD3", "HD2", "HD1", "Control"))
DimPlot(Paryani_MSN_Accumbens, group.by = "Grade", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Accumbens_GRADE.pdf", width = 8, height = 8)
```


```R
#extract phase C genes in the model that are also expressed in the dataset
phaseC_genes_Paryani_MSN_Accumbens <- intersect(colnames(phase_all), rownames(Paryani_MSN_Accumbens))
phaseC_genes_notExpParyani_MSN_Accumbens_names <- setdiff(colnames(phase_all), rownames(Paryani_MSN_Accumbens))
phaseC_genes_notExpParyani_MSN_Accumbens <- matrix(data = 0, nrow = length(phaseC_genes_notExpParyani_MSN_Accumbens_names), ncol = dim(Paryani_MSN_Accumbens@assays$SCT$data)[2])
rownames(phaseC_genes_notExpParyani_MSN_Accumbens) <- phaseC_genes_notExpParyani_MSN_Accumbens_names
colnames(phaseC_genes_notExpParyani_MSN_Accumbens) <- colnames(Paryani_MSN_Accumbens@assays$SCT$data)
phase_counts_Paryani_MSN_Accumbens <- Paryani_MSN_Accumbens@assays$SCT$data[phaseC_genes_Paryani_MSN_Accumbens, ]
phase_counts_Paryani_MSN_Accumbens <- rbind(phase_counts_Paryani_MSN_Accumbens, phaseC_genes_notExpParyani_MSN_Accumbens)
```


```R
phase_counts_Paryani_MSN_Accumbens <- phase_counts_Paryani_MSN_Accumbens[colnames(phase_all), ]
phase_test_Paryani_MSN_Accumbens <- data.frame(t(phase_counts_Paryani_MSN_Accumbens))
#head(phase_test_Paryani_HD_MSN_Accumbens)
```


```R
#predict phase
class_thr <- perf_training_test_noNA$thr
predicted_phase_test_Paryani_MSN_Accumbens_wprob <- predict(model_phase, newx = as.matrix(phase_test_Paryani_MSN_Accumbens), s = lambda_val, type = "response") #glmnet
predicted_phase_test_Paryani_MSN_Accumbens <- rep("A-B", dim(predicted_phase_test_Paryani_MSN_Accumbens_wprob)[1])  
predicted_phase_test_Paryani_MSN_Accumbens[which(predicted_phase_test_Paryani_MSN_Accumbens_wprob > class_thr)] <- "C-D-E" #glmnet
Paryani_MSN_Accumbens[["PREDICTED_PHASE"]] <- NA
Paryani_MSN_Accumbens@meta.data[rownames(predicted_phase_test_Paryani_MSN_Accumbens_wprob), "PREDICTED_PHASE"] <- predicted_phase_test_Paryani_MSN_Accumbens
Paryani_MSN_Accumbens[["PROB_PHASE_CDE"]] <- NA
Paryani_MSN_Accumbens@meta.data[rownames(predicted_phase_test_Paryani_MSN_Accumbens_wprob), "PROB_PHASE_CDE"] <- predicted_phase_test_Paryani_MSN_Accumbens_wprob
```


```R
table(predicted_phase_test_Paryani_MSN_Accumbens)
summary(predicted_phase_test_Paryani_MSN_Accumbens_wprob)
```


```R
#do plots
DimPlot(Paryani_MSN_Accumbens, group.by = "PREDICTED_PHASE", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Accumbens_PREDPHASE.pdf", width = 8, height = 8)
DimPlot(Paryani_MSN_Accumbens, group.by = "PREDICTED_PHASE", reduction = "harmony", split.by = "Condition", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Accumbens_PREDPHASE_splitByCondition.pdf", width = 8, height = 8)
FeaturePlot(Paryani_MSN_Accumbens, features = "PROB_PHASE_CDE", cols = c("blue", "red"), pt.size = 1, reduction = "harmony")
ggsave("UMAP_Paryani_MSN_Accumbens_PROBPREDPHASE.pdf", width = 8, height = 8)
FeaturePlot(Paryani_MSN_Accumbens, features = "PROB_PHASE_CDE", cols = c("blue", "red"), pt.size = 1, split.by = "Condition", reduction = "harmony")
ggsave("UMAP_Paryani_MSN_Accumbens_PROBPREDPHASE_splitByCondition.pdf", width = 8, height = 8)
DimPlot(Paryani_MSN_Accumbens, group.by = "Donor", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Accumbens_DONOR.pdf", width = 8, height = 8)
Paryani_MSN_Accumbens@meta.data$Grade <- factor(Paryani_MSN_Accumbens@meta.data$Grade, levels = c("Control", "HD1", "HD2", "HD3", "HD4", "HDJ"))
DimPlot(Paryani_MSN_Accumbens, group.by = "Grade", reduction = "harmony", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Accumbens_GRADE.pdf", width = 8, height = 8)
DimPlot(Paryani_MSN_Accumbens, group.by = "Grade", reduction = "harmony", split.by = "Condition", pt.size = 1)
ggsave("UMAP_Paryani_MSN_Accumbens_GRADE_splitByCondition.pdf", width = 8, height = 8)
```


```R
#plot the fraction of SPN in C-D-E phase vs the number of CAG in the germline
fract_CDE_tmp_Accumbens <- lapply(split(Paryani_MSN_Accumbens@meta.data, Paryani_MSN_Accumbens@meta.data$CAG), function(x) {
    num_SPN <- dim(x)[1]
    num_CDE <- length(which(x$PREDICTED_PHASE == "C-D-E"))
    fraction_CDE <- num_CDE/num_SPN
    num_CAG <- gsub(x = x$CAG[1], pattern = "-.*", replacement = "")
    #cat(sprintf("Fraction MSNs in C-D-E phase: %.2f\n", fraction_CDE*100))
    return(rbind(fraction_CDE, num_SPN, num_CAG))
})
fract_CDE_Accumbens <- t(do.call(cbind, fract_CDE_tmp_Accumbens))
fract_CDE_Paryani_Accumbens <- data.frame(num_CAG_germ = as.numeric(fract_CDE_Accumbens[, 3]), fract_CDE = as.numeric(fract_CDE_Accumbens[, 1]), num_SPN = as.numeric(fract_CDE_Accumbens[, 2]))
#cor_numCAGGerm_fractCDE_Accumbens <- cor(x = fract_CDE_Accumbens$num_CAG_germ, y = fract_CDE_Accumbens$fract_CDE, method = "pearson", use="complete.obs")

cor_numCAGGerm_fractCDE_Paryani_Accumbens <- weighted.cor(x = fract_CDE_Paryani_Accumbens$num_CAG_germ, y = fract_CDE_Paryani_Accumbens$fract_CDE, weights = fract_CDE_Paryani_Accumbens$num_SPN, method = "pearson", na.rm = TRUE)
print(cor_numCAGGerm_fractCDE_Paryani_Accumbens)

#plot(x = fract_CDE_Accumbens$num_CAG_germ, y = fract_CDE_Accumbens$fract_CDE)

ggplot(fract_CDE_Paryani_Accumbens, aes(x=num_CAG_germ, y=100*fract_CDE)) + 
  geom_point(col = "blue", aes(size = num_SPN))+
  geom_smooth(method=lm, se = TRUE, mapping = aes(weight = num_SPN)) +
  theme_classic() +
  xlab("Num. CAG germline") +
  ylab("Fraction cells in C-D-E phase (%)")
  #+ ggtitle(sprintf("Sp. corr = %.2f", cor_numCAGGerm_fractCDE_Accumbens))
  #+ xlim(c(36, 70)) + ylim(c(0, 1))
ggsave("Paryani_NumCAGGermline_FractCDE_Accumbens.pdf", width = 8, height = 8)
```


```R
wtd.cor(x = fract_CDE_Paryani_Accumbens$num_CAG_germ, y = fract_CDE_Paryani_Accumbens$fract_CDE, weight = fract_CDE_Paryani_Accumbens$num_SPN)
```


```R
Calculate_confInt(r = wtd.cor(x = fract_CDE_Paryani_Accumbens$num_CAG_germ, y = fract_CDE_Paryani_Accumbens$fract_CDE, weight = fract_CDE_Paryani_Accumbens$num_SPN)[, "correlation"], w = fract_CDE_Paryani_Accumbens$num_SPN)
```


```R
fract_CDE_tmp_Accumbens
```


```R
#create metadata for plotting
Paryani_Accumbens_Grade_CDE_tmp <- lapply(split(Paryani_MSN_Accumbens@meta.data, Paryani_MSN_Accumbens@meta.data$Donor), function(x) {
    GRADE <- factor(unique(x$Grade), levels = c("Control", "HD1", "HD2", "HD3", "HD4", "HDJ"))
    SN <- unique(x$Donor)
    #NUM_MSN <- length(x$PREDICTED_PHASE)
    tmp <- lapply(split(x, x$sub_type_4), function(y) {
        NUM_MSN_CT <- length(y$PREDICTED_PHASE)
        NUM_MSN_CT_SQ <- NUM_MSN_CT**2
        CELL_TYPE <- y$sub_type_4
        if (length(unique(y$PREDICTED_PHASE)) > 1) {
          val <- 100*table(y$PREDICTED_PHASE)[2]/(table(y$PREDICTED_PHASE)[1] + table(y$PREDICTED_PHASE)[2])
        } else {
          val <- 0
        }
        names(val) <- "Fract_CDE";
        Fract_CDE <- rep(val, NUM_MSN_CT_SQ);
    return(data.frame(NUM_MSN_CT_SQ, CELL_TYPE, Fract_CDE))})
    df <- do.call(rbind, tmp)
    #print(df)
    #NUM_MSN <- length(x$PREDICTED_PHASE)
    #Fract_CDE <- 100*length(which(x$PREDICTED_PHASE == "C-D-E"))/length(x$PREDICTED_PHASE)
    #CELL_TYPE <- x$sub_type_4
    #return(df)
    #return(data.frame(SAMPLE = rep(SN, df$NUM_MSN_CT_SQ), GRADE = rep(GRADE, df$NUM_MSN_CT_SQ), CELL_TYPE = df$CELL_TYPE, Fract_CDE = df$Fract_CDE, NUM_MSN_SQ=df$NUM_MSN_CT_SQ))
    return(data.frame(SAMPLE = SN, GRADE = GRADE, CELL_TYPE = df$CELL_TYPE, Fract_CDE = df$Fract_CDE, NUM_MSN_SQ=df$NUM_MSN_CT_SQ))
})
Paryani_Accumbens_Grade_CDE <- do.call(rbind, Paryani_Accumbens_Grade_CDE_tmp)
Paryani_Accumbens_Grade_CDE <- Paryani_Accumbens_Grade_CDE[order(Paryani_Accumbens_Grade_CDE$GRADE), ]
Paryani_Accumbens_Grade_CDE$SAMPLE <- factor(Paryani_Accumbens_Grade_CDE$SAMPLE, levels = unique(Paryani_Accumbens_Grade_CDE$SAMPLE))
#head(Paryani_Accumbens_Grade_CDE)
```


```R
uni_Paryani_Accumbens_Grade_CDE <- unique(Paryani_Accumbens_Grade_CDE)
uni_Paryani_Accumbens_Grade_CDE$NUM_SPN <- sqrt(uni_Paryani_Accumbens_Grade_CDE$NUM_MSN_SQ)
#uni_Paryani_Accumbens_Grade_CDE
#uni_Paryani_Caudate_Grade_CDE[which(uni_Paryani_Accumbens_Grade_CDE$CELL_TYPE == "iSPN_2"), ]
```


```R
lapply(split(Paryani_MSN_Accumbens@meta.data[rownames(predicted_phase_test_Paryani_MSN_Accumbens_wprob), "PREDICTED_PHASE"], Paryani_MSN_Accumbens@meta.data[rownames(predicted_phase_test_Paryani_MSN_Accumbens_wprob), "Donor"]), function(x) table(x))
lapply(split(Paryani_MSN_Accumbens@meta.data[rownames(predicted_phase_test_Paryani_MSN_Accumbens_wprob), "PREDICTED_PHASE"], Paryani_MSN_Accumbens@meta.data[rownames(predicted_phase_test_Paryani_MSN_Accumbens_wprob), "Grade"]), function(x) table(x))
```


```R
#add metadata
Fract_SPN_Accumbens <- c(num_cells_Paryani["T-4915", "num_dSPN1_accumbens"]/num_cells_Paryani["T-4915", "num_cells_accumbens"],
                       num_cells_Paryani["T-4915", "num_dSPN2_accumbens"]/num_cells_Paryani["T-4915", "num_cells_accumbens"],
                       num_cells_Paryani["T-4915", "num_iSPN1_accumbens"]/num_cells_Paryani["T-4915", "num_cells_accumbens"],
                       num_cells_Paryani["T-4915", "num_iSPN2_accumbens"]/num_cells_Paryani["T-4915", "num_cells_accumbens"],
                         
                       num_cells_Paryani["T-5227", "num_dSPN1_accumbens"]/num_cells_Paryani["T-5227", "num_cells_accumbens"],
                       num_cells_Paryani["T-5227", "num_dSPN2_accumbens"]/num_cells_Paryani["T-5227", "num_cells_accumbens"],
                       num_cells_Paryani["T-5227", "num_iSPN1_accumbens"]/num_cells_Paryani["T-5227", "num_cells_accumbens"],
                       num_cells_Paryani["T-5227", "num_iSPN2_accumbens"]/num_cells_Paryani["T-5227", "num_cells_accumbens"],
                        
                       num_cells_Paryani["T-5596", "num_dSPN1_accumbens"]/num_cells_Paryani["T-5596", "num_cells_accumbens"],
                       num_cells_Paryani["T-5596", "num_dSPN2_accumbens"]/num_cells_Paryani["T-5596", "num_cells_accumbens"],
                       num_cells_Paryani["T-5596", "num_iSPN1_accumbens"]/num_cells_Paryani["T-5596", "num_cells_accumbens"],
                       num_cells_Paryani["T-5596", "num_iSPN2_accumbens"]/num_cells_Paryani["T-5596", "num_cells_accumbens"],
                         
                       num_cells_Paryani["T-5700", "num_dSPN1_accumbens"]/num_cells_Paryani["T-5700", "num_cells_accumbens"],
                       num_cells_Paryani["T-5700", "num_dSPN2_accumbens"]/num_cells_Paryani["T-5700", "num_cells_accumbens"],
                       num_cells_Paryani["T-5700", "num_iSPN1_accumbens"]/num_cells_Paryani["T-5700", "num_cells_accumbens"],
                       num_cells_Paryani["T-5700", "num_iSPN2_accumbens"]/num_cells_Paryani["T-5700", "num_cells_accumbens"],
                    
                       num_cells_Paryani["T-4285", "num_dSPN1_accumbens"]/num_cells_Paryani["T-4285", "num_cells_accumbens"],
                       num_cells_Paryani["T-4285", "num_dSPN2_accumbens"]/num_cells_Paryani["T-4285", "num_cells_accumbens"],
                       num_cells_Paryani["T-4285", "num_iSPN1_accumbens"]/num_cells_Paryani["T-4285", "num_cells_accumbens"],
                       num_cells_Paryani["T-4285", "num_iSPN2_accumbens"]/num_cells_Paryani["T-4285", "num_cells_accumbens"],
                       
                       num_cells_Paryani["T-5266", "num_dSPN1_accumbens"]/num_cells_Paryani["T-5266", "num_cells_accumbens"],
                       num_cells_Paryani["T-5266", "num_dSPN2_accumbens"]/num_cells_Paryani["T-5266", "num_cells_accumbens"],
                       num_cells_Paryani["T-5266", "num_iSPN1_accumbens"]/num_cells_Paryani["T-5266", "num_cells_accumbens"],
                       num_cells_Paryani["T-5266", "num_iSPN2_accumbens"]/num_cells_Paryani["T-5266", "num_cells_accumbens"],
                      
                       num_cells_Paryani["T-5263", "num_dSPN1_accumbens"]/num_cells_Paryani["T-5263", "num_cells_accumbens"],
                       num_cells_Paryani["T-5263", "num_dSPN2_accumbens"]/num_cells_Paryani["T-5263", "num_cells_accumbens"],
                       num_cells_Paryani["T-5263", "num_iSPN1_accumbens"]/num_cells_Paryani["T-5263", "num_cells_accumbens"],
                       num_cells_Paryani["T-5263", "num_iSPN2_accumbens"]/num_cells_Paryani["T-5263", "num_cells_accumbens"],
                            
                       num_cells_Paryani["T-5608", "num_dSPN1_accumbens"]/num_cells_Paryani["T-5608", "num_cells_accumbens"],
                       num_cells_Paryani["T-5608", "num_dSPN2_accumbens"]/num_cells_Paryani["T-5608", "num_cells_accumbens"],
                       num_cells_Paryani["T-5608", "num_iSPN1_accumbens"]/num_cells_Paryani["T-5608", "num_cells_accumbens"],
                       num_cells_Paryani["T-5608", "num_iSPN2_accumbens"]/num_cells_Paryani["T-5608", "num_cells_accumbens"], 
                          
                       num_cells_Paryani["T-5693", "num_dSPN1_accumbens"]/num_cells_Paryani["T-5693", "num_cells_accumbens"],
                       num_cells_Paryani["T-5693", "num_dSPN2_accumbens"]/num_cells_Paryani["T-5693", "num_cells_accumbens"],
                       num_cells_Paryani["T-5693", "num_iSPN1_accumbens"]/num_cells_Paryani["T-5693", "num_cells_accumbens"],
                       num_cells_Paryani["T-5693", "num_iSPN2_accumbens"]/num_cells_Paryani["T-5693", "num_cells_accumbens"],
                      
                       num_cells_Paryani["T-4161", "num_dSPN1_accumbens"]/num_cells_Paryani["T-4161", "num_cells_accumbens"],
                       num_cells_Paryani["T-4161", "num_dSPN2_accumbens"]/num_cells_Paryani["T-4161", "num_cells_accumbens"],
                       num_cells_Paryani["T-4161", "num_iSPN1_accumbens"]/num_cells_Paryani["T-4161", "num_cells_accumbens"],
                       num_cells_Paryani["T-4161", "num_iSPN2_accumbens"]/num_cells_Paryani["T-4161", "num_cells_accumbens"],
                         
                       num_cells_Paryani["T-5109", "num_dSPN1_accumbens"]/num_cells_Paryani["T-5109", "num_cells_accumbens"],
                       num_cells_Paryani["T-5109", "num_dSPN2_accumbens"]/num_cells_Paryani["T-5109", "num_cells_accumbens"],
                       num_cells_Paryani["T-5109", "num_iSPN1_accumbens"]/num_cells_Paryani["T-5109", "num_cells_accumbens"],
                       num_cells_Paryani["T-5109", "num_iSPN2_accumbens"]/num_cells_Paryani["T-5109", "num_cells_accumbens"],
                      
                       num_cells_Paryani["T-5493", "num_dSPN1_accumbens"]/num_cells_Paryani["T-5493", "num_cells_accumbens"],
                       num_cells_Paryani["T-5493", "num_dSPN2_accumbens"]/num_cells_Paryani["T-5493", "num_cells_accumbens"],
                       num_cells_Paryani["T-5493", "num_iSPN1_accumbens"]/num_cells_Paryani["T-5493", "num_cells_accumbens"],
                       num_cells_Paryani["T-5493", "num_iSPN2_accumbens"]/num_cells_Paryani["T-5493", "num_cells_accumbens"],
                      
                       num_cells_Paryani["T-5575", "num_dSPN1_accumbens"]/num_cells_Paryani["T-5575", "num_cells_accumbens"],
                       #num_cells_Paryani["T-5575", "num_dSPN2_accumbens"]/num_cells_Paryani["T-5575", "num_cells_accumbens"],
                       num_cells_Paryani["T-5575", "num_iSPN1_accumbens"]/num_cells_Paryani["T-5575", "num_cells_accumbens"],
                       num_cells_Paryani["T-5575", "num_iSPN2_accumbens"]/num_cells_Paryani["T-5575", "num_cells_accumbens"],
                      
                       num_cells_Paryani["T-4638", "num_dSPN1_accumbens"]/num_cells_Paryani["T-4638", "num_cells_accumbens"],
                       num_cells_Paryani["T-4638", "num_dSPN2_accumbens"]/num_cells_Paryani["T-4638", "num_cells_accumbens"],
                       num_cells_Paryani["T-4638", "num_iSPN1_accumbens"]/num_cells_Paryani["T-4638", "num_cells_accumbens"],
                       num_cells_Paryani["T-4638", "num_iSPN2_accumbens"]/num_cells_Paryani["T-4638", "num_cells_accumbens"],
                      
                       num_cells_Paryani["T-5534", "num_dSPN1_accumbens"]/num_cells_Paryani["T-5534", "num_cells_accumbens"],
                       num_cells_Paryani["T-5534", "num_dSPN2_accumbens"]/num_cells_Paryani["T-5534", "num_cells_accumbens"],
                       num_cells_Paryani["T-5534", "num_iSPN1_accumbens"]/num_cells_Paryani["T-5534", "num_cells_accumbens"],
                       num_cells_Paryani["T-5534", "num_iSPN2_accumbens"]/num_cells_Paryani["T-5534", "num_cells_accumbens"],
                      
                       num_cells_Paryani["T-5714", "num_dSPN1_accumbens"]/num_cells_Paryani["T-5714", "num_cells_accumbens"],
                       num_cells_Paryani["T-5714", "num_dSPN2_accumbens"]/num_cells_Paryani["T-5714", "num_cells_accumbens"],
                       num_cells_Paryani["T-5714", "num_iSPN1_accumbens"]/num_cells_Paryani["T-5714", "num_cells_accumbens"],
                       num_cells_Paryani["T-5714", "num_iSPN2_accumbens"]/num_cells_Paryani["T-5714", "num_cells_accumbens"])
#Fract_SPN_accumbens
uni_Paryani_Accumbens_Grade_CDE$Fract_SPN <- Fract_SPN_Accumbens
uni_Paryani_Accumbens_Grade_CDE
```


```R
unique(uni_Paryani_Accumbens_Grade_CDE[, c("SAMPLE", "GRADE")])
```


```R
# ggplot(Paryani_Accumbens_Grade_CDE, aes(x = SAMPLE, y = Fract_CDE, color = CELL_TYPE)) +
# geom_boxplot(varwidth = TRUE, lwd = 5, position = position_dodge2(preserve = "single")) +
# scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
# theme_classic() +
# xlab("SAMPLE") +
# ylim(c(0, 100)) +
# ylab("Fraction cells in C-D-E phase (%)") 
# ggsave("Paryani_Accumbens_CELLTYPE_vs_fraction_SPN_CDE.pdf", width = 8, height = 8)

uni_Paryani_Accumbens_Grade_CDE_HD <- uni_Paryani_Accumbens_Grade_CDE[which(uni_Paryani_Accumbens_Grade_CDE$GRADE != "Control"), ]
#uni_Paryani_Accumbens_Grade_CDE_CTRL <- uni_Paryani_Accumbens_Grade_CDE[which(uni_Paryani_Accumbens_Grade_CDE$GRADE == "Control"), ]

#plot fraction of SPNs in C-D-E phase by HD grade
ggplot(uni_Paryani_Accumbens_Grade_CDE_HD, aes(x = GRADE, y = Fract_CDE, color = CELL_TYPE)) +
geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.1, height = 0)) +
scale_color_manual(values =c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
theme_classic() +
xlab("HD GRADE") +
ylim(c(0, 100)) +
ylab("Fraction cells in C-D-E phase (%)") 
ggsave("Paryani_Accumbens_HD_GRADE_vs_fraction_SPN_CDE.pdf", width = 6, height = 6)
```


```R
# ggplot(uni_Paryani_Accumbens_Grade_CDE[which(uni_Paryani_Accumbens_Grade_CDE$NUM_SPN > 30), ], aes(x = GRADE, y = Fract_CDE, color = CELL_TYPE)) +
# geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.1, height = 0)) +
# scale_color_manual(values =c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
# theme_classic() +
# xlab("HD GRADE") +
# ylim(c(0, 100)) +
# ylab("Fraction cells in C-D-E phase (%)") 
# ggsave("Paryani_Accumbens_GRADE_vs_fraction_SPN_CDE_gt30.pdf", width = 6, height = 6)
```


```R
#subset HD
Paryani_HD_MSN_Accumbens <- subset(x = Paryani_MSN_Accumbens, subset = Condition == "HD")
```

# Plot_controls


```R
Handsaker_donor_metadata_CTRL <- Handsaker_donor_metadata_full[which(Handsaker_donor_metadata_full$Status == "Control"), ]
Handsaker_fract_CDE_tmp_Caudate_CTRL <- Handsaker_donor_metadata_CTRL
Handsaker_fract_CDE_tmp_Caudate_CTRL$SAMPLE <- rownames(Handsaker_donor_metadata_CTRL)
Handsaker_fract_CDE_tmp_Caudate_CTRL$DATASET <- "Handsaker"
Handsaker_fract_CDE_tmp_Caudate_CTRL$TISSUE = "Caudate"
Handsaker_fract_CDE_tmp_Caudate_CTRL <- Handsaker_fract_CDE_tmp_Caudate_CTRL[, c("DATASET", "SAMPLE", "PRED_FRACT_CDE", "NUM_SPN", "TISSUE")]
Handsaker_fract_CDE_tmp_Caudate_CTRL$Fract_SPN <- num_cells_Handsaker["S04577", "num_spns"]/num_cells_Handsaker["S04577", "num_cells"]
```


```R
Lee_CTRL_MSN_Caudate <- subset(x = Lee_MSN_Caudate, subset = Condition=="Control")
Lee_fract_CDE_tmp_Caudate_CTRL <- lapply(split(Lee_CTRL_MSN_Caudate@meta.data, Lee_CTRL_MSN_Caudate@meta.data$NBB_ID), function(x) {
    SAMPLE=unique(x$NBB_ID)
    num_SPN <- dim(x)[1]
    num_CDE <- length(which(x$PREDICTED_PHASE == "C-D-E"))
    fraction_CDE <- num_CDE/num_SPN
    num_CAG <- gsub(x = x$CAG[1], pattern = "-.*", replacement = "")
    #cat(sprintf("Fraction MSNs in C-D-E phase: %.2f\n", fraction_CDE*100))
    df <- data.frame(DATASET = "Lee", SAMPLE, PRED_FRACT_CDE=fraction_CDE, NUM_SPN = num_SPN, TISSUE = "Caudate")
    return(df)
})
#Lee_fract_CDE_tmp_Caudate_CTRL
Lee_fract_CDE_Caudate_CTRL <- do.call(rbind, Lee_fract_CDE_tmp_Caudate_CTRL)
```


```R
Fract_SPN_Caudate <- c((num_cells_Lee["3345", "num_D1_SPN_caudate"] + num_cells_Lee["3345", "num_D2_SPN_caudate"])/num_cells_Lee["3345", "num_cells_caudate"],
                       (num_cells_Lee["4294", "num_D1_SPN_caudate"] + num_cells_Lee["4294", "num_D2_SPN_caudate"])/num_cells_Lee["4294", "num_cells_caudate"],
                       (num_cells_Lee["4308", "num_D1_SPN_caudate"] + num_cells_Lee["4308", "num_D2_SPN_caudate"])/num_cells_Lee["4308", "num_cells_caudate"],                  
                       (num_cells_Lee["4494", "num_D1_SPN_caudate"] + num_cells_Lee["4494", "num_D2_SPN_caudate"])/num_cells_Lee["4494", "num_cells_caudate"],
                       (num_cells_Lee["4621", "num_D1_SPN_caudate"] + num_cells_Lee["4621", "num_D2_SPN_caudate"])/num_cells_Lee["4621", "num_cells_caudate"],
                       (num_cells_Lee["A39R", "num_D1_SPN_caudate"] + num_cells_Lee["A39R", "num_D2_SPN_caudate"])/num_cells_Lee["A39R", "num_cells_caudate"],
                       (num_cells_Lee["A47L", "num_D1_SPN_caudate"] + num_cells_Lee["A47L", "num_D2_SPN_caudate"])/num_cells_Lee["A47L", "num_cells_caudate"])
Lee_fract_CDE_Caudate_CTRL$Fract_SPN <- Fract_SPN_Caudate
#Lee_fract_CDE_Caudate_CTRL
```


```R
Lee_CTRL_MSN_Putamen <- subset(x = Lee_MSN_Putamen, subset = Condition=="Control")
Lee_fract_CDE_tmp_Putamen_CTRL <- lapply(split(Lee_CTRL_MSN_Putamen@meta.data, Lee_CTRL_MSN_Putamen@meta.data$NBB_ID), function(x) {
    SAMPLE=unique(x$NBB_ID)
    num_SPN <- dim(x)[1]
    num_CDE <- length(which(x$PREDICTED_PHASE == "C-D-E"))
    fraction_CDE <- num_CDE/num_SPN
    num_CAG <- gsub(x = x$CAG[1], pattern = "-.*", replacement = "")
    #cat(sprintf("Fraction MSNs in C-D-E phase: %.2f\n", fraction_CDE*100))
    df <- data.frame(DATASET = "Lee", SAMPLE, PRED_FRACT_CDE=fraction_CDE, NUM_SPN = num_SPN, TISSUE = "Putamen")
    return(df)
})
#Lee_fract_CDE_tmp_Caudate_CTRL
Lee_fract_CDE_Putamen_CTRL <- do.call(rbind, Lee_fract_CDE_tmp_Putamen_CTRL)
```


```R
Fract_SPN_Putamen <- c((num_cells_Lee["4294", "num_D1_SPN_putamen"] + num_cells_Lee["4294", "num_D2_SPN_putamen"])/num_cells_Lee["4294", "num_cells_putamen"],
                       (num_cells_Lee["4308", "num_D1_SPN_putamen"] + num_cells_Lee["4308", "num_D2_SPN_putamen"])/num_cells_Lee["4308", "num_cells_putamen"],                  
                       (num_cells_Lee["4494", "num_D1_SPN_putamen"] + num_cells_Lee["4494", "num_D2_SPN_putamen"])/num_cells_Lee["4494", "num_cells_putamen"],
                       (num_cells_Lee["A39R", "num_D1_SPN_putamen"] + num_cells_Lee["A39R", "num_D2_SPN_putamen"])/num_cells_Lee["A39R", "num_cells_putamen"],
                       (num_cells_Lee["A47L", "num_D1_SPN_putamen"] + num_cells_Lee["A47L", "num_D2_SPN_putamen"])/num_cells_Lee["A47L", "num_cells_putamen"])
Lee_fract_CDE_Putamen_CTRL$Fract_SPN <- Fract_SPN_Putamen
#Lee_fract_CDE_Putamen_CTRL
```


```R
Paryani_CTRL_MSN_Caudate <- subset(x = Paryani_MSN_Caudate, subset = Condition=="Con")
Paryani_fract_CDE_tmp_Caudate_CTRL <- lapply(split(Paryani_CTRL_MSN_Caudate@meta.data, Paryani_CTRL_MSN_Caudate@meta.data$Donor), function(x) {
    SAMPLE=unique(x$Donor)
    num_SPN <- dim(x)[1]
    num_CDE <- length(which(x$PREDICTED_PHASE == "C-D-E"))
    fraction_CDE <- num_CDE/num_SPN
    num_CAG <- gsub(x = x$CAG[1], pattern = "-.*", replacement = "")
    #cat(sprintf("Fraction MSNs in C-D-E phase: %.2f\n", fraction_CDE*100))
    df <- data.frame(DATASET = "Paryani", SAMPLE, PRED_FRACT_CDE=fraction_CDE, NUM_SPN = num_SPN, TISSUE = "Caudate")
    return(df)
})
#Paryani_fract_CDE_tmp_Caudate_CTRL
Paryani_fract_CDE_Caudate_CTRL <- do.call(rbind, Paryani_fract_CDE_tmp_Caudate_CTRL)
#Paryani_fract_CDE_Caudate_CTRL
```


```R
Fract_SPN_Caudate <- c((num_cells_Paryani["T-4812", "num_dSPN1_caudate"] +  num_cells_Paryani["T-4812", "num_dSPN2_caudate"] + num_cells_Paryani["T-4812", "num_iSPN1_caudate"] + num_cells_Paryani["T-4812", "num_iSPN2_caudate"])/num_cells_Paryani["T-4812", "num_cells_caudate"],
                       (num_cells_Paryani["T-4915", "num_dSPN1_caudate"] +  num_cells_Paryani["T-4915", "num_dSPN2_caudate"] + num_cells_Paryani["T-4915", "num_iSPN1_caudate"] + num_cells_Paryani["T-4915", "num_iSPN2_caudate"])/num_cells_Paryani["T-4915", "num_cells_caudate"],
                       (num_cells_Paryani["T-5227", "num_dSPN1_caudate"] +  num_cells_Paryani["T-5227", "num_dSPN2_caudate"] + num_cells_Paryani["T-5227", "num_iSPN1_caudate"] + num_cells_Paryani["T-5227", "num_iSPN2_caudate"])/num_cells_Paryani["T-5227", "num_cells_caudate"],
                       (num_cells_Paryani["T-5596", "num_dSPN1_caudate"] +  num_cells_Paryani["T-5596", "num_dSPN2_caudate"] + num_cells_Paryani["T-5596", "num_iSPN1_caudate"] + num_cells_Paryani["T-5596", "num_iSPN2_caudate"])/num_cells_Paryani["T-5596", "num_cells_caudate"],
                       (num_cells_Paryani["T-5700", "num_dSPN1_caudate"] +  num_cells_Paryani["T-5700", "num_dSPN2_caudate"] + num_cells_Paryani["T-5700", "num_iSPN1_caudate"] + num_cells_Paryani["T-5700", "num_iSPN2_caudate"])/num_cells_Paryani["T-5700", "num_cells_caudate"])

#Fract_SPN_Caudate
Paryani_fract_CDE_Caudate_CTRL$Fract_SPN <- Fract_SPN_Caudate
#Paryani_fract_CDE_Caudate_CTRL
```


```R
Paryani_CTRL_MSN_Accumbens <- subset(x = Paryani_MSN_Accumbens, subset = Condition=="Con")
Paryani_fract_CDE_tmp_Accumbens_CTRL <- lapply(split(Paryani_CTRL_MSN_Accumbens@meta.data, Paryani_CTRL_MSN_Accumbens@meta.data$Donor), function(x) {
    SAMPLE=unique(x$Donor)
    num_SPN <- dim(x)[1]
    num_CDE <- length(which(x$PREDICTED_PHASE == "C-D-E"))
    fraction_CDE <- num_CDE/num_SPN
    num_CAG <- gsub(x = x$CAG[1], pattern = "-.*", replacement = "")
    #cat(sprintf("Fraction MSNs in C-D-E phase: %.2f\n", fraction_CDE*100))
    df <- data.frame(DATASET = "Paryani", SAMPLE, PRED_FRACT_CDE=fraction_CDE, NUM_SPN = num_SPN, TISSUE = "Accumbens")
    return(df)
})
#Paryani_fract_CDE_tmp_Accumbens_CTRL
Paryani_fract_CDE_Accumbens_CTRL <- do.call(rbind, Paryani_fract_CDE_tmp_Accumbens_CTRL)
#Paryani_fract_CDE_Accumbens_CTRL
```


```R
Fract_SPN_Accumbens <- c((num_cells_Paryani["T-4915", "num_dSPN1_accumbens"] +  num_cells_Paryani["T-4915", "num_dSPN2_accumbens"] + num_cells_Paryani["T-4915", "num_iSPN1_accumbens"] + num_cells_Paryani["T-4915", "num_iSPN2_accumbens"])/num_cells_Paryani["T-4915", "num_cells_accumbens"],
                       (num_cells_Paryani["T-5227", "num_dSPN1_accumbens"] +  num_cells_Paryani["T-5227", "num_dSPN2_accumbens"] + num_cells_Paryani["T-5227", "num_iSPN1_accumbens"] + num_cells_Paryani["T-5227", "num_iSPN2_accumbens"])/num_cells_Paryani["T-5227", "num_cells_accumbens"],
                       (num_cells_Paryani["T-5596", "num_dSPN1_accumbens"] +  num_cells_Paryani["T-5596", "num_dSPN2_accumbens"] + num_cells_Paryani["T-5596", "num_iSPN1_accumbens"] + num_cells_Paryani["T-5596", "num_iSPN2_accumbens"])/num_cells_Paryani["T-5596", "num_cells_accumbens"],
                       (num_cells_Paryani["T-5700", "num_dSPN1_accumbens"] +  num_cells_Paryani["T-5700", "num_dSPN2_accumbens"] + num_cells_Paryani["T-5700", "num_iSPN1_accumbens"] + num_cells_Paryani["T-5700", "num_iSPN2_accumbens"])/num_cells_Paryani["T-5700", "num_cells_accumbens"])

#Fract_SPN_accumbens
Paryani_fract_CDE_Accumbens_CTRL$Fract_SPN <- Fract_SPN_Accumbens
#Paryani_fract_CDE_Accumbens_CTRL
```


```R
CTRL_all_datasets <- rbind(Handsaker_fract_CDE_tmp_Caudate_CTRL, Lee_fract_CDE_Caudate_CTRL, Lee_fract_CDE_Putamen_CTRL, Paryani_fract_CDE_Caudate_CTRL, Paryani_fract_CDE_Accumbens_CTRL)
CTRL_all_datasets
```


```R
#plot results for CTRL donors
ggplot(CTRL_all_datasets, aes(x = TISSUE, color = DATASET, y = 100*PRED_FRACT_CDE, size = Fract_SPN)) +
geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.3, height = 0)) +
#scale_color_manual(values =c("#F8766D")) +
theme_classic() +
xlab("Tissue") +
ylim(c(0, 100)) +
ylab("Non-HD donors - Fraction cells in C-D-E phase (%)")
ggsave("CTRL_all_datasets_fraction_CDE.pdf", width = 8, height = 8)

ggplot(CTRL_all_datasets[which(CTRL_all_datasets$TISSUE == "Caudate"), ], aes(x = TISSUE, color = DATASET, y = 100*PRED_FRACT_CDE, size = Fract_SPN)) +
geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.5, height = 0)) +
#scale_color_manual(values =c("#F8766D")) +
theme_classic() +
xlab("Tissue") +
ylim(c(0, 100)) +
ylab("Non-HD donors - Fraction cells in C-D-E phase (%)")
ggsave("CTRL_all_datasets_Caudate_fraction_CDE.pdf", width = 8, height = 8)

ggplot(CTRL_all_datasets[which(CTRL_all_datasets$TISSUE != "Caudate"), ], aes(x = TISSUE, color = DATASET, y = 100*PRED_FRACT_CDE, size = Fract_SPN)) +
geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.3, height = 0)) +
#scale_color_manual(values =c("#F8766D")) +
theme_classic() +
xlab("Tissue") +
ylim(c(0, 100)) +
ylab("Non-HD donors - Fraction cells in C-D-E phase (%)")
ggsave("CTRL_all_datasets_nonCaudate_fraction_CDE.pdf", width = 8, height = 8)
#"#7CAE00", "#00BFC4", "#C77CFF"
```


```R
length(unique(uni_Paryani_Caudate_Grade_CDE$SAMPLE))
length(unique(uni_Lee_Caudate_Grade_CDE$SAMPLE))

length(unique(uni_Paryani_Accumbens_Grade_CDE$SAMPLE))
length(unique(uni_Lee_Putamen_Grade_CDE$SAMPLE))
```

# Import_Xu_dataset


```R
Read_10X_matrixes = function(folder, project = "") {
    #extract counts
    counts <- Read10X(data.dir = folder)
    #create a single cell experiment object
    sce <- SingleCellExperiment(list(counts = counts))
    #create a seurat object
    seuratObject <- CreateSeuratObject(counts(sce), project = project, min.cells = 0, min.features = 0)
    cat(sprintf("%s; Num. features = %d; Num. cells = %d\n", basename(folder), dim(seuratObject)[1], dim(seuratObject)[2]))
    return(seuratObject)
}
```


```R
matrix_folders <- list.dirs(path = "/mnt/projects/labs/CLAB/PROJECT_ELongATE/Xu/", full.names = TRUE)
matrix_folders <- matrix_folders[grep(x = basename(matrix_folders), pattern = "GSM")]
names(matrix_folders) <- basename(matrix_folders)
```


```R
Xu_counts <- list()
for (i in 1:length(matrix_folders)) {
    Xu_counts[[i]] <- Read_10X_matrixes(folder = matrix_folders[i], project = "Xu")
}
```


```R
names(Xu_counts) <- names(matrix_folders)
```


```R
Xu_all <- merge(x = Xu_counts[[1]],
                         y = Xu_counts[-1],
                         add.cell.ids = names(Xu_counts))
```


```R
head(Xu_all@meta.data)
```


```R
Xu_all[["SAMPLE"]] <- factor(gsub(x = rownames(Xu_all@meta.data), pattern = "_GSM.*", replacement = ""))
Xu_all[["condition"]] <- factor(gsub(x = gsub(x = gsub(x = rownames(Xu_all@meta.data), pattern = "_GSM.*", replacement = ""), pattern = ".*_", replacement = ""), pattern = "\\d", replacement = ""))
Xu_all[["percent.mt"]] <- PercentageFeatureSet(Xu_all, pattern = "^MT-")
```


```R
Xu_all_unfiltered <- Xu_all
```


```R
options(repr.plot.width=20, repr.plot.height=15)   #for plotting

nFeature_RNA_min <- 500
nFeature_RNA_max <- 8000
nCount_RNA_min <- 5000
nCount_RNA_max <- 30000
percent.mt_max <- 10

#filter cells
  # Visualize QC metrics as a violin plot
  #plot(VlnPlot(Xu_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  p1 <- VlnPlot(Xu_all_unfiltered, features = c("nFeature_RNA"), group.by = "orig.ident") + geom_hline(yintercept = nFeature_RNA_min, linetype = "dashed") + geom_hline(yintercept = nFeature_RNA_max, linetype = "dashed")
  p2 <- VlnPlot(Xu_all_unfiltered, features = c("nCount_RNA"), group.by = "orig.ident") + geom_hline(yintercept = nCount_RNA_min, linetype = "dashed") + geom_hline(yintercept = nCount_RNA_max, linetype = "dashed")
  #+ scale_y_continuous(limits = c(0, 2000)) 
  p3 <- VlnPlot(Xu_all_unfiltered, features = c("percent.mt"), group.by = "orig.ident") + geom_hline(yintercept = percent.mt_max, linetype = "dashed")
  plot(ggarrange(p1, p2, p3, ncol = 3))
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  plot1 <- FeatureScatter(Xu_all_unfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
  plot2 <- FeatureScatter(Xu_all_unfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
  plot(plot1 + plot2)
  cat(sprintf("Number of cells before QC: %d\n", length(Cells(x = Xu_all_unfiltered))))
  Xu_all <- subset(Xu_all, subset = nFeature_RNA > nFeature_RNA_min & nFeature_RNA < nFeature_RNA_max & percent.mt < percent.mt_max & nCount_RNA > nCount_RNA_min & nCount_RNA < nCount_RNA_max)
  cat(sprintf("Number of cells after QC: %d\n", length(Cells(x = Xu_all))))
  p1 <- VlnPlot(Xu_all, features = c("nFeature_RNA"), group.by = "orig.ident") + geom_hline(yintercept = nFeature_RNA_min, linetype = "dashed") + geom_hline(yintercept = nFeature_RNA_max, linetype = "dashed") + ggtitle("Post-filter nFeature_RNA")
  #+ scale_y_continuous(limits = c(0, 2000)) 
  p2 <- VlnPlot(Xu_all, features = c("nCount_RNA"), group.by = "orig.ident") + geom_hline(yintercept = nCount_RNA_min, linetype = "dashed") + geom_hline(yintercept = nCount_RNA_max, linetype = "dashed") + ggtitle("Post-filter nCount_RNA")
  p3 <- VlnPlot(Xu_all, features = c("percent.mt"), group.by = "orig.ident") + geom_hline(yintercept = percent.mt_max, linetype = "dashed") + ggtitle("Post-filter percent.mt")
  plot(ggarrange(p1, p2, p3, ncol = 3))
  plot3 <- FeatureScatter(Xu_all, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident") + ggtitle("Post-filter")
  plot4 <- FeatureScatter(Xu_all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident") + ggtitle("Post-filter")
  plot(plot3 + plot4)
```


```R
# #Run SCTransform
Xu_all <- SCTransform(Xu_all, variable.features.n = 3000,  ncells = 5000, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
#Run PCA
Xu_all <- RunPCA(Xu_all, features = VariableFeatures(object = Xu_all), npcs = 50, reduction.name = "pca")
ElbowPlot(Xu_all, reduction = "pca", ndims = 50)
num_dims <- 40
#Run UMAP
Xu_all <- RunUMAP(Xu_all, dims = 1:num_dims, reduction = "pca", reduction.name = "umap", reduction.key = "umap")
```


```R
DimPlot(Xu_all, group.by = "SAMPLE", reduction = "umap", pt.size = 1)
Xu_all <- RunHarmony(object = Xu_all, reduction.use = "pca", group.by.vars = "SAMPLE", 
                           reduction.save = "harmonyPca", assay = "SCT", verbose = FALSE, normalization.method = "SCT")
Xu_all <- RunUMAP(Xu_all, dims = 1:num_dims, reduction = "harmonyPca", reduction.name = "harmony", reduction.key = "harmony")
DimPlot(Xu_all, group.by = "SAMPLE", reduction = "harmony", pt.size = 1)
```


```R
FeaturePlot(Xu_all, reduction = "harmony", features = c("DRD1", "DRD2"), cols = c("blue", "red"), pt.size = 1)
FeaturePlot(Xu_all, reduction = "harmony", features = "DRD1", cols = c("blue", "red"), pt.size = 1)
ggsave("harmony_Xu_ALL_DRD1.pdf", width = 8, height = 8)
FeaturePlot(Xu_all, reduction = "harmony", features = "DRD2", cols = c("blue", "red"), pt.size = 1)
ggsave("harmony_Xu_ALL_DRD2.pdf", width = 8, height = 8)
#FeaturePlot(Xu_all, reduction = "harmony", features = c("DRD1", "DRD2"), blend = TRUE, pt.size = 1)
Xu_all <- AddModuleScore(object = Xu_all, features = list("DRD1", "DRD2"), name = "DRD1_DRD2")
FeaturePlot(Xu_all, reduction = "harmony", features = "DRD1_DRD21", cols = c("blue", "red"), pt.size = 1) 
FeaturePlot(Xu_all, reduction = "harmony", features = "nCount_RNA", cols = c("blue", "red"), pt.size = 1)
FeaturePlot(Xu_all, reduction = "harmony", features = "nFeature_RNA", cols = c("blue", "red"), pt.size = 1)
FeaturePlot(Xu_all, reduction = "harmony", features = "percent.mt", cols = c("blue", "red"), pt.size = 1)
```

# Identify_SPNs_in_Xu_dataset


```R
#do clustering
Xu_all <- FindNeighbors(Xu_all, reduction = "pca", dims = 1:num_dims)
Xu_all <- FindClusters(Xu_all, resolution = 0.1, method = 4)
DimPlot(Xu_all, reduction = "harmony", group.by = "seurat_clusters", pt.size = 1)
Xu_all_SPN <- subset(Xu_all, subset = seurat_clusters %in% c("0", "4"))
```


```R
FeaturePlot(Xu_all_SPN, reduction = "harmony", features = c("DRD1", "DRD2"), cols = c("blue", "red"), pt.size = 1)
DimPlot(Xu_all_SPN, reduction = "harmony", group.by = "SAMPLE", pt.size = 1)
ggsave("UMAP_Xu_MSN_ALL_SAMPLE.pdf", width = 8, height = 8)
DimPlot(Xu_all_SPN, reduction = "harmony", group.by = "condition", pt.size = 1)
ggsave("UMAP_Xu_MSN_ALL_SAMPLE.pdf", width = 8, height = 8)
#FeaturePlot(Xu_all_SPN, reduction = "umap", features = "nCount_RNA", cols = c("blue", "red"), pt.size = 1)
#FeaturePlot(Xu_all_SPN, reduction = "umap", features = "nFeature_RNA", cols = c("blue", "red"), pt.size = 1)
#FeaturePlot(Xu_all_SPN, reduction = "umap", features = "percent.mt", cols = c("blue", "red"), pt.size = 1)
```


```R
# Xu_CTRL_SPN <- subset(Xu_all_SPN, subset = condition == "CTRL")
# Xu_AD_SPN <- subset(Xu_all_SPN, subset = condition == "AD")
# Xu_PD_SPN <- subset(Xu_all_SPN, subset = condition == "PD")
```


```R
# FeaturePlot(Xu_CTRL_SPN, reduction = "harmony", features = c("DRD1", "DRD2"), cols = c("blue", "red"), pt.size = 1)
# FeaturePlot(Xu_CTRL_SPN, reduction = "harmony", features = c("DRD1", "DRD2"), cols = c("blue", "red"), pt.size = 1)
# FeaturePlot(Xu_CTRL_SPN, reduction = "harmony", features = c("DRD1", "DRD2"), cols = c("blue", "red"), pt.size = 1)

# DimPlot(Xu_CTRL_SPN, reduction = "harmony", group.by = "SAMPLE", pt.size = 1)
# DimPlot(Xu_CTRL_SPN, reduction = "harmony", group.by = "condition", pt.size = 1)

# FeaturePlot(Xu_AD_SPN, reduction = "harmony", features = c("DRD1", "DRD2"), cols = c("blue", "red"), pt.size = 1)
# DimPlot(Xu_AD_SPN, reduction = "harmony", group.by = "SAMPLE", pt.size = 1)
# DimPlot(Xu_AD_SPN, reduction = "harmony", group.by = "condition", pt.size = 1)

# FeaturePlot(Xu_PD_SPN, reduction = "harmony", features = c("DRD1", "DRD2"), cols = c("blue", "red"), pt.size = 1)
# DimPlot(Xu_PD_SPN, reduction = "harmony", group.by = "SAMPLE", pt.size = 1)
# DimPlot(Xu_PD_SPN, reduction = "harmony", group.by = "condition", pt.size = 1)
```

# Run_model_on_Xu_dataset


```R
phaseC_genes_Xu_all_SPN <- intersect(colnames(phase_all), rownames(Xu_all_SPN))
phaseC_genes_notExpXu_all_SPN_names <- setdiff(colnames(phase_all), rownames(Xu_all_SPN))
phaseC_genes_notExpXu_all_SPN <- matrix(data = 0, nrow = length(phaseC_genes_notExpXu_all_SPN_names), ncol = dim(Xu_all_SPN@assays$SCT$data)[2])
rownames(phaseC_genes_notExpXu_all_SPN) <- phaseC_genes_notExpXu_all_SPN_names
colnames(phaseC_genes_notExpXu_all_SPN) <- colnames(Xu_all_SPN@assays$SCT$data)
phase_counts_Xu_all_SPN <- Xu_all_SPN@assays$SCT$data[phaseC_genes_Xu_all_SPN, ]
phase_counts_Xu_all_SPN <- rbind(phase_counts_Xu_all_SPN, phaseC_genes_notExpXu_all_SPN)
```


```R
phase_counts_Xu_all_SPN <- phase_counts_Xu_all_SPN[colnames(phase_all), ]
phase_test_Xu_all_SPN <- data.frame(t(phase_counts_Xu_all_SPN))
head(phase_test_Xu_all_SPN)
```


```R
class_thr <- perf_training_test_noNA$thr
predicted_phase_test_Xu_all_SPN_wprob <- predict(model_phase,
                                                 newx = as.matrix(phase_test_Xu_all_SPN),
                                                 s = lambda_val,
                                                 type = "response")
predicted_phase_test_Xu_all_SPN <- rep("A-B", dim(predicted_phase_test_Xu_all_SPN_wprob)[1])  
predicted_phase_test_Xu_all_SPN[which(predicted_phase_test_Xu_all_SPN_wprob > class_thr)] <- "C-D-E"
Xu_all_SPN[["PREDICTED_PHASE"]] <- NA
Xu_all_SPN@meta.data[rownames(predicted_phase_test_Xu_all_SPN_wprob), "PREDICTED_PHASE"] <- predicted_phase_test_Xu_all_SPN
Xu_all_SPN[["PROB_PHASE_CDE"]] <- NA
Xu_all_SPN@meta.data[rownames(predicted_phase_test_Xu_all_SPN_wprob), "PROB_PHASE_CDE"] <- predicted_phase_test_Xu_all_SPN_wprob
```


```R
Xu_all_SPN$condition <- factor(Xu_all_SPN$condition, levels = c("CTRL", "AD", "PD"))
```


```R
Xu_CTRL_SPN <- subset(Xu_all_SPN, subset = condition == "CTRL")
Xu_AD_SPN <- subset(Xu_all_SPN, subset = condition == "AD")
Xu_PD_SPN <- subset(Xu_all_SPN, subset = condition == "PD")
```


```R
table(Xu_CTRL_SPN[["PREDICTED_PHASE"]] )
table(Xu_AD_SPN[["PREDICTED_PHASE"]] )
table(Xu_PD_SPN[["PREDICTED_PHASE"]] )
```


```R
DimPlot(Xu_CTRL_SPN, reduction = "harmony", group.by = "PREDICTED_PHASE", pt.size = 1)
ggsave("UMAP_Xu_CTRL_MSN_PREDPHASE.pdf", width = 8, height = 8)
DimPlot(Xu_AD_SPN, reduction = "harmony", group.by = "PREDICTED_PHASE", pt.size = 1)
ggsave("UMAP_Xu_AD_MSN_PREDPHASE.pdf", width = 8, height = 8)
DimPlot(Xu_PD_SPN, reduction = "harmony", group.by = "PREDICTED_PHASE", pt.size = 1)
ggsave("UMAP_Xu_PD_MSN_PREDPHASE.pdf", width = 8, height = 8)
```

# All_datasets


```R
DefaultAssay(SPN_HD) <- "RNA"
DefaultAssay(Paryani_HD_MSN_Caudate) <- "RNA"
DefaultAssay(Paryani_HD_MSN_Accumbens) <- "RNA"
DefaultAssay(Lee_HD_MSN_Caudate) <- "RNA"
DefaultAssay(Lee_HD_MSN_Putamen) <- "RNA"
```


```R
objs <- list(SPN_HD_noSCT, 
             Paryani_HD_MSN_Caudate, 
             Paryani_HD_MSN_Accumbens, 
             Lee_HD_MSN_Caudate, 
             Lee_HD_MSN_Putamen)

# Rimuovo eventuali SCT
objs <- lapply(objs, function(x) {
  if ("SCT" %in% names(x@assays)) {
    x[["SCT"]] <- NULL
  }
  DefaultAssay(x) <- "RNA"
  return(x)
})

ALL_DATASETS <- merge(x = objs[[1]], y = objs[-1], add.cell.ids = c("HandsakerCaudate", "ParyaniCaudate", "ParyaniAccumbens", "LeeCaudate", "LeePutamen"))
```


```R
ALL_DATASETS
```


```R
ALL_DATASETS@meta.data$Dataset <- factor(gsub(x = rownames(ALL_DATASETS@meta.data), pattern = "_.*", replacement = ""))
```


```R
ALL_DATASETS@meta.data
```


```R
#Run SCTransform
ALL_DATASETS <- SCTransform(ALL_DATASETS, vars.to.regress = c("nCount_RNA", "nFeature_RNA"))
cat(sprintf("Num. features = %d; Num. SPN = %d\n", dim(ALL_DATASETS)[1], dim(ALL_DATASETS)[2]))
```


```R
#Run PCA and UMAP for HD donors
ALL_DATASETS <- RunPCA(ALL_DATASETS, features = VariableFeatures(object = ALL_DATASETS), npcs = 50, reduction.name = "pca")
ElbowPlot(ALL_DATASETS, reduction = "pca", ndims = 50)
num_dims <- 40
ALL_DATASETS <- RunUMAP(ALL_DATASETS, dims = 1:num_dims, reduction = "pca", reduction.name = "umap", reduction.key = "umap")
DimPlot(ALL_DATASETS, group.by = "Dataset", reduction = "umap", pt.size = 1)

ALL_DATASETS <- RunHarmony(object = ALL_DATASETS, reduction.use = "pca", group.by.vars = "Dataset", 
                           reduction.save = "harmonyPca", assay = "SCT", verbose = FALSE, normalization.method = "SCT")
ALL_DATASETS <- RunUMAP(ALL_DATASETS, dims = 1:num_dims, reduction = "harmonyPca", reduction.name = "harmony")
```


```R
p1 <- DimPlot(ALL_DATASETS, group.by = "Dataset", reduction = "harmony", pt.size = 1)
plot(p1)
ggsave("All_datasets_origin.pdf", width = 8, height = 8)
```


```R
DimPlot(ALL_DATASETS, group.by = "PHASE", reduction = "harmony", pt.size = 1)
```


```R
p2 <- DimPlot(ALL_DATASETS, group.by = "PREDICTED_PHASE", reduction = "harmony", pt.size = 1)
plot(p2)
ggsave("All_datasets_predPhase.pdf", width = 8, height = 8)
```

# Mounted Figures


```R
par(mar = c(3, 5, 1, 2) + 0.1)
options(repr.plot.width=22, repr.plot.height=15) 
p1_B <- ggplot(SPN_HD@meta.data, aes(x = EXP_CAGLENGTH)) +
#p1_B <- ggplot(SPN_filt@meta.data, aes(x = EXP_CAGLENGTH)) +
geom_histogram(binwidth = 5, alpha = 0.5, position = "identity", fill = "#F8766D") +
geom_vline(aes(xintercept = 150), colour="black", linetype = 2) +
theme_classic()

# p1_B
# ggsave("Figure_1B.pdf", height = 8, width = 8)

p1_B <- p1_B + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

#p1_C <- FeaturePlot(SPN_HD, reduction = "umap", features = "EXP_CAGLENGTH", cols = c("blue", "red"), pt.size = 1) + 
p1_C <- FeaturePlot(SPN_filt, reduction = "umap", features = "EXP_CAGLENGTH", cols = c("blue", "red"), pt.size = 1) + 
ggtitle("Handsaker Caudate - Expanded CAG size") + theme(plot.title = element_text(face = "plain", size = 15))

# p1_C
# ggsave("Figure_1C.pdf", height = 8, width = 8)

p1_C <- p1_C + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

#p1_D <- DimPlot(SPN_HD, reduction = "umap", group.by = "PHASE", pt.size = 1) + 
p1_D <- DimPlot(SPN_filt, reduction = "umap", group.by = "PHASE", pt.size = 1) + 
ggtitle("Handsaker Caudate - Phase") + theme(plot.title = element_text(face = "plain", size = 15))

# p1_D
# ggsave("Figure_1D.pdf", height = 8, width = 8)

p1_D <- p1_D + labs(tag = "D") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p1_E <- ggplot(Handsaker_fract_CDE_HD, aes(x = GRADE, y = FRACT_CDE, size = FRACT_SPN)) +
#p1_E <- ggplot(Handsaker_fract_CDE, aes(x = GRADE, y = FRACT_CDE, size = FRACT_SPN)) +
geom_point(aes(size = FRACT_SPN), color = "#F8766D", alpha = 0.8, position = position_jitter(width = 0, height = 0)) +
theme_classic() +
xlab("HD Grade") +
ylim(c(0, 25)) +
ylab("Fraction cells in C-D-E phase (%)")

# p1_E
# ggsave("Figure_1E.pdf", height = 8, width = 8)

p1_E <- p1_E + labs(tag = "E") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

Handsaker_donor_metadata_HD <- Handsaker_donor_metadata_full[which(Handsaker_donor_metadata_full$Status == "Case"), ]
p1_F <- ggplot(Handsaker_donor_metadata_HD, aes(x=CAP, y=100*FRACT_CDE, size = NUM_SPN)) + 
#p1_F <- ggplot(Handsaker_donor_metadata_full, aes(x=CAP, y=100*FRACT_CDE, size = NUM_SPN)) + 
  geom_point(color = "blue") +
  geom_smooth(method=lm, se=TRUE, mapping = aes(weight = NUM_SPN)) +
  theme_classic() +
  xlab("CAP score") +
  ylab("Fraction cells in C-D-E phase (%)")

# p1_F
# ggsave("Figure_1F.pdf", height = 8, width = 8)

p1_F <- p1_F + labs(tag = "F") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p1_G <- ggplot(Handsaker_donor_metadata_HD, aes(x=GERM_CAGEXP, y=100*FRACT_CDE, size = NUM_SPN)) + 
#p1_G <- ggplot(Handsaker_donor_metadata_full, aes(x=GERM_CAGEXP, y=100*FRACT_CDE, size = NUM_SPN)) + 
  geom_point(color = "blue") +
  geom_smooth(method=lm, se=TRUE, mapping = aes(weight = NUM_SPN)) +
  theme_classic() +
   xlab("Num. CAG germline") +
  ylab("Fraction cells in C-D-E phase (%)")

# p1_G
# ggsave("Figure_1G.pdf", height = 8, width = 8)

p1_G <- p1_G + labs(tag = "G") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p1 <- (p1_B | p1_C | p1_D) / (p1_E | p1_F | p1_G)
p1
ggsave("Figure_1_noA.pdf", width = 18, height = 14)
ggsave("Figure_1_noA.tiff", width = 18, height = 14)
```


```R
par(mar = c(3, 5, 1, 2) + 0.1)
options(repr.plot.width=22, repr.plot.height=15) 


# Data_filt2_tmp <- data.frame(CAG_measured = res_CAGsizing_test$S05202[which(res_CAGsizing_test$S05202$CAG_measured > 37), "CAG_measured"],
#                          CAG_predicted = res_CAGsizing_test$S05202[which(res_CAGsizing_test$S05202$CAG_measured > 37), "CAG_predicted"],
#                          PROB_PHASE_CDE = res_phase_test_noNA$S05202[, "lambda.min"],
#                          PRED_PHASE =  res_phase_test_noNA$S05202[, "predicted_phase_test_noNA"])
# Data_filt2 <- Data_filt2_tmp[order(Data_filt2_tmp$CAG_measured), ]
# Data_filt2$index = sort(Data_filt2$CAG_measured, index.return = TRUE, decreasing = FALSE)$ix

# ind_150 <- Data_filt2[which(Data_filt2$CAG_measured > 150), "index"][1]
# p2_A <- ggplot(Data_filt2, aes(x = index)) +
# geom_point(aes(y = CAG_predicted, color = PRED_PHASE, alpha = PROB_PHASE_CDE), size = 1) +
# geom_point(aes(y = CAG_measured), size = 0.1) +
# theme_classic() +
# geom_vline(aes(xintercept = ind_150), colour="black", linetype = 2) +
# xlab("SPNs sorted by CAG size") +
# ylab("Number of CAG repeats")

df_pr <- data.frame(
  recall = PR_curves_test_noNA$S05202$curve[, 1],
  precision = PR_curves_test_noNA$S05202$curve[, 2],
  cutoff = PR_curves_test_noNA$S05202$curve[, 3])

target_precision <- 0.99
idx <- which.min(abs(df_pr$precision - target_precision))
point_99 <- df_pr[idx, ]

baseline <- PR_curves_test_noNA$S05202$pos / (PR_curves_test_noNA$S05202$pos + PR_curves_test_noNA$S05202$neg)

p2_A <- ggplot(df_pr, aes(recall, precision, color = cutoff)) +
  geom_path(linewidth = 1) +
  geom_hline(
    yintercept = baseline,
    color = "blue",
    linewidth = 0.8
  ) +
  scale_color_gradientn(
    colours = rainbow(nrow(df_pr)),
    name = "Cutoff"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "Recall",
    y = "Precision"
  ) +
#    geom_point(
#     data = point_99,
#     aes(x = recall, y = precision),
#     color = "black",
#     size = 5
#   ) +
  theme_classic() +
ggtitle("Test set for held-out donor: S02205") +  theme(plot.title = element_text(face = "plain", size = 15))

# p2_A
# ggsave("Figure_2A.pdf", width = 8, height = 8)

p2_A <- p2_A + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1)) 

p2_B <- DimPlot(SPN_filt, reduction = "umap", group.by = "PREDICTED_PHASE", pt.size = 1) +
ggtitle("Handsaker Caudate - Predicted Phase") + theme(plot.title = element_text(face = "plain", size = 15))
#p2_B <- DimPlot(SPN_HD, reduction = "umap", group.by = "PREDICTED_PHASE", pt.size = 1) +

# p2_B
# ggsave("Figure_2B.pdf", height = 8, width = 8)

p2_B <- p2_B + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

#p2_C <- ggplot(Handsaker_donor_metadata_full, aes(x=100*FRACT_CDE, y=100*PRED_FRACT_CDE, size = NUM_SPN)) + 
p2_C <- ggplot(Handsaker_donor_metadata_HD, aes(x=100*FRACT_CDE, y=100*PRED_FRACT_CDE, size = NUM_SPN)) +
  geom_point(color = "blue") +
  geom_smooth(method=lm, se=TRUE, mapping = aes(weight = NUM_SPN)) +
  theme_classic() +
  xlab("Fraction cells in C-D-E phase (%)") +
  ylab("Predicted fraction cells in C-D-E phase (%)")

# p2_C
# ggsave("Figure_2C.pdf", height = 8, width = 8)

p2_C <- p2_C + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

# p2_C <- ggplot(df_allPhaseCgenes_sorted, aes(x = index, y = avg_expr_PhaseC, color = PHASE)) +
# geom_point(aes(x = index, y = avg_expr_PhaseC, color = PHASE), size = 1) +
# theme_classic() +
# geom_vline(aes(xintercept = ind_150_all), colour="black", linetype = 2) +
# xlab("SPNs sorted by CAG length") +
# ylab("Average expression Phase C genes") +
# labs(tag = "(C)") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))


#Data_filt_tmp <- data.frame(CAG_measured = training_test_data[, "CAGLENGTH_SPN_sized"], CAG_predicted = pred_all[, 1])
#Data_filt <- Data_filt_tmp[which(Data_filt_tmp$CAG_measured >= 150), ]

Data_filt_tmp <- data.frame(CAG_measured = res_CAGsizing_test_expMeas$S05202[, "CAG_measured"],
                            CAG_predicted = res_CAGsizing_test_expMeas$S05202[, "CAG_predicted"])
Data_filt <- Data_filt_tmp[which(Data_filt_tmp$CAG_measured >= 150), ]
             
# ind_noNA <- which(apply(results_test_phaseC_exp, 1, function(x) !any(is.na(x))))
# Data_filt <- results_test_phaseC_exp[ind_noNA, ]
df <- data.frame(CAG = c(Data_filt[, "CAG_measured"], Data_filt[, "CAG_predicted"]), method = c(rep("Measured", dim(Data_filt)[1]), rep("Predicted", dim(Data_filt)[1])))
p2_D <- ggplot(df, aes(x = CAG, fill = method, color = method)) +
  geom_density(alpha = 0.5, adjust = 1) +
  #geom_histogram(binwidth = 5, alpha = 0.5, position = "identity") +
  scale_color_manual(values = c("purple", "green")) +
  scale_fill_manual(values = c("purple", "green")) +
  xlab("Num. CAG") + ylab("Probability Density") +
  theme_classic() + 
ggtitle("Test set for held-out donor: S02205") +  theme(plot.title = element_text(face = "plain", size = 15))

# p2_D
# ggsave("Figure_2D.pdf", width = 8, height = 8)

p2_D <- p2_D + labs(tag = "D") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

quant_thr <- 0.999
p2_E <- as.ggplot(~ {
                par(mar = c(3, 5, 1, 2) + 0.1)
                smoothScatter(Data_filt[, "CAG_measured"], Data_filt[, "CAG_predicted"], xlab = "CAG measured", ylab = "CAG predicted", cex.main = 0.6, cex.axis = 0.8,
                xlim = c(0, quantile(c(Data_filt[, "CAG_measured"], Data_filt[, "CAG_predicted"]), quant_thr, na.rm = TRUE)), 
                ylim = c(0, quantile(c(Data_filt[, "CAG_measured"], Data_filt[, "CAG_predicted"]), quant_thr, na.rm = TRUE)))}) +
geom_abline(intercept = 0, slope = 1, color = "green") + 
ggtitle("Test set for held-out donor: S02205") +  theme(plot.title = element_text(face = "plain", size = 15))

# p2_E
# ggsave("Figure_2E.pdf", width = 8, height = 8)

p2_E <- p2_E + labs(tag = "E") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

#abline(0, 1, col = "green") 

# ind_150 <- df_sorted[which(df_sorted$CAGLENGTH_SPN_sized > 150), "index"][1]
# p2_F <- ggplot(df_sorted, aes(x = index)) +
# geom_point(aes(y = score, color = PRED_PHASE, alpha = PROB_PHASE_CDE), size = 1) +
# geom_point(aes(y = CAGLENGTH_SPN_sized), size = 0.1) +
# theme_classic() +
# geom_vline(aes(xintercept = ind_150), colour="black", linetype = 2) +
# xlab("SPNs sorted by CAG length") +
# ylab("Number of CAG repeats") + 
# labs(tag = "(F)") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

Data_filt2_tmp <- data.frame(CAG_measured = res_CAGsizing_test$S05202[which(res_CAGsizing_test$S05202$CAG_measured > 37), "CAG_measured"],
                         CAG_predicted = res_CAGsizing_test$S05202[which(res_CAGsizing_test$S05202$CAG_measured > 37), "CAG_predicted"],
                         PROB_PHASE_CDE = res_phase_test_noNA$S05202[, "lambda.min"],
                         PRED_PHASE =  res_phase_test_noNA$S05202[, "predicted_phase_test_noNA"])
Data_filt2 <- Data_filt2_tmp[order(Data_filt2_tmp$CAG_measured), ]
Data_filt2$index = sort(Data_filt2$CAG_measured, index.return = TRUE, decreasing = FALSE)$ix

ind_150 <- Data_filt2[which(Data_filt2$CAG_measured > 150), "index"][1]
p2_F <- ggplot(Data_filt2, aes(x = index)) +
geom_point(aes(y = CAG_predicted, color = PRED_PHASE, alpha = PROB_PHASE_CDE), size = 1) +
geom_point(aes(y = CAG_measured), size = 0.1) +
theme_classic() +
geom_vline(aes(xintercept = ind_150), colour="black", linetype = 2) +
xlab("SPNs sorted by CAG size") +
ylab("Number of CAG repeats") +
ggtitle("Test set for held-out donor: S02205") +  theme(plot.title = element_text(face = "plain", size = 15))

# p2_F
# ggsave("Figure_2F.pdf", width = 8, height = 8)

p2_F <- p2_F + labs(tag = "F") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p2 <- (p2_A | p2_B | p2_C) / (p2_D | p2_E | p2_F)
p2
ggsave("Figure_2.pdf", width = 18, height = 14)
ggsave("Figure_2.tiff", width = 18, height = 14)
```


```R
par(mar = c(3, 5, 1, 2) + 0.1)
options(repr.plot.width=22, repr.plot.height=15)

p3_A <- DimPlot(Paryani_HD_MSN_Caudate, group.by = "PREDICTED_PHASE", reduction = "harmony", pt.size = 1) + 
ggtitle("Paryani HD Caudate - Predicted Phase") +  theme(plot.title = element_text(face = "plain", size = 15))

# p3_A
# ggsave("Figure_3A.pdf", width = 8, height = 8)

p3_A <- p3_A + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p3_B <- ggplot(uni_Paryani_Caudate_Grade_CDE_HD, aes(x = GRADE, y = Fract_CDE, color = CELL_TYPE)) +
geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.1, height = 0)) +
scale_color_manual(values =c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
theme_classic() +
xlab("HD GRADE") +
#ylim(c(0, 100)) +
ylab("Fraction cells in C-D-E phase (%)") + 
ggtitle("Paryani HD Caudate") +  theme(plot.title = element_text(face = "plain", size = 15))

# p3_B
# ggsave("Figure_3B.pdf", width = 8, height = 8)

p3_B <- p3_B + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p3_C <- ggplot(fract_CDE_Paryani_Caudate, aes(x=num_CAG_germ, y=100*fract_CDE)) + 
  geom_point(col = "blue", aes(size = num_SPN))+
  geom_smooth(method=lm, se = TRUE, mapping = aes(weight = num_SPN)) +
  theme_classic() +
  xlab("Num. CAG germline") +
  ylab("Fraction cells in C-D-E phase (%)") +
ggtitle("Paryani HD Caudate") +  theme(plot.title = element_text(face = "plain", size = 15))

# p3_C
# ggsave("Figure_3C.pdf", width = 8, height = 8)

p3_C <- p3_C + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p3_D <- DimPlot(Lee_HD_MSN_Caudate, group.by = "PREDICTED_PHASE", reduction = "harmony", pt.size = 1) +
ggtitle("Lee HD Caudate - Predicted Phase") + theme(plot.title = element_text(face = "plain", size = 15))

# p3_D
# ggsave("Figure_3D.pdf", width = 8, height = 8)

p3_D <- p3_D + labs(tag = "D") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p3_E <- ggplot(uni_Lee_Caudate_Grade_CDE_HD, aes(x = GRADE, y = Fract_CDE, color = CELL_TYPE)) +
geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.1, height = 0)) +
scale_color_manual(values = c("#F8766D", "#619CFF")) +
theme_classic() +
xlab("HD GRADE") +
ylim(c(0, 100)) +
ylab("Fraction cells in C-D-E phase (%)")  + 
ggtitle("Lee HD Caudate") +  theme(plot.title = element_text(face = "plain", size = 15))

# p3_E
# ggsave("Figure_3E.pdf", width = 8, height = 8)

p3_E <- p3_E + labs(tag = "E") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p3_F <- ggplot(CTRL_all_datasets[which(CTRL_all_datasets$TISSUE == "Caudate"), ], aes(x = TISSUE, color = DATASET, y = 100*PRED_FRACT_CDE, size = Fract_SPN)) +
geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.3, height = 0)) +
#scale_color_manual(values =c("#F8766D")) +
theme_classic() +
xlab("Tissue") +
ylim(c(0, 100)) +
ylab("Non-HD donors - Fraction cells in C-D-E phase (%)") + 
ggtitle("Non-HD donors - Caudate") + theme(plot.title = element_text(face = "plain", size = 15))

# p3_F
# ggsave("Figure_3F.pdf", width = 8, height = 8)

p3_F <- p3_F + labs(tag = "F") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))                        

p3 <- (p3_A | p3_B | p3_C) / (p3_D | p3_E | p3_F)
p3
ggsave("Figure_3.pdf", width = 18, height = 14)
ggsave("Figure_3.tiff", width = 18, height = 14)
```


```R
par(mar = c(3, 5, 1, 2) + 0.1)
options(repr.plot.width=22, repr.plot.height=15) 
p4_A <- DimPlot(Paryani_HD_MSN_Accumbens, group.by = "PREDICTED_PHASE", reduction = "harmony", pt.size = 1) + 
ggtitle("Paryani HD Accumbens - Predicted Phase") +  theme(plot.title = element_text(face = "plain", size = 15))

# p4_A
# ggsave("Figure_4A.pdf", width = 8, height = 8)

p4_A <- p4_A + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p4_B <- ggplot(uni_Paryani_Accumbens_Grade_CDE_HD, aes(x = GRADE, y = Fract_CDE, color = CELL_TYPE)) +
geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.1, height = 0)) +
scale_color_manual(values =c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
theme_classic() +
xlab("HD GRADE") +
#ylim(c(0, 100)) +
ylab("Fraction cells in C-D-E phase (%)")

# p4_B
# ggsave("Figure_4B.pdf", width = 8, height = 8)

p4_B <- p4_B + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p4_C <- ggplot(fract_CDE_Paryani_Accumbens, aes(x=num_CAG_germ, y=100*fract_CDE)) + 
  geom_point(col = "blue", aes(size = num_SPN))+
  geom_smooth(method=lm, se = TRUE, mapping = aes(weight = num_SPN)) +
  theme_classic() +
  xlab("Num. CAG germline") +
  ylab("Fraction cells in C-D-E phase (%)")

# p4_C
# ggsave("Figure_4C.pdf", width = 8, height = 8)

p4_C <- p4_C + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p4_D <- DimPlot(Lee_HD_MSN_Putamen, group.by = "PREDICTED_PHASE", reduction = "harmony", pt.size = 1) + 
ggtitle("Lee HD Putamen - Predicted Phase") +  theme(plot.title = element_text(face = "plain", size = 15))

# p4_D
# ggsave("Figure_4D.pdf", width = 8, height = 8)

p4_D <- p4_D + labs(tag = "D") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p4_E <- ggplot(uni_Lee_Putamen_Grade_CDE_HD, aes(x = GRADE, y = Fract_CDE, color = CELL_TYPE)) +
geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.1, height = 0)) +
scale_color_manual(values = c("#F8766D", "#619CFF")) +
theme_classic() +
xlab("HD GRADE") +
ylim(c(0, 100)) +
ylab("Fraction cells in C-D-E phase (%)")

# p4_E
# ggsave("Figure_4E.pdf", width = 8, height = 8)

p4_E <- p4_E + labs(tag = "E") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p4_F <- ggplot(CTRL_all_datasets[which(CTRL_all_datasets$TISSUE != "Caudate"), ], aes(x = TISSUE, color = DATASET, y = 100*PRED_FRACT_CDE, size = Fract_SPN)) +
geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.3, height = 0)) +
#scale_color_manual(values =c("#F8766D")) +
theme_classic() +
xlab("Tissue") +
ylim(c(0, 100)) +
ylab("Non-HD donors - Fraction cells in C-D-E phase (%)") + 
ggtitle("Non-HD donors - Accumbens and Putamen") + theme(plot.title = element_text(face = "plain", size = 15))

# p4_F
# ggsave("Figure_4F.pdf", width = 8, height = 8)

p4_F <- p4_F + labs(tag = "F") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))
                        
p4 <- (p4_A | p4_B | p4_C) / (p4_D | p4_E | p4_F)
p4
ggsave("Figure_4.pdf", width = 18, height = 14)
ggsave("Figure_4.tiff", width = 18, height = 14)
```


```R
par(mar = c(3, 5, 1, 2) + 0.1)
options(repr.plot.width=22, repr.plot.height=15)

p5_A <- FeaturePlot(Xu_all, reduction = "harmony", features = "DRD1", cols = c("blue", "red"), pt.size = 1) + 
ggtitle("DRD1 expression") + theme(plot.title = element_text(face = "plain", size = 15))

# p5_A
# ggsave("Figure_5A.pdf", width = 8, height = 8)
p5_A <- p5_A + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p5_B <- FeaturePlot(Xu_all, reduction = "harmony", features = "DRD2", cols = c("blue", "red"), pt.size = 1) +
ggtitle("DRD2 expression") + theme(plot.title = element_text(face = "plain", size = 15))

# p5_B
# ggsave("Figure_5B.pdf", width = 8, height = 8)
p5_B <- p5_B + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p5_C <- DimPlot(Xu_all_SPN, reduction = "harmony", group.by = "SAMPLE", pt.size = 1) + 
ggtitle("Donor") + theme(plot.title = element_text(face = "plain", size = 15))

# p5_C
# ggsave("Figure_5C.pdf", width = 8, height = 8)
p5_C <- p5_C + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p5_D <- DimPlot(Xu_CTRL_SPN, reduction = "harmony", group.by = "PREDICTED_PHASE", pt.size = 1) +
ggtitle("CTRL - Predicted Phase") + theme(plot.title = element_text(face = "plain", size = 15))

# p5_D
# ggsave("Figure_5D.pdf", width = 8, height = 8)
p5_D <- p5_D + labs(tag = "D") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p5_E <- DimPlot(Xu_AD_SPN, reduction = "harmony", group.by = "PREDICTED_PHASE", pt.size = 1) +
ggtitle("AD - Predicted Phase") + theme(plot.title = element_text(face = "plain", size = 15))

# p5_E
# ggsave("Figure_5E.pdf", width = 8, height = 8)
p5_E <- p5_E + labs(tag = "E") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p5_F <- DimPlot(Xu_PD_SPN, reduction = "harmony", group.by = "PREDICTED_PHASE", pt.size = 1) +
ggtitle("PD - Predicted Phase") + theme(plot.title = element_text(face = "plain", size = 15))

# p5_F
# ggsave("Figure_5F.pdf", width = 8, height = 8)

p5_F <- p5_F + labs(tag = "F") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p5 <- (p5_A | p5_B | p5_C) / (p5_D | p5_E | p5_F)
p5
ggsave("Figure_5.pdf", width = 18, height = 14)
ggsave("Figure_5.tiff", width = 18, height = 14)
```


```R
par(mar = c(3, 5, 1, 2) + 0.1)
options(repr.plot.width=22, repr.plot.height=15) 
pS1_A <- DimPlot(SPN_all, group.by = "SAMPLE", reduction = "umap", pt.size = 1) +
ggtitle("Donor") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS1_A
# ggsave("Figure_S1A.pdf", width = 8, height = 8)
pS1_A <- pS1_A + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS1_B <- FeaturePlot(SPN_all, reduction = "umap", features = "nFeature_RNA", cols = c("blue", "red"), pt.size = 1) +
ggtitle("nFeature_RNA") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS1_B
# ggsave("Figure_S1B.pdf", width = 8, height = 8)
pS1_B <- pS1_B + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS1_C <- FeaturePlot(SPN_all, reduction = "umap", features = "nCount_RNA", cols = c("blue", "red"), pt.size = 1) + 
ggtitle("nCount_RNA") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS1_C
# ggsave("Figure_S1C.pdf", width = 8, height = 8)
pS1_C <- pS1_C + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS1_D <- FeaturePlot(SPN_all, reduction = "umap", features = "percent.mt", cols = c("blue", "red"), pt.size = 1) + 
ggtitle("percent.mt") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS1_D
# ggsave("Figure_S1D.pdf", width = 8, height = 8)
pS1_D <- pS1_D + labs(tag = "D") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS1_E <- FeaturePlot(SPN_all, reduction = "umap", features = "CAGLENGTH", cols = c("blue", "red"), pt.size = 1) + 
ggtitle("CAG size") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS1_E
# ggsave("Figure_S1E.pdf", width = 8, height = 8)
pS1_E <- pS1_E + labs(tag = "E") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS1_F <- DimPlot(SPN_all, reduction = "umap", group.by = "seurat_clusters", pt.size = 1) + 
ggtitle("Clusters") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS1_F
# ggsave("Figure_S1F.pdf", width = 8, height = 8)
pS1_F <- pS1_F + labs(tag = "F") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS1 <- (pS1_A | pS1_B | pS1_C) / (pS1_D | pS1_E | pS1_F)
pS1
ggsave("Figure_S1.pdf", width = 18, height = 14)
ggsave("Figure_S1.tiff", width = 18, height = 14)
```


```R
par(mar = c(3, 5, 1, 2) + 0.1)
options(repr.plot.width=22, repr.plot.height=15) 
pS3_A <- DimPlot(Paryani_MSN_Caudate, group.by = "Batch", reduction = "harmony", pt.size = 1) + 
ggtitle("Paryani Caudate - Batch") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS3_A
# ggsave("Figure_S3A.pdf", width = 8, height = 8)
pS3_A <- pS3_A + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS3_B <- DimPlot(Paryani_MSN_Caudate, group.by = "sub_type_4", reduction = "harmony", pt.size = 1) + 
ggtitle("Paryani Caudate - sub_type_4") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS3_B
# ggsave("Figure_S3B.pdf", width = 8, height = 8)
pS3_B <- pS3_B + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

Paryani_MSN_Caudate@meta.data$Grade <- factor(Paryani_MSN_Caudate@meta.data$Grade, levels = c("HDJ", "HD4", "HD3", "HD2", "HD1", "Control"))
pS3_C <- DimPlot(Paryani_MSN_Caudate, group.by = "Grade", reduction = "harmony", pt.size = 1) +
ggtitle("Paryani Caudate - Grade") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS3_C
# ggsave("Figure_S3C.pdf", width = 8, height = 8)
pS3_C <- pS3_C + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS3_D <- DimPlot(Lee_MSN_Caudate, group.by = "Batch", reduction = "harmony", pt.size = 1) +
ggtitle("Lee Caudate - Batch") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS3_D
# ggsave("Figure_S3D.pdf", width = 8, height = 8)
pS3_D <- pS3_D + labs(tag = "D") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS3_E <- DimPlot(Lee_MSN_Caudate, group.by = "CellType", reduction = "harmony", pt.size = 1) + 
ggtitle("Lee Caudate - Cell type") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS3_E
# ggsave("Figure_S3E.pdf", width = 8, height = 8)
pS3_E <- pS3_E + labs(tag = "E") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS3_F <- DimPlot(Lee_MSN_Caudate, group.by = "Grade", reduction = "harmony", pt.size = 1) + 
ggtitle("Lee Caudate - Grade") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS3_F
# ggsave("Figure_S3F.pdf", width = 8, height = 8)
pS3_F <- pS3_F + labs(tag = "F") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))
                        
pS3 <- (pS3_A | pS3_B | pS3_C) / (pS3_D | pS3_E | pS3_F)
pS3
ggsave("Figure_S3.pdf", width = 18, height = 14)
ggsave("Figure_S3.tiff", width = 18, height = 14)
```


```R
par(mar = c(3, 5, 1, 2) + 0.1)
options(repr.plot.width=22, repr.plot.height=15) 
pS4_A <- DimPlot(SPN_CTRL, reduction = "umap", group.by = "PREDICTED_PHASE", pt.size = 1) + 
ggtitle("Handsaker CTRL Caudate - Predicted Phase") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS4_A
# ggsave("Figure_S4A.pdf", width = 8, height = 8)
pS4_A <- pS4_A + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS4_B <- DimPlot(Paryani_CTRL_MSN_Caudate, group.by = "PREDICTED_PHASE", reduction = "harmony", pt.size = 1) + 
ggtitle("Paryani CTRL Caudate - Predicted Phase") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS4_B
# ggsave("Figure_S4B.pdf", width = 8, height = 8)
pS4_B <- pS4_B + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS4_C <- DimPlot(Lee_CTRL_MSN_Caudate, group.by = "PREDICTED_PHASE", reduction = "harmony", pt.size = 1) +
ggtitle("Lee CTRL Caudate - Predicted Phase") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS4_C
# ggsave("Figure_S4C.pdf", width = 8, height = 8)
pS4_C <- pS4_C + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS4_D <- DimPlot(Paryani_CTRL_MSN_Accumbens, group.by = "PREDICTED_PHASE", reduction = "harmony", pt.size = 1) +
ggtitle("Paryani CTRL Accumbens - Predicted Phase") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS4_D
# ggsave("Figure_S4D.pdf", width = 8, height = 8)
pS4_D <- pS4_D + labs(tag = "D") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS4_E <- DimPlot(Lee_CTRL_MSN_Putamen, group.by = "PREDICTED_PHASE", reduction = "harmony", pt.size = 1) + 
ggtitle("Lee CTRL Putamen - Predicted Phase") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS4_E
# ggsave("Figure_S4E.pdf", width = 8, height = 8)
pS4_E <- pS4_E + labs(tag = "E") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS4_F <- ggplot(CTRL_all_datasets, aes(x = TISSUE, color = DATASET, y = 100*PRED_FRACT_CDE, size = Fract_SPN)) +
geom_point(aes(size = Fract_SPN), alpha = 0.8, position = position_jitter(width = 0.3, height = 0)) +
#scale_color_manual(values =c("#F8766D")) +
theme_classic() +
xlab("Tissue") +
ylim(c(0, 100)) +
ylab("Non-HD donors - Fraction cells in C-D-E phase (%)") + 
ggtitle("Non-HD donors") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS4_F
# ggsave("Figure_S4F.pdf", width = 8, height = 8)
pS4_F <- pS4_F + labs(tag = "F") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS4 <- (pS4_A | pS4_B | pS4_C) / (pS4_D | pS4_E | pS4_F)
pS4
ggsave("Figure_S4.pdf", width = 18, height = 14)
ggsave("Figure_S4.tiff", width = 18, height = 14)
```


```R
par(mar = c(3, 5, 1, 2) + 0.1)
options(repr.plot.width=15, repr.plot.height=15) 

pS5_A <- FeaturePlot(SPN_HD, reduction = "umap", features = "PROB_PHASE_CDE", pt.size = 1) +
scale_color_gradient(limits = c(0, 1), low = "blue", high = "red", oob = scales::squish) +
ggtitle("Handaker et al. - Caudate HD") + theme(plot.title = element_text(face = "plain", size = 15))

# pS5_A
# ggsave("Figure_S5A.pdf", width = 8, height = 8)
pS5_A <- pS5_A + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS5_B <- FeaturePlot(SPN_CTRL, reduction = "umap", features = "PROB_PHASE_CDE", pt.size = 1) + 
scale_color_gradient(limits = c(0, 1), low = "blue", high = "red", oob = scales::squish) +
ggtitle("Handaker et al. - Caudate CTRL") + theme(plot.title = element_text(face = "plain", size = 15))

# pS5_B
# ggsave("Figure_S5B.pdf", width = 8, height = 8)
pS5_B <- pS5_B + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS5_C <- FeaturePlot(Paryani_HD_MSN_Caudate, reduction = "harmony", features = "PROB_PHASE_CDE", pt.size = 1) +
scale_color_gradient(limits = c(0, 1), low = "blue", high = "red", oob = scales::squish) +
ggtitle("Paryani et al. - Caudate HD") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS5_C
# ggsave("Figure_S5C.pdf", width = 8, height = 8)
pS5_C <- pS5_C + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS5_D <- FeaturePlot(Paryani_CTRL_MSN_Caudate, reduction = "harmony", features = "PROB_PHASE_CDE") +
scale_color_gradient(limits = c(0, 1), low = "blue", high = "red", oob = scales::squish) +
ggtitle("Paryani et al. - Caudate CTRL") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS5_D
# ggsave("Figure_S5D.pdf", width = 8, height = 8)
pS5_D <- pS5_D + labs(tag = "D") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS5_E <- FeaturePlot(Lee_HD_MSN_Caudate, reduction = "harmony", features = "PROB_PHASE_CDE", pt.size = 1) + 
scale_color_gradient(limits = c(0, 1), low = "blue", high = "red", oob = scales::squish) +
ggtitle("Lee et al. - Caudate HD") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS5_E
# ggsave("Figure_S5E.pdf", width = 8, height = 8)
pS5_E <- pS5_E + labs(tag = "E") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS5_F <- FeaturePlot(Lee_CTRL_MSN_Caudate, reduction = "harmony", features = "PROB_PHASE_CDE", pt.size = 1) + 
scale_color_gradient(limits = c(0, 1), low = "blue", high = "red", oob = scales::squish) +
ggtitle("Lee et al. - Caudate CTRL") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS5_F
# ggsave("Figure_S5F.pdf", width = 8, height = 8)
pS5_F <- pS5_F + labs(tag = "F") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS5 <- (pS5_A | pS5_B) / (pS5_C | pS5_D) / (pS5_E | pS5_F)

pS5
ggsave("Figure_S5.pdf", width = 18, height = 14)
ggsave("Figure_S5.tiff", width = 18, height = 14)
```


```R
par(mar = c(3, 5, 1, 2) + 0.1)
options(repr.plot.width=22, repr.plot.height=15) 
pS6_A <- DimPlot(Paryani_MSN_Accumbens, group.by = "Batch", reduction = "harmony", pt.size = 1) + 
ggtitle("Paryani Accumbens - Batch") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS6_A
# ggsave("Figure_S6A.pdf", width = 8, height = 8)
pS6_A <- pS6_A + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS6_B <- DimPlot(Paryani_MSN_Accumbens, group.by = "sub_type_4", reduction = "harmony", pt.size = 1) + 
ggtitle("Paryani Accumbens - sub_type_4") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS6_B
# ggsave("Figure_S6B.pdf", width = 8, height = 8)
pS6_B <- pS6_B + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

Paryani_MSN_Accumbens@meta.data$Grade <- factor(Paryani_MSN_Accumbens@meta.data$Grade, levels = c("HDJ", "HD4", "HD3", "HD2", "HD1", "Control"))
pS6_C <- DimPlot(Paryani_MSN_Accumbens, group.by = "Grade", reduction = "harmony", pt.size = 1) +
ggtitle("Paryani Accumbens - Grade") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS6_C
# ggsave("Figure_S6C.pdf", width = 8, height = 8)
pS6_C <- pS6_C + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS6_D <- DimPlot(Lee_MSN_Putamen, group.by = "Batch", reduction = "harmony", pt.size = 1) +
ggtitle("Lee Putamen - Batch") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS6_D
# ggsave("Figure_S6D.pdf", width = 8, height = 8)
pS6_D <- pS6_D + labs(tag = "D") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS6_E <- DimPlot(Lee_MSN_Putamen, group.by = "CellType", reduction = "harmony", pt.size = 1) + 
ggtitle("Lee Putamen - Cell type") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS6_E
# ggsave("Figure_S6E.pdf", width = 8, height = 8)
pS6_E <- pS6_E + labs(tag = "E") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))


pS6_F <- DimPlot(Lee_MSN_Putamen, group.by = "Grade", reduction = "harmony", pt.size = 1) + 
ggtitle("Lee Putamen - Grade") +  theme(plot.title = element_text(face = "plain", size = 15))
# pS6_F
# ggsave("Figure_S6F.pdf", width = 8, height = 8)

pS6_F <- pS6_F + labs(tag = "F") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))
                        
pS6 <- (pS6_A | pS6_B | pS6_C) / (pS6_D | pS6_E | pS6_F)
pS6
ggsave("Figure_S6.pdf", width = 18, height = 14)
ggsave("Figure_S6.tiff", width = 18, height = 14)
```


```R
#CTRL_all_datasets[which(CTRL_all_datasets$TISSUE == "Caudate"), ]
#1 - sum(CTRL_all_datasets[which(CTRL_all_datasets$TISSUE == "Caudate"), "NUM_SPN"]*CTRL_all_datasets[which(CTRL_all_datasets$TISSUE == "Caudate"), "PRED_FRACT_CDE"])/sum(CTRL_all_datasets[which(CTRL_all_datasets$TISSUE == "Caudate"), "NUM_SPN"])
```


```R
par(mar = c(3, 5, 1, 2) + 0.1)
options(repr.plot.width=15, repr.plot.height=15) 

pS7_A <- FeaturePlot(Paryani_HD_MSN_Accumbens, reduction = "harmony", features = "PROB_PHASE_CDE", pt.size = 1) +
scale_color_gradient(limits = c(0, 1), low = "blue", high = "red", oob = scales::squish) +
ggtitle("Paryani et al. - Accumbens HD") + theme(plot.title = element_text(face = "plain", size = 15))

# pS7_A
# ggsave("Figure_S7A.pdf", width = 8, height = 8)
pS7_A <- pS7_A + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS7_B <- FeaturePlot(Paryani_CTRL_MSN_Accumbens, reduction = "harmony", features = "PROB_PHASE_CDE") +
scale_color_gradient(limits = c(0, 1), low = "blue", high = "red", oob = scales::squish) +
ggtitle("Paryani et al. - Accumbens CTRL") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS7_B
# ggsave("Figure_S7B.pdf", width = 8, height = 8)
pS7_B <- pS7_B + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS7_C <- FeaturePlot(Lee_HD_MSN_Putamen, reduction = "harmony", features = "PROB_PHASE_CDE", pt.size = 1) + 
scale_color_gradient(limits = c(0, 1), low = "blue", high = "red", oob = scales::squish) +
ggtitle("Lee et al. - Putamen HD") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS7_C
# ggsave("Figure_S7C.pdf", width = 8, height = 8)
pS7_C <- pS7_C + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS7_D <- FeaturePlot(Lee_CTRL_MSN_Putamen, reduction = "harmony", features = "PROB_PHASE_CDE", pt.size = 1) + 
scale_color_gradient(limits = c(0, 1), low = "blue", high = "red", oob = scales::squish) +
ggtitle("Lee et al. - Putamen CTRL") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS7_D
# ggsave("Figure_S7D.pdf", width = 8, height = 8)
pS7_D <- pS7_D + labs(tag = "D") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS7_E <- FeaturePlot(Xu_CTRL_SPN, reduction = "harmony", features = "PROB_PHASE_CDE", pt.size = 1) + 
scale_color_gradient(limits = c(0, 1), low = "blue", high = "red", oob = scales::squish) +
ggtitle("Xu et al. - Putamen CTRL") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS7_E
# ggsave("Figure_S7E.pdf", width = 8, height = 8)
pS7_E <- pS7_E + labs(tag = "E") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS7_F <- FeaturePlot(Xu_AD_SPN, reduction = "harmony", features = "PROB_PHASE_CDE", pt.size = 1) + 
scale_color_gradient(limits = c(0, 1), low = "blue", high = "red", oob = scales::squish) +
ggtitle("Xu et al. - Putamen AD") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS7_F
# ggsave("Figure_S7F.pdf", width = 8, height = 8)
pS7_F <- pS7_F + labs(tag = "F") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS7_G <- FeaturePlot(Xu_PD_SPN, reduction = "harmony", features = "PROB_PHASE_CDE", pt.size = 1) + 
scale_color_gradient(limits = c(0, 1), low = "blue", high = "red", oob = scales::squish) +
ggtitle("Xu et al. - Putamen PD") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS7_G
# ggsave("Figure_S7G.pdf", width = 8, height = 8)
pS7_G <- pS7_G + labs(tag = "G") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS7 <- (pS7_A | pS7_B) / (pS7_C | pS7_D) / (pS7_E | pS7_F | pS7_G)

pS7
ggsave("Figure_S7.pdf", width = 18, height = 14)
ggsave("Figure_S7.tiff", width = 18, height = 14)
```


```R
probs_split_Lee_Putamen <- split(Lee_MSN_Putamen@meta.data[, c("PROB_PHASE_CDE", "Condition")], Lee_MSN_Putamen@meta.data$Condition)
probs_split_Lee_Putamen_df <- do.call(rbind, probs_split_Lee_Putamen)
probs_split_Lee_Putamen_df$Condition <- factor(probs_split_Lee_Putamen_df$Condition, levels = c("HD", "Control"))

p_Lee_Putamen <- ggplot(probs_split_Lee_Putamen_df, aes(x=Condition, y=PROB_PHASE_CDE, color = Condition, fill = Condition)) + 
  geom_violin() +
theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank()) +
xlab("Condition") +
ylab("Prob. Phase C-D-E") +
stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = list(c("HD", "Control")))

probs_split_Paryani_Accumbens <- split(Paryani_MSN_Accumbens@meta.data[, c("PROB_PHASE_CDE", "Condition")], Paryani_MSN_Accumbens@meta.data$Condition)
probs_split_Paryani_Accumbens_df <- do.call(rbind, probs_split_Paryani_Accumbens)
probs_split_Paryani_Accumbens_df[which(probs_split_Paryani_Accumbens_df$Condition == "Con"), "Condition"] <- "Control"
probs_split_Paryani_Accumbens_df$Condition <- factor(probs_split_Paryani_Accumbens_df$Condition, levels = c("HD", "Control"))

p_Paryani_Accumbens <- ggplot(probs_split_Paryani_Accumbens_df, aes(x=Condition, y=PROB_PHASE_CDE, color = Condition, fill = Condition)) + 
  geom_violin() +
theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank()) +
xlab("Condition") +
ylab("Prob. Phase C-D-E") +
stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = list(c("HD", "Control")))


probs_split_Xu <- split(Xu_all_SPN@meta.data[, c("PROB_PHASE_CDE", "condition")], Xu_all_SPN@meta.data$condition)
probs_split_Xu_df <- do.call(rbind, probs_split_Xu)
colnames(probs_split_Xu_df)[2] <- "Condition"
levels(probs_split_Xu_df$Condition)[which(levels(probs_split_Xu_df$Condition) == "CTRL")] <- "Control"
probs_split_Xu_df$Condition <- factor(probs_split_Xu_df$Condition, levels = c("AD", "PD", "Control"))

p_Xu <- ggplot(probs_split_Xu_df, aes(x=Condition, y=PROB_PHASE_CDE, color = Condition, fill = Condition)) + 
  geom_violin() +
theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank()) +
xlab("Condition") +
ylab("Prob. Phase C-D-E") +
stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = list(c("PD", "Control"), c("AD", "Control")))
```


```R
pS8_A <- p_Paryani_Accumbens + ggtitle("Paryani et al. - Accumbens") +  theme(plot.title = element_text(face = "plain", size = 15))
pS8_A
ggsave("Figure_S8A.pdf", width = 8, height = 8)
pS8_A <- pS8_A + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

pS8_B <- p_Lee_Putamen + ggtitle("Lee et al. - Putamen") +  theme(plot.title = element_text(face = "plain", size = 15))
pS8_B
ggsave("Figure_S8B.pdf", width = 8, height = 8)
pS8_B <- pS8_B + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))


pS8_C <- p_Xu + ggtitle("Xu et al. - Putamen ") +  theme(plot.title = element_text(face = "plain", size = 15))
pS8_C
ggsave("Figure_S8C.pdf", width = 8, height = 8)
pS8_C <- pS8_C + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))

p_S8 <- pS8_A | pS8_B | pS8_C

p_S8
ggsave("Figure_S8.pdf", width = 12, height = 8)
ggsave("Figure_S8.tiff", width = 12, height = 8)
```


```R
wilcox.test(x = probs_split_Lee_Putamen_df[which(probs_split_Lee_Putamen_df$Condition == "Control"), "PROB_PHASE_CDE"],
       y = probs_split_Lee_Putamen_df[which(probs_split_Lee_Putamen_df$Condition == "HD"), "PROB_PHASE_CDE"])

wilcox.test(x = probs_split_Paryani_Accumbens_df[which(probs_split_Paryani_Accumbens_df$Condition == "Control"), "PROB_PHASE_CDE"],
       y = probs_split_Paryani_Accumbens_df[which(probs_split_Paryani_Accumbens_df$Condition == "HD"), "PROB_PHASE_CDE"])

wilcox.test(x = probs_split_Xu_df[which(probs_split_Xu_df$Condition == "Control"), "PROB_PHASE_CDE"],
       y = probs_split_Xu_df[which(probs_split_Xu_df$Condition == "AD"), "PROB_PHASE_CDE"])

wilcox.test(x = probs_split_Xu_df[which(probs_split_Xu_df$Condition == "Control"), "PROB_PHASE_CDE"],
       y = probs_split_Xu_df[which(probs_split_Xu_df$Condition == "PD"), "PROB_PHASE_CDE"])
```


```R
par(mar = c(3, 5, 1, 2) + 0.1)
options(repr.plot.width=15, repr.plot.height=10)
pS9_A <- DimPlot(ALL_DATASETS, group.by = "Dataset", reduction = "harmony", pt.size = 1) +
ggtitle("All datasets") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS9_A
# ggsave("Figure_S9A.pdf", width = 8, height = 8)
pS9_A <- pS9_A + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))


pS9_B <- DimPlot(ALL_DATASETS, group.by = "PREDICTED_PHASE", reduction = "harmony", pt.size = 1) +
ggtitle("All datasets - Predicted Phase") +  theme(plot.title = element_text(face = "plain", size = 15))

# pS9_B
# ggsave("Figure_S9B.pdf", width = 8, height = 8)
pS9_B <- pS9_B + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold", vjust = 1))
                        
p_S9 <- pS9_A | pS9_B
p_S9
ggsave("Figure_S9.pdf", width = 12, height = 8)
ggsave("Figure_S9.tiff", width = 12, height = 8)
```
