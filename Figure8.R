library(readxl)
library(cowplot)
library(dplyr)
library(ggplot2)
library(magrittr)
library(MAST)
library(Matrix)
library(monocle3)
library(patchwork)
library(plyr)
library(RCurl)
library(reshape2)
library(scales)
library(scater)
library(Seurat)
library(SingleCellExperiment)
library(stringi)
library(stringr)
library(tidyverse)

#read the raw data
mat <- readMM(file = 'path/filtered_feature_bc_matrix/matrix.mtx.gz')
feature.names = read.delim('path/filtered_feature_bc_matrix/features.tsv.gz', header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim('/path/barcodes.tsv', header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2
df <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 100)

#assign cells base on their sex, condition
metadata <- df@meta.data
metadata$Sample.ID <- NA
metadata$Sample.ID[which(str_detect(row.names(metadata), "-1$"))] <- "WT_Female_1"
metadata$Sample.ID[which(str_detect(row.names(metadata), "-2$"))] <- "WT_Female_2"
metadata$Sample.ID[which(str_detect(row.names(metadata), "-3$"))] <- "WT_Female_3"
metadata$Sample.ID[which(str_detect(row.names(metadata), "-4$"))] <- "KO_Female_1"
metadata$Sample.ID[which(str_detect(row.names(metadata), "-5$"))] <- "KO_Female_2"
metadata$Sample.ID[which(str_detect(row.names(metadata), "-6$"))] <- "KO_Female_3"
metadata$Sample.ID[which(str_detect(row.names(metadata), "-7$"))] <- "WT_Male_1"
metadata$Sample.ID[which(str_detect(row.names(metadata), "-8$"))] <- "WT_Male_2"
metadata$Sample.ID[which(str_detect(row.names(metadata), "-9$"))] <- "WT_Male_3"
metadata$Sample.ID[which(str_detect(row.names(metadata), "-10$"))] <- "KO_Male_1"
metadata$Sample.ID[which(str_detect(row.names(metadata), "-11$"))] <- "KO_Male_2"
metadata$Sample.ID[which(str_detect(row.names(metadata), "-12$"))] <- "KO_Male_3"
metadata$Genotype <- NA
metadata$Genotype[which(str_detect(row.names(metadata), "-1$"))] <- "WT_Female"
metadata$Genotype[which(str_detect(row.names(metadata), "-2$"))] <- "WT_Female"
metadata$Genotype[which(str_detect(row.names(metadata), "-3$"))] <- "WT_Female"
metadata$Genotype[which(str_detect(row.names(metadata), "-4$"))] <- "KO_Female"
metadata$Genotype[which(str_detect(row.names(metadata), "-5$"))] <- "KO_Female"
metadata$Genotype[which(str_detect(row.names(metadata), "-6$"))] <- "KO_Female"
metadata$Genotype[which(str_detect(row.names(metadata), "-7$"))] <- "WT_Male"
metadata$Genotype[which(str_detect(row.names(metadata), "-8$"))] <- "WT_Male"
metadata$Genotype[which(str_detect(row.names(metadata), "-9$"))] <- "WT_Male"
metadata$Genotype[which(str_detect(row.names(metadata), "-10$"))] <- "KO_Male"
metadata$Genotype[which(str_detect(row.names(metadata), "-11$"))] <- "KO_Male"
metadata$Genotype[which(str_detect(row.names(metadata), "-12$"))] <- "KO_Male"
metadata$Condition <- NA
metadata$Condition[which(str_detect(row.names(metadata), "-1$"))] <- "WT"
metadata$Condition[which(str_detect(row.names(metadata), "-2$"))] <- "WT"
metadata$Condition[which(str_detect(row.names(metadata), "-3$"))] <- "WT"
metadata$Condition[which(str_detect(row.names(metadata), "-4$"))] <- "KO"
metadata$Condition[which(str_detect(row.names(metadata), "-5$"))] <- "KO"
metadata$Condition[which(str_detect(row.names(metadata), "-6$"))] <- "KO"
metadata$Condition[which(str_detect(row.names(metadata), "-7$"))] <- "WT"
metadata$Condition[which(str_detect(row.names(metadata), "-8$"))] <- "WT"
metadata$Condition[which(str_detect(row.names(metadata), "-9$"))] <- "WT"
metadata$Condition[which(str_detect(row.names(metadata), "-10$"))] <- "KO"
metadata$Condition[which(str_detect(row.names(metadata), "-11$"))] <- "KO"
metadata$Condition[which(str_detect(row.names(metadata), "-12$"))] <- "KO"
metadata$Sex <- NA
metadata$Sex[which(str_detect(row.names(metadata), "-1$"))] <- "Female"
metadata$Sex[which(str_detect(row.names(metadata), "-2$"))] <- "Female"
metadata$Sex[which(str_detect(row.names(metadata), "-3$"))] <- "Female"
metadata$Sex[which(str_detect(row.names(metadata), "-4$"))] <- "Female"
metadata$Sex[which(str_detect(row.names(metadata), "-5$"))] <- "Female"
metadata$Sex[which(str_detect(row.names(metadata), "-6$"))] <- "Female"
metadata$Sex[which(str_detect(row.names(metadata), "-7$"))] <- "Male"
metadata$Sex[which(str_detect(row.names(metadata), "-8$"))] <- "Male"
metadata$Sex[which(str_detect(row.names(metadata), "-9$"))] <- "Male"
metadata$Sex[which(str_detect(row.names(metadata), "-10$"))] <- "Male"
metadata$Sex[which(str_detect(row.names(metadata), "-11$"))] <- "Male"
metadata$Sex[which(str_detect(row.names(metadata), "-12$"))] <- "Male"
df@meta.data <- metadata
df$log10GenesPerUMI <- log10(df$nFeature_RNA) / log10(df$nCount_RNA)

# Add mitochondrial percenatge
mouse_mito = read_excel("Mouse.MitoCarta2.0.xls", sheet = 2)
mouse_mito = mouse_mito %>% select(c(Symbol, MCARTA2.0_score)) 
mito.genes = as.character(mouse_mito$Symbol)
mito.genes = mito.genes[mito.genes %in% rownames(df@assays$RNA)]
df$percent.mt <- Matrix::colSums(df@assays$RNA[mito.genes,]) / Matrix::colSums(df@assays$RNA)
# Add ribosomal RNA content
ribo.genes <- grep(pattern = "^Rp[l|s]", x = rownames(df), value = TRUE)
df$percent.rp <- PercentageFeatureSet(df, features = ribo.genes)
metadata <- df@meta.data
# filter for qualified cells 
df <- subset(df, subset = nFeature_RNA > 200 & 
               nFeature_RNA < 6000&
               percent.mt < 0.2 &
               nCount_RNA > 1000)

# Save the data
saveRDS(df, file="project.rds")

# Read in the cell barcodes annotated from the 10x Genomics Loupe Browser 
dc <- read.csv('Dendritic Cells.csv')
mac1 <- read.csv('Mac 1.csv')
mac2 <- read.csv('Mac 2.csv')
mac3 <- read.csv('Mac 3.csv')
mac4 <- read.csv('Mac 4.csv')
mono <- read.csv('Monocytes.csv')
cellType <- rbind(dc, mac1, mac2, mac3, mac4, mono)

df <- readRDS("project.rds.rds")
DefaultAssay(object = dat0) <- "RNA"

### add cluster info to the single cell data
meta <- df@meta.data
meta$Barcode <- row.names(meta)
meta <- left_join(meta, cellType, by = "Barcode")
row.names(meta) <- meta$Barcode
df@meta.data <- meta
dat1 <- subset(df, CellType != c(''))

# WT_Female_Monocle3_Traject
dat2 <- subset(dat1, Sex == "Female")
dat2 <- subset(dat2, Condition == "WT")

exp_mtx <- as.matrix(GetAssayData(dat2[["RNA"]], slot = "counts"))
cell_meta <- dat2@meta.data
gene_ano <- data.frame(gene_short_name=row.names(exp_mtx))
rownames(gene_ano) <- gene_ano$gene_short_name
cds <- new_cell_data_set(exp_mtx,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_ano)

cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "Sample.ID")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

dpi = 300
png(file = "Figure8a_WT_Female.png", width = dpi * 6,height = dpi * 4,units = "px",res = dpi,type = 'cairo')
plot_cells(cds,
           color_cells_by = "CellType",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           label_principal_points = TRUE,
           graph_label_size=1.5) 
dev.off()

# KO_Female_Monocle3_Traject
dat2 <- subset(dat1, Sex == "Female")
dat2 <- subset(dat2, Condition == "KO")
exp_mtx <- as.matrix(GetAssayData(dat2[["RNA"]], slot = "counts"))
cell_meta <- dat2@meta.data
gene_ano <- data.frame(gene_short_name=row.names(exp_mtx))
rownames(gene_ano) <- gene_ano$gene_short_name
cds <- new_cell_data_set(exp_mtx,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_ano)

cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "Sample.ID")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

dpi = 300
png(file = "Figure8a_KO_Female.png", width = dpi * 9,height = dpi * 6,units = "px",res = dpi,type = 'cairo')
plot_cells(cds,
           color_cells_by = "CellType",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           label_principal_points = TRUE,
           graph_label_size=1.5) +
  scale_color_manual(values = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f'))
dev.off()

# WT_Male_Monocle3_Traject
dat2 <- subset(dat1, Sex == "Male")
dat2 <- subset(dat2, Condition == "WT")

exp_mtx <- as.matrix(GetAssayData(dat2[["RNA"]], slot = "counts"))
cell_meta <- dat2@meta.data
gene_ano <- data.frame(gene_short_name=row.names(exp_mtx))
rownames(gene_ano) <- gene_ano$gene_short_name
cds <- new_cell_data_set(exp_mtx,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_ano)

cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "Sample.ID")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

dpi = 300
png(file = "Figure8a_WT_Male.png", width = dpi * 6,height = dpi * 4,units = "px",res = dpi,type = 'cairo')
plot_cells(cds,
           color_cells_by = "CellType",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           label_principal_points = TRUE,
           graph_label_size=1.5) 
dev.off()

# WT_Male_Monocle3_Traject
dat2 <- subset(dat1, Sex == "Male")
dat2 <- subset(dat2, Condition == "KO")

exp_mtx <- as.matrix(GetAssayData(dat2[["RNA"]], slot = "counts"))
cell_meta <- dat2@meta.data
gene_ano <- data.frame(gene_short_name=row.names(exp_mtx))
rownames(gene_ano) <- gene_ano$gene_short_name
cds <- new_cell_data_set(exp_mtx,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_ano)

cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "Sample.ID")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

dpi = 300
png(file = "Figure8a_KO_Male.png", width = dpi * 6,height = dpi * 4,units = "px",res = dpi,type = 'cairo')
plot_cells(cds,
           color_cells_by = "CellType",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           label_principal_points = TRUE,
           graph_label_size=1.5) 
dev.off()




