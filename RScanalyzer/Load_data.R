# Load_data.R

# Loading single-cell RNA-seq data and creating a Seurat object
#
# Supports:
#  - 10X (directory produced by CellRanger / 10X Genomics)
#  - MTX format (matrix.mtx + features.tsv + barcodes.tsv)
#  - CSV/TSV count matrix (genes as rows, cells as columns)
#
# Features:
#  - optional metadata file merge (CSV/TSV where first column is cell barcodes)
#  - calculates percent mitochondrial reads (robust to MT- and mt- prefixes)
#  - optional filtering by min.cells, min.features, max.percent.mt
#  - runs NormalizeData, FindVariableFeatures, ScaleData, PCA and optional UMAP + clustering
#  - saves Seurat object to an RDS file

# Usage examples:
# Rscript Load_data.R --input /path/to/10x_dir --format 10x --project MySample --out sample.rds
# Rscript Load_data.R --input /path/to/mtx_dir --format mtx --out sample.rds
# Rscript Load_data.R --input /path/to/counts.csv --format csv --metadata meta.csv --out sample.rds

library(Seurat)
library(Matrix)
library(readr)
library(dplyr)

# Simple argument parser (no external packages required)
args <- commandArgs(trailingOnly = TRUE)
arg_list <- list()
for (i in seq(1, length(args), by = 2)) {
	name <- sub('^--', '', args[i])
	value <- ifelse(i+1 <= length(args), args[i+1], TRUE)
	arg_list[[name]] <- value
}

get_arg <- function(name, default = NULL) {
	if (!is.null(arg_list[[name]])) return(arg_list[[name]])
	return(default)
}

input_path <- get_arg('input', NULL)
fmt <- tolower(get_arg('format', '10x'))
project_name <- get_arg('project', 'scProject')
metadata_file <- get_arg('metadata', NULL)
out_file <- get_arg('out', file.path(getwd(), paste0(project_name, '.rds')))
min_cells <- as.integer(get_arg('min.cells', '3'))
min_features <- as.integer(get_arg('min.features', '200'))
max_percent_mt <- as.numeric(get_arg('max.percent.mt', '20'))
run_umap <- as.logical(get_arg('umap', 'TRUE'))
run_clustering <- as.logical(get_arg('cluster', 'TRUE'))

if (is.null(input_path)) {
	stop('Please provide --input <path> and --format <10x|mtx|csv>')
}

message('Loading data...')
counts <- NULL
if (fmt == '10x') {
	# Read10X will detect features/barcodes in the input directory
	counts <- Read10X(data.dir = input_path)
} else if (fmt %in% c('mtx', 'matrix')) {
	# Expect matrix.mtx, features.tsv (or genes.tsv), barcodes.tsv in input_path
	mtx <- file.path(input_path, 'matrix.mtx')
	features <- file.path(input_path, 'features.tsv')
	barcodes <- file.path(input_path, 'barcodes.tsv')
	if (!file.exists(mtx)) stop('matrix.mtx not found in the provided directory')
	if (!file.exists(features)) stop('features.tsv not found in the provided directory')
	if (!file.exists(barcodes)) stop('barcodes.tsv not found in the provided directory')
	counts <- ReadMtx(mtx = mtx, features = features, cells = barcodes)
} else if (fmt == 'csv' || fmt == 'tsv' || fmt == 'txt') {
	# csv/tsv: first column genes, header contains cell barcodes
	if (!file.exists(input_path)) stop('Count file not found: ', input_path)
	# attempt to detect separator
	if (fmt == 'csv' || grepl('\\.csv$', input_path, ignore.case = TRUE)) {
		mat <- read_csv(input_path)
	} else {
		mat <- read_delim(input_path, delim = '\\t')
	}
	# assume first column is gene names
	gene_col <- colnames(mat)[1]
	rownames(mat) <- mat[[gene_col]]
	mat[[gene_col]] <- NULL
	counts <- as.matrix(mat)
} else {
	stop('Unsupported format: ', fmt, '. Use 10x, mtx, or csv')
}

message('Creating Seurat object...')
seurat_obj <- CreateSeuratObject(counts = counts, project = project_name,
																 min.cells = min_cells, min.features = min_features)

# Merge metadata if provided
if (!is.null(metadata_file)) {
	if (!file.exists(metadata_file)) stop('Metadata file not found: ', metadata_file)
	meta <- read_csv(metadata_file)
	# If metadata first column is not named, try to detect it
	if (ncol(meta) >= 1) {
		# assume first column contains cell barcodes (or rownames)
		meta_col1 <- colnames(meta)[1]
		rownames(meta) <- as.character(meta[[meta_col1]])
		meta[[meta_col1]] <- NULL
		# subset to cells present in Seurat object
		common <- intersect(colnames(seurat_obj), rownames(meta))
		if (length(common) == 0) stop('No matching cell barcodes between counts and metadata')
		meta_subset <- meta[common, , drop = FALSE]
		seurat_obj <- AddMetaData(seurat_obj, metadata = meta_subset)
	}
}

# Calculate percent mitochondrial reads (robust to MT-/mt- prefixes)
seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = '^MT-|^mt-')

message('Initial cell count: ', ncol(seurat_obj), ' ; feature count: ', nrow(seurat_obj))

# Apply filtering based on features and percent.mt
if (!is.na(max_percent_mt) && max_percent_mt > 0) {
	seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= min_features & percent.mt <= max_percent_mt)
} else {
	seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= min_features)
}

message('Post-filter cell count: ', ncol(seurat_obj))

# Basic preprocessing workflow
seurat_obj <- NormalizeData(seurat_obj, normalization.method = 'LogNormalize', scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = 'vst', nfeatures = 2000)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

if (isTRUE(run_umap) || isTRUE(run_clustering)) {
	seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
}
if (isTRUE(run_clustering)) {
	seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
}
if (isTRUE(run_umap)) {
	seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
}

message('Saving Seurat object to: ', out_file)
saveRDS(seurat_obj, file = out_file)

message('Done.')


