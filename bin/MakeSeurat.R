#!/opt/conda/envs/dolphinnext/bin/Rscript

# ╔══════════════════════════════════════════════════════════════════════════════════════════════╗
# ╠═                                     Title: MakeSeurat.R                                    ═╣
# ╠═                                     Updated: May 31 2024                                   ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═                                       nf-core/scscape                                      ═╣
# ╠═                                  MDI Biological Laboratory                                 ═╣
# ╠═                          Comparative Genomics and Data Science Core                        ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═ Description:                                                                               ═╣
# ╠═     Takes the raw count data (features, barcodes, matrix) and generates a Seurat           ═╣
# ╠═     Object. Options to subset features are avalible.                                       ═╣
# ╚══════════════════════════════════════════════════════════════════════════════════════════════╝
# ╔══════════════════╗
# ╠═ Load Libraries ═╣
# ╚══════════════════╝
library(stringr)
library(dplyr)
library(Matrix)
library(Seurat)
library(SeuratObject)


# ╔══════════════════════╗
# ╠═ Read in Parameters ═╣
# ╚══════════════════════╝
message("Reading in Parameters")
# Read in trailing arguments
args <- commandArgs(trailingOnly = TRUE)

# File path to features, barcodes, mtx directory
params.data_directory <- args[1]

# Which column to use for Seurat Obj. instantiation
params.gene_identifier <- args[2]
if (toupper(params.gene_identifier) == "GENE_ID"){
    params.gene_column <- 1
}else if (toupper(params.gene_identifier) == "GENE_NAME" || toupper(params.gene_identifier) == "COMBINED"){
    params.gene_column <- 2
}else {
    params.gene_column <- 2
}

# File path for target gene list
message1 <- NULL
params.genes_2_rm <- NULL
tryCatch({
    params.genes_2_rm <- read.csv(args[3])$RMgenes
    params.genes_2_rm <- params.genes_2_rm[!(params.genes_2_rm == "" | is.na(params.genes_2_rm))]
}, error = function(e) {
    message1 <- "NOTE:MAKE_SEURAT:Gene List File not Provided -- Skipping Custom Gene Removal"
    message(message1)
})

# Sample Name
params.sample_name <- args[4]

# Min Cells
params.min_cells <- as.integer(args[5])

# Min Features
params.min_features <- as.integer(args[6])

# Meta Data
params.meta_data <- args[7]

message(params.meta_data)

# ╔═══════════════╗
# ╠═ META Object ═╣
# ╚═══════════════╝

toDict <- function(input_string) {
    # Remove the outer brackets
    input_string <- gsub("^\\[|\\]$", "", input_string)

    # Split the string by commas, but not commas inside brackets
    parts <- strsplit(input_string, ",\\s*(?![^\\[\\]]*\\])", perl = TRUE)[[1]]

    result <- list()

    for (part in parts) {
        # Check if it contains a nested list
        if (grepl("\\[.*\\]", part)) {
            # Extract key and value separately
            key <- gsub("^(.*?):\\[.*$", "\\1", part)
            nested_value <- gsub("^.*?:\\[(.*?)\\]$", "\\1", part)

            # Split nested values
            nested_list <- strsplit(nested_value, ",\\s*")[[1]]
            result[[key]] <- nested_list
        } else {
            # Regular key-value pair
            key_value <- strsplit(part, ":")[[1]]
            if (length(key_value) == 2) {
                key <- trimws(key_value[1])
                value <- trimws(key_value[2])
                result[[key]] <- value
            }
        }
    }

    return(result)
}

meta <- toDict(params.meta_data)

# ╔═══════════════════╗
# ╠═ Subset Features ═╣
# ╚═══════════════════╝
message("Subsetting Features")
if (length(params.genes_2_rm) == 0){
  if (substr(params.data_directory, nchar(params.data_directory), nchar(params.data_directory)) == "/"){
    feature_list  <- read.csv(paste0(params.data_directory, "features.tsv.gz"), sep = "\t", header = F)
    new_feat_list <- feature_list
  }else {
    feature_list  <- read.csv(paste0(params.data_directory, "/features.tsv.gz"), sep = "\t", header = F)
    new_feat_list <- feature_list
  }
}else{
  if (substr(params.data_directory, nchar(params.data_directory), nchar(params.data_directory)) == "/"){
    gene_list     <- params.genes_2_rm
    feature_list  <- read.csv(paste0(params.data_directory, "features.tsv.gz"), sep = "\t", header = F)
    if (params.gene_column == 1){
        new_feat_list <- feature_list[-c(which(feature_list$V1 %in% gene_list)),]
    }else{
        new_feat_list <- feature_list[-c(which(feature_list$V2 %in% gene_list)),]
    }
  }else{
    gene_list     <- params.genes_2_rm
    feature_list  <- read.csv(paste0(params.data_directory, "/features.tsv.gz"), sep = "\t", header = F)
    if (params.gene_column == 1){
        new_feat_list <- feature_list[-c(which(feature_list$V1 %in% gene_list)),]
    }else{
        new_feat_list <- feature_list[-c(which(feature_list$V2 %in% gene_list)),]
    }
  }
}

# ╔════════════════════════════════╗
# ╠═ Create 10X Object and Subset ═╣
# ╚════════════════════════════════╝
message("Creating 10X Object and Subsetting")
# Create raw 10X object
Name10X <- params.sample_name
assign(Name10X, Read10X(data.dir = params.data_directory, strip.suffix = T, gene.column = params.gene_column))

# Subset 10X Object
Name10XAnnotated <- paste0(params.sample_name,"_Ann")
assign(Name10XAnnotated, get(params.sample_name)[which(rownames(get(params.sample_name)) %in% new_feat_list[,params.gene_column]),])

# ╔════════════════════════╗
# ╠═ Create Seurat Object ═╣
# ╚════════════════════════╝
message("Creating Seurat Object")
NameSO <- paste0("SO_",params.sample_name)
assign(NameSO, CreateSeuratObject(counts = get(Name10XAnnotated), project = params.sample_name, min.cells = params.min_cells, min.features = params.min_features))

# ╔═════════════════╗
# ╠═ Add Meta Data ═╣
# ╚═════════════════╝
message("Adding Meta Data")
SO <- get(NameSO)
for (name in names(meta)){
    SO[[name]] <- meta[[name]]
}
assign(NameSO, SO)

# ╔══════════════════════╗
# ╠═ Save Seurat Object ═╣
# ╚══════════════════════╝
message("Saving Seurat Object")
SaveSeuratRds(get(NameSO), file = paste0("00_",params.sample_name, "_InitalSO.rds"))


# ╔═════════════════╗
# ╠═ Save Log File ═╣
# ╚═════════════════╝
sink(paste0("00_",params.sample_name,"_InitialValidation.log"))
cat("╔══════════════════════════════════════════════════════════════════════════════════════════════╗\n")
cat("╠  MakeSeurat.R Validation log\n")
cat(paste0("╠  Sample: ", params.sample_name,"\n"))
cat("╚══════════════════════════════════════════════════════════════════════════════════════════════╝\n")
cat(paste0("Gene Identifier: ", params.gene_identifier,"\n"))
cat(paste0("Minimum Cells: ", params.min_cells,"\n"))
cat(paste0("Min Features: ", params.min_features,"\n"))
cat("\n")
cat("Seurat Object Status:\n")
print(get(NameSO))
cat("\n")
if(length(message1)!= 0){
    cat(message1)
    cat("\n")
}
sink()

sink(paste0("00_",params.sample_name,"_InitialVersions.log"))
cat("╔══════════════════════════════════════════════════════════════════════════════════════════════╗\n")
cat("╠  MakeSeurat.R Versions\n")
cat(paste0("╠  Sample: ", params.sample_name,"\n"))
cat("╚══════════════════════════════════════════════════════════════════════════════════════════════╝\n")
sessionInfo()
sink()
