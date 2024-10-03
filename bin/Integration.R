#!/opt/conda/envs/dolphinnext/bin/Rscript

# ╔══════════════════════════════════════════════════════════════════════════════════════════════╗
# ╠═                                     Title: Integration.R                                   ═╣
# ╠═                                     Updated: May 31 2024                                   ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═                                       nf-core/scscape                                      ═╣
# ╠═                                  MDI Biological Laboratory                                 ═╣
# ╠═                          Comparative Genomics and Data Science Core                        ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═ Description:                                                                               ═╣
# ╠═     Takes in the Seurat Object generated from RunPCA.R. Runs integration on the Seurat     ═╣
# ╠═     Object if the Integration Method parameter is not set to NULL. Otherwise, this script  ═╣
# ╠═     is skipped.                                                                            ═╣
# ╚══════════════════════════════════════════════════════════════════════════════════════════════╝
# ╔══════════════════╗
# ╠═ Load Libraries ═╣
# ╚══════════════════╝
library(dplyr)
library(stringr)
library(Matrix)
library(Seurat)
library(SeuratObject)

# ╔══════════════════════╗
# ╠═ Read in Parameters ═╣
# ╚══════════════════════╝
message("Reading in Parameters")
args <- commandArgs(trailingOnly = TRUE)

# RDS file from QC
params.SeuratObject <- args[1]

# Integration Method: Options( CCA, RPCA, Harmony, FastMNN, NULL) where NULL is to not run
params.IntegrationMethod <- args[2]

# Project Name
params.ProjectName <- args[3]

# Scale Options ( SD or default SCT)
params.scaleMethod <- args[4]

# ╔════════════════════╗
# ╠═ Load Seurat .rds ═╣
# ╚════════════════════╝
message("Loading Seurat Object")
MergedSO <- readRDS(params.SeuratObject)

# ╔═══════════════════╗
# ╠═ Run Integration ═╣
# ╚═══════════════════╝
message("Running Integration")
if(params.IntegrationMethod == "FastMNN"){
    if(params.scaleMethod == "SD"){
        MergedSO <- IntegrateLayers(object = MergedSO, method = paste0(params.IntegrationMethod,"Integration"), new.reduction = paste0("integrated.",params.IntegrationMethod))
    }else{
        MergedSO <- IntegrateLayers(object = MergedSO, method = paste0(params.IntegrationMethod,"Integration"), new.reduction = paste0("integrated.",params.IntegrationMethod), normalization.method = "SCT")
    }
}else{
    if(params.scaleMethod == "SD"){
        MergedSO <- IntegrateLayers(object = MergedSO, method = paste0(params.IntegrationMethod,"Integration"), orig.reduction = "pca", new.reduction = paste0("integrated.",params.IntegrationMethod))
    }else{
        MergedSO <- IntegrateLayers(object = MergedSO, method = paste0(params.IntegrationMethod,"Integration"), orig.reduction = "pca", new.reduction = paste0("integrated.",params.IntegrationMethod), normalization.method = "SCT")
    }
}

# ╔══════════════════════╗
# ╠═ Save Seurat Object ═╣
# ╚══════════════════════╝
message("Saving Seurat Object")
SaveSeuratRds(MergedSO, file = paste0("05_",params.ProjectName, "_IntegrateSO.rds"))

# ╔═════════════════╗
# ╠═ Save Log File ═╣
# ╚═════════════════╝
sink(paste0("05_",params.ProjectName,"_IntegrateValidation.log"))
cat("╔══════════════════════════════════════════════════════════════════════════════════════════════╗\n")
cat("╠  Integration.R Validation log\n")
cat(paste0("╠  Analysis Group: ", params.ProjectName,"\n"))
cat("╚══════════════════════════════════════════════════════════════════════════════════════════════╝\n")
cat("Seurat Object Status:\n")
print(MergedSO)
sink()

sink(paste0("05_",params.ProjectName,"_IntegrateVersions.log"))
cat("╔══════════════════════════════════════════════════════════════════════════════════════════════╗\n")
cat("╠  Integration.R Versions\n")
cat(paste0("╠  Analysis Group: ", params.ProjectName,"\n"))
cat("╚══════════════════════════════════════════════════════════════════════════════════════════════╝\n")
sessionInfo()
sink()
