#!/opt/conda/envs/dolphinnext/bin/Rscript

# ╔══════════════════════════════════════════════════════════════════════════════════════════════╗
# ╠═                                      Title: Plotting.R                                     ═╣
# ╠═                                     Updated: May 31 2024                                   ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═                                       nf-core/scscape                                      ═╣
# ╠═                                  MDI Biological Laboratory                                 ═╣
# ╠═                          Comparative Genomics and Data Science Core                        ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═ Description:                                                                               ═╣
# ╠═     Takes in the Seurat Object generated from FindNeighborsClustersMarkers.R. Runs UMAP    ═╣
# ╠═     and TSNE, makes/saves Preliminary Plots, and creates a Loupe File for interaciive      ═╣
# ╠═     exploration.                                                                           ═╣
# ╚══════════════════════════════════════════════════════════════════════════════════════════════╝
# ╔══════════════════╗
# ╠═ Load Libraries ═╣
# ╚══════════════════╝
library(dplyr)
library(stringr)
library(Matrix)
library(viridis)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)

# ╔══════════════════════╗
# ╠═ Read in Parameters ═╣
# ╚══════════════════════╝
message("Reading in Parameters")
args <- commandArgs(trailingOnly = TRUE)

# RDS file from QC
params.SeuratObject <- args[1]

# Resolutions
params.Resolutions <- args[2]
params.Resolutions <- as.numeric(unlist(strsplit(params.Resolutions, ",")))
print(params.Resolutions)

# PC Max
params.pcMax <- as.integer(args[3])

# Integration Method: Options( CCA, RPCA, Harmony, FastMNN, NULL) where NULL is to not run
params.IntegrationMethod <- args[4]

# Project Name
params.ProjectName <- args[5]

# Make Loupe File T/F
params.MakeLoupe <- args[6]

# Loupe EULA
params.10xEULA <- args[7]

# ╔════════════════════╗
# ╠═ Load Seurat .rds ═╣
# ╚════════════════════╝
message("Loading in Seurat Object")
MergedSO <- readRDS(params.SeuratObject)

# ╔══════════════════════╗
# ╠═ Make UMAP and TSNE ═╣
# ╚══════════════════════╝
message("Running UMAP and TSNE Reductions")
MergedSO <- RunUMAP(MergedSO, dims = 1:params.pcMax, reduction = "pca", reduction.name = "umap.unintegrated")
if(params.IntegrationMethod != "NULL" ){
    MergedSO <- RunUMAP(MergedSO, dims = 1:params.pcMax, reduction = paste0("integrated.",params.IntegrationMethod), reduction.name = paste0("umap.",params.IntegrationMethod))
}
MergedSO <- RunTSNE(MergedSO, dims = 1:params.pcMax, reduction = "pca", reduction.name = "tsne.unintegrated")
if(params.IntegrationMethod != "NULL" ){
    MergedSO <- RunTSNE(MergedSO, dims = 1:params.pcMax, reduction = paste0("integrated.",params.IntegrationMethod), reduction.name = paste0("tsne.",params.IntegrationMethod))
}

# Set Front Page Layout
Page1layout <- "
AAAAAA
AAAAAA
AAAAAA
BBCCDD
"

# ╔══════════════════════════╗
# ╠═ Plot Unintegrated UMAP ═╣
# ╚══════════════════════════╝
message("Generating Unintegrated UMAP Plots")
p1 <- DimPlot(object = MergedSO, reduction = 'umap.unintegrated', pt.size =(-0.00007653*length(MergedSO$orig.ident))+4, label = T, group.by = "orig.ident", shuffle = T)
p2 <- FeaturePlot(MergedSO, reduction = "umap.unintegrated", pt.size = (-0.00001837*length(MergedSO$orig.ident))+1, features = "nFeature_RNA", order = T) + scale_color_viridis(limits =c(min(MergedSO$nFeature_RNA),max(MergedSO$nFeature_RNA)), direction = -1)
p3 <- FeaturePlot(MergedSO, reduction = "umap.unintegrated", pt.size = (-0.00001837*length(MergedSO$orig.ident))+1, features = "nCount_RNA", order = T) + scale_color_viridis(limits =c(min(MergedSO$nCount_RNA),max(MergedSO$nCount_RNA)), direction = -1)
p4 <- FeaturePlot(MergedSO, reduction = "umap.unintegrated", pt.size = (-0.00001837*length(MergedSO$orig.ident))+1, features = "percent.mt", order = T) + scale_color_viridis(limits =c(min(MergedSO$percent.mt),max(MergedSO$percent.mt)), direction = -1)

pdf(paste0(params.ProjectName,"UnintegratedUMAP.pdf"),width = 20, height = 15)
p1 + p2 + p3 + p4 + plot_layout(design = Page1layout)
try(expr = {DimPlot(object = MergedSO, reduction = 'umap.unintegrated', pt.size =(-0.00007653*length(MergedSO$orig.ident))+4, label = T, group.by = "CIscCATCH")})
for (i in params.Resolutions){
    print(DimPlot(object = MergedSO, reduction = 'umap.unintegrated', pt.size =(-0.00007653*length(MergedSO$orig.ident))+4, label = T, group.by = paste0("unintegratedRes.",i), shuffle = T))
}
dev.off()

# ╔══════════════════════════╗
# ╠═ Plot Unintegrated TSNE ═╣
# ╚══════════════════════════╝
message("Generating Unintegrated TSNE Plots")
p1 <- DimPlot(object = MergedSO, reduction = 'tsne.unintegrated', pt.size =(-0.00007653*length(MergedSO$orig.ident))+4, label = T, group.by = "orig.ident", shuffle = T)
p2 <- FeaturePlot(MergedSO, reduction = "tsne.unintegrated", pt.size = (-0.00001837*length(MergedSO$orig.ident))+1, features = "nFeature_RNA", order = T) + scale_color_viridis(limits =c(min(MergedSO$nFeature_RNA),max(MergedSO$nFeature_RNA)), direction = -1)
p3 <- FeaturePlot(MergedSO, reduction = "tsne.unintegrated", pt.size = (-0.00001837*length(MergedSO$orig.ident))+1, features = "nCount_RNA", order = T) + scale_color_viridis(limits =c(min(MergedSO$nCount_RNA),max(MergedSO$nCount_RNA)), direction = -1)
p4 <- FeaturePlot(MergedSO, reduction = "tsne.unintegrated", pt.size = (-0.00001837*length(MergedSO$orig.ident))+1, features = "percent.mt", order = T) + scale_color_viridis(limits =c(min(MergedSO$percent.mt),max(MergedSO$percent.mt)), direction = -1)

pdf(paste0(params.ProjectName,"UnintegratedTSNE.pdf"),width = 20, height = 15)
p1 + p2 + p3 + p4 + plot_layout(design = Page1layout)
try(expr = {DimPlot(object = MergedSO, reduction = 'tsne.unintegrated', pt.size =(-0.00007653*length(MergedSO$orig.ident))+4, label = T, group.by = "CIscCATCH")})
for (i in params.Resolutions){
    print(DimPlot(object = MergedSO, reduction = 'tsne.unintegrated', pt.size =(-0.00007653*length(MergedSO$orig.ident))+4, label = T, group.by = paste0("unintegratedRes.",i), shuffle = T))
}
dev.off()

# ╔═════════════════════════════════╗
# ╠═ Plot Integrated UMAP (if run) ═╣
# ╚═════════════════════════════════╝
pdf(paste0(params.ProjectName,params.IntegrationMethod,"IntegratedUMAP.pdf"),width = 20, height = 15)
if (params.IntegrationMethod != "NULL"){
    message("Generating Integrated UMAP Plots")
    p1 <- DimPlot(object = MergedSO, reduction = paste0("umap.",params.IntegrationMethod), pt.size =(-0.00007653*length(MergedSO$orig.ident))+4, label = T, group.by = "orig.ident", shuffle = T)
    p2 <- FeaturePlot(MergedSO, reduction = paste0("umap.",params.IntegrationMethod), pt.size = (-0.00001837*length(MergedSO$orig.ident))+1, features = "nFeature_RNA", order = T) + scale_color_viridis(limits =c(min(MergedSO$nFeature_RNA),max(MergedSO$nFeature_RNA)), direction = -1)
    p3 <- FeaturePlot(MergedSO, reduction = paste0("umap.",params.IntegrationMethod), pt.size = (-0.00001837*length(MergedSO$orig.ident))+1, features = "nCount_RNA", order = T) + scale_color_viridis(limits =c(min(MergedSO$nCount_RNA),max(MergedSO$nCount_RNA)), direction = -1)
    p4 <- FeaturePlot(MergedSO, reduction = paste0("umap.",params.IntegrationMethod), pt.size = (-0.00001837*length(MergedSO$orig.ident))+1, features = "percent.mt", order = T) + scale_color_viridis(limits =c(min(MergedSO$percent.mt),max(MergedSO$percent.mt)), direction = -1)

    print(p1 + p2 + p3 + p4 + plot_layout(design = Page1layout))
    try(expr = {print(DimPlot(object = MergedSO, reduction = paste0("umap.",params.IntegrationMethod), pt.size =(-0.00007653*length(MergedSO$orig.ident))+4, label = T, group.by = "CIscCATCH"))})
    for (i in params.Resolutions){
        print(DimPlot(object = MergedSO, reduction = paste0("umap.",params.IntegrationMethod), pt.size =(-0.00007653*length(MergedSO$orig.ident))+4, label = T, group.by = paste0(params.IntegrationMethod,"Res.",i)))
    }
}
dev.off()

if(params.IntegrationMethod == "NULL"){
    file.remove(paste0(params.ProjectName,params.IntegrationMethod,"IntegratedUMAP.pdf"))
}

# ╔═════════════════════════════════╗
# ╠═ Plot Integrated TSNE (if run) ═╣
# ╚═════════════════════════════════╝
pdf(paste0(params.ProjectName,params.IntegrationMethod,"IntegratedTSNE.pdf"),width = 20, height = 15)
if (params.IntegrationMethod != "NULL"){
    message("Generating Integrated TSNE Plots")
    p1 <- DimPlot(object = MergedSO, reduction = paste0("tsne.",params.IntegrationMethod), pt.size =(-0.00007653*length(MergedSO$orig.ident))+4, label = T, group.by = "orig.ident", shuffle = T)
    p2 <- FeaturePlot(MergedSO, reduction = paste0("tsne.",params.IntegrationMethod), pt.size = (-0.00001837*length(MergedSO$orig.ident))+1, features = "nFeature_RNA", order = T) + scale_color_viridis(limits =c(min(MergedSO$nFeature_RNA),max(MergedSO$nFeature_RNA)), direction = -1)
    p3 <- FeaturePlot(MergedSO, reduction = paste0("tsne.",params.IntegrationMethod), pt.size = (-0.00001837*length(MergedSO$orig.ident))+1, features = "nCount_RNA", order = T) + scale_color_viridis(limits =c(min(MergedSO$nCount_RNA),max(MergedSO$nCount_RNA)), direction = -1)
    p4 <- FeaturePlot(MergedSO, reduction = paste0("tsne.",params.IntegrationMethod), pt.size = (-0.00001837*length(MergedSO$orig.ident))+1, features = "percent.mt", order = T) + scale_color_viridis(limits =c(min(MergedSO$percent.mt),max(MergedSO$percent.mt)), direction = -1)

    print(p1 + p2 + p3 + p4 + plot_layout(design = Page1layout))
    try(expr = {print(DimPlot(object = MergedSO, reduction = paste0("tsne.",params.IntegrationMethod), pt.size =(-0.00007653*length(MergedSO$orig.ident))+4, label = T, group.by = "CIscCATCH"))})
    for (i in params.Resolutions){
        print(DimPlot(object = MergedSO, reduction = paste0("tsne.",params.IntegrationMethod), pt.size =(-0.00007653*length(MergedSO$orig.ident))+4, label = T, group.by = paste0(params.IntegrationMethod,"Res.",i)))
    }
}
dev.off()

if(params.IntegrationMethod == "NULL"){
    file.remove(paste0(params.ProjectName,params.IntegrationMethod,"IntegratedTSNE.pdf"))
}

# ╔═══════════════════╗
# ╠═ Make Loupe File ═╣
# ╚═══════════════════╝
message("Generating Loupe File")
EULAmessage <- NULL
generate_loupe_file <- function(MergedSO, params.ProjectName) {
    library(loupeR)
    loupeR::setup()
    create_loupe(count_mat = MergedSO@assays$RNA$counts,
                 clusters = select_clusters(MergedSO),
                 projections = select_projections(MergedSO),
                 output_name = params.ProjectName
    )
}

try_generate_loupe_file <- function(MergedSO, params.ProjectName, max_tries = 5) {
    for (i in 1:max_tries) {
        tryCatch({
            generate_loupe_file(MergedSO, params.ProjectName)
            message("Loupe file generated successfully.")
            return(TRUE)
        }, error = function(e) {
            message(paste("Attempt", i, "failed:", e$message))
            if (i < max_tries) {
                Sys.sleep(60)  # Wait for 1 minute before retrying
            } else {
                message(paste0("Failed to generate Loupe file after ", max_tries, " attempts."))
                message("Loupe browser will need to be generated externally.")
                return(FALSE)
            }
        })
    }
}

if(toupper(params.MakeLoupe) == "TRUE"){
    if (toupper(params.10xEULA) == "AGREE"){
        try_generate_loupe_file(MergedSO, params.ProjectName)
    } else {
        EULAmessage <- "WARNING: Loupe File set to TRUE but you have not agreed to the 10x EULA -- Set params.10xEULA to Agree to create Loupe File."
        message(EULAmessage)
    }
}else {
    if (toupper(params.10xEULA) == "AGREE"){
        try_generate_loupe_file(MergedSO, params.ProjectName)
        EULAmessage <- "NOTE: Loupe File set to FALSE but you have agreed to the 10x EULA -- Loupe File has been made."
        message(EULAmessage)
    }
}

# ╔══════════════════════╗
# ╠═ Save Seurat Object ═╣
# ╚══════════════════════╝
message("Saving Seurat Object")
SaveSeuratRds(MergedSO, file = paste0("08",params.ProjectName, "_FinalSO.rds"))

# ╔════════════════════════╗
# ╠═ Save Meta Data Table ═╣
# ╚════════════════════════╝
meta_table <- MergedSO@meta.data
meta_table <- cbind(meta_table, as.data.table(MergedSO@reductions$pca@cell.embeddings))
if (params.IntegrationMethod != "NULL"){
    meta_table <- cbind(meta_table, MergedSO@reductions[[paste0("tsne.",params.IntegrationMethod)]]@cell.embeddings)
    meta_table <- cbind(meta_table, MergedSO@reductions[[paste0("umap.",params.IntegrationMethod)]]@cell.embeddings)
} else {
    meta_table <- cbind(meta_table, MergedSO@reductions[["tsne.unintegrated"]]@cell.embeddings)
    meta_table <- cbind(meta_table, MergedSO@reductions[["umap.unintegrated"]]@cell.embeddings)
}
write.csv(file = paste0("08_",params.ProjectName,"_MetaTable.tsv"), meta_table, col.names = T)



# ╔═════════════════╗
# ╠═ Save Log File ═╣
# ╚═════════════════╝
sink(paste0("08_",params.ProjectName,"_PlotValidation.log"))
cat("╔══════════════════════════════════════════════════════════════════════════════════════════════\n╗")
cat("╠  Plotting.R Validation log\n")
cat(paste0("╠  Analysis Group: ", params.ProjectName,"\n"))
cat("╚══════════════════════════════════════════════════════════════════════════════════════════════╝\n")
cat("Seurat Object Status:\n")
print(MergedSO)
if (length(EULAmessage) != 0){
    cat(EULAmessage)
    cat("\n")
}
sink()

sink(paste0("08_",params.ProjectName,"_PlotVersions.log"))
cat("╔══════════════════════════════════════════════════════════════════════════════════════════════╗\n")
cat("╠  Plotting.R Versions\n")
cat(paste0("╠  Analysis Group: ", params.ProjectName,"\n"))
cat("╚══════════════════════════════════════════════════════════════════════════════════════════════╝\n")
sessionInfo()
sink()
