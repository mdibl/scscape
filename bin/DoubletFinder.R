#!/opt/conda/envs/dolphinnext/bin/Rscript

# ╔══════════════════════════════════════════════════════════════════════════════════════════════╗
# ╠═                                    Title: DoubletFinder.R                                  ═╣
# ╠═                                     Updated: May 31 2024                                   ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═                                       nf-core/scscape                                      ═╣
# ╠═                                  MDI Biological Laboratory                                 ═╣
# ╠═                          Comparative Genomics and Data Science Core                        ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═ Description:                                                                               ═╣
# ╠═     Takes in the Seurat Object generated from Normalize_QC.R and runs through a process    ═╣
# ╠═     to identify Doublets within the sample. This is run and is reliant on the number of    ═╣
# ╠═     cells recovered and the 10X published doublet rate.                                    ═╣
# ╚══════════════════════════════════════════════════════════════════════════════════════════════╝
# ╔══════════════════╗
# ╠═ Load Libraries ═╣
# ╚══════════════════╝
library(dplyr)
library(stringr)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(DoubletFinder)
library(findPC)

# ╔══════════════════════╗
# ╠═ Read in Parameters ═╣
# ╚══════════════════════╝
message("Reading in Parameters")
args <- commandArgs(trailingOnly = TRUE)

# RDS file from QC
params.SeuratObject <- args[1]

# Variables To Regress in the Scaling
params.VarsToRegress <- args[2]
params.VarsToRegress <- unlist(strsplit(params.VarsToRegress, ","))
print(params.VarsToRegress)

# Data Dir
params.DataDir <- args[3]

# SampleName
params.sampleName <- args[4]

# Variable Features (VF) or ALL (default VF)
params.scaleFeatures <- args[5]

# Scale Options ( SD or default SCT)
params.scaleMethod <- args[6]

# ╔═════════════════════════════════╗
# ╠═ Preparation for Find Doublets ═╣
# ╚═════════════════════════════════╝

# set future max
options(future.globals.maxSize = 4 * 1024^3)  # 4 GiB

message("Preparing for Find Doublets")
# read in .rds
samp <- readRDS(params.SeuratObject)

# create a list of all genes
all.genes <- rownames(samp)

# Scale Data and PCA
if (params.scaleMethod == "SD"){
    if (params.scaleFeatures == "ALL"){
        samp <- ScaleData(samp, features = all.genes ,vars.to.regress = params.VarsToRegress)
        samp <- RunPCA(samp ,npcs = 100)
    }else{
        samp <- ScaleData(samp, vars.to.regress = params.VarsToRegress)
        samp <- RunPCA(samp ,npcs = 100)
    }
}else{
    if (params.scaleFeatures == "ALL"){
        samp <- SCTransform(samp, vars.to.regress = params.VarsToRegress, return.only.var.genes = F )
        samp <- RunPCA(samp ,npcs = 100)
    }else{
        samp <- SCTransform(samp, vars.to.regress = params.VarsToRegress)
        samp <- RunPCA(samp,  npcs = 100)
    }
}

# Generate Elbow Plots and Calculate PCs to use
Elbow <- ElbowPlot(samp,  ndims = 100, reduction = "pca")
ElbowPoints <- Elbow$data

for (i in seq(0.1, 0.9, 0.1)){
    loess <- loess(stdev ~ dims,data=ElbowPoints, span = i)
    ElbowPoints[paste0("loessS", i)] <- loess$fitted
}

for (i in seq(0.1, 0.9, 0.1)){
    first_deriv <- diff(ElbowPoints[,paste0("loessS", i)])/diff(ElbowPoints$dims)
    ElbowPoints[3:100,paste0("stddev_Der2_", i)]<- (diff(first_deriv)/diff(ElbowPoints$dims))[1:length(diff(first_deriv)/diff(ElbowPoints$dims))-1]
}

columns <- colnames(ElbowPoints)
idents <- seq(0.1, 0.9, 0.1)
SSE <- c()
deriv2_var <- c()
for (item in columns){
    if (length(grep("Der2", item)) > 0){
        variance <- var(ElbowPoints[item], na.rm = T)
        deriv2_var <- append(deriv2_var, variance)
    }else if (length(grep("loess", item)) > 0 & !(length(grep("Der", item)) > 0)){
        err <- ElbowPoints[item] - ElbowPoints$stdev
        norms <- sum(((err - min(err))**2)/((max(err) - min(err))**2))
        SSE <- append(SSE, norms)
    }
}

rm_ind <- which(deriv2_var > 0.010)
if (length(rm_ind) > 0){
    SSE <- SSE[-rm_ind]
    deriv2_var <- deriv2_var[-rm_ind]
    idents <- idents[-rm_ind]
    index <- which(SSE == min(SSE))
    ident <- paste0('loessS', idents[index])
}else {
    index <- 1
    ident <- paste0('loessS', idents[index])
}

x = ElbowPoints$dims
y = ElbowPoints$stdev
df <- data.frame(x,y)


for (i in 1:length(ElbowPoints[[ident]])) {
    tmp <- ElbowPoints[[ident]][i]
    if (i == length(ElbowPoints[[ident]])){
        break
    }
    if (ElbowPoints[[ident]][i+1] >= tmp){
        ElbowPoints[[ident]][i+1] <- tmp - 0.00001
    }
}

pc_tbl <- findPC(sdev = ElbowPoints[[ident]], number = 100, method = "all", figure = F)
params.pcMax <- mean(x = c(pc_tbl[1,2], pc_tbl[1,3], pc_tbl[1,4]))
params.pcMax <- ceiling(params.pcMax)

pdf(paste0(params.sampleName,"_ElbowPlot.pdf"), width = 20, height = 15)
ggplot(df, aes(x, y)) +
    geom_point() +
    geom_line(aes(x = seq(1,100), y = ElbowPoints[[ident]]), color = "green") +
    xlab("PC") +
    ylab("Std Dev") +
    geom_vline(xintercept = params.pcMax, linetype="dotted",
               color = "red", size=1.5) +
    ggtitle(paste0("Loess Regression of Std Dev ~ PC  :  PC Chosen = ", params.pcMax))
dev.off()

# Find Neighbors
samp <- FindNeighbors(samp, dims = 1:params.pcMax, reduction = "pca")

# Find Clusters
samp <- FindClusters(samp)

# Build UMAP
samp <- RunUMAP(samp, dims = 1:params.pcMax, reduction = "pca")

# Get the number of recovered cells
params.CellsRecovered <- length(readLines((paste0(params.DataDir,"/",list.files(path = params.DataDir, pattern = "barcodes")))))

# ╔═════════════════════╗
# ╠═ Identify Doublets ═╣
# ╚═════════════════════╝
message("Identifying Doublets")
#    (Percent Doublet is calculated based on a 0.8% increase per 1000 cells recovered as per the 10X website.)
sweep.res.list_meta <- paramSweep(samp, PCs = 1:params.pcMax, sct = ifelse(toupper(params.scaleMethod) == "SD", F,T) )
sweep.stats_meta <- summarizeSweep(sweep.res.list_meta, GT = FALSE)
bcmvn_meta <- find.pK(sweep.stats_meta)
pK <- bcmvn_meta %>%
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
pK <- as.numeric(as.character(pK[[1]]))
annotations <- samp$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(((params.CellsRecovered/1000)*0.008)*nrow(samp@meta.data))
nExp_poi.adj <- round(nExp_poi*(1 - homotypic.prop))
samp <- doubletFinder(samp,
                      PCs = 1:params.pcMax,
                      pN = 0.25,
                      pK = pK,
                      nExp = nExp_poi.adj,
                      reuse.pANN = F,
                      sct = ifelse(toupper(params.scaleMethod) == "SD", F,T))


# Visualize and Count Doublets
doubletmeta <- as.character(colnames(samp[[grep("DF",(colnames(samp[[]])))]]))
pdf(paste0(params.sampleName,"VisualizeDoublets.pdf"), width = 20, height = 15)
DimPlot(samp, group.by = doubletmeta)+ ggtitle(params.sampleName)
dev.off()
DoubletCount <- table(samp[[doubletmeta]])

# ╔═══════════════════════╗
# ╠═ Subset out Doublets ═╣
# ╚═══════════════════════╝
message("Subsetting out Doublets")
samp <- subset(samp, subset = !!as.name(doubletmeta) == "Singlet"  )

# ╔════════════════════╗
# ╠═ Metadata Cleanup ═╣
# ╚════════════════════╝
message("Cleaning up Metadata in Seurat Object")
try(expr = {samp[[doubletmeta]] <- NULL})
try(expr = {samp[["RNA_snn_res.0.8"]] <- NULL})
try(expr = {samp[["SCT_snn_res.0.8"]] <- NULL})


# ╔══════════════════════╗
# ╠═ Save Seurat Object ═╣
# ╚══════════════════════╝
message("Save Seurat Object")
SaveSeuratRds(samp, file = paste0("02_", params.sampleName, "_DoubletsRmSO.rds"))

# ╔═════════════════╗
# ╠═ Save Log File ═╣
# ╚═════════════════╝
sink(paste0("02_", params.sampleName,"_DoubletsRmValidation.log"))
cat("╔══════════════════════════════════════════════════════════════════════════════════════════════╗\n")
cat("╠  DoubletFinder.R Validation log\n")
cat(paste0("╠  Sample: ", params.sampleName,"\n"))
cat("╚══════════════════════════════════════════════════════════════════════════════════════════════╝")
cat(paste0("PCs used: 1 - ", params.pcMax,"\n"))
cat("pN (proportion of artifical doublets) set to: 0.25\n")
cat(paste0("pK (PC neighborhood size used to compute pANN) set to: ", pK,"\n"))
cat(paste0("nExp (Expected Doublets) set to: ", nExp_poi,"\n"))
cat(paste0("nExp Adjusted (Expected Doublets adjusted to proportion of homotypic doublets) set to: ", nExp_poi.adj,"\n"))
cat("Singlet-Doublet Split: \n")
print(DoubletCount)
cat("Seurat Object Status:\n")
print(samp)
sink()

sink(paste0("02_",params.sampleName,"_DoubletsRmVersions.log"))
cat("╔══════════════════════════════════════════════════════════════════════════════════════════════╗\n")
cat("╠  DoubletFinder.R Versions\n")
cat(paste0("╠  Sample: ", params.sampleName,"\n"))
cat("╚══════════════════════════════════════════════════════════════════════════════════════════════╝\n")
sessionInfo()
sink()
