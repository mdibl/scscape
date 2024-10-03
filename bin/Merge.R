#!/opt/conda/envs/dolphinnext/bin/Rscript

# ╔══════════════════════════════════════════════════════════════════════════════════════════════╗
# ╠═                                        Title: Merge.R                                      ═╣
# ╠═                                     Updated: May 31 2024                                   ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═                                       nf-core/scscape                                      ═╣
# ╠═                                  MDI Biological Laboratory                                 ═╣
# ╠═                          Comparative Genomics and Data Science Core                        ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═ Description:                                                                               ═╣
# ╠═     Takes all of the Seurat Objects post DoubletFinder.R that belong to a specific         ═╣
# ╠═     analysis group and run them through the following steps: Merge, Find Variable Features ═╣
# ╠═     Generate merged Quality Control Plots, and Scaling.                                    ═╣
# ╚══════════════════════════════════════════════════════════════════════════════════════════════╝
# ╔══════════════════╗
# ╠═ Load Libraries ═╣
# ╚══════════════════╝
library(dplyr)
library(stringr)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(patchwork)


# ╔══════════════════════╗
# ╠═ Read in Parameters ═╣
# ╚══════════════════════╝
message("Reading in Parameters")
args <- commandArgs(trailingOnly = TRUE)

# Project Name for this analysis group
params.ProjectName <- args[1]

# Variable To Regress For Scaling
params.VarsToRegress <- args[2]
params.VarsToRegress <- unlist(strsplit(params.VarsToRegress, ","))
print(params.VarsToRegress)

# Variable Features (VF) or ALL (default VF)
params.scaleFeatures <- args[3]

# Scale Options ( SD or default SCT)
params.scaleMethod <- args[4]


# ╔════════════════════╗
# ╠═ Load Seurat .rds ═╣
# ╚════════════════════╝
message("Loading in Seurat Objects")
SeuratRDSfiles <- list.files(pattern = "\\.rds$", ignore.case = T)
SampleNames <- c()

counter <- 1 
for( file in SeuratRDSfiles){
    ObjName <- gsub("^02_","",gsub("_Doublets.*","", file))
    SampleNames <- append(SampleNames,ObjName)
    assign(ObjName, readRDS(file))
}

count <- 2
countMax <- length(SampleNames)
SOlist2up <- c()
while (count <= countMax){
    SOlist2up <- append(SOlist2up, get(SampleNames[count]))
    count <- count + 1
}

MergedSO <- merge(get(SampleNames[1]), y= SOlist2up, add.cell.ids = SampleNames, project = params.ProjectName)
MergedSO@assays$SCT <- NULL
DefaultAssay(object = MergedSO) <- "RNA"


# ╔══════════════════╗
# ╠═ Normalize Data ═╣
# ╚══════════════════╝
message("Normalizing Data")
MergedSO <- NormalizeData(MergedSO)

# ╔═════════════════════════════════════╗
# ╠═ Generate Merged Post-Filter Plots ═╣
# ╚═════════════════════════════════════╝
message("Generating Merged Post-Filter Plots")
plot1 <- FeatureScatter(MergedSO, feature1 = "nCount_RNA", feature2 = "percent.mt", shuffle = T, group.by = "orig.ident")
plot2 <- FeatureScatter(MergedSO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", shuffle = T, group.by = "orig.ident")
vplot1 <- VlnPlot(MergedSO, features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = -1, group.by = "orig.ident")

pdf(file = paste0(params.ProjectName,'postMergeQC.pdf'), title = paste0(params.ProjectName,' Post Merge QC'), width = 11, height = 8.5)
plot1 + plot2
vplot1
dev.off()

# ╔══════════════════════════════╗
# ╠═ Scale Merged Seurat Object ═╣
# ╚══════════════════════════════╝
message("Scaling Merged Seurat Object")
# create a list of all genes
all.genes <- rownames(MergedSO)

# Scale Data and PCA
if (params.scaleMethod == "SD"){
    if (params.scaleFeatures == "ALL"){
        MergedSO <- FindVariableFeatures(MergedSO)
        MergedSO <- ScaleData(MergedSO, features = all.genes ,vars.to.regress = params.VarsToRegress)
    }else{
        MergedSO <- FindVariableFeatures(MergedSO)
        MergedSO <- ScaleData(MergedSO, vars.to.regress = params.VarsToRegress)
    }
}else{
    if (params.scaleFeatures == "ALL"){
        MergedSO <- SCTransform(MergedSO, vars.to.regress = params.VarsToRegress, return.only.var.genes = F )
    }else{
        MergedSO <- SCTransform(MergedSO, vars.to.regress = params.VarsToRegress)
    }
}

# ╔══════════════════════╗
# ╠═ Save Seurat Object ═╣
# ╚══════════════════════╝
message("Saving Seurat Object")
saveRDS(MergedSO, paste0("03_", params.ProjectName,"_MergedSO.rds"))


# ╔═════════════════╗
# ╠═ Save Log File ═╣
# ╚═════════════════╝
sink(paste0("03_",params.ProjectName,"_MergeValidation.log"))
cat("╔══════════════════════════════════════════════════════════════════════════════════════════════╗\n")
cat("╠  Merge.R Validation log\n")
cat(paste0("╠  Analysis Group: ", params.ProjectName,"\n"))
cat("╚══════════════════════════════════════════════════════════════════════════════════════════════╝\n")
cat("Seurat Object Status:\n")
print(MergedSO)
sink()

sink(paste0("03_",params.ProjectName,"_MergeVersions.log"))
cat("╔══════════════════════════════════════════════════════════════════════════════════════════════╗\n")
cat("╠  Merge.R Versions\n")
cat(paste0("╠  Analysis Group: ", params.ProjectName,"\n"))
cat("╚══════════════════════════════════════════════════════════════════════════════════════════════╝\n")
sessionInfo()
sink()



