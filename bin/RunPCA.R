#!/opt/conda/envs/dolphinnext/bin/Rscript

# ╔══════════════════════════════════════════════════════════════════════════════════════════════╗
# ╠═                                       Title: RunPCA.R                                      ═╣
# ╠═                                     Updated: May 31 2024                                   ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═                                       nf-core/scscape                                      ═╣
# ╠═                                  MDI Biological Laboratory                                 ═╣
# ╠═                          Comparative Genomics and Data Science Core                        ═╣
# ╠══════════════════════════════════════════════════════════════════════════════════════════════╣
# ╠═ Description:                                                                               ═╣
# ╠═     Takes in the Merged Seurat Object generated from Merge.R. PCA is run on the Seurat     ═╣
# ╠═     Object and then a max PC value is selected for follow on steps. The max PC value is    ═╣
# ╠═     selected by smoothing the elbow plot and then averaging the first derivative, second   ═╣
# ╠═     derivative and preceding residual selection from FindPC.                               ═╣
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
library(findPC)

# ╔══════════════════════╗
# ╠═ Read in Parameters ═╣
# ╚══════════════════════╝
message("Reading in Parameters")
args <- commandArgs(trailingOnly = TRUE)

# Load RDS
MergedSO <- LoadSeuratRds(args[1])

# Maximum PC to be used in dimensional reduciton/ find neighbors
params.pcMax <- args[2]

# Sample Name
params.ProjectName <- args[3]

# ╔═══════════╗
# ╠═ Run PCA ═╣
# ╚═══════════╝
message("Running PCA")
MergedSO <- RunPCA(MergedSO, npcs = 100)

# ╔═══════════════╗
# ╠═ Find PC Max ═╣
# ╚═══════════════╝
message("Finding PC Max")
Elbow <- ElbowPlot(MergedSO,  ndims = 100, reduction = "pca")
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

if (toupper(params.pcMax) == "NULL"){
    for (i in 1:length(ElbowPoints[[ident]])) {
        tmp <- ElbowPoints[[ident]][i]
        if (i == length(ElbowPoints[[ident]])){
            break
        }
        if (ElbowPoints[[ident]][i+1] >= tmp){
            ElbowPoints[[ident]][i+1] <- tmp - 0.00001
        }
    }
    
    pc_tbl <- findPC(sdev = ElbowPoints[[ident]], number = 100, method = "all", figure = T)
    params.pcMax <- mean(x = c(pc_tbl[1,2], pc_tbl[1,3], pc_tbl[1,4]))
    params.pcMax <- ceiling(params.pcMax)
}else{
  params.pcMax
}

x = ElbowPoints$dims
y = ElbowPoints$stdev
df <- data.frame(x,y)

pdf(paste0(params.ProjectName,"_Merged_ElbowPlot.pdf"), width = 20, height = 15)
ggplot(df, aes(x, y)) +
  geom_point() +
  geom_line(aes(x = seq(1,100), y = ElbowPoints[[ident]]), color = "green") +
  xlab("PC") +
  ylab("Std Dev") +
  geom_vline(xintercept = params.pcMax, linetype="dotted", 
             color = "red", size=1.5) +
  ggtitle(paste0("Loess Regression of Std Dev ~ PC  :  PC Chosen = ", params.pcMax))
dev.off()

# ╔══════════════════════╗
# ╠═ Save Seurat Object ═╣
# ╚══════════════════════╝
message("Saving Seurat Object")
SaveSeuratRds(MergedSO, file = paste0("04_", params.ProjectName, "_PCASO.rds"))

# ╔═════════════════╗
# ╠═ Save Log File ═╣
# ╚═════════════════╝
loess <- loess(stdev ~ dims,data=ElbowPoints, span = idents[index])
my_summary <- summary(loess)

sink(paste0("04_", params.ProjectName,"_PCAValidation.log"))
cat("╔══════════════════════════════════════════════════════════════════════════════════════════════╗\n")
cat("╠  RunPCA.R Validation log\n")
cat(paste0("╠  Analysis Group: ", params.ProjectName,"\n"))
cat("╚══════════════════════════════════════════════════════════════════════════════════════════════╝\n")
cat("Loess Summary: \n")
print(my_summary)
cat("\n")
cat("PC cutoff calculations:\n")
cat(paste0("First Derivative : ", pc_tbl[1,2], "\n"))
cat(paste0("Second Derivative : ", pc_tbl[1,3], "\n"))
cat(paste0("Preceding Residual : ", pc_tbl[1,4], "\n"))
cat("\n")
cat(paste0("PCs used: 1 - ", params.pcMax,"\n"))
sink()

sink(paste0("04_",params.ProjectName,"_PCAVersions.log"))
cat("╔══════════════════════════════════════════════════════════════════════════════════════════════╗\n")
cat("╠  RunPCA.R Versions\n")
cat(paste0("╠  Analysis Group: ", params.ProjectName,"\n"))
cat("╚══════════════════════════════════════════════════════════════════════════════════════════════\n╝")
sessionInfo()
sink()
