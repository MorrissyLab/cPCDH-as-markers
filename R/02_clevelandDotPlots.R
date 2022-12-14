##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Updated on November 30, 2022
# Written by Heewon Seo (Heewon.Seo@UCalgary.ca)
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
setwd(file.path(getwd(), "R"))

dataDir <- file.path("../data")
resultsDir <- file.path("../results")

profileFileName <- "GeoMean_SD_Max.rds"

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(ggplot2)
library(stringr)

dataL <- readRDS(file.path(dataDir, profileFileName))
pcdhGenes <- rownames(dataL[[1]])

fullMat <- Reduce(function(...) merge(..., by = "Gene", all = TRUE), dataL)
rownames(fullMat) <- fullMat[,1]; fullMat <- fullMat[pcdhGenes,]

mData <- as.data.frame(fullMat[, seq(2, ncol(fullMat), 3)])
sdData <- as.data.frame(fullMat[, seq(3, ncol(fullMat), 3)])
maxData <- as.data.frame(fullMat[, seq(4, ncol(fullMat), 3)])

colnames(mData) <- names(dataL)
rownames(mData) <- fullMat$Gene
colnames(sdData) <- names(dataL)
rownames(sdData) <- fullMat$Gene
colnames(maxData) <- names(dataL)
rownames(maxData) <- fullMat$Gene

for (idx in seq_along(fullMat$Gene)) { # foreach gene
        geneName <- fullMat$Gene[idx]
        print(geneName)

        df <- data.frame(Cohort = names(dataL), GeoMean = unlist(mData[idx, ]), SD = unlist(sdData[idx, ]), MAX = unlist(maxData[idx, ]))
        df$Cohort <- factor(df$Cohort,
                levels = rev(c("GTEx", "CPTAC-GTEx", "BrainSpanI", "BrainSpanII", "BrainSpanIII", # non-malignant samples
                        "TCGA-GBM", "CPTAC-GBM", "TFRI-GBM", "TCGA-LGG", # adult brain tumors
                        "CBTN-ATRT", "CBTN-CP", "CBTN-CPP", "CBTN-DIPG", "CBTN-DNET", "CBTN-EPD", "CBTN-GG", "CBTN-HGG", "CBTN-LGG", "CBTN-MB", "CBTN-MNG", "CBTN-PNF", "CBTN-SWM", # Ped-brain tumors
                        "Mayo-PDX", "TFRI-Xeno", "TFRI-BTIC")), # preclinical models
                labels = rev(c("GTEx (230)", "CPTAC-GTEx (9)", "BrainSpan I (231)", "BrainSpan II (43)", "BrainSpan III (227)", # non-malignant samples
                        "TCGA-GBM (152)", "CPTAC-GBM (82)", "TFRI-GBM (44)", "TCGA-LGG (530)", # adult brain tumors
                        "CBTN-ATRT (25)", "CBTN-CP (10)", "CBTN-CPP (21)", "CBTN-DIPG (15)", "CBTN-DNET (21)", "CBTN-EPD (59)", "CBTN-GG (33)", "CBTN-HGG (83)", "CBTN-LGG (212)", "CBTN-MB (95)", "CBTN-MNG (22)", "CBTN-PNF (13)", "CBTN-SWM (14)", # Ped-brain tumors
                        "Mayo-PDX (51)", "TFRI-Xeno (13)", "TFRI-BTIC (61)")) # preclinical models
        )

        clevelandDotPlot <- ggplot(df) +
                geom_segment(aes(x = Cohort, xend = Cohort, y = GeoMean, yend = SD), color = "grey80") +
                geom_segment(aes(x = Cohort, xend = Cohort, y = GeoMean, yend = MAX), color = "grey80") +
                geom_point(aes(x = Cohort, y = MAX), color = rgb(0, 0, 0, 0.5), size = 3) +
                geom_point(aes(x = Cohort, y = SD), color = rgb(0.74, 0.23, 0.02, 0.5), size = 3) +
                geom_point(aes(x = Cohort, y = GeoMean), color = rgb(0.55, 0.44, 0, 0.5), size = 3) +
                geom_vline(xintercept = c(3.5, 16.5, 20.5), col = "grey40", lty = 2) +
                geom_hline(yintercept = 0, col = "grey80", lty = 1, lwd = 0.5) +
                coord_flip() +
                scale_y_continuous(
                        breaks = c(0, 0.1, 0.2, 0.3), 
                        labels = c("0", "0.1", "0.2", "0.3"), 
                        limits = c(0, 0.375),
                        sec.axis = sec_axis(trans=~.*3, name="Max(black)")
                ) +
                labs(
                        title = str_replace(geneName, "PCDH", ""),
                        x = "", y = "GeoMean(gold)/SD(red)"
                ) +
                theme_bw() +
                theme(
                        axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
                        legend.position = "none",
                        text = element_text(size = 12)
                )
        pdf(file.path(resultsDir, paste0("ClevelandDotPlot_", geneName, ".pdf")), width = 7, height = 7)
        print(clevelandDotPlot)
        dev.off()
}

q("no")