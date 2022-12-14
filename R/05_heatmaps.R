##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Updated on December 14, 2022
# Written by Heewon Seo (Heewon.Seo@UCalgary.ca)
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
setwd(file.path(getwd(), "R"))

dataDir <- file.path("../data")
resultsDir <- file.path("../results")

profileFileName <- "tripletStats.rds"

selectCohorts <- c("TCGA-LGG-IDHwt", "TCGA-LGG-IDHmut-codel", "TCGA-LGG-IDHmut-non-codel", "TCGA-GBM", "CPTAC-GBM", "CBTN-ATRT", "CBTN-CP", "CBTN-CPP", "CBTN-DIPG", "CBTN-DNET", "CBTN-EPD", "CBTN-GG", "CBTN-HGG", "CBTN-LGG", "CBTN-MB", "CBTN-MNG", "CBTN-PNF", "CBTN-SWM")

#### #### #### #### #### #### #### #### #### #### #### ####
library(stringr)
library(ComplexHeatmap)

triplets <- readRDS(file.path(dataDir, profileFileName))

# Overrepresented triplets (top 10) in heatmap
freq <- triplets[["Triplet"]][["Freq"]]
freq <- freq[, which(colnames(freq) %in% selectCohorts)]

topTenTriplets <- c()
for (idx in c(1:ncol(freq))) {
        sub <- freq[, idx]
        names(sub) <- rownames(freq)
        sub <- sub[order(sub, decreasing = T)]
        topTenTriplets <- c(topTenTriplets, names(sub)[which(sub >= sub[10])])
}
topTenTriplets <- unique(topTenTriplets)

prop <- triplets[["Triplet"]][["Prop"]]
prop <- prop[which(rownames(prop) %in% topTenTriplets), which(colnames(prop) %in% selectCohorts)]

topAnnotation <- HeatmapAnnotation(
        df = data.frame(
                A4 = factor(
                        sapply(rownames(prop), function(x) {
                                y <- str_replace_all(x, "PCDH", "")
                                buff <- unlist(str_split(y, "\\$"))
                                hasIt <- "N"
                                if (any(buff %in% "A4")) {
                                        hasIt <- "Y"
                                }
                                return(hasIt)
                        }),
                        levels = c("Y", "N")
                ),
                GB7 = factor(
                        sapply(rownames(prop), function(x) {
                                y <- str_replace_all(x, "PCDH", "")
                                buff <- unlist(str_split(y, "\\$"))
                                hasIt <- "N"
                                if (any(buff %in% "GB7")) {
                                        hasIt <- "Y"
                                }
                                return(hasIt)
                        }),
                        levels = c("Y", "N")
                )
        ),
        col = list(
                A4 = c("Y" = "indianred2", "N" = "grey60"),
                GB7 = c("Y" = "skyblue", "N" = "grey60")
        )
)
labels <- triplets[["Label"]]
labels <- labels[which(labels$cohort %in% selectCohorts), ]
rownames(labels) <- labels$cohort
labels <- labels[colnames(prop), ]
newRowNames <- labels$labels
newColNames <- str_replace_all(str_replace_all(rownames(prop), "PCDH", ""), "\\$", "-")

heatmap1 <- Heatmap(
        t(prop),
        name = "Proportion",
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_title_side = "top",
        column_title = "",
        row_title_side = "right",
        row_title = "",
        show_heatmap_legend = TRUE,
        col = circlize::colorRamp2(c(0, 0.5), c("white", "#4062A1")),
        top_annotation = topAnnotation,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_max_width = unit(12, "cm"),
        row_labels = newRowNames,
        column_labels = newColNames
)
pdf(file.path(resultsDir, "TopTenTriplets.pdf"), width=20, height=10)
draw(heatmap1)
dev.off()

# Isoform usage from triplets above in heatmap
prop <- triplets[["Isoform"]][["Prop"]]
prop <- prop[, which(colnames(prop) %in% selectCohorts)]
propSum <- apply(prop, 1, sum)
prop <- prop[which(propSum > 0), ]

labels <- triplets[["Label"]]
labels <- labels[which(labels$cohort %in% selectCohorts), ]
rownames(labels) <- labels$cohort
labels <- labels[colnames(prop), ]
newRowNames <- labels$labels

# heatmap
topAnnotation <- HeatmapAnnotation(
        df = data.frame(
                Cluster = factor(sapply(rownames(prop), function(x) {
                        return(paste(unlist(str_extract_all(x, pattern = regex("[^0-9]"))), collapse = ""))
                }), levels = c("A", "B", "GA", "GB"), labels = c("Alpha", "Beta", "Gamma A", "Gamma B"))
        ),
        col = list(
                Cluster = c(
                        "Alpha" = "#4062A1",
                        "Beta" = "#BC312A",
                        "Gamma A" = "#149103",
                        "Gamma B" = "#0C5C01"
                )
        )
)

clusterNames <- c("Alpha", "Beta", "Gamma")
names(clusterNames) <- c("A", "B", "G")

heatmap2 <- Heatmap(
        t(prop),
        name = "Proportion",
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_title_side = "top",
        column_title = "",
        row_title_side = "right",
        row_title = "",
        show_heatmap_legend = TRUE,
        col = circlize::colorRamp2(c(0, 0.5), c("white", "#4062A1")),
        top_annotation = topAnnotation,
        column_split = clusterNames[substr(rownames(prop), 1, 1)],
        row_labels = newRowNames
)

pdf(file.path(resultsDir, "TopTenTriplets_isoforms.pdf"), width=10, height=10)
draw(heatmap2)
dev.off()

q("no")