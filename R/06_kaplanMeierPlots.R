##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Updated on December 14, 2022
# Written by Heewon Seo (Heewon.Seo@UCalgary.ca)
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
setwd(file.path(getwd(), "R"))

dataDir <- file.path("../data")
resultsDir <- file.path("../results")

profileFileName <- "isoformUsages.rds"
osFileName <- "overallSurvival.rds"

cohort <- "TCGA-LGG"
geneOne <- "PCDHA4"
geneTwo <- "PCDHGB7"

#### #### #### #### #### #### #### #### #### #### #### ####
library(stringr)
library(survival)
library(survminer)
library(RColorBrewer)

# User defined functions
.reform <- function(x) { # TCGA-GBM
        buff <- unlist(stringr::str_split(x, "\\."))
        return(paste0(buff[1], "-", buff[2], "-", buff[3], "-", buff[4]))
}
.reform2 <- function(x) { # TCGA-LGG
        buff <- unlist(stringr::str_split(x, "\\."))
        return(paste0(buff[1], "-", buff[2], "-", buff[3], "-", str_sub(buff[4], 1, 2)))
}

isoformUsages <- readRDS(file.path(dataDir, profileFileName))
overallSurvival <- readRDS(file.path(dataDir, osFileName))

if (geneOne == geneTwo) {
        geneTwo = "NULL"
}

group <- unlist(str_split(cohort, "-"))[1]

if (group == "CBTN") { # Childhood brain tumors
        dat <- isoformUsages[[cohort]]
        mat <- data.frame(Sample = colnames(dat), Group = "None")
        osInfo <- overallSurvival[[group]]
        colnames(osInfo)[3] <- c("Sample")
} else {
        dat <- isoformUsages[[cohort]]
        mat <- data.frame(Sample = colnames(dat), Group = "None")
        osInfo <- overallSurvival[[cohort]]
}

for (idx in c(1:nrow(mat))) {
        x <- dat[, idx]
        if (sum(x) > 2) { # triplet
                gene1 <- x[which(names(x) == geneOne)]
                if (geneTwo != "NULL") {
                        gene2 <- x[which(names(x) == geneTwo)]

                        if (gene1 == 1 && gene2 == 1) {
                                mat$Group[idx] <- paste0(str_replace(geneOne, "PCDH", ""), "+", str_replace(geneTwo, "PCDH", ""))
                        } else if (gene1 == 1 && gene2 != 1) {
                                mat$Group[idx] <- str_replace(geneOne, "PCDH", "")
                        } else if (gene1 != 1 && gene2 == 1) {
                                mat$Group[idx] <- str_replace(geneTwo, "PCDH", "")
                        } else if (gene1 != 1 && gene2 != 1) {
                                mat$Group[idx] <- "Neither"
                        }
                } else {
                        if (gene1 == 1) {
                                mat$Group[idx] <- paste0(str_replace(geneOne, "PCDH", ""), "+")
                        } else {
                                mat$Group[idx] <- paste0(str_replace(geneOne, "PCDH", ""), "-")
                        }
                }
        }
}

if (geneTwo != "NULL") {
        mat$Group <- factor(mat$Group, levels = c(paste0(str_replace(geneOne, "PCDH", ""), "+", str_replace(geneTwo, "PCDH", "")), str_replace(geneOne, "PCDH", ""), str_replace(geneTwo, "PCDH", ""), "Neither"))
} else {
        mat$Group <- factor(mat$Group, levels = c(paste0(str_replace(geneOne, "PCDH", ""), "+"), paste0(str_replace(geneOne, "PCDH", ""), "-")))
}

if (cohort == "TCGA-GBM") {
        mat$Sample <- sapply(mat$Sample, .reform)
} else if (cohort == "TCGA-LGG") {
        mat$Sample <- sapply(mat$Sample, .reform2)
}

survival <- merge(mat, osInfo, by = "Sample", all = FALSE)
if (length(levels(mat$Group)) != length(unique(survival$Group))) {
        survival$Group <- factor(survival$Group)
}

os <- survfit(Surv(OS.days, Deceased) ~ Group, data = survival)

km <- survminer::ggsurvplot(
        os,
        data = survival,
        conf.int = FALSE,
        pval.method = TRUE, pval = TRUE,
        legend.title = "",
        xlab = "OS, days",
        title = cohort,
        palette = brewer.pal(n = length(levels(survival$Group)), name = "Dark2"),
        risk.table = TRUE,
        legend.labs = levels(survival$Group)
)

pdf(file.path(resultsDir, paste0(cohort, "_", geneOne, "_", geneTwo, ".pdf")), width=6, height=6, onefile = FALSE)
print(km)
dev.off()

q("no")