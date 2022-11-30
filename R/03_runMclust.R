##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Updated on November 30, 2022
# Written by Heewon Seo (Heewon.Seo@UCalgary.ca)
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
setwd(file.path(getwd(), "R"))

dataDir <- file.path("../data")
resultsDir <- file.path("../results")

profileFileName <- "cPCDH_expr_25_cohorts.rds"

betaIdx <- c(33:48)

#### #### #### #### #### #### #### #### #### #### #### ####
# User defined functions
.getG <- function(x) {
    if (length(which(is.na(x))) > 0) {
        x <- x[-which(is.na(x))]
    }
    if (which(x == max(x)) <= 3) {
        idx <- which(x == max(x))
    } else {
        idx <- 2
        while (idx <= length(x)) {
            x0 <- x[idx-1]
            x1 <- x[idx]
            if (x0 * 1.05 < x1) {
                idx <- idx + 1
            } else {
                break
            }
        }
        idx <- idx - 1
    }
    return(idx)
}

#### #### #### #### #### #### #### #### #### #### #### ####
library(mclust)

dataL <- readRDS(file.path(dataDir, profileFileName))
geneSymbols <- rownames(dataL[[1]])

binaryExpressionL <- lapply(seq_along(dataL), function(idx) {
        cohort <- names(dataL)[idx]
        dat <- dataL[[idx]]
        message(cohort)

        thresholds <- c()
        buffBinaryExpr <- c()
        
        for (jdx in seq_along(geneSymbols)) {
                geneSymbol <- geneSymbols[jdx]
                message(geneSymbol)

                x <- dat[jdx, ] # epxression

                if (sum(x) > 0 && jdx %in% betaIdx) { # beta genes
                        mod <- densityMclust(x, modelNames="E", plot = FALSE)
                        gdx <- .getG(as.vector(mod$BIC[,"E"]))
                        newMod <- densityMclust(x, modelNames="E", G = gdx, plot = FALSE)
                        df <- data.frame(Sample = rownames(newMod$data), ExonProp = newMod$data, Classification = newMod$classification)

                        if (gdx > 1) {
                                threshold <- max(df$ExonProp[which(df$Classification == gdx-1)])
                                if (threshold < 0.1) {
                                        threshold <- 0.1
                                }
                        } else {
                                threshold <- 1
                        }
                } else { # alpha/gamma genes
                        threshold <- 0.1
                        df <- data.frame(Sample = names(x), ExonProp = x, Classification = 1)
                }
                thresholds <- rbind(thresholds, c(geneSymbol, threshold))

                # Expression in binary
                df$ExonProp[which(is.na(df$ExonProp))] <- 0
                df$ExonProp[which(df$ExonProp <= threshold)] <- 0
                df$ExonProp[which(df$ExonProp >  threshold)] <- 1

                buffBinaryExpr <- cbind(buffBinaryExpr, df$ExonProp)
        }

        thresholdDf <- data.frame(Gene = thresholds[,1], Threshold = as.numeric(thresholds[, 2]))
        write.table(thresholdDf, file.path(resultsDir, paste0("EJPthreshold_", cohort, ".txt")), row.names = F, col.names = T, quote = F, sep = "\t")

        colnames(buffBinaryExpr) <- geneSymbols
        rownames(buffBinaryExpr) <- colnames(dat)
        
        return(buffBinaryExpr)
})

names(binaryExpressionL) <- names(dataL)
saveRDS(binaryExpressionL, file.path(dataDir, paste0("mclust_binaryExpr.RDS")))

q("no")