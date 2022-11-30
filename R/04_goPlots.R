##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Updated on November 30, 2022
# Written by Heewon Seo (Heewon.Seo@UCalgary.ca)
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
setwd(file.path(getwd(), "R"))

dataDir <- file.path("../data")
resultsDir <- file.path("../results")

profileFileName <- "mclust_binaryExpr.rds"

alphaIdx <- c(1:13)
gammaIdx <- c(14:32)
betaIdx <- c(33:48)

howManyPt <- 3

nodeCol <- c(
    rep("#4062A1", length(alphaIdx)), 
    "#149103", "#149103", "#149103", 
    "#0C5C01", "#149103", "#0C5C01", 
    "#149103", "#0C5C01", "#149103", 
    "#149103", "#0C5C01", "#149103", 
    "#0C5C01", "#149103", "#0C5C01", 
    "#149103", "#0C5C01", "#149103", 
    "#149103",
    rep("#BC312A", length(betaIdx))
)

#### #### #### #### #### #### #### #### #### #### #### ####
library(stringr)
library(ggplot2)

dataL <- readRDS(file.path(dataDir, profileFileName))
pcdhGenes <- colnames(dataL[[1]])

names(nodeCol) <- substr(pcdhGenes, 5, nchar(pcdhGenes))

pcdhGeneIdx <- c(length(pcdhGenes):1)
names(pcdhGeneIdx) <- substr(pcdhGenes, 5, nchar(pcdhGenes))

returnNoMeaning <- lapply(seq_along(dataL), function(idx) {
        cohort <- names(dataL)[idx]
        dat <- dataL[[idx]]
        samples <- rownames(dat)
        message(cohort)

        results <- list()
        for (jdx in seq_along(samples)) { # foreach sample
                tmp <- dat[jdx,]
                hits <- tmp[tmp >= 1]
                if (length(hits) > 2) {
                        geneTriplets <- apply(combn(names(hits),3), 2, paste, collapse='$')
                for (geneTriplet in geneTriplets) {
                        if (is.null(results[[geneTriplet]])) {
                                results[[geneTriplet]] <- 1
                        } else {
                                results[[geneTriplet]] <- results[[geneTriplet]] + 1
                        }
                }
            }
        }

        if (length(results) > 1) {
                resTmp <- unlist(results)
                df <- data.frame(
                        t1 = stringr::str_replace(sapply(stringr::str_split(names(resTmp), "\\$"), "[[", 1), "PCDH", ""),
                        t2 = stringr::str_replace(sapply(stringr::str_split(names(resTmp), "\\$"), "[[", 2), "PCDH", ""),
                        t3 = stringr::str_replace(sapply(stringr::str_split(names(resTmp), "\\$"), "[[", 3), "PCDH", ""),
                        frequency = as.numeric(resTmp)
                )                
                df <- df[which(df$frequency >= howManyPt),]

                if (nrow(df) > 0) { 
                        tmpDf <- c() # gene indices instead of gene symbols
                        for (jdx in c(1:nrow(df))) {
                                x <- as.matrix(df[jdx,])
                                tmpDf <- rbind(tmpDf, c(pcdhGeneIdx[as.matrix(x[1])], pcdhGeneIdx[as.matrix(x[2])], pcdhGeneIdx[as.matrix(x[3])], as.numeric(x[4])))
                        }

                        for (jdx in pcdhGeneIdx) {
                                res <- c()
                                for (kdx in c(1:nrow(tmpDf))) {
                                        x <- as.matrix(tmpDf[kdx,])
                                        if (x[1] == jdx || x[2] == jdx || x[3] == jdx) {
                                                res <- rbind(res, t(x))
                                        }
                                }

                                if (!is.null(res)) {
                                        if (nrow(res) > 1) {
                                                data <- as.data.frame(res[order(res[,1], res[,2], res[,3], decreasing=c(T,T,T)),])
                                                colnames(data) <- c("T1", "T2", "T3", "Size")
                                        } else {
                                                data <- data.frame(T1 = res[1], T2 = res[2], T3 = res[3], Size = res[4])
                                        }
                                        
                                        data$OPT <- factor(paste0("OPT_", c(1:nrow(res))), levels=paste0("OPT_", c(1:nrow(res))))
                                        data$Size_cut <- cut(data$Size, breaks = c(3,5,10,15,50), include.lowest = TRUE, right = TRUE)
                                        levels_count <- length(levels(data$Size_cut))

                                        manual_size <- seq(1, by = 1, length.out = levels_count)
                                        names(manual_size) <- levels(data$Size_cut)
                                        
                                        goPlot <- ggplot(data) +
                                                geom_hline(yintercept = c(48:1), col=nodeCol, lty=1, lwd=0.5) + 
                                                geom_segment( aes(x=OPT, xend=OPT, y=T1, yend=T2), color="grey20", lwd=0.3) +
                                                geom_segment( aes(x=OPT, xend=OPT, y=T2, yend=T3), color="grey20", lwd=0.3) +
                                                geom_point(data=data, aes(x=OPT, y=T1, size=Size_cut), color=rgb(0, 0, 0, 0.7)) +
                                                geom_point(data=data, aes(x=OPT, y=T2, size=Size_cut), color=rgb(0, 0, 0, 0.7)) +
                                                geom_point(data=data, aes(x=OPT, y=T3, size=Size_cut), color=rgb(0, 0, 0, 0.7)) +
                                                scale_y_continuous(breaks=c(length(pcdhGeneIdx):1), labels=names(pcdhGeneIdx)) +
                                                scale_size_manual(name="No. of Samples", values = manual_size * 2, labels=c("3-5","6-10","11-15",">15")) +
                                                labs(title=paste0("PCDH", names(pcdhGeneIdx)[which(pcdhGeneIdx == jdx)], " in ", cohort),
                                                        x="OPTs", y = "cPCDH") + 
                                                        theme_bw() +
                                                        theme(
                                                        axis.line = element_line(colour = "black"),
                                                        panel.grid.major = element_blank(),
                                                        panel.grid.minor = element_blank(),
                                                        panel.border = element_blank(),
                                                        panel.background = element_blank(),
                                                        axis.text.x=element_blank(), 
                                                        text = element_text(size = 16)
                                                        )

                                        pdf(file.path(resultsDir, paste0("GOplot_", cohort, "_PCDH", names(pcdhGeneIdx)[which(pcdhGeneIdx == jdx)], ".pdf")), width=9, height=9)
                                        print(goPlot)
                                        dev.off()
                                }
                        }
                }
        }
        return(TRUE)
})

q("no")