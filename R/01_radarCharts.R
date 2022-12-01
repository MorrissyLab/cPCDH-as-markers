##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Updated on November 30, 2022
# Written by Heewon Seo (Heewon.Seo@UCalgary.ca)
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
setwd(file.path(getwd(), "R"))

dataDir <- file.path("../data")
resultsDir <- file.path("../results")

profileFileName <- "cPCDH_expr_25_cohorts.rds"

col <- c("#E7B800", "#FC4E07", "black") # GeoMean, SD, MAX

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Wrapper
.plotRadarChart <- function(data, title) {
        fmsb::radarchart(
                data,
                axistype = 1,
                seg = 3,
                title = title,
                centerzero = F,
        # Polygon options
                pcol = col, # line color
                pfcol = c(scales::alpha(col[1], 0.5), scales::alpha(col[1], 0.3), NA), # fill color
                plwd = 1, # line width
                plty = c(1, 1, 2), # line type
        # Grid options:
                cglcol = "grey", # line color
                cglty = 1, # line type
                cglwd = 0.5, # line width
                axislabcol = "grey80", # axis label color
                caxislabels = c("0 (0)", "0.1 (0.33)", "0.2 (0.66)", "0.3 (1.00)"), # axis label
        # Variable options
                vlcex = 0.8, # controls the font size of variable labels
                vlabels = stringr::str_replace(colnames(data), "PCDH", "") # variable labels
        )
}

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Load library
library(EnvStats)
library(scales)
library(fmsb)

dataL <- readRDS(file.path(dataDir, profileFileName))

summaryMatL <- lapply(seq_along(dataL), function(idx) {
        cohort <- names(dataL)[idx]
        dat <- dataL[[idx]]
        message(cohort)

        sdData <- apply(dat, 1, sd, na.rm=T)
        maxData <- apply(dat, 1, max, na.rm=T)
        dat[dat == 0] <- 1.0E-08
        mData <- apply(dat, 1, EnvStats::geoMean, na.rm=T)
        
        cpcdh <- data.frame(Gene = names(mData), mData = mData, sdData = sdData, maxData = scales::rescale(maxData, to=c(0, 0.3), from=c(0, 1)))

        alpha_max_min <- data.frame(matrix(rep(c(0.3, 0), 13), ncol=13, byrow=F)); alphaIdx <- c(1:13)
        beta_max_min <- data.frame(matrix(rep(c(0.3, 0), 16), ncol=16, byrow=F)); betaIdx <- c(33:48)
        gamma_max_min <- data.frame(matrix(rep(c(0.3, 0), 19), ncol=19, byrow=F)); gammaIdx <- c(14:32)

        alpha <- as.data.frame(rbind(alpha_max_min, as.numeric(t(cpcdh$mData[alphaIdx])), as.numeric(t(cpcdh$sdData[alphaIdx])), as.numeric(t(cpcdh$maxData[alphaIdx]))))
        rownames(alpha) <- c("LimitMax", "LimitMin", "Mean", "SD", "Max")
        colnames(alpha) <- c(rownames(cpcdh)[alphaIdx])

        beta <- as.data.frame(rbind(beta_max_min, as.numeric(t(cpcdh$mData[betaIdx])), as.numeric(t(cpcdh$sdData[betaIdx])), as.numeric(t(cpcdh$maxData[betaIdx]))))
        rownames(beta) <- c("LimitMax", "LimitMin", "Mean", "SD", "Max")
        colnames(beta) <- c(rownames(cpcdh)[betaIdx])

        gamma <- as.data.frame(rbind(gamma_max_min, as.numeric(t(cpcdh$mData[gammaIdx])), as.numeric(t(cpcdh$sdData[gammaIdx])), as.numeric(t(cpcdh$maxData[gammaIdx]))))
        rownames(gamma) <- c("LimitMax", "LimitMin", "Mean", "SD", "Max")
        colnames(gamma) <- c(rownames(cpcdh)[gammaIdx])

        pdf(file.path(resultsDir, paste0("radarChart_", cohort, "_", ncol(dat), "_geoMean.pdf")), width=12, height=4)
        par(mfrow=c(1,3), oma=c(0,0,0,0))

                .plotRadarChart(data = alpha, title = "Alpha")
                .plotRadarChart(data = beta, title = "Beta")
                mtext(side = 3, line = 2.5, at = 0, cex = 1.1, paste0(cohort, ", N=", ncol(dat)), font = 2) # Title
                .plotRadarChart(data = gamma, title = "Gamma")
        dev.off()

        return(cpcdh)
})

names(summaryMatL) <- names(dataL)
saveRDS(summaryMatL, file.path(dataDir, "GeoMean_SD_Max.rds"))

q("no")