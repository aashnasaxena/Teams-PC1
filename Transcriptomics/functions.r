### Functions to calculate PCA metrics including gene coefficients if desired,
## generating housekeeping genes as random genes in the network,
## and replacing genes with housekeeping genes to see how the PCA metrics change
pc1StabilityAnalysis <- function(df) {
    # browser()
    dSummarised <- df %>%
        group_by(nGenes) %>%
        summarise(PC1 = mean(pc1Var), numPC = mean(numPC), 
            # PC1SD = sd(pc1Var), numPCSD = sd(numPC),
            PC1median = median(pc1Var), numPCmedian = median(numPC)) %>%
        arrange(nGenes) %>%
        mutate(nGenes = nGenes/max(nGenes))
    halMax <- sapply(colnames(dSummarised)[-1], function(x) {
        d <- dSummarised[[x]]
        nGenes <- dSummarised$nGenes
        M <- d[1]
        m <- d[length(d)]
        hm <- m + (M-m)/2
        idless <- which(d < hm)
        if (m > M)
            idless <- which(d > hm)
        idless <- idless[1]-1
        halfMax <- nGenes[idless]/max(nGenes)
        if (m > M) {
            halfMax <- -1*halfMax
        }
        slopeMore <- (d[idless] - d[1])/nGenes[idless]
        slopeLess <- (d[length(d)] - d[idless])/nGenes[length(d)] - nGenes[idless]
        res <- c(halfMax = halfMax, slopeMore = slopeMore, slopeLess = slopeLess, gsVal = M, hkVal = m)
        # print(res)
        return(res)
    }) %>% t %>% as.data.frame %>%
        mutate(metrics = colnames(dSummarised)[-1])
    return(halMax)
}


conditionalProb <- function(df, ref, node, replaced = F) {
    if (replaced)
        df <- df %>% select(repGenes, pc1Var) %>% set_names(c("Genes", "PC1"))
    else
        df <- df %>% select(leftGenes, pc1Var) %>% set_names(c("Genes", "PC1"))
    prob <- df %>% filter(str_detect(Genes, node)) %>%
        mutate(Increase = PC1 > ref) %>%
        pull(Increase)
    prob <- sum(prob)/length(prob)
    return(prob)
}

dTrans <- function(df, gs) {
    d <- df %>% filter(Genes %in% gs)
    v <- d$Genes
    d <- d %>% select(-Genes) %>% t
    d <- data.frame(d)
    colnames(d) <- v
    d
}

getPCA <- function(df, gs = NULL, getTeams = F) {
    if (is.null(gs))
        gs <- colnames(df)
    df <- df %>% select(all_of(gs))
    # df <- df %>% mutate_all(function(x){log2(x + 1)})
    pca <- prcomp(df)
    varDat <- pca$sdev^2
    varDat <- varDat / sum(varDat)
    pc1Var <- varDat[1]
    pc2Var <- varDat[2]
    numPC <- sum(cumsum(varDat) < 0.9) + 1
    if (!getTeams)
        return(c(pc1Var, pc2Var, numPC))
    pc1 <- pca$rotation[, 1]
    d <- data.frame(Gene = names(pc1), Coeff = pc1) %>% arrange(Coeff)
    return(list(numbers = c(pc1Var, pc2Var, numPC), Coefficients = d))
}

generateHouseKeeping <- function(df, gs = NULL, nGenes = 100, 
    random = T) {
    if (is.null(gs))
        gs <- colnames(df)
    dVec <- df %>% select(all_of(gs)) %>% unlist %>% 
        as.numeric %>% sort(decreasing = T)
    eliquote <- ceiling(0.1*length(dVec))
    dVec <- dVec[1:eliquote]
    if (random) dVec <- runif(100)
    hkGenes <- sample(dVec, nrow(df)*nGenes, replace = T) %>% 
        matrix(ncol = nGenes, byrow = T) %>%
        as.data.frame %>%
        set_names(paste0("gene", 1:nGenes))
    return(hkGenes)
}

geneReplace <- function(df, geneSet, hkGenes = NULL, nIters = 100) {
    if (is.null(hkGenes)) {
        hkDat <- generateHouseKeeping(df, gs = geneSet, nGenes = length(geneSet))
        hkGenes <- colnames(hkDat)
    }
    else {
        hkGenes <- hkGenes[hkGenes %in% colnames(df)]
        hkDat <- df %>% select(any_of(hkGenes))
    }
    df <- df %>% select(all_of(geneSet))
    
    repSet <- 1:length(geneSet)
    if (length(geneSet)>50)
        repSet <- c(seq(1, length(geneSet), by = 5), length(geneSet))
    d <- lapply(repSet, function(x) {
        d1 <- sapply(1:nIters, function(it) {
            repCor <- hkCor <- leftCor <- hkPC1 <- repPC1 <- leftPC1 <- hkCoeffs <- leftCoeffs <- 0
            replaceSet <- sample(geneSet, x, replace=F)
            leftSet <- geneSet[!(geneSet %in% replaceSet)]
            hkSet <- sample(hkGenes, x, replace = T)
            dNew <- df %>% select(-all_of(replaceSet)) %>%
                bind_cols(hkDat %>% select(all_of(hkSet)))
            pcaDat <- getPCA(dNew, getTeams = F)
            # if (length(replaceSet)) {
            #     repCor <- df %>% select(all_of(replaceSet)) %>% cor %>% abs %>% sum
            #     repCor <- repCor/(x*x)
            #     repPC1 <- getPCA(df %>% select(all_of(replaceSet)), getTeams = F)[1]
            # }
            
            # if (length(hkSet)) {
            #     hkCor <- hkDat %>% select(all_of(hkSet)) %>% cor %>% abs %>% sum
            #     hkPC1 <- getPCA(hkDat %>% select(all_of(hkSet)), getTeams = F)[1]
            #     hkCor <- hkCor/(x*x)
            #     hkCoeffs <- pcaDat$Coefficients %>% filter(Gene %in% hkSet) %>% 
            #     pull(Coeff) %>% sum
            # }
            
            # if (x != ncol(df)) {
            #     leftCor <- df %>% select(-all_of(replaceSet)) %>% cor %>% abs %>% sum
            #     leftCor <- leftCor/((ncol(df) - x)*(ncol(df) - x))
            #     leftPC1 <- getPCA(df %>% select(-all_of(replaceSet)), getTeams = F)[1]
            #     leftCoeffs <- pcaDat$Coefficients %>% filter(Gene %in% geneSet) %>% 
            #         pull(Coeff) %>% sum
            # }
            # pcaDat <- pcaDat$numbers
            if (is.character(pcaDat)) {
                browser()
            }
            return(c(pc1Var = pcaDat[1], pc2Var = pcaDat[2], numPC = pcaDat[3],
                hkCor = hkCor, repCor = repCor, repPC1 = repPC1, leftCor = leftCor,
                hkPC1 = hkPC1, leftPC1 = leftPC1, hkCoeffs = hkCoeffs,
                leftCoeffs = leftCoeffs, repGenes = paste0(replaceSet, collapse = "_"),
                leftGenes = paste0(leftSet, collapse = "_")))
        }) %>% t %>% as.data.frame %>%
            mutate(nGenes = x)
        return(d1)
    }) %>% bind_rows
    sapply(1:(ncol(d) - 3), function(x) {
        d[[x]] <<- as.numeric(d[[x]])
    })
    pcWild <- df %>% select(all_of(geneSet)) %>% getPCA()
    if (is.character(pcWild)) {
        browser()
    }
    d1 <- data.frame(pc1Var = pcWild[1], pc2Var = pcWild[2], numPC = pcWild[3], nGenes = 0)
    if (class(d$pc1Var) == "character") {
        browser()
    }
    if (class(d1$pc1Var) == "character") {
        browser()
    }
    return(d %>% bind_rows(d1))
}

plotData <- function(df, dfHk,key, fld = ".", w = 7.5, reduceGenes = T) {
    # browser()
    wd <- getwd()
    if (!dir.exists(fld))
        dir.create(fld, recursive = T)
    setwd(fld)
    hkGenes <- colnames(dfHk)
    
    pca <- getPCA(df, getTeams = T)
    pc1Var <- pca$numbers[1]
    numPC <- pca$numbers[3]
    pcaDat <- pca$Coefficients
    if (reduceGenes && ncol(df) > 60) {
        # browser()
        pcaDat <- pcaDat %>% filter (abs(Coeff) < 0.9) %>% 
            mutate(AbsCoeff = abs(Coeff)) %>%
            arrange(desc(AbsCoeff)) %>% head(60)
        GenestoUse <- pcaDat$Gene
        df <- df %>% select(all_of(GenestoUse))
        pca <- getPCA(df, getTeams = T)
        pc1Var <- pca$numbers[1]
        numPC <- pca$numbers[3]
        pcaDat <- pca$Coefficients
    }
    pcaDatPlot <- pcaDat %>% 
        arrange(desc(Coeff))
    nodeReplacementDf <- geneReplace(cbind.data.frame(df, dfHk), 
        colnames(df), hkGenes = hkGenes) %>%
        mutate(numNodes = factor(nGenes, levels = nGenes %>% unique %>% sort))
    ### Standrard plots
    textSize <- min(0.9, 5*(w-2.5)/ncol(df))
    p3 <- ggplot(pcaDatPlot, aes(y = reorder(Gene, -Coeff), x = Coeff)) +
        geom_bar(stat = "identity") +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = rel(textSize)), 
            axis.title.x = element_text(size = rel(0.9)),
            legend.position = "none") +
        labs(y = "", x = "PC1 Coefficient") +
        scale_fill_viridis_d()
    textSize <- min(0.9, 5*(w-2.5)/nlevels(nodeReplacementDf$numNodes))
    p4 <- ggplot(nodeReplacementDf, aes(x = numNodes, y = pc1Var)) +
        geom_boxplot() +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = rel(textSize)), 
            axis.title.x = element_text(size = rel(0.9))) +
        labs(y = "PC1 Variance", x = "# Nodes replaced")
    p5 <- ggplot(nodeReplacementDf, aes(x = numNodes, y = numPC)) +
        geom_boxplot() +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = rel(textSize)), 
            axis.title.x = element_text(size = rel(0.9))) +
        labs(y = "# PC Axes", x = "# Nodes replaced")

    png(paste0(key, "_nodeReplacePlots.png"), 
        width = 2*w + (w-1.5)/2 , height = w-0.5, units = "in", res = 300)
    gridExtra::grid.arrange(p3, p4, p5, nrow = 1, widths = c(0.9, 0.9, 0.9),
        top = textGrob(paste0(key, "; PC1 : ", round(pc1Var, 2),
        "; Number of PC axes : ", numPC), 
            gp = gpar(fontsize = 24, fontface = "bold")))
    dev.off()

    ### PC1 stability analysis -------------------------
    dSummarised <- pc1StabilityAnalysis(nodeReplacementDf)


    ### diagnostic plots--------------------------------
    # dir.create("Diagnostics", showWarnings = F)
    # nodeRepMeanSD <- nodeReplacementDf %>% group_by(numNodes) %>%
    #     summarise(across(everything(), list(mean = mean, sd = sd)))
    # # hk coeff vs left coeff
    # p1 <- ggplot(nodeRepMeanSD, aes(x = hkCoeffs_mean, y = leftCoeffs_mean, color = nGenes_mean)) +
    #     geom_point() +
    #     geom_errorbar(aes(ymin = leftCoeffs_mean - leftCoeffs_sd, 
    #         ymax = leftCoeffs_mean + leftCoeffs_sd), width = 0.1) +
    #     geom_errorbarh(aes(xmin = hkCoeffs_mean - hkCoeffs_sd,
    #         xmax = hkCoeffs_mean + hkCoeffs_sd), height = 0.1) +
    #         scale_color_viridis_c()+
    #     theme_Publication() +
    #     theme(legend.position = "none")+
    #     labs(x = "HK Coefficients", y = "Left Coefficients")
    # # hk PC1 vs left PC1
    # p2 <- ggplot(nodeRepMeanSD, aes(x = hkPC1_mean, y = leftPC1_mean, color = nGenes_mean)) +
    #     geom_point() +
    #     geom_errorbar(aes(ymin = leftPC1_mean - leftPC1_sd, 
    #         ymax = leftPC1_mean + leftPC1_sd), width = 0.1) +
    #     geom_errorbarh(aes(xmin = hkPC1_mean - hkPC1_sd,
    #         xmax = hkPC1_mean + hkPC1_sd), height = 0.1) +
    #     scale_color_viridis_c()+
    #     theme_Publication() +
    #     theme(legend.position = "none") +
    #     xlim(c(0, 1)) + ylim(c(0, 1)) +
    #     labs(x = "HK PC1", y = "Left PC1")
    # # hk cor vs left cor
    # p6 <- ggplot(nodeRepMeanSD, aes(x = hkCor_mean, y = leftCor_mean, color = nGenes_mean)) +
    #     geom_point() +
    #     geom_errorbar(aes(ymin = leftCor_mean - leftCor_sd, 
    #         ymax = leftCor_mean + leftCor_sd), width = 0.1) +
    #     geom_errorbarh(aes(xmin = hkCor_mean - hkCor_sd,
    #         xmax = hkCor_mean + hkCor_sd), height = 0.1) +
    #     scale_color_viridis_c() +
    #     theme_Publication() +
    #     theme(legend.position = "none") +
    #     xlim(c(0, 1)) + ylim(c(0, 1)) +
    #     labs(x = "HK Correlation", y = "Left Correlation")
    
    # # hk pc1 across numNodes

    # p7 <- ggplot(nodeReplacementDf, aes(x = numNodes, y = hkPC1)) +
    #     geom_boxplot() +
    #     theme_Publication() +
    #     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    #         axis.text = element_text(size = rel(textSize)), 
    #         axis.title.x = element_text(size = rel(0.9))) +
    #     ylim(c(0, 1)) +
    #     labs(y = "HK PC1", x = "# Nodes replaced")
    # # left pc1 across numNodes
    # p8 <- ggplot(nodeReplacementDf, aes(x = numNodes, y = leftPC1)) +
    #     geom_boxplot() +
    #     theme_Publication() +
    #     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    #         axis.text = element_text(size = rel(textSize)), 
    #         axis.title.x = element_text(size = rel(0.9))) +
    #     ylim(c(0, 1)) +
    #     labs(y = "Left PC1", x = "# Nodes replaced")
    # # replaced pc1 across numNodes
    # p9 <- ggplot(nodeReplacementDf, aes(x = numNodes, y = repPC1)) +
    #     geom_boxplot() +
    #     theme_Publication() +
    #     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    #         axis.text = element_text(size = rel(textSize)), 
    #         axis.title.x = element_text(size = rel(0.9))) +
    #     ylim(c(0, 1)) +
    #     labs(y = "Replaced PC1", x = "# Nodes replaced")
    # # replaced cor across numNodes
    # p10 <- ggplot(nodeReplacementDf, aes(x = numNodes, y = repCor)) +
    #     geom_boxplot() +
    #     theme_Publication() +
    #     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    #         axis.text = element_text(size = rel(textSize)), 
    #         axis.title.x = element_text(size = rel(0.9))) +
    #     ylim(c(0, 1)) +
    #     labs(y = "Replaced Correlation", x = "# Nodes replaced")
    
    # png(paste0("Diagnostics/", key, "_diagnosticPlots.png"),
    #     width = 2*w + (w-1.5)/2 , height = 2*(w-0.5), units = "in", res = 300)
    # gridExtra::grid.arrange(p1, p2, p6, p7, p8, p9, p10, nrow = 3, widths = c(0.9, 0.9, 0.9),
    #     top = textGrob(paste0(key, "; PC1 : ", round(pc1Var, 2),
    #     "; Number of PC axes : ", numPC), 
    #         gp = gpar(fontsize = 24, fontface = "bold")))
    # dev.off()

    # ### Which genes increse PC1 the most?
    # # browser()
    # basePC1 <- pc1Var
    # nodes <- pcaDat$Gene
    # repSet <- nodeReplacementDf %>% pull(nGenes) %>% unique %>% sort
    # probabilities <- lapply(nodes, function(node) {
    #     dReplace <- sapply(repSet, function(r) {
    #         conditionalProb(nodeReplacementDf %>% 
    #             filter(nGenes == r), basePC1, node, replaced = T)
    #     })
    #     dKeep <- sapply(repSet, function(r) {
    #         conditionalProb(nodeReplacementDf %>% 
    #             filter(nGenes == r), basePC1, node, replaced = F)
    #     })
    #     dAll <- data.frame(Probabilities = c(dReplace, dKeep),
    #         nGenes = rep(repSet, 2),
    #         Condition = rep(c("Replaced", "Kept"), each = length(repSet)),
    #         Gene = node)
    #     return(dAll)
    # }) %>% bind_rows
    # probabilities <- probabilities %>%
    #     full_join(pcaDatPlot, by = "Gene")
    # ggplot(probabilities, aes(x = Coeff, y = Probabilities, color = Condition)) +
    #     geom_point() + geom_line() + 
    #     facet_wrap(~nGenes, scales ="fixed") +
    #     theme_Publication() +
    #     xlim(c(-1, 1)) + ylim(c(0, 1)) +
    #     labs(x = "PC1 Coefficient", y = "Probability of increase in PC1")
    # ggsave(paste0("Diagnostics/", key, "_probabilities.png"), 
    #     width = 15, height = 15)
    ###### diagnosis end -------------------------------
    setwd(wd)
    return(dSummarised)
}

plotSummaries <- function(df, figDir) {
    # browser()
    if (!dir.exists(figDir))
        dir.create(figDir, recursive = T)
    ggplot(df %>% filter(metrics == "PC1"), aes(x = halfMax, y = (slopeLess - slopeMore), color = gsVal)) +
            geom_point() + #geom_line(aes(group = gsVal)) + 
            scale_color_viridis_c() +
            theme_Publication() +
            labs(x = "Half Maximum", y = "Slope Difference", title = "Mean PC1",
                color = "Gs PC1 Variance") +
            theme(legend.position = "bottom", legend.key.width = unit(0.5, "cm"),
            legend.text = element_text(angle = 90))
    ggsave(paste0(figDir, "/HalfMaxSlopeDiffPC1mean.png"), width = 5.5, height = 6)
    ggplot(df %>% filter(metrics == "numPC"), 
        aes(x = halfMax, y = (slopeLess - slopeMore), color = gsVal)) +
        geom_point() + #geom_line(aes(group = gsVal)) + 
        scale_color_viridis_c() +
        theme_Publication() +
        labs(x = "Half Maximum", y = "Slope Difference", 
        title = "Mean # PC axes", color = "Gs # PC axes") +
        theme(legend.position = "bottom", legend.key.width = unit(0.5, "cm"),
            legend.text = element_text(angle = 90))
    ggsave(paste0(figDir, "/HalfMaxSlopeDiffNumPCmean.png"), width = 5.5, height = 6)
    ggplot(df %>% filter(metrics == "PC1median"), aes(x = halfMax, y = (slopeLess - slopeMore), color = gsVal)) +
        geom_point() + geom_line(aes(group = gsVal)) + 
        scale_color_viridis_c() +
        theme_Publication() +
        labs(x = "Half Maximum", y = "Slope Difference", title = "Median PC1") +
        theme(legend.position = "bottom", legend.key.width = unit(0.5, "cm"))
    ggsave(paste0(figDir, "/HalfMaxSlopeDiffPC1median.png"), width = 5.5, height = 6)
    ggplot(df %>% filter(metrics == "numPCmedian"), aes(x = halfMax, y = (slopeLess - slopeMore), color = gsVal)) +
        geom_point() + geom_line(aes(group = gsVal)) + 
        scale_color_viridis_c() +
        theme_Publication() +
        labs(x = "Half Maximum", y = "Slope Difference", title = "Median Number of PC axes") +
        theme(legend.position = "bottom", legend.key.width = unit(0.5, "cm"))
    ggsave(paste0(figDir, "/HalfMaxSlopeDiffNumPCmedian.png"), width = 5.5, height = 6)

    ggplot(df %>% filter(metrics == "PC1"), aes(x = gsVal, y = halfMax)) +
        geom_point() + geom_smooth() + 
        theme_Publication() +
        labs(x = "Gene set PC1", y = "Half Maximum", title = "Mean PC1")
    ggsave(paste0(figDir, "/PC1vsHalfMaxMean.png"), width = 5.5, height = 5.3)


    ggplot(df %>% filter(metrics == "numPC"), aes(x = gsVal, y = halfMax)) +
        geom_point() + geom_smooth() + 
        theme_Publication() +
        labs(x = "Gene set PC1", y = "Half Maximum", title = "Mean Number of PC axes")
    ggsave(paste0(figDir, "/NumPCvsHalfMaxMean.png"), width = 5.5, height = 5.3)
}

plotSummariesDataSets <- function(df, figDir) {
    # browser()
    if (!dir.exists(figDir))
        dir.create(figDir, recursive = T)
    ### heatmaps for each metric for geneset-dataset combinations
    metrics <- df %>% pull(metrics) %>% unique
    sapply(metrics, function(met) {
        d <- df %>% filter(metrics == met)
        p <- ggplot(d, aes(x = dataSet, y = geneSet, fill = halfMax)) +
            geom_tile() +
            scale_fill_viridis_c() +
            theme_Publication() +
            labs(x = "Data Set", y = "Gene Set", fill = "Half Maximum") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
                legend.position = "right", legend.direction = "vertical", 
                legend.key.height = unit(0.7, "cm"))
        ggsave(paste0(figDir, "/", met, "HeatmaphalfMax.png"), width = 5.5, height = 7)
        p <- ggplot(d, aes(x = dataSet, y = geneSet, fill = slopeLess - slopeMore)) +
            geom_tile() +
            scale_fill_gradient2(midpoint = 0) +
            theme_Publication() +
            labs(x = "Data Set", y = "Gene Set", fill = "Slope Difference") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
                legend.position = "right", legend.direction = "vertical", 
                legend.key.height = unit(0.7, "cm"))
        ggsave(paste0(figDir, "/", met, "HeatmapSlopeDiff.png"), width = 5.5, height = 7)
        p <- ggplot(d, aes(x = dataSet, y = geneSet, fill = gsVal)) +
            geom_tile() +
            scale_fill_viridis_c() +
            theme_Publication() +
            labs(x = "Data Set", y = "Gene Set", fill = "Gene Set Value") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
                legend.position = "right", legend.direction = "vertical", 
                legend.key.height = unit(0.7, "cm"))
        ggsave(paste0(figDir, "/", met, "HeatmapgsVal.png"), width = 5.5, height = 7)
    })
}

dataSetAnalysis <- function(datasetList, genesetList, 
    dataKey = NULL, figureDir = ".", hkGenes = NULL, 
    HKthreshold = 0.5, reduceGenes = F) {
    if (!dir.exists(figureDir))
        dir.create(figureDir, recursive = T)
    if (!is.list(datasetList))
        {
            if (!is.null(dataKey))
                datasetList <- list(dataKey = datasetList)
            else {
                stop("Please provide a key for the dataset")
            }
        }
    dataSets <- names(datasetList)
    geneSets <- names(genesetList)
    dAll <- lapply(dataSets, function(data) {
        print(paste0("Working on ", data))
        figDir <- paste0(figureDir, "/", data)
        dir.create(figDir, showWarnings = F)
        dSummaries <- lapply(geneSets, function(geneKey) {
            print(paste0("  Working on ", geneKey))
            df <- datasetList[[data]]
            dfFull <- df
            gs <- genesetList[[geneKey]]
            hkGenesUpd <- hkGenes[!(hkGenes %in% gs)]
            if ("Genes" %in% colnames(df)) {
                gs <- gs[gs %in% df$Genes]
                df  <- df %>% filter(Genes %in% gs)
                v <- df$Genes
                df <- df %>% select(-Genes) %>% t %>% data.frame %>% set_names(v)
                if (is.null(hkGenes))
                    dfHk <- generateHouseKeeping(df, nGenes = max(100, length(gs)))
                else {
                    dfHk <- dfFull %>% filter(Genes %in% hkGenesUpd)
                    v <- dfHk$Genes
                    dfHk <- dfHk %>% select(-Genes) %>% t %>% data.frame %>% set_names(v)
                }
            }
            else {
                gs <- gs[gs %in% colnames(df)]
                df <- df %>% select(all_of(gs))
                if (is.null(hkGenes))
                    dfHk <- generateHouseKeeping(df, nGenes = max(100, length(gs)))
                else
                    dfHk <- dfFull %>% select(any_of(hkGenes))
            }
            
            if (reduceGenes && ncol(df) > 60) {
                pcaDat <- getPCA(df, getTeams = T)
                pcaDat <- pcaDat$Coefficients
                pcaDat <- pcaDat %>% filter (abs(Coeff) < 0.9) %>% 
                    mutate(AbsCoeff = abs(Coeff)) %>%
                    arrange(desc(AbsCoeff)) %>% head(60)
                GenestoUse <- pcaDat$Gene
                df <- df %>% select(all_of(GenestoUse))
                gs <- GenestoUse
            }

            pc1Gs <- getPCA(df, gs, getTeams = F)[1]
            pc1Hk <- getPCA(dfHk, getTeams = F)[1]
            iter <- 0
            hkSample <- colnames(dfHk)
            while(pc1Hk > HKthreshold*pc1Gs && iter < 1000) {
                # browser()
                nSample <- length(gs)*2
                hkSample <- sample(colnames(dfHk), nSample, replace = F)
                dfH <- dfHk %>% select(all_of(hkSample))
                pc1Hk <- getPCA(dfH, getTeams = F)[1]
                if (iter > 900) {
                    if (HKthreshold < 0.9) {
                        HKthreshold <- 0.9
                        iter <- 0
                    }
                }
                if (iter == 999) {
                    print(paste0(data, "_", geneKey, ": Could not find housekeeping genes"))
                }
                iter <- iter + 1
            }
            # print(pc1Gs)
            # print(pc1Hk)
            # print(length(hkSample))
            # print(HKthreshold)
            dfH <- dfHk %>% select(all_of(hkSample))
            plotData(df, dfH,
                key = paste0(data, "_", geneKey), fld = figDir) %>% 
                mutate(geneSet = geneKey, dataSet = data)
        }) %>% bind_rows()
        # write.csv(dSummaries, paste0(figDir, "/", paste0(data, "_", geneKey), "Summary.csv"))
        # plotSummaries(dSummaries, figDir)
        dSummaries
    }) %>% bind_rows()
    write_csv(dAll, paste0(figureDir, "/Summary.csv"))
    summaryDir <- paste0(figureDir, "/SummaryHeatmaps")
    plotSummariesDataSets(dAll, paste0(figureDir, "/SummaryHeatmaps"))
    sapply(dataSets, function(data) {
        plotSummaries(dAll %>% filter(dataSet == data), 
            paste0(summaryDir, "/", data))
    })
}

