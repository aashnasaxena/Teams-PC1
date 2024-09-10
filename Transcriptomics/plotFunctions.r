### Functions to calculate PCA metrics including gene coefficients if desired,
## generating housekeeping genes as random genes in the network,
## and replacing genes with housekeeping genes to see how the PCA metrics change

geneSetDiagnostics <- function(nodeReplacementDf, pc1Var, 
        pcaDat, pcaDatPlot) {
    dir.create("Diagnostics", showWarnings = F)
    nodeRepMeanSD <- nodeReplacementDf %>% group_by(numNodes) %>%
        summarise(across(everything(), list(mean = mean, sd = sd)))
    # hk coeff vs left coeff
    p1 <- ggplot(nodeRepMeanSD, aes(x = hkCoeffs_mean, y = leftCoeffs_mean, color = nGenes_mean)) +
        geom_point() +
        geom_errorbar(aes(ymin = leftCoeffs_mean - leftCoeffs_sd, 
            ymax = leftCoeffs_mean + leftCoeffs_sd), width = 0.1) +
        geom_errorbarh(aes(xmin = hkCoeffs_mean - hkCoeffs_sd,
            xmax = hkCoeffs_mean + hkCoeffs_sd), height = 0.1) +
            scale_color_viridis_c()+
        theme_Publication() +
        theme(legend.position = "none")+
        labs(x = "HK Coefficients", y = "Left Coefficients")
    # hk PC1 vs left PC1
    p2 <- ggplot(nodeRepMeanSD, aes(x = hkPC1_mean, y = leftPC1_mean, color = nGenes_mean)) +
        geom_point() +
        geom_errorbar(aes(ymin = leftPC1_mean - leftPC1_sd, 
            ymax = leftPC1_mean + leftPC1_sd), width = 0.1) +
        geom_errorbarh(aes(xmin = hkPC1_mean - hkPC1_sd,
            xmax = hkPC1_mean + hkPC1_sd), height = 0.1) +
        scale_color_viridis_c()+
        theme_Publication() +
        theme(legend.position = "none") +
        xlim(c(0, 1)) + ylim(c(0, 1)) +
        labs(x = "HK PC1", y = "Left PC1")
    # hk cor vs left cor
    p6 <- ggplot(nodeRepMeanSD, aes(x = hkCor_mean, y = leftCor_mean, color = nGenes_mean)) +
        geom_point() +
        geom_errorbar(aes(ymin = leftCor_mean - leftCor_sd, 
            ymax = leftCor_mean + leftCor_sd), width = 0.1) +
        geom_errorbarh(aes(xmin = hkCor_mean - hkCor_sd,
            xmax = hkCor_mean + hkCor_sd), height = 0.1) +
        scale_color_viridis_c() +
        theme_Publication() +
        theme(legend.position = "none") +
        xlim(c(0, 1)) + ylim(c(0, 1)) +
        labs(x = "HK Correlation", y = "Left Correlation")
    
    # hk pc1 across numNodes

    p7 <- ggplot(nodeReplacementDf, aes(x = numNodes, y = hkPC1)) +
        geom_boxplot() +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = rel(textSize)), 
            axis.title.x = element_text(size = rel(0.9))) +
        ylim(c(0, 1)) +
        labs(y = "HK PC1", x = "# Nodes replaced")
    # left pc1 across numNodes
    p8 <- ggplot(nodeReplacementDf, aes(x = numNodes, y = leftPC1)) +
        geom_boxplot() +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = rel(textSize)), 
            axis.title.x = element_text(size = rel(0.9))) +
        ylim(c(0, 1)) +
        labs(y = "Left PC1", x = "# Nodes replaced")
    # replaced pc1 across numNodes
    p9 <- ggplot(nodeReplacementDf, aes(x = numNodes, y = repPC1)) +
        geom_boxplot() +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = rel(textSize)), 
            axis.title.x = element_text(size = rel(0.9))) +
        ylim(c(0, 1)) +
        labs(y = "Replaced PC1", x = "# Nodes replaced")
    # replaced cor across numNodes
    p10 <- ggplot(nodeReplacementDf, aes(x = numNodes, y = repCor)) +
        geom_boxplot() +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = rel(textSize)), 
            axis.title.x = element_text(size = rel(0.9))) +
        ylim(c(0, 1)) +
        labs(y = "Replaced Correlation", x = "# Nodes replaced")
    
    png(paste0("Diagnostics/", key, "_diagnosticPlots.png"),
        width = 2*w + (w-1.5)/2 , height = 2*(w-0.5), units = "in", res = 300)
    gridExtra::grid.arrange(p1, p2, p6, p7, p8, p9, p10, nrow = 3, widths = c(0.9, 0.9, 0.9),
        top = textGrob(paste0(key, "; PC1 : ", round(pc1Var, 2),
        "; Number of PC axes : ", numPC), 
            gp = gpar(fontsize = 24, fontface = "bold")))
    dev.off()

    ### Which genes increse PC1 the most?
    # browser()
    basePC1 <- pc1Var
    nodes <- pcaDat$Gene
    repSet <- nodeReplacementDf %>% pull(nGenes) %>% unique %>% sort
    probabilities <- lapply(nodes, function(node) {
        dReplace <- sapply(repSet, function(r) {
            conditionalProb(nodeReplacementDf %>% 
                filter(nGenes == r), basePC1, node, replaced = T)
        })
        dKeep <- sapply(repSet, function(r) {
            conditionalProb(nodeReplacementDf %>% 
                filter(nGenes == r), basePC1, node, replaced = F)
        })
        dAll <- data.frame(Probabilities = c(dReplace, dKeep),
            nGenes = rep(repSet, 2),
            Condition = rep(c("Replaced", "Kept"), each = length(repSet)),
            Gene = node)
        return(dAll)
    }) %>% bind_rows
    probabilities <- probabilities %>%
        full_join(pcaDatPlot, by = "Gene")
    ggplot(probabilities, aes(x = Coeff, y = Probabilities, color = Condition)) +
        geom_point() + geom_line() + 
        facet_wrap(~nGenes, scales ="fixed") +
        theme_Publication() +
        xlim(c(-1, 1)) + ylim(c(0, 1)) +
        labs(x = "PC1 Coefficient", y = "Probability of increase in PC1")
    ggsave(paste0("Diagnostics/", key, "_probabilities.png"), 
        width = 15, height = 15)
}


plotData <- function(df, dfHk,key, fld = ".", w = 7.5, 
    reduceGenes = T, diagnostics = F) {
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
        colnames(df), hkGenes = hkGenes, diagnostics = diagnostics) %>%
        mutate(numNodes = factor(nGenes, levels = nGenes %>% unique %>% sort))
    ### Standrard plots
    # textSize <- min(0.9, 5*(w-2.5)/ncol(df))
    # p3 <- ggplot(pcaDatPlot, aes(y = reorder(Gene, -Coeff), x = Coeff)) +
    #     geom_bar(stat = "identity") +
    #     theme_Publication() +
    #     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    #         axis.text = element_text(size = rel(textSize)), 
    #         axis.title.x = element_text(size = rel(0.9)),
    #         legend.position = "none") +
    #     labs(y = "", x = "PC1 Coefficient") +
    #     scale_fill_viridis_d()
    # textSize <- min(0.9, 5*(w-2.5)/nlevels(nodeReplacementDf$numNodes))
    # p4 <- ggplot(nodeReplacementDf, aes(x = numNodes, y = pc1Var)) +
    #     geom_boxplot() +
    #     theme_Publication() +
    #     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    #         axis.text = element_text(size = rel(textSize)), 
    #         axis.title.x = element_text(size = rel(0.9))) +
    #     labs(y = "PC1 Variance", x = "# Nodes replaced")
    # p5 <- ggplot(nodeReplacementDf, aes(x = numNodes, y = numPC)) +
    #     geom_boxplot() +
    #     theme_Publication() +
    #     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    #         axis.text = element_text(size = rel(textSize)), 
    #         axis.title.x = element_text(size = rel(0.9))) +
    #     labs(y = "# PC Axes", x = "# Nodes replaced")

    # png(paste0(key, "_nodeReplacePlots.png"), 
    #     width = 2*w + (w-1.5)/2 , height = w-0.5, units = "in", res = 300)
    # gridExtra::grid.arrange(p3, p4, p5, nrow = 1, widths = c(0.9, 0.9, 0.9),
    #     top = textGrob(paste0(key, "; PC1 : ", round(pc1Var, 2),
    #     "; Number of PC axes : ", numPC), 
    #         gp = gpar(fontsize = 24, fontface = "bold")))
    # dev.off()
    write_csv(nodeReplacementDf, paste0(key, "_nodeReplaceDf.csv"))

    ### PC1 stability analysis -------------------------
    dSummarised <- pc1StabilityAnalysis(nodeReplacementDf)


    ### diagnostic plots--------------------------------
    if (diagnostics) {
        geneSetDiagnostics(nodeReplacementDf, pc1Var, pcaDat, pcaDatPlot)
    }
    setwd(wd)
    return(dSummarised)
}

plotSummaries <- function(df, figDir) {
    # browser()
    if (!dir.exists(figDir))
        dir.create(figDir, recursive = T)
    df <- df %>% mutate(slopeDiff = abs(slopeLess) - abs(slopeMore))
    ggplot(df %>% filter(metrics == "PC1"), 
        aes(x = halfMax, y = (slopeDiff), color = gsVal)) +
        geom_point() + #geom_line(aes(group = gsVal)) + 
        scale_color_viridis_c() +
        theme_Publication() +
        labs(x = "Half Maximum", y = "Slope Difference", title = "Mean PC1",
            color = "Gs PC1 Variance") +
        theme(legend.position = "bottom", legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(angle = 90))
    ggsave(paste0(figDir, "/HalfMaxSlopeDiffPC1mean.png"), width = 5.5, height = 6)
    ggplot(df %>% filter(metrics == "numPC"), 
        aes(x = halfMax, y = (slopeDiff), color = gsVal)) +
        geom_point() + #geom_line(aes(group = gsVal)) + 
        scale_color_viridis_c() +
        theme_Publication() +
        labs(x = "Half Maximum", y = "Slope Difference", 
        title = "Mean # PC axes", color = "Gs # PC axes") +
        theme(legend.position = "bottom", legend.key.width = unit(0.5, "cm"),
            legend.text = element_text(angle = 90))
    ggsave(paste0(figDir, "/HalfMaxSlopeDiffNumPCmean.png"), width = 5.5, height = 6)
    ggplot(df %>% filter(metrics == "PC1median"), aes(x = halfMax, y = (slopeDiff), color = gsVal)) +
        geom_point() + geom_line(aes(group = gsVal)) + 
        scale_color_viridis_c() +
        theme_Publication() +
        labs(x = "Half Maximum", y = "Slope Difference", title = "Median PC1") +
        theme(legend.position = "bottom", legend.key.width = unit(0.5, "cm"))
    ggsave(paste0(figDir, "/HalfMaxSlopeDiffPC1median.png"), width = 5.5, height = 6)
    ggplot(df %>% filter(metrics == "numPCmedian"), aes(x = halfMax, y = (slopeDiff), color = gsVal)) +
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
    df <- df %>% mutate(slopeDiff = abs(slopeLess) - abs(slopeMore))
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
        p <- ggplot(d, aes(x = dataSet, y = geneSet, fill = slopeDiff)) +
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


