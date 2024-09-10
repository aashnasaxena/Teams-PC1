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
        slopeLess <- (d[length(d)] - d[idless])/(nGenes[length(d)] - nGenes[idless])
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
    d <- df %>% filter(Genes %in% gs) %>%
        filter(!is.na(Genes)) %>%
        filter(!duplicated(Genes))
    v <- d$Genes
    v <- v[!is.na(v)] %>% unique
    d <- d %>% select(-Genes) %>% t
    d <- data.frame(d)
    colnames(d) <- v
    d
}

getPCA <- function(df, gs = NULL, getTeams = F, write = F) {
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
    if (write) {
        write_csv(data.frame(Axis = paste0("PC", 1:numPC),
            Variance = varDat[1:numPC]), "PCVariance.csv")
    }
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

geneReplace <- function(df, geneSet, hkGenes = NULL, nIters = 100, diagnostics = F) {
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
    if (length(geneSet)>=50)
        repSet <- c(seq(1, length(geneSet), length.out = 20), length(geneSet))
    # plan(multisession, workers = 10)
    d <- lapply(repSet, function(x) {
        plan(multisession, workers = 10)
        d1 <- future_sapply(1:nIters, function(it) {
            repCor <- hkCor <- leftCor <- hkPC1 <- repPC1 <- leftPC1 <- hkCoeffs <- leftCoeffs <- 0
            replaceSet <- sample(geneSet, x, replace=F)
            leftSet <- geneSet[!(geneSet %in% replaceSet)]
            hkSet <- sample(hkGenes, x, replace = T)
            dNew <- df %>% select(-all_of(replaceSet)) %>%
                bind_cols(hkDat %>% select(all_of(hkSet)))
            pcaDat <- getPCA(dNew, getTeams = F)
            if (diagnostics) {
                pcaDat <- getPCA(dNew, getTeams = T)
                if (length(replaceSet)) {
                    repCor <- df %>% select(all_of(replaceSet)) %>% cor %>% abs %>% sum
                    repCor <- repCor/(x*x)
                    repPC1 <- getPCA(df %>% select(all_of(replaceSet)), getTeams = F)[1]
                }
                
                if (length(hkSet)) {
                    hkCor <- hkDat %>% select(all_of(hkSet)) %>% cor %>% abs %>% sum
                    hkPC1 <- getPCA(hkDat %>% select(all_of(hkSet)), getTeams = F)[1]
                    hkCor <- hkCor/(x*x)
                    hkCoeffs <- pcaDat$Coefficients %>% filter(Gene %in% hkSet) %>% 
                    pull(Coeff) %>% sum
                }
                
                if (x != ncol(df)) {
                    leftCor <- df %>% select(-all_of(replaceSet)) %>% cor %>% abs %>% sum
                    leftCor <- leftCor/((ncol(df) - x)*(ncol(df) - x))
                    leftPC1 <- getPCA(df %>% select(-all_of(replaceSet)), getTeams = F)[1]
                    leftCoeffs <- pcaDat$Coefficients %>% filter(Gene %in% geneSet) %>% 
                        pull(Coeff) %>% sum
                }
                pcaDat <- pcaDat$numbers
            }
            
            return(c(pc1Var = pcaDat[1], pc2Var = pcaDat[2], numPC = pcaDat[3],
                hkCor = hkCor, repCor = repCor, repPC1 = repPC1, leftCor = leftCor,
                hkPC1 = hkPC1, leftPC1 = leftPC1, hkCoeffs = hkCoeffs,
                leftCoeffs = leftCoeffs, repGenes = paste0(replaceSet, collapse = "_"),
                leftGenes = paste0(leftSet, collapse = "_")))
        }, future.seed = T) %>% t %>% as.data.frame %>%
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

getHkGenes <- function(df, gs, dfHk, HKthreshold, nGenes = 100) {
    pc1Gs <- getPCA(df, gs, getTeams = F)[1]
    pc1Hk <- getPCA(dfHk, getTeams = F)[1]
    iter <- 0
    hkSample <- colnames(dfHk)
    sampleAlt <- hkSample
    pc1Hk <- 1
    while(pc1Hk > HKthreshold*pc1Gs && iter < 1000) {
        # browser()
        nSample <- ifelse(is.null(nGenes), length(gs)*2, nGenes)
        plan(multisession, workers = 10)
        hkSampleList <- lapply(1:100, function(it) {
            hkSample <- sample(colnames(dfHk), nSample, replace = F)
            return(hkSample)
        })
        pc1Vals <- future_lapply(1:100, function(it) {
            hkSample <- hkSampleList[[it]]
            dfH <- dfHk %>% select(all_of(hkSample))
            pc1Hk <- getPCA(dfH, getTeams = F)[1]
            return(pc1Hk)
        }) %>% unlist
        future:::ClusterRegistry("stop")
        pc1Temp <- min(c(pc1Vals))
        if (pc1Temp < pc1Hk) {
            hkSample <- hkSampleList[[which.min(pc1Vals)]]
            pc1Hk <- pc1Temp
        }
        iter <- iter + 100
        if (iter > 900) {
            if (HKthreshold < 0.9) {
                HKthreshold <- 0.9
                iter <- 0
                print("Resetting threshold")
            }
        }
        if (iter >= 999) {
            print(paste0(": Could not find housekeeping genes"))
        }
        iter <- iter + 1
    }
    return(hkSample)
}

randomGeneSetAnalysis <- function(df, hkGenes, HKthreshold = 0.5, nGs = 50, nIters = 100,
    highPC = F, pcThresh = 0.4) {
        # browser()
    hkGenes <- hkGenes[hkGenes %in% df$Genes]
    dfHk <- df %>% filter(Genes %in% hkGenes) %>% dTrans(gs = hkGenes)
    hkGenes <- colnames(dfHk)
    geneLabels <- read_csv(geneLabels)
    useful <- geneLabels %>% filter(Gene_type %in% c("protein_coding", "miRNA")) %>% 
        pull(Gene_name)
    useful <- useful[useful %in% df$Genes]
    useful <- useful[!(useful %in% hkGenes)]
    df <- df %>% filter(Genes %in% useful) %>% dTrans(gs = useful)
    dfRandom <- lapply(1:nIters, function(it) {
        d <- df %>% select(sample(colnames(df), nGs, replace = F))
        hkSample <- getHkGenes(d, colnames(d), dfHk, HKthreshold = HKthreshold, nGenes = max(100, 2*nGs))
        nodeReplacementDf <- geneReplace(cbind.data.frame(d, dfHk), colnames(d), hkGenes = hkGenes)
        dSummarised <- pc1StabilityAnalysis(nodeReplacementDf) %>% 
            mutate(geneSet = paste0("Random_", it))
    }) %>% bind_rows()
    return(dfRandom)
}

dataSetAnalysis <- function(datasetList, genesetList, 
    dataKey = NULL, figureDir = ".", hkGenes = NULL, 
    HKthreshold = 0.5, reduceGenes = F, 
    randomGenes = T, diagnostics = F, nRand = 100, nThreads = 10) {
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
                if (is.null(gs)) {
                    geneLabels <- read_csv(geneLabels)
                    useful <- geneLabels %>% filter(Gene_type %in% c("protein_coding", "miRNA")) %>% 
                        pull(Gene_name)
                    useful <- useful[useful %in% df$Genes]
                    useful <- useful[!(useful %in% hkGenes)]
                    gs <- sample(useful, 50, replace = F)
                }
                gs <- gs[gs %in% df$Genes]
                df  <- df %>% filter(Genes %in% gs)
                v <- df$Genes
                df <- df %>% select(-Genes) %>% t %>% data.frame %>% set_names(v)
                
                # if (is.null(hkGenes))
                #     dfHk <- generateHouseKeeping(df, nGenes = max(100, length(gs)))
                # else {
                #     dfHk <- dfFull %>% filter(Genes %in% hkGenesUpd)
                #     v <- dfHk$Genes
                #     dfHk <- dfHk %>% select(-Genes) %>% t %>% data.frame %>% set_names(v)
                # }
            }
            else {
                gs <- gs[gs %in% colnames(df)]
                df <- df %>% select(all_of(gs))
                # if (is.null(hkGenes))
                #     dfHk <- generateHouseKeeping(df, nGenes = max(100, length(gs)))
                # else
                #     dfHk <- dfFull %>% select(any_of(hkGenes))
            }
            
            if (reduceGenes && ncol(df) > 100) {
                pcaDat <- getPCA(df, getTeams = T)
                pcaDat <- pcaDat$Coefficients
                pcaDat <- pcaDat %>% filter (abs(Coeff) < 0.9) %>% 
                    mutate(AbsCoeff = abs(Coeff)) %>%
                    arrange(desc(AbsCoeff)) %>% head(100)
                GenestoUse <- pcaDat$Gene
                df <- df %>% select(all_of(GenestoUse))
                gs <- GenestoUse
            }
            # browser()
            d <- getPCA(df, write = T)
            d <- read_csv("PCVariance.csv")
            write_csv(d, paste0(figDir, "/", data, "_", geneKey, "_PCVariance.csv"))

            # hkSample <- getHkGenes(df, gs, dfHk, HKthreshold = HKthreshold, nGenes = max(100, 2*length(gs)))
            # dfH <- dfHk %>% select(all_of(hkSample))
            # plotData(df, dfH,
            #     key = paste0(data, "_", geneKey), fld = figDir,
            #     reduceGenes = reduceGenes, diagnostics = diagnostics) %>% 
            #     mutate(geneSet = geneKey, dataSet = data)
            return(NULL)
        }) #%>% bind_rows()
        # if (randomGenes) {
        #     dfRandom <- randomGeneSetAnalysis(datasetList[[data]], hkGenes, nGs = 50, nIters = nRand) %>% mutate(dataSet = data)
        #     dSummaries <- bind_rows(dSummaries, dfRandom)
        # }
        # dSummaries
    }) #%>% bind_rows()
    # write_csv(dAll, paste0(figureDir, "/Summary.csv"))
    # summaryDir <- paste0(figureDir, "/SummaryHeatmaps")
    # plotSummariesDataSets(dAll, paste0(figureDir, "/SummaryHeatmaps"))
    # sapply(dataSets, function(data) {
    #     plotSummaries(dAll %>% filter(dataSet == data), 
    #         paste0(summaryDir, "/", data))
    # })
}

pcPlotsFull <- function(df, key, figureDir = ".", geneListDf = NULL) {
    if(!dir.exists(figureDir))
        dir.create(figureDir, recursive = T)
    pcDat <- prcomp(df)

    pcvar <- pcDat$sdev^2
    pcV <- data.frame(PCaxes = 1:length(pcvar), Variance = log10(pcvar), 
        Percentage = pcvar*100/sum(pcvar), 
        CumulativeVar = cumsum(pcvar/sum(pcvar))) %>%
        filter(PCaxes < 101)
    pcVar <- pcV %>% gather(key = "Metric", value = "Value", -PCaxes)
    Key  <- c("Log10(Variance)", "Percentage Variance", "Cumulative Variance")
    names(Key) <- c("Variance", "Percentage", "CumulativeVar")
    pcVar$Metric <- Key[pcVar$Metric] %>% factor(levels = Key)
    ggplot(pcVar, aes(x = PCaxes, y = Value)) + 
        geom_point() +
        geom_line() + 
        facet_wrap(~Metric, scales = "free_y") + 
        theme_Publication()
    ggsave(paste0(figureDir, "/", key, "_all.png"), width = 12, height = 5)

    ### PC1 coefficients
    pc1 <- pcDat$rotation[, 1]
    d <- data.frame(Genes = names(pc1), Coefficient = pc1) %>% arrange(Coefficient)
    ggplot(d, aes(x = Coefficient)) +
        geom_histogram() +
        theme_Publication() +
        scale_y_log10() +
        labs(x = "PC1 Coefficient", y = "Frequency")
    ggsave(paste0(figureDir, "/", key, "_PC1Coefficients.png"), width = 5.5, height = 5)

    if (!is.null(geneListDf)) {
        numLists <- length(unique(geneListDf$Gs))
        rows <- numLists %/% 3
        geneListDf <- merge(geneListDf, d, by = "Genes", all.x  = T)
        geneListDf <- geneListDf[complete.cases(geneListDf),]
        ggplot(geneListDf, aes(x = Coefficient)) +
            geom_histogram() +
            facet_wrap(~Gs, scales = "free_y", ncol = 4) +
            theme_Publication() +
            scale_y_log10() +
            labs(x = "PC1 Coefficient", y = "Frequency")
        ggsave(paste0(figureDir, "/", key, "_PC1Coefficients_GeneSets.png"), 
            width = 15, height = 5*rows)
        gs <- unique(geneListDf$Gs)
        gsVar <- lapply(gs, function(x) {
            g <- geneListDf %>% filter(Gs == x) %>% pull(Genes)
            d <- df %>% select(all_of(g))
            pca <- prcomp(d)
            pcvar <- pca$sdev^2
            pVar <- data.frame(PCaxes = 1:length(pcvar), Variance = log10(pcvar), 
                Percentage = pcvar*100/sum(pcvar), 
                CumulativeVar = cumsum(pcvar/sum(pcvar)), 
                Geneset = x %>% str_replace(" ", "\n")) %>%
                filter(PCaxes < 101) %>%
                mutate(PCaxes = paste0("PC", PCaxes))
            # axis <- pca$rotation[, 1]
            # get the variance explained by the first axis
            pVar
        }) %>% bind_rows()
        gsVar <- gsVar %>% 
            rbind.data.frame(pcV %>% mutate(Geneset = "Whole\nGenome",
            PCaxes = paste0("PC", PCaxes))) %>%
            filter(PCaxes %in% c("PC1", "PC2"))
        ggplot(gsVar, aes(x = Geneset, y = Variance, fill = PCaxes)) +
            geom_bar(stat = "identity", position = "dodge") +
            theme_Publication() +
            theme(legend.position = "top") +
            labs(x = "", y = "Log10(Variance)") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        ggsave(paste0(figureDir, "/", key, "_GeneSetVariance.png"),
            width = 6.3, height = 6.3)
        
        ggplot(gsVar, aes(x = Geneset, y = Percentage, fill = PCaxes)) +
            geom_bar(stat = "identity", position = "dodge") +
            theme_Publication() +
            theme(legend.position = "top") +
            labs(x = "", y = "Percentage Variance") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        ggsave(paste0(figureDir, "/", key, "_GeneSetPercentage.png"), 
            width = 6.3, height = 6.3)
    }
}
