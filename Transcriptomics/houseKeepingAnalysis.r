source("setupScript.r")
source("functions.r")

dTrans <- function(df, gs) {
    d <- df %>% filter(Genes %in% gs)
    v <- d$Genes
    d <- d %>% select(-Genes) %>% t
    d <- data.frame(d)
    colnames(d) <- v
    d
}

sclc <- readRDS(SCLCgenesets)
sclcNet <- sclc$network
sclcDat <- sclc$fiftyGenes

emt <- readRDS(EMTgenesets)
emtAiello <- emt$Aiello_et_al.2018

CCLE <- readRDS(ccleCounts)
CCLETPM <- readRDS(ccletpm)
ccleOldnonSCLC <- readRDS(ccleOld)

GTEX <- readRDS(gtexFile)
GTExTPM <- readRDS(GTExtpm)

hkGenes <- readRDS(houseKeepingGenes)

# cclenonSCLC <- ccleOld %>%
#     select(-id) %>%
#     select(Genes, !contains("SCLC"))

houseKeepingAnalysis <- function(df, gs, hkGenes, key, fld, nReplace = round(length(gs)/2)) {
    if (!dir.exists(fld))
        dir.create(fld, recursive = T)
    ccleOldnonSCLC <- df
    sclcDat <- unique(gs)
    hkGenes <- unique(hkGenes)
    n <- length(gs)
    pc1HK <- sapply(1:100, function(i) {
        d <- ccleOldnonSCLC
        hk <- sample(hkGenes, n)
        d <- dTrans(d, hk)
        getPCA(d)
    }) %>% t %>% data.frame %>%
        set_names(c("PC1", "PC2", "numPC"))
    if(n>(length(hkGenes)/2)) {
        nNew <- min(length(hkGenes)/2, 100)
        pc1Gs <- dTrans(ccleOldnonSCLC, sclcDat) %>% getPCA(getTeams = T)
        coeffs <- pc1Gs$Coefficients %>% mutate(Coeff = abs(Coeff)) %>% arrange(desc(Coeff)) %>% pull(Gene)
        sclcDat <- coeffs[1:nNew]
        n <- nNew
    }
    pc1Gs <- dTrans(ccleOldnonSCLC, sclcDat) %>% getPCA

    ggplot(pc1HK, aes(x = PC1)) +
        geom_histogram(bins = 30) +
        geom_vline(xintercept = pc1Gs[1], color = "red") +
        theme_Publication() + 
        labs(title = "Housekeeping genes PC1 variance",
            x = "PC1 variance",
            y = "Frequency")
    ggsave(paste0(fld, "/PC1variance.png"), width = 5.5, height = 5)

    ggplot(pc1HK, aes(x = numPC)) +
        geom_histogram(bins = 30) +
        geom_vline(xintercept = pc1Gs[3], color = "red") +
        theme_Publication() + 
        labs(title = "Housekeeping genes # PC axes",
            x = "# PC axes",
            y = "Frequency")
    ggsave(paste0(fld, "/numPC.png"), width = 5.5, height = 5)


    # How correlated are the housekeeping genes

    corHK <- sapply(1:100, function(i) {
        d <- ccleOldnonSCLC
        hk <- sample(hkGenes, n)
        d <- dTrans(d, hk)
        res <- cor(d) %>% abs
        sum(res)/(n*n)
    }) 
    corHK <- data.frame(Correlation = corHK)

    corGs <-  dTrans(ccleOldnonSCLC, sclcDat) %>% cor %>% abs
    corGs <- data.frame(Correlation = sum(corGs)/(n*n))

    ggplot(corHK, aes(x = Correlation)) +
        geom_histogram(bins = 30) +
        geom_vline(xintercept = corGs$Correlation, color = "red") +
        theme_Publication() + 
        labs(title = "Housekeeping genes correlation strength",
            x = "Correlation",
            y = "Frequency")
    ggsave(paste0(fld, "/correlation.png"), width = 5.5, height = 5)

    # Are any of the house keeping genes correlated with the sclc genes?

    # nReplace <- 30
    hkGenes <- hkGenes[!hkGenes %in% sclcDat]
    sclcDat <- sclcDat[sclcDat %in% ccleOldnonSCLC$Genes]
    hkGenes <- hkGenes[hkGenes %in% ccleOldnonSCLC$Genes]
    pc1Gs <- dTrans(ccleOldnonSCLC, sclcDat) %>% getPCA(getTeams = T)
    coeffs <- pc1Gs$Coefficients
    t1 <- coeffs %>% filter(Coeff > 0) %>% pull(Gene)
    t2 <- coeffs %>% filter(Coeff < 0) %>% pull(Gene)
    replacementData <- sapply(1:100, function(i) {
        d <- ccleOldnonSCLC
        gsSample <- sample(sclcDat, n-nReplace)
        gsLeft <- sclcDat[!sclcDat %in% gsSample]
        hkSample <- sample(hkGenes, nReplace)
        dGs <- dTrans(d, c(gsSample, hkSample))
        dLeft <- dTrans(d, gsLeft)
        pc <- getPCA(dGs, getTeams = T)
        coeffs <- pc$Coefficients
        T1 <- coeffs %>% filter(Coeff > 0) %>% pull(Gene)
        T2 <- coeffs %>% filter(Coeff < 0) %>% pull(Gene)
        hkCoeff <- coeffs %>% filter(Gene %in% hkSample) %>% pull(Coeff) %>% abs %>% sum
        gsCoeff <- coeffs %>% filter(Gene %in% gsSample) %>% pull(Coeff) %>% abs %>% sum
        pc1 <- pc$numbers
        corDat <- cor(dGs) %>% abs %>% as.data.frame
        browser()
        rownames(corDat) <- colnames(corDat) <- colnames(dGs)
        hkCor <- corDat[hkSample, hkSample] %>% sum
        gsCor <- corDat[gsSample, gsSample] %>% sum
        hkGs <- corDat[hkSample, gsSample] %>% sum
        gsLeftCor <- cor(dLeft) %>% abs %>% sum
        c(hkCoeff = hkCoeff/length(hkSample), gsCoeff = gsCoeff/length(gsSample), 
            hkCor = hkCor/(length(hkSample)^2), 
            gsCor = gsCor/(length(gsSample)^2),
            hkGs = hkGs/(length(gsSample)*length(hkSample)),
            gsLeftCor = gsLeftCor/(length(gsLeft)^2), 
            pc1 = pc1[1], numPC = pc1[3])
    }) %>% t %>% data.frame

    ggplot(replacementData, aes(x = hkCoeff, y = gsCoeff, color = pc1)) +
        geom_point() +
        theme_Publication() + 
        scale_color_viridis_c() +
        theme(legend.key.width = unit(0.8, "cm")) +
        labs(x = "HK genes PC1 Coefficients",
            y = "Geneset PC1 Coefficients")
    ggsave(paste0(fld, "/", nReplace, "hkCoeffVsGsCoeff.png"), width = 5.5, height = 5.5)

    ggplot(replacementData, aes(x = hkCor, y = gsCor, color = pc1)) +
        geom_point() +
        theme_Publication() + 
        scale_color_viridis_c() +
        theme(legend.key.width = unit(0.8, "cm")) +
        labs(x = "HK genes correlation",
            y = "Geneset correlation")
    ggsave(paste0(fld, "/", nReplace, "hkCorVsGsCor.png"), width = 5.5, height = 5.5)

    ggplot(replacementData, aes(x = hkGs)) +
        geom_histogram(bins = 30) +
        theme_Publication() + 
        labs(x = "HK genes vs Geneset correlation",
            y = "Frequency")
    ggsave(paste0(fld, "/", nReplace, "hkGsCor.png"), width = 5.5, height = 5.5)

    ggplot(replacementData, aes(x = gsLeftCor)) +
        geom_histogram(bins = 30) +
        theme_Publication() + 
        labs(x = "Remove genes correlation",
            y = "Frequency")
    ggsave(paste0(fld, "/", nReplace, "gsLeftCor.png"), width = 5.5, height = 5.5)
}


houseKeepingAnalysis(ccleOldnonSCLC, emtAiello, hkGenes, "CCLEOld", "Figures/HK/CCLEOld_EMT/", 18)

gtexLung <- GTExTPM$LUNG %>% bind_rows
gs <- emtAiello
nReplace <- 18
houseKeepingAnalysis(gtexLung, gs, hkGenes, "GTExLung", "Figures/HK/GTExLung_EMT/", nReplace)

houseKeepingAnalysis(gtexLung, eVm, hkGenes, "GTExLung", "Figures/HK/GTExLung_eVm/")

houseKeepingAnalysis(ccleLung, at1, hkGenes, "CCLELung", "Figures/HK/CCLELung_AT1/", 40)
# What happens at the humps?
# What happens in cases where we see an increase to 1 (that must be an error)?