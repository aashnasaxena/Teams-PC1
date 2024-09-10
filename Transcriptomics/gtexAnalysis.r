source("setupScript.r")
source("functions.r")

logMutate <- function(df) {
    df %>% mutate(across(where(is.numeric), function(x){log2(x+1)}))
}

sclc <- readRDS(SCLCgenesets)
sclcNet <- sclc$network
sclcDat <- sclc$fiftyGenes

emt <- readRDS(EMTgenesets)
emtAiello <- emt$Aiello_et_al.2018

CCLE <- readRDS(ccleCounts)
CCLETPM <- readRDS(ccletpm)
CCLEOLD <- readRDS(ccleOld)
CCLERPKM <- readRDS(ccleRPKM)

GTEX <- readRDS(gtexFile)
GTExTPM <- readRDS(GTExtpm)

# Tissues in functional genesets
#  [1] "Adipose Tissue"  "Adrenal Gland"   "Bladder"         "Blood"          
#  [5] "Blood Vessel"    "Brain"           "Breast"          "Cervix Uteri"   
#  [9] "Colon"           "Esophagus"       "Fallopian Tube"  "Heart"          
# [13] "Kidney"          "Liver"           "Lung"            "Muscle"         
# [17] "Nerve"           "Ovary"           "Pancreas"        "Pituitary"      
# [21] "Prostate"        "Salivary Gland"  "Skin"            "Small Intestine"
# [25] "Spleen"          "Stomach"         "Testis"          "Thyroid"        
# [29] "Uterus"          "Vagina" 

funcGenes <- function(tis, thresh) {
    readRDS(functionalGenesets) %>% filter(str_detect(toupper(Tissue), tis)) %>%
        filter(Label %in% c("P", "N")) %>%
        filter(Probability > thresh)
} 

getUniqueGs <- function(g1, g2){
    v1 <- setdiff(g1, g2)
    v2 <- setdiff(g2, g1)
    c(v1, v2)
}

# "TRAVAGLINI_LUNG_ALVEOLAR_EPITHELIAL_TYPE_1_CELL"    - Epithelial           
# "TRAVAGLINI_LUNG_ALVEOLAR_EPITHELIAL_TYPE_2_CELL"   - Epithelial            
# "TRAVAGLINI_LUNG_ALVEOLAR_FIBROBLAST_CELL" - Mes (maybe)
# "TRAVAGLINI_LUNG_AIRWAY_SMOOTH_MUSCLE_CELL" -Mesenchymal



#### lung vs sclc --------------------------------------------------------------
msigDBgs <- readRDS(msigDBGTEx) %>%
    filter(gs_name %in% c("TRAVAGLINI_LUNG_ALVEOLAR_EPITHELIAL_TYPE_1_CELL",
        "TRAVAGLINI_LUNG_ALVEOLAR_EPITHELIAL_TYPE_2_CELL",
        "TRAVAGLINI_LUNG_ALVEOLAR_FIBROBLAST_CELL",
        "TRAVAGLINI_LUNG_AIRWAY_SMOOTH_MUSCLE_CELL")) %>%
    mutate(gs_name = str_remove(gs_name, "TRAVAGLINI_LUNG_")) %>%
    mutate(gs_name = str_remove(gs_name, "_CELL")) %>%
    mutate(gs_name = str_replace(gs_name, "ALVEOLAR_EPITHELIAL_TYPE_1", "AT1")) %>%
    mutate(gs_name = str_replace(gs_name, "ALVEOLAR_EPITHELIAL_TYPE_2", "AT2")) %>%
    mutate(gs_name = str_replace(gs_name, "AIRWAY_SMOOTH_MUSCLE", "AirSM")) %>%
    mutate(gs_name = str_replace(gs_name, "ALVEOLAR_FIBROBLAST", "AF")) %>%
    select(gs_name, gene_symbol)
at1 <- msigDBgs %>% filter(gs_name == "AT1") %>% pull(gene_symbol)
at2 <- msigDBgs %>% filter(gs_name == "AT2") %>% pull(gene_symbol)
af <- msigDBgs %>% filter(gs_name == "AF") %>% pull(gene_symbol)
airSM <- msigDBgs %>% filter(gs_name == "AirSM") %>% pull(gene_symbol)

at1at2 <- getUniqueGs(at1, at2)
afat <- getUniqueGs(af, c(at1, at2))
airSMat <- getUniqueGs(airSM, c(at1, at2))
eVm <- getUniqueGs(c(at1, at2), c(af, airSM))

functionalGenes <- readRDS(functionalGenesets) %>% filter(str_detect(toupper(Tissue), "LUNG")) %>%
    filter(Label %in% c("P", "N")) %>%
    filter(Probability > 0.8)
gsFunctional <- functionalGenes %>% pull(Gene)

geneSets <- list(sclc = sclcNet, sclcfifty = sclcDat, emt = emtAiello, 
    at1 = at1, at2 = at2, af = af, airSM = airSM, 
    at1at2 = at1at2, afat = afat, airSMat = airSMat, eVm = eVm, 
    functional = gsFunctional, Random = NULL)
gtexLung <- GTExTPM$LUNG %>% bind_rows
    # mutate(Genes = Description) %>%
    # select(-Name, -Description)

# ccleDat <- CCLE %>% select(-id)
ccleDat <- CCLETPM %>% select(-id)
ccleLung <- ccleDat %>% 
    select(Genes, contains("LUNG"))
ccleSCLC <- ccleDat %>%
    select(Genes, contains("SCLC"))
cclenonSCLC <- ccleDat %>%
    select(Genes, !contains("SCLC"))

# dataSets <- list(#CCLELung = ccleLung, CCLESCLC = ccleSCLC, 
# GTEx = gtexLung)

# dataSetAnalysis(dataSets, geneSets, dataKey = NULL, figureDir = "Figures/LungTPM", 
#     hkGenes = readRDS(houseKeepingGenes))

# dataSetAnalysis(list(CCLELung = ccleLung, GTEx = gtexLung, 
#     GTExLog = gtexLung %>% logMutate,
#     ccleLungLog = ccleLung %>% logMutate,
#     cclenonSCLC = ), geneSets, dataKey = NULL, 
#     figureDir = "Figures/LungTPM_noreduction", hkGenes = readRDS(houseKeepingGenes))

# dataSetAnalysis(list(CCLELung = ccleLung, GTEx = gtexLung, 
#     GTExLog = gtexLung %>% 
#     mutate(across(where(is.numeric), function(x){log2(x+1)})),
#     ccleLungLog = ccleLung %>% 
#     mutate(across(where(is.numeric), function(x){log2(x+1)}))), geneSets, dataKey = NULL, 
#     reduceGenes = T, 
#     figureDir = "Figures/LungTPM_reduction", hkGenes = readRDS(houseKeepingGenes))

dataSetAnalysis(list(#CCLELung = ccleLung, GTEx = gtexLung, 
    GTExLog = gtexLung %>% logMutate,
    ccleLungLog = ccleLung %>% logMutate
    cclenonSCLCLog = cclenonSCLC %>% logMutate,
    ccleSCLClog = ccleSCLC %>% logMutate
    ), geneSets, dataKey = NULL, 
    reduceGenes = T, randomGenes = T, #HKthreshold = 0.1,
    figureDir = "Figures/LungTPM_reduction_random", hkGenes = readRDS(houseKeepingGenes))
dataSetAnalysis(list(#CCLELung = ccleLung, GTEx = gtexLung, 
    GTExLog = gtexLung %>% logMutate,
    ccleLungLog = ccleLung %>% logMutate,
    cclenonSCLCLog = cclenonSCLC %>% logMutate,
    ccleSCLClog = ccleSCLC %>% logMutate
    ), geneSets, dataKey = NULL, 
    reduceGenes = F, randomGenes = T, #HKthreshold = 0.1,
    figureDir = "Figures/LungTPM_noreduction_random", hkGenes = readRDS(houseKeepingGenes))

correlationNPC1Plots(list(CCLELungLog = ccleLung %>% logMutate, 
        GTEx = gtexLung %>% logMutate,
        ccleSCLCLog = ccleSCLC %>% logMutate,
        cclenonSCLCLog = cclenonSCLC %>% logMutate), geneSets[-length(geneSets)], reduceGenes = F,
    figureDir = "Figures/LungTPM_noreduction_rand")
correlationNPC1Plots(list(CCLELungLog = ccleLung %>% logMutate, 
        GTEx = gtexLung %>% logMutate,
        ccleSCLCLog = ccleSCLC %>% logMutate,
        cclenonSCLCLog = cclenonSCLC %>% logMutate), geneSets[-length(geneSets)], reduceGenes = T,
    figureDir = "Figures/LungTPM_reduction_rand")

onlyPlots(figureDir = "Figures/LungTPM_noreduction_random")
onlyPlots(figureDir = "Figures/LungTPM_reduction_random")
#### breast --------------------------------------------------------------    


# functionalGenesBR <- funcGenes("BREAST", 0.8)

# geneSets <- list(sclc = sclcNet, sclcfifty = sclcDat, emt = emtAiello, 
#     functional = functionalGenes %>% pull(Gene))

# CCLEbreast <- CCLETPM %>% select(-id) %>%
#     select(Genes, contains("BREAST"))


# gtexBreast <- GTExTPM$BREAST %>% bind_rows
#     # mutate(Genes = Description) %>%
#     # select(-Name, -Description)

# dataSetAnalysis(list(#CCLELung = ccleLung, GTEx = gtexLung, 
#     GTExLog = gtexBreast %>% 
#     mutate(across(where(is.numeric), function(x){log2(x+1)})),
#     ccleLungLog = CCLEbreast %>% 
#     mutate(across(where(is.numeric), function(x){log2(x+1)}))), geneSets, dataKey = NULL, 
#     reduceGenes = T, randomGenes = T, #HKthreshold = 0.1,
#     figureDir = "Figures/BreastTPM_reduction_random", hkGenes = readRDS(houseKeepingGenes))

# #### kidney --------------------------------------------------------------



# functionalGenesBL <- funcGenes("BLOOD", 0.8)

# geneSets <- list(sclc = sclcNet, sclcfifty = sclcDat, emt = emtAiello, 
#     functional = functionalGenes %>% pull(Gene))

# CCLEblood <- CCLETPM %>% select(-id) %>%
#     select(Genes, contains("BLOOD"))