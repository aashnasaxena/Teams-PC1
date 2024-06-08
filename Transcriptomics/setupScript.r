library(funcsKishore)
library(msigdbr)

pc1Data <- "/Users/kishorehari/Desktop/PostDoc/TeamsPC1Nets/DataAnalysis"

transcriptomicsFolder <- "/Users/kishorehari/Desktop/PostDoc/Transcriptomics"
gtexFile <- paste0(transcriptomicsFolder, "/GTEx/GTEx_reads.rds")
gtexAnnotations <- paste0(transcriptomicsFolder, "/GTEx/GTEx_SampleAnnotations.rds")
categorisedGTEx <- paste0(transcriptomicsFolder, "/GTEx/GTExCategorized.rds")
GTExtpm <- paste0(transcriptomicsFolder, "/GTEx/GTEx_tpm.RDS")

ccletpm <- paste0(transcriptomicsFolder, "/CCLE/allCCLE.rds")
ccleCounts <- paste0(transcriptomicsFolder, "/CCLE/allCCLE_counts.rds")
ccleOld <- paste0(transcriptomicsFolder, "/CCLE/EMT_CCLE_all_nonSCLC.RDS")
ccleRPKM <- paste0(transcriptomicsFolder, "/CCLE/CCLE_RPKM.RDS")

msigDBGTEx <- paste0(transcriptomicsFolder, "/GTEx/msigDB.rds")

EMTgenesets <- paste0(transcriptomicsFolder,"/EMT/EMTGenes.RDS")
SCLCgenesets <- paste0(transcriptomicsFolder, "/SCLC/SCLCGenes.RDS")
houseKeepingGenes <- paste0(transcriptomicsFolder, "/GeneLists/HouseKeepingGenes.RDS")

functionalGenesets <- paste0(transcriptomicsFolder, "/GeneLists/FunctionalGenes.RDS")
