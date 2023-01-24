library(funcsKishore)

EMTWTData <- "D:\\Github\\Projects\\Finished\\Mubasher\\RACIPE\\RawData\\EMT_RACIPE"

df  <- read_csv(paste0(EMTWTData, "/wild_EMT_solution_1.csv"))
colnames(df) <- colnames(df) %>% str_remove_all("-")

pcDat <- prcomp(df %>% select(-(1:3)))
pc1 <- pcDat$rotation[,1]
pcRotation <- pcDat$rotation %>% data.frame %>% 
    mutate(Nodes = rownames(.)) %>% arrange(PC1)

nodes <- pcRotation$Nodes

pcDatAll <- pcDat$x

zscore <- function(x) {
    (x-mean(x))/sd(x)
}

clusteredData <- hclust(dist(pcDatAll))
hc <- cutree(clusteredData, 3)
df2 <- df %>% select(all_of(nodes)) %>% mutate(across(.fns = zscore))
sapply(unique(hc), function(x) {
    df1 <- df2[which(hc == x),] %>%
        mutate(ID = 1:nrow(.)) %>%
        gather(key = "Nodes", value = "Expression", -ID) %>%
        mutate(Nodes = factor(Nodes, levels = nodes))
    ggplot(df1, aes(x = Nodes, y = ID, fill = Expression)) +
        geom_tile() +
        theme_Publication() +
        theme(legend.key.height = unit(0.8, "cm"),
        legend.direction = "vertical",
        legend.position = "right",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, family = "Palatino")) +
        labs(y = "", x = "") +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradient2(low = "red", high = "blue", limits = c(-5,5))
    ggsave(paste0("EMT_RACIPE_", x, ".png"), width = 6, height = 5)
})

pc1 <- data.frame(Nodes = names(pc1), Loading = pc1, Class = ifelse(Loading > 0, "S", "N"))

ggplot(pc1, aes(x = reorder(Nodes, -Loading), y = Loading, fill = Class)) +
    geom_bar(stat = "identity") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = NULL) +
    labs(x = "", y = "PC1 Loading")

ggsave("PCALoading.png", width = 9, height = 5.5)


# multicolor matrix
setwd("D:\\Github\\Projects\\Finished\\Mubasher\\FinalResults\\TopoFiles")

topoFile <- "EMT_RACIPE.topo"

hex2dec <- function(x) {
    x <- str_split(x, "") %>% rev
    hexcode <- 0:15
    hexKey <- c(0:9, letters[1:6])
    names(hexcode) <- hexKey
    s <- sapply(1:length(x), function(i) {
        hexCode[x[i]]*(16^(i-1))
    }) %>% sum
    return(s)
}

dec2hex <- function(x) {
    hexcode <- 0:15
    hexKey <- c(0:9, letters[1:6])

}

colorGen <- function(high = "#ff0000", low = "#ffffff", 
                        limits = c(0,1), breakSize = 0.05) {
    
}


plotCustomColors <- function(topoFiles = NULL, negCol = "#696969", 
    T1 = "#ff0000", T2 = "#0000ff", zeroColr = "#ffffff",
    alpha = T) {
    # DirectoryNav("influenceTest")
    # file.copy(paste0("../", topoFile), ".")
    if (is.null(topoFiles)) 
        topoFiles <- list.files(".", ".topo$")
    if (length(topoFiles) == 0 || !all(topoFiles %>% sapply(file.exists)))
    {
        message("Invalid input topofiles. Do topofiles exist in your folder?")
        return()

    }
    getGsVec(method = "Cluster", topoFiles = topoFiles)
    sapply(topoFiles, function(topoFile) {
        inflMat <- read_csv(paste0("Influence/", 
            str_replace(topoFile, ".topo", "_reducedInfl.csv")))
        teams <- readLines(str_replace(topoFile, ".topo", ".teams")) %>%
            str_split(",")
        teamsOrder <- teams %>% unlist
        df <- inflMat %>% 
            gather(key = "Target", value = "Influence", -Source) %>%
            mutate(colorVal = negCol, alphaVal = 1)
        if (alpha) {
            df <- df %>%
                mutate(alphaVal = abs(Influence),
                    colorVal = ifelse(Source %in% teams[[1]] & 
                        Target %in% teams[[1]], T1, colorVal)) %>%
                mutate(colorVal = ifelse(Source %in% teams[[2]] & 
                    Target %in% teams[[2]], T2, colorVal))
        }
        df <- df  %>%
            mutate(Target = factor(Target, levels = teamsOrder),
                Source = factor(Source, levels = teamsOrder))
        ggplot(df, aes(x = Target, y = Source, fill = colorVal, 
            alpha = alphaVal)) +
            geom_tile() + 
            theme_Publication() + 
            scale_fill_identity() + 
            scale_alpha_identity() +
            theme(axis.text.x = element_text(angle = 90, 
                hjust = 1, vjust = 0.5))

        DirectoryNav("MatrixPlots")
        ggsave(str_replace(topoFile, 
            ".topo", "_customInfluence.png"), width = 6.5, height = 6)
        setwd("..")
    })
    
        
}
