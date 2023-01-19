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
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        labs(y = "", x = "") +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradient2(low = "red", high = "blue", limits = c(-5,5))
    ggsave(paste0("EMT_RACIPE_", x, ".png"), width = 6, height = 5)
})
