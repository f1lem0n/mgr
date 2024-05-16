library(RNAsmc)
library(RRNA)
library(pheatmap)
library(reshape2)
library(tidyverse)
library(dplyr)

genes <- list.dirs("output", full.names = FALSE, recursive = FALSE)
genes

invivo_vs_insilico <- c()
invivo_vs_consensus <- c()

for (gene in genes) {
  insilico <- loadCt(paste0("output/", gene, "/CT/insilico.ct"))
  invivo <- loadCt(paste0("output/", gene, "/CT/invivo.ct"))
  consensus <- loadCt(paste0("output/", gene, "/CT/consensus.ct"))
  substructures <- list(
    insilico = getSubStr(insilico),
    invivo = getSubStr(invivo),
    consensus = getSubStr(consensus)
  )
  score_matrix <- getCompare(substructures)
  invivo_vs_insilico <- append(invivo_vs_insilico, score_matrix[2, 1])
  invivo_vs_consensus <- append(invivo_vs_consensus, score_matrix[2, 3])
  labs <- c("in silico", "in vivo", "consensus")
  colnames(score_matrix) <- labs
  rownames(score_matrix) <- labs
  score_df <- melt(score_matrix)
  colnames(score_df) <- c("x", "y", "value")
  plot <- score_df %>%
    ggplot(aes(x = x, y = y, fill = value)) +
    geom_tile(
      color = "black",
      lwd = 0.5, linetype = 1
    ) +
    scale_fill_gradient(
      low = "white", high = "red",
      limits = c(0, 10)
    ) +
    geom_text(
      aes(label = round(value, digits = 2)),
      color = "black", size = 4
    ) +
    guides(
      fill = guide_colourbar(
        barwidth = 1,
        barheight = 20,
        frame.colour = "black",
        ticks.colour = "black",
        title = "Similarity score"
      )
    ) +
    coord_fixed() +
    ggtitle(paste0("RNA structure similarity of ", gene)) +
    theme_minimal() +
    theme(
      title = element_text(size = 16),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_text(size = 14),
      plot.background = element_rect(fill = "white")
    )
  ggsave(
    paste0("output/", gene, "/similarity_heatmap.png"),
    plot = plot,
    width = 10, height = 10, dpi = 300
  )
}

df <- data.frame(
  gene = genes,
  VS = invivo_vs_insilico,
  VC = invivo_vs_consensus
)

write.table(
  df, "output/similarity_scores.tsv",
  sep = "\t",
  row.names = FALSE, quote = FALSE
)

df <- read.table("output/similarity_scores.tsv", header = TRUE, sep = "\t")
df <- melt(df)
colnames(df) <- c("gene", "type", "value")
df <- df %>%
  mutate(type = ifelse(type == "VS", "in silico", "consensus"))
plot <- df %>%
  ggplot(aes(x = type, y = value)) +
  geom_boxplot(staplewidth = 0.1, lwd = 1) +
  geom_dotplot(
    aes(fill = gene),
    binaxis = "y", stackdir = "center", dotsize = 0.5,
    alpha = 0.8, col = NA
  ) +
  ylab("Similarity score") +
  theme_classic() +
  theme(
    title = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 14),
    plot.background = element_rect(fill = "white")
  )
ggsave(
  "output/similarity_scores_boxplot.png",
  plot = plot,
  width = 10, height = 10, dpi = 300
)
