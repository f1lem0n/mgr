library(RNAsmc)
library(RRNA)
library(pheatmap)
library(reshape2)
library(tidyverse)
library(dplyr)

types <- list.dirs("output/structures/", full.names = FALSE, recursive = FALSE)
programs <- c("vienna")

invivo_vs_insilico <- c()
invivo_vs_consensus <- c()
invivo_vs_guided <- c()

for (program in programs) {
  for (type in types) {
    genes <- list.dirs(
      paste0("output/structures/", type),
      full.names = FALSE, recursive = FALSE
    )
    for (gene in genes) {
      insilico <- loadCt(
        paste0(
          "output/structures/", type, "/",
          gene, "/", program, "/CT/insilico.ct"
        )
      )
      invivo <- loadCt(
        paste0(
          "output/structures/", type, "/",
          gene, "/", program, "/CT/invivo.ct"
        )
      )
      consensus <- loadCt(
        paste0(
          "output/structures/", type, "/",
          gene, "/", program, "/CT/consensus.ct"
        )
      )
      if (type == "mRNA") {
        consensus_guided <- loadCt(
          paste0(
            "output/structures/", type, "/",
            gene, "/", program, "/CT/consensus_guided.ct"
          )
        )
        substructures <- list(
          insilico = getSubStr(insilico),
          invivo = getSubStr(invivo),
          consensus = getSubStr(consensus),
          consensus_guided = getSubStr(consensus_guided)
        )
        score_matrix <- getCompare(substructures)
        invivo_vs_insilico <- append(invivo_vs_insilico, score_matrix[2, 1])
        invivo_vs_consensus <- append(invivo_vs_consensus, score_matrix[2, 3])
        invivo_vs_guided <- append(invivo_vs_guided, score_matrix[2, 4])
        labs <- c("in silico", "in vivo", "consensus", "guided")
        colnames(score_matrix) <- labs
        rownames(score_matrix) <- labs
        score_df <- melt(score_matrix)
        colnames(score_df) <- c("x", "y", "value")
      } else {
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
      }
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
        paste0(
          "output/structures/", type, "/",
          gene, "/", program, "/similarity_heatmap.png"
        ),
        plot = plot,
        width = 10, height = 10, dpi = 300
      )
    }
    if (type == "mRNA") {
      df <- data.frame(
        gene = genes,
        VS = invivo_vs_insilico,
        VC = invivo_vs_consensus,
        VG = invivo_vs_guided
      )
    } else {
      df <- data.frame(
        gene = genes,
        VS = invivo_vs_insilico,
        VC = invivo_vs_consensus
      )
    }

    write.table(
      df,
      paste0(
        "output/structures/", type,
        "/similarity_scores_", program, ".tsv"
      ),
      sep = "\t",
      row.names = FALSE, quote = FALSE
    )

    df <- read.table(
      paste0(
        "output/structures/", type,
        "/similarity_scores_", program, ".tsv"
      ),
      header = TRUE, sep = "\t"
    )
    df <- melt(df)
    colnames(df) <- c("gene", "type", "value")
    replacements <- c("VS" = "in silico", "VC" = "consensus", "VG" = "guided")
    df$type <- str_replace_all(df$type, replacements)
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
      paste0(
        "output/structures/", type,
        "/similarity_scores_boxplot_", program, ".png"
      ),
      plot = plot,
      width = 10, height = 10, dpi = 300
    )
  }
}
