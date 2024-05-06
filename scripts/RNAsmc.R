library(RNAsmc)
library(RRNA)
library(pheatmap)

gene <- "ASH1"

insilico <- loadCt(paste0("output/", gene, "/CT/insilico.ct"))
invivo <- loadCt(paste0("output/", gene, "/CT/invivo.ct"))
consensus <- loadCt(paste0("output/", gene, "/CT/consensus.ct"))

subStrList <- list(
  insilico = getSubStr(insilico),
  invivo = getSubStr(invivo),
  consensus = getSubStr(consensus)
)
score_matrix <- getCompare(subStrList)
score_matrix
labs <- c("in silico", "in vivo", "consensus")
colnames(score_matrix) <- labs
rownames(score_matrix) <- labs
pheatmap(score_matrix, angle_col = 0, display_numbers = TRUE, color = blues9)
