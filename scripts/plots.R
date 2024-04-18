peaks124 <- read.table(
  "output/peaks/DMSext124_peaks.narrowPeak",
  header = FALSE
)
peaks200 <- read.table(
  "output/peaks/DMS_peaks.narrowPeak",
  header = FALSE
)
peaks300 <- read.table(
  "output/peaks/DMSext300_peaks.narrowPeak",
  header = FALSE
)
cols <- c(
  "chr",
  "start",
  "end",
  "name",
  "score",
  "strand",
  "signalValue",
  "pValue",
  "qValue",
  "peak"
)

colnames(peaks124) <- cols
colnames(peaks200) <- cols
colnames(peaks300) <- cols

plot(
  peaks124$start, peaks124$pValue,
  type = "h", xlab = "Genomic position", ylab = "-log10(p)"
)
plot(
  peaks200$start, peaks200$pValue,
  type = "h", xlab = "Genomic position", ylab = "-log10(p)"
)
plot(
  peaks300$start, peaks300$pValue,
  type = "h", xlab = "Genomic position", ylab = "-log10(p)"
)

plot(peaks124$start, peaks124$end)
abline(a = 0, b = 1, col = "red")
