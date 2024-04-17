peaks  <- read.table("output/peaks/macs2.out_summits.bed", header=FALSE)
head(peaks)

plot(peaks$V5, type="l")
