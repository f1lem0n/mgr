library(reshape2)

type <- "mRNA"

df <- read.table(
  paste0("output/structures/", type, "/similarity_scores_vienna.tsv"),
  header = TRUE, sep = "\t"
)

df <- melt(df)
df.anova <- aov(value ~ variable, data = df)
# summary(df.anova)

TukeyHSD(df.anova)

library(tidyverse)
library(ggplot2)

df <- read.csv("output/structures/combined.csv")
colnames(df) <- c("gen", "struktura", "SS", "typ", "IA")
unique(df$struktura)
strukt <- "VG"
df <- df %>%
  filter(struktura == strukt)

m.df <- df[df$typ == "mRNA", ]
t.df <- df[df$typ == "tRNA", ]
r.df <- df[df$typ == "rRNA", ]


m.lm <- lm(SS ~ IA, data = m.df)
summary(m.lm)
cor.test(m.df$IA, m.df$SS)

t.lm <- lm(SS ~ IA, data = t.df)
summary(t.lm)
cor.test(t.df$IA, t.df$SS)

r.lm <- lm(SS ~ IA, data = r.df)
summary(r.lm)
cor.test(r.df$IA, r.df$SS)

plot <- ggplot(df, aes(x = IA, y = SS, colour = typ)) +
  geom_point(size = 5) +
  geom_abline(slope = coef(m.lm)[2], intercept = coef(m.lm)[1], colour = "red", size = 1.5) +
  geom_abline(slope = coef(r.lm)[2], intercept = coef(r.lm)[1], colour = "green", size = 1.5) +
  geom_abline(slope = coef(t.lm)[2], intercept = coef(t.lm)[1], colour = "blue", size = 1.5) +
  guides(colour = guide_legend(title = "")) +
  guides(colour = "none") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "italic"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 24),
    axis.text = element_text(size = 24),
  )
plot
ggsave(
  paste0(
    "output/structures/correlation/", strukt, ".png"
  ),
  plot = plot,
  width = 10, height = 10, dpi = 300
)
