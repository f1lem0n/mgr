library(reshape2)
library(tidyverse)
library(ggplot2)

df <- read.csv("output_yeast_all/structures/combined.csv")
colnames(df) <- c("gen", "struktura", "SS", "typ", "IA", "CS")
head(df)

m_df <- df[df$typ == "mRNA", ]
t_df <- df[df$typ == "tRNA", ]
r_df <- df[df$typ == "rRNA", ]

## SS CS corr
cor.test(df$SS, df$CS)
ss_cs_lm <- lm(SS ~ CS, data = df)
summary(ss_cs_lm)
plot <- ggplot(df, aes(x = SS, y = CS)) +
  geom_point(size = 2.5, col = "blue") +
  geom_smooth(method = "lm") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "italic"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 24),
    axis.text = element_text(size = 24),
  )

plot

ggsave(
  "output/structures/SS_CS_corr.png",
  plot = plot, width = 10,
  height = 10, dpi = 300
)


# Linear models
unique(df$struktura)

strukt <- "VS"
df <- df %>%
  filter(struktura == strukt)

# korelacja IA i SS
m_lm <- lm(SS ~ IA, data = m_df)
summary(m_lm)
cor.test(m_df$IA, m_df$SS)

t_lm <- lm(SS ~ IA, data = t_df)
summary(t_lm)
cor.test(t_df$IA, t_df$SS)

r_lm <- lm(SS ~ IA, data = r_df)
summary(r_lm)
cor.test(r_df$IA, r_df$SS)

plot <- ggplot(df, aes(x = IA, y = SS, colour = typ)) +
  geom_point(size = 5) +
  geom_abline(
    slope = coef(m_lm)[2],
    intercept = coef(m_lm)[1],
    colour = "red", size = 1.5
  ) +
  geom_abline(
    slope = coef(r_lm)[2],
    intercept = coef(r_lm)[1],
    colour = "green", size = 1.5
  ) +
  geom_abline(
    slope = coef(t_lm)[2],
    intercept = coef(t_lm)[1],
    colour = "blue", size = 1.5
  ) +
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
    "output/structures/", strukt, "_IA_SS_corr.png"
  ),
  plot = plot,
  width = 10, height = 10, dpi = 300
)

# korelacja IA i CS
m_lm <- lm(SS ~ IA, data = m_df)
summary(m_lm)
cor.test(m_df$IA, m_df$SS)

t_lm <- lm(SS ~ IA, data = t_df)
summary(t_lm)
cor.test(t_df$IA, t_df$SS)

r_lm <- lm(SS ~ IA, data = r_df)
summary(r_lm)
cor.test(r_df$IA, r_df$SS)

plot <- ggplot(df, aes(x = IA, y = CS, colour = typ)) +
  geom_point(size = 5) +
  geom_abline(
    slope = coef(m_lm)[2],
    intercept = coef(m_lm)[1],
    colour = "red", size = 1.5
  ) +
  geom_abline(
    slope = coef(r_lm)[2],
    intercept = coef(r_lm)[1],
    colour = "green", size = 1.5
  ) +
  geom_abline(
    slope = coef(t_lm)[2],
    intercept = coef(t_lm)[1],
    colour = "blue", size = 1.5
  ) +
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
    "output/structures/", strukt, "_IA_CS_corr.png"
  ),
  plot = plot,
  width = 10, height = 10, dpi = 300
)

#####################################################################

# ANOVA (SS)
m_anova <- aov(SS ~ struktura, data = m_df)
summary(m_anova)
TukeyHSD(m_anova)

t_anova <- aov(SS ~ struktura, data = t_df)
summary(t_anova)
TukeyHSD(t_anova)

r_anova <- aov(SS ~ struktura, data = r_df)
summary(r_anova)
TukeyHSD(r_anova)


# ANOVA (CS)
m_anova <- aov(CS ~ struktura, data = m_df)
summary(m_anova)
TukeyHSD(m_anova)

t_anova <- aov(CS ~ struktura, data = t_df)
summary(t_anova)
TukeyHSD(t_anova)

r_anova <- aov(CS ~ struktura, data = r_df)
summary(r_anova)
TukeyHSD(r_anova)

# boxploty
types <- list.dirs("output/structures/", full.names = FALSE, recursive = FALSE)
programs <- c("vienna")

for (program in programs) {
  for (type in types) {
    df <- read.table(
      paste0(
        "output/structures/", type,
        "/similarity_scores_", program, ".tsv"
      ),
      header = TRUE, sep = "\t"
    )
    df <- melt(df)
    colnames(df) <- c("gen", "type", "value")
    replacements <- c("VS" = "in silico", "VC" = "consensus", "VG" = "guided")
    df$type <- str_replace_all(df$type, replacements)
    plot <- df %>%
      ggplot(aes(x = type, y = value)) +
      geom_boxplot(staplewidth = 0.1, lwd = 1, outliers = FALSE) +
      geom_dotplot(
        aes(fill = gen),
        binaxis = "y", stackdir = "center", dotsize = 0.5,
        alpha = 0.8, col = NA
      ) +
      ylab("SS") +
      theme_classic() +
      theme(
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 24, face = "italic"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 24, face = "italic"),
        axis.text.y = element_text(size = 24),
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
