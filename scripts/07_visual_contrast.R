###############################################################################
## 07_visual_contrast.R
## Chromatic and achromatic contrast (JND) between taxidermied ginger and grey
## rabbits and desert substrates, as seen by eagles, dingoes and cats.
## Contrasts were calculated with the receptor noise-limited model in
## micaToolbox (ImageJ); this script summarises and plots them.
## Supports: Fig. 5b
## Input:  data/jnd_taxidermied_rabbits.csv (Supplementary Data 7)
## Output: results/07_visual_contrast/
## Run from the repository root: Rscript scripts/07_visual_contrast.R
###############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

out_dir <- "results/07_visual_contrast"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

contrasts <- read.csv("data/jnd_taxidermied_rabbits.csv", check.names = FALSE)

# One value per rabbit: mean contrast across the 25 substrate photographs
per_rabbit <- contrasts %>%
  rename(rabbit = `Rabbit image`, metric = `Visual metric`) %>%
  pivot_longer(starts_with("Contrast"), names_to = "viewer",
               values_to = "contrast") %>%
  mutate(
    viewer = factor(sub("Contrast \\((.*)\\)", "\\1", viewer),
                    levels = c("eagle", "dingo", "cat"),
                    labels = c("Eagles", "Dingoes", "Cats")),
    metric = factor(metric, levels = c("chromatic", "achromatic"),
                    labels = c("Chromatic", "Achromatic")),
    Colour = factor(Colour, levels = c("Grey", "Ginger"))
  ) %>%
  group_by(rabbit, Colour, metric, viewer) %>%
  summarise(contrast = mean(contrast), .groups = "drop")

summary_table <- per_rabbit %>%
  group_by(viewer, metric, Colour) %>%
  summarise(
    n      = n(),
    q1     = quantile(contrast, 0.25),
    median = median(contrast),
    q3     = quantile(contrast, 0.75),
    .groups = "drop"
  )
print(summary_table, n = Inf)
write.csv(summary_table, file.path(out_dir, "visual_contrast_summary.csv"),
          row.names = FALSE)

set.seed(1)
plot <- ggplot(per_rabbit, aes(x = viewer, y = contrast, fill = Colour)) +
  geom_boxplot(outlier.shape = NA, width = 0.6,
               position = position_dodge(width = 0.7)) +
  geom_point(aes(colour = Colour), size = 0.9,
             position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.7)) +
  facet_wrap(~ metric) +
  scale_fill_manual(values = c(Grey = "#E6E6E6", Ginger = "#F9D9BF")) +
  scale_colour_manual(values = c(Grey = "#656565", Ginger = "#E67F52")) +
  labs(x = NULL, y = "Just noticeable difference (JND)",
       fill = NULL, colour = NULL) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank())

ggsave(file.path(out_dir, "visual_contrast.pdf"), plot, width = 7, height = 3.5)
