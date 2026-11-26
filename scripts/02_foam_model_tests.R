###############################################################################
## 02_foam_model_tests.R
## Fisher's exact tests and relative risk for the foam-model predation
## experiment in the Strzelecki Desert.
## Supports: Fig. 5d
## Input:  attack counts entered below
## Output: results/02_foam_model/
## Run from the repository root: Rscript scripts/02_foam_model_tests.R
###############################################################################
library(ggplot2)
library(dplyr)
library(grid)


out_dir <- "results/02_foam_model"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Attacked models by predator guild (ginger = orange fur)
bop_orange   <- 19
bop_grey     <- 27

canid_orange <- 26
canid_grey   <- 13

cat_orange   <- 4
cat_grey     <- 6

# Colour composition of attacked models: canids vs birds of prey
ft1 <- fisher.test(
  matrix(
    c(
      canid_grey,  canid_orange,
      bop_grey,    bop_orange
    ),
    nrow = 2,
    byrow = TRUE
  )
)

# Colour composition of attacked models: cats vs birds of prey
ft2 <- fisher.test(
  matrix(
    c(
      cat_grey,   cat_orange,
      bop_grey,   bop_orange
    ),
    nrow = 2,
    byrow = TRUE
  )
)

cat("\n=== FISHER'S EXACT TESTS — colour composition of attacked models ===\n\n")

cat(sprintf(
  "1. Canid vs Avian: OR = %.3f, p = %.4f\n",
  ft1$estimate,
  ft1$p.value
))

cat(sprintf(
  "2. Cat vs Avian:   OR = %.3f, p = %.4f\n",
  ft2$estimate,
  ft2$p.value
))

# Within guild: attacked vs deployed models of each colour, and relative risk
n_ginger_deployed <- 346
n_grey_deployed   <- 329

within_guild <- function(orange, grey, label) {
  m <- matrix(
    c(
      orange, n_ginger_deployed - orange,
      grey,   n_grey_deployed   - grey
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Model  = c("Ginger", "Grey"),
      Status = c("Attacked", "Not attacked")
    )
  )

  ft <- fisher.test(m)

  rr <- (orange / n_ginger_deployed) / (grey / n_grey_deployed)

  data.frame(
    Guild            = label,
    Ginger_attacked  = orange,
    Ginger_deployed  = n_ginger_deployed,
    Grey_attacked    = grey,
    Grey_deployed    = n_grey_deployed,
    RR               = rr,
    p_value          = ft$p.value,
    stringsAsFactors = FALSE
  )
}

within_guild_results <- rbind(
  within_guild(canid_orange, canid_grey, "Canid"),
  within_guild(bop_orange,   bop_grey,   "Avian")
)

cat("\n=== FISHER'S EXACT TESTS - attack rate vs models deployed ===\n\n")
print(within_guild_results, row.names = FALSE, digits = 4)

write.csv(
  within_guild_results,
  file.path(out_dir, "foam_model_within_guild_tests.csv"),
  row.names = FALSE
)

plot_data <- data.frame(
  Predator = factor(
    rep(c("Avian", "Canid"), each = 2),
    levels = c("Avian", "Canid")
  ),
  
  Model = factor(
    rep(c("Grey", "Ginger"), times = 2),
    levels = c("Grey", "Ginger")
  ),
  
  Attacks = c(
    bop_grey,
    bop_orange,
    canid_grey,
    canid_orange
  ),
  
  Percent = c(
    bop_grey   / (bop_grey + bop_orange) * 100,
    bop_orange / (bop_grey + bop_orange) * 100,
    canid_grey   / (canid_grey + canid_orange) * 100,
    canid_orange / (canid_grey + canid_orange) * 100
  )
)

p_bracket <- ifelse(
  ft1$p.value < 0.001,
  "Fisher's Exact Test, p < 0.001",
  sprintf("Fisher's Exact Test, p = %.3f", ft1$p.value)
)

y_max     <- max(plot_data$Percent)
y_bracket <- y_max * 1.12
y_top     <- y_max * 1.25

p <- ggplot(
  plot_data,
  aes(x = Predator, y = Percent, fill = Model)
) +
  
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.6,
    colour = "black",
    linewidth = 0.1,
    alpha = 0.4
  ) +
  
  geom_text(
    aes(label = sprintf("%.0f%%", Percent)),
    position = position_dodge(width = 0.7),
    vjust = -0.5,
    size = 4,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c(
      "Grey"   = "#BDBDBD",
      "Ginger" = "#EDA975"
    )
  ) +
  
  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = y_bracket,
    yend = y_bracket,
    linewidth = 0.1
  ) +
  
  annotate(
    "segment",
    x = 1,
    xend = 1,
    y = y_bracket,
    yend = y_bracket * 0.97,
    linewidth = 0.1
  ) +
  
  annotate(
    "segment",
    x = 2,
    xend = 2,
    y = y_bracket,
    yend = y_bracket * 0.97,
    linewidth = 0.1
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = y_top,
    label = p_bracket,
    size = 3.8
  ) +
  
  scale_y_continuous(
    limits = c(0, y_top * 1.10),
    expand = c(0, 0),
    labels = function(x) round(x)
  ) +
  
  labs(
    x = NULL,
    y = "Proportion of attacked models (%)"
  ) +
  
  theme_classic(base_family = "Helvetica") +
  
  theme(
    axis.text.x = element_text(
      size = 14,
      face = "bold",
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 11,
      face = "plain",
      colour = "black"
    ),
    
    axis.title.y = element_text(
      size = 13,
      face = "plain",
      margin = margin(r = 18)
    ),
    
    axis.line = element_line(
      linewidth = 0.1
    ),
    
    axis.ticks = element_line(
      linewidth = 0.1
    ),
    
    legend.position = "none"
  )

print(p)

ggsave(
  file.path(out_dir, "foam_model_colour_composition.pdf"),
  p,
  width = 4,
  height = 3.5,
  dpi = 300
)
