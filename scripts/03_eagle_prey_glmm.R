###############################################################################
## 03_eagle_prey_glmm.R
## Proportion of ginger rabbits among wedge-tailed eagle prey compared with
## field surveys in the same area.
## Supports: Fig. 6c
## Input:  data/eagle_nest_phenotypes.tsv (202 genotyped prey individuals;
##         Supplementary Data 8)
##         data/field_survey_counts.csv (Supplementary Data 2)
## Output: results/03_eagle_prey_glmm/
## Run from the repository root: Rscript scripts/03_eagle_prey_glmm.R
###############################################################################
library(ggplot2)
library(dplyr)
library(glmmTMB)


out_dir <- "results/03_eagle_prey_glmm"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.table(
  "data/eagle_nest_phenotypes.tsv",
  sep = "\t",
  header = TRUE
)

transects <- read.csv(
  "data/field_survey_counts.csv"
)

# Genotyped eagle-prey individuals, counted per nest
nest_data <- data

nest_data$nest <- substr(nest_data$Sample, 4, 7)

nest_counts <- table(nest_data$nest, nest_data$Phenotype)
nest_df <- as.data.frame.matrix(nest_counts)

nest_df$survey_type <- "nest"
nest_df$survey_type2 <- "nest"

nest_df <- nest_df[, c(
  "Ginger",
  "Grey",
  "survey_type",
  "survey_type2"
)]

transects <- transects[, c(
  "Ginger",
  "Grey",
  "survey_type",
  "survey_type2"
)]

# Combine nests and field survey units
all_data <- rbind(nest_df, transects)
all_data$sample <- as.factor(1:nrow(all_data))

# Binomial GLMM: eagle prey vs field surveys, with separate residual
# variances for nests and survey units
model4 <- glmmTMB(
  cbind(Ginger, Grey) ~ survey_type2 +
    diag(0 + survey_type2 | sample),
  family = "binomial",
  data = all_data
)

summary_m4 <- summary(model4)

print(
  summary_m4$coefficients$cond[
    "survey_type2nest",
    "Pr(>|z|)"
  ]
)

p_val <- summary_m4$coefficients$cond[
  "survey_type2nest",
  "Pr(>|z|)"
]

# Plot (Fig. 6c)
plot_summary <- all_data %>%
  group_by(survey_type) %>%
  summarise(
    Ginger = sum(Ginger),
    Grey = sum(Grey),
    Total_Rabbits = sum(Ginger + Grey),
    N_Units = n()
  ) %>%
  mutate(
    Freq = (Ginger / Total_Rabbits) * 100,
    
    x_pos = case_when(
      survey_type == "nest" ~ 1,
      survey_type == "spotlight transect" ~ 2.5,
      survey_type == "hunter harvest" ~ 3.5
    ),
    
    Label = case_when(
      survey_type == "nest" ~
        paste0(
          "Eagle\nNest\n(n=",
          Total_Rabbits,
          ", N=",
          N_Units,
          ")"
        ),
      
      survey_type == "spotlight transect" ~
        paste0(
          "Spotlight\nTransect\n(n=",
          Total_Rabbits,
          ", N=",
          N_Units,
          ")"
        ),
      
      survey_type == "hunter harvest" ~
        paste0(
          "Hunter\nHarvest\n(n=",
          Total_Rabbits,
          ", N=",
          N_Units,
          ")"
        )
    ),
    
    Pct_Label = sprintf("%.1f%%", Freq)
  )

y_max <- 40

group_bracket_y <- 31
stat_bracket_y <- 36

if (p_val < 0.001) {
  stat_label <- "'GLMM, ' ~ italic(p) < 0.001"
} else {
  stat_label <- paste0(
    "'GLMM, ' ~ italic(p) == ",
    round(p_val, 4)
  )
}

a <- ggplot(
  plot_summary,
  aes(
    x = x_pos,
    y = Freq,
    fill = survey_type
  )
) +
  
  geom_col(
    width = 0.8,
    color = "black",
    linewidth = 0.1,
    show.legend = FALSE
  ) +
  
  geom_text(
    aes(label = Pct_Label),
    vjust = -0.8,
    size = 5,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c(
      "nest" = "#f8d194",
      "spotlight transect" = "#ce7334",
      "hunter harvest" = "#835844"
    )
  ) +
  
  annotate(
    "segment",
    x = 2.5,
    xend = 3.5,
    y = group_bracket_y,
    yend = group_bracket_y,
    size = 0.1
  ) +
  
  annotate(
    "segment",
    x = 2.5,
    xend = 2.5,
    y = group_bracket_y,
    yend = group_bracket_y - 1,
    size = 0.1
  ) +
  
  annotate(
    "segment",
    x = 3.5,
    xend = 3.5,
    y = group_bracket_y,
    yend = group_bracket_y - 1,
    size = 0.1
  ) +
  
  annotate(
    "text",
    x = 3,
    y = group_bracket_y + 2,
    label = "General Population",
    size = 5.5
  ) +
  
  annotate(
    "segment",
    x = 1,
    xend = 3,
    y = stat_bracket_y,
    yend = stat_bracket_y,
    size = 0.1
  ) +
  
  annotate(
    "segment",
    x = 1,
    xend = 1,
    y = stat_bracket_y,
    yend = group_bracket_y - 10,
    size = 0.1
  ) +
  
  annotate(
    "segment",
    x = 3,
    xend = 3,
    y = stat_bracket_y,
    yend = group_bracket_y + 3.5,
    size = 0.1
  ) +
  
  annotate(
    "text",
    x = 2,
    y = stat_bracket_y + 2.5,
    label = stat_label,
    size = 5.5,
    parse = TRUE
  ) +
  
  scale_x_continuous(
    breaks = plot_summary$x_pos,
    labels = plot_summary$Label,
    limits = c(0.5, 4)
  ) +
  
  scale_y_continuous(
    limits = c(0, y_max + 5),
    expand = c(0, 0)
  ) +
  
  theme_classic(base_family = "Helvetica") +
  
  labs(
    y = "Proportion of ginger rabbits (%)",
    x = NULL
  ) +
  
  theme(
    axis.text.x = element_text(
      size = 13,
      color = "black",
      margin = margin(t = 10),
      lineheight = 0.9
    ),
    
    axis.text.y = element_text(
      size = 13,
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    
    axis.line = element_line(linewidth = 0.1),
    axis.ticks.x = element_blank()
  )

print(a)

ggsave(
  file.path(out_dir, "eagle_prey_barplot.pdf"),
  a,
  width = 5,
  height = 5,
  dpi = 300
)