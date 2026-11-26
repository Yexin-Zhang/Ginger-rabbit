###############################################################################
## 06_upgma_tree.R
## UPGMA tree of 129 rabbits homozygous for the ginger or grey MC1R allele,
## built from identity-by-state distances (PLINK v1.9, --distance square 1-ibs)
## across a 50-kb region flanking MC1R (NC_091449.1:69,400,000-69,450,000).
## Supports: Fig. 2c,d
## Input:  data/tree/ (IBS distance matrix, sample IDs and sample information)
## Output: results/06_upgma_tree/
## Run from the repository root: Rscript scripts/06_upgma_tree.R
###############################################################################

library(ape)
library(phangorn)

out_dir <- "results/06_upgma_tree"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ids <- read.table("data/tree/mc1r_region_ibs_distances.mdist.id")$V1
distances <- as.matrix(read.table("data/tree/mc1r_region_ibs_distances.mdist"))
dimnames(distances) <- list(ids, ids)

samples <- read.delim("data/tree/tree_samples.txt")
samples <- samples[match(ids, samples$Sample), ]

tree <- upgma(as.dist(distances))
write.tree(tree, file.path(out_dir, "mc1r_upgma_tree.nwk"))

# Plot: tips coloured by MC1R allele; wild rabbits as circles, domestic as
# triangles
tip_info   <- samples[match(tree$tip.label, samples$Sample), ]
tip_colour <- ifelse(tip_info$MC1R == "Ginger", "#E67F52", "#9E9E9E")
tip_shape  <- ifelse(tip_info$Population == "Domestic", 17, 16)

pdf(file.path(out_dir, "mc1r_upgma_tree.pdf"), width = 8, height = 8)
plot(tree, type = "fan", show.tip.label = FALSE, edge.width = 0.5)
tiplabels(pch = tip_shape, col = tip_colour, cex = 0.8)
legend("topleft", bty = "n", cex = 0.8,
       legend = c("Ginger allele", "Grey allele", "Wild", "Domestic"),
       pch = c(15, 15, 16, 17),
       col = c("#E67F52", "#9E9E9E", "black", "black"))
dev.off()

cat("Tree with", length(tree$tip.label), "tips written to", out_dir, "\n")
print(table(samples$Population, samples$MC1R))
