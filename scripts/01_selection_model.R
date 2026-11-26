###############################################################################
## 01_selection_model.R
## Eagle predation rate required for the ginger allele to reach the frequency
## observed in the Strzelecki Desert, as a function of starting frequency.
## Supports: Extended Data Fig. 6; Supplementary Table 6
## Input:  none (parameters set below)
## Output: results/01_selection_model/
## Run from the repository root: Rscript scripts/01_selection_model.R
###############################################################################
# Parameters (Methods; SI Section 16)
out_dir <- "results/01_selection_model"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Target: ginger allele frequency at Kain Bore, Strzelecki Desert (42.9%)
targets <- list(
  genetic    = list(p = 0.429,
                    file = "genetic")
)

n_gen_seq <- c(100, 130, 150, 200)     # generations (one curve each)

p0_seq <- seq(0.01, 0.20, by = 0.0025)    # starting allele frequency, 1-20%
k_seq  <- seq(0.01, 1.00, by = 0.001)     # % rabbits killed by eagles, 1-100%

# Proportional reduction in eagle predation of ginger rabbits (from nest data)
d_seq  <- 0.41

# Recursion for a recessive allele A: w_AA = 1, w_Aa = w_aa = 1 - s
# p' = (p^2 + p*q*(1 - s)) / (1 - s*(1 - p^2))
sim_p <- function(p0, s, n) {
  p <- p0
  for (i in seq_len(n)) {
    q <- 1 - p
    p <- (p * p + p * q * (1 - s)) / (1 - s * (1 - p * p))
  }
  p
}

# Evaluate all parameter combinations
grid <- expand.grid(p0 = p0_seq, k = k_seq, d = d_seq, n_gen = n_gen_seq)
grid$s     <- grid$k * grid$d
grid$p_end <- mapply(sim_p, grid$p0, grid$s, grid$n_gen)

ceiling_p <- with(grid[grid$k == max(k_seq), ],
                  setNames(p_end, paste(p0, n_gen)))

# For each starting frequency, the predation rate whose final frequency is
# closest to the target; "reached" is FALSE if k = 1 cannot attain the target
solve_for <- function(obs_p) {
  g <- grid
  g$dist <- abs(g$p_end - obs_p)
  b <- do.call(rbind, lapply(split(g, list(g$p0, g$n_gen), drop = TRUE),
                             function(x) x[which.min(x$dist), ]))
  b <- b[order(b$n_gen, b$p0), ]
  rownames(b) <- NULL
  b$reached <- ceiling_p[paste(b$p0, b$n_gen)] >= obs_p
  b
}

# Plot (Extended Data Fig. 6)
cols <- c("#B3541E", "#D98324", "#2E6F73", "#41436A")

make_plot <- function(best, obs_p, sub) {
  plot(NA, xlim = c(0, 21), ylim = c(0, 100), las = 1, bty = "l",
       xlab = "Starting ginger allele frequency (%)",
       ylab = "Rabbits killed by eagles required (%)",
       cex.main = 0.95)
  mtext(sub, side = 3, line = 0.15, cex = 0.72, col = "grey35")
  grid(col = "grey90", lty = 1)
  for (i in seq_along(n_gen_seq)) {
    b  <- best[best$n_gen == n_gen_seq[i], ]
    ok <- b$reached
    lines(b$p0[ok] * 100, b$k[ok] * 100, col = cols[i], lwd = 2.5)
    points(b$p0[ok] * 100, b$k[ok] * 100, col = cols[i], pch = 16, cex = 0.55)
  }
  abline(h = 100, lty = 3, col = "grey40")
  legend("topright", bty = "n", cex = 0.85, lwd = 2.5, pch = 16,
         col = cols[seq_along(n_gen_seq)],
         legend = sprintf("%d generations", n_gen_seq))
}

show_p0 <- c(0.01, 0.02, 0.04, 0.05, 0.08, 0.10, 0.15, 0.20)
results <- list()

for (nm in names(targets)) {
  tg   <- targets[[nm]]
  best <- solve_for(tg$p)
  results[[nm]] <- best

  write.csv(best, file.path(out_dir, sprintf("required_predation_rate_%s.csv", tg$file)),
            row.names = FALSE)

  png(file.path(out_dir, sprintf("required_predation_rate_%s.png", tg$file)),
      width = 1900, height = 1450, res = 220)
  par(mar = c(4.5, 4.5, 3.2, 1)); make_plot(best, tg$p, tg$lab); dev.off()

  pdf(file.path(out_dir, sprintf("required_predation_rate_%s.pdf", tg$file)),
      width = 7.5, height = 5.5)
  par(mar = c(4.5, 4.5, 3.2, 1)); make_plot(best, tg$p, tg$lab); dev.off()
}

# Summary table (Supplementary Table 6)
cell <- function(best, p0, n) {
  r <- best[best$n_gen == n & abs(best$p0 - p0) < 1e-9, ]
  if (!r$reached) "not reachable" else sprintf("%.0f%%  (s=%.2f)", r$k * 100, r$s)
}

for (nm in names(targets)) {
  tg   <- targets[[nm]]
  best <- results[[nm]]

  summ <- data.frame(`starting freq` = sprintf("%.0f%%", show_p0 * 100),
                     check.names = FALSE, stringsAsFactors = FALSE)
  for (n in n_gen_seq) {
    summ[[sprintf("%d gen", n)]] <- sapply(show_p0, cell, best = best, n = n)
  }
  write.csv(summ, file.path(out_dir, sprintf("summary_table_%s.csv", tg$file)), row.names = FALSE)

  cat(sprintf("\n\n%s\n", toupper(nm)))
  cat(sprintf("%s\n", tg$lab))
  cat(sprintf("d = %.2f (ginger %.0f%% less likely to be taken; RR = %.2f)\n\n",
              d_seq, d_seq * 100, 1 - d_seq))
  cat(formatC("starting freq", width = 16, flag = "-"))
  for (n in n_gen_seq) cat(formatC(sprintf("%d gen", n), width = 16, flag = "-"))
  cat("\n")
  for (i in seq_along(show_p0)) {
    cat(formatC(summ[[1]][i], width = 16, flag = "-"))
    for (n in n_gen_seq) {
      cat(formatC(summ[[sprintf("%d gen", n)]][i], width = 16, flag = "-"))
    }
    cat("\n")
  }
}
cat("\n\"not reachable\" = the target is not attained even if eagles killed\n")
cat("every rabbit. The recursion sees k and d only through their product, so\n")
cat("every figure rescales by 1/d if the predation difference is revised.\n")
