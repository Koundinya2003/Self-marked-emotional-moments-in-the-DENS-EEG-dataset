# 05: confirmatory analyses H1, H2, H3a, H3b (preregistration section 5).
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2); library(patchwork) })
source("code/functions/hypotheses.R")

ep <- read_csv("results/features_epoch.csv", show_col_types = FALSE)
dir.create("results/models", showWarnings = FALSE)
dir.create("figures/diagnostics", recursive = TRUE, showWarnings = FALSE)

res <- lapply(names(H), function(k) {
  x <- run_hyp(H[[k]], ep, ep$rejected)
  saveRDS(x$model, file.path("results/models", paste0(k, ".rds")))
  x$row <- cbind(hypothesis = k, label = H[[k]]$label, x$row)
  x$log <- cbind(hypothesis = k, x$log)
  x
})
names(res) <- names(H)

conf <- bind_rows(lapply(res, `[[`, "row"))
conf$p_fdr <- p.adjust(conf$p, method = "BH")
conf$predicted_sign <- c(H1 = -1, H2 = 1, H3a = 1, H3b = 1)[conf$hypothesis]
conf$supported <- conf$p_fdr < 0.05 & sign(conf$estimate) == conf$predicted_sign
write_csv(conf, "results/confirmatory.csv")
write_csv(bind_rows(lapply(res, `[[`, "log")), "results/model_selection.csv")
print(conf |> select(hypothesis, estimate, ci_low, ci_high, df, p, p_fdr, std_effect, n_obs, n_subjects, structure, supported))

# descriptive condition means per subject for H1 (used in figures and text)
h1d <- res$H1$data
write_csv(h1d |> group_by(subject, condition) |> summarise(ln_alpha_posterior = mean(ln_alpha_posterior), n = n(), .groups = "drop"),
          "results/h1_subject_means.csv")

# ------------------------------------------------------------ diagnostics
diag_plots <- function(m, title) {
  df <- data.frame(fitted = fitted(m), resid = resid(m))
  re <- ranef(m)
  ri <- bind_rows(lapply(names(re), function(g) data.frame(group = g, value = re[[g]][, "(Intercept)"])))
  p1 <- ggplot(df, aes(sample = resid)) + stat_qq(size = 0.6, alpha = 0.6) + stat_qq_line(colour = "#D55E00") +
    labs(title = "Residual QQ", x = "Theoretical", y = "Sample")
  p2 <- ggplot(df, aes(fitted, resid)) + geom_point(size = 0.6, alpha = 0.5) +
    geom_hline(yintercept = 0, colour = "#D55E00") + geom_smooth(se = FALSE, colour = "#0072B2", method = "loess", formula = y ~ x) +
    labs(title = "Residuals vs fitted", x = "Fitted", y = "Residual")
  p3 <- ggplot(ri, aes(sample = value)) + stat_qq(size = 0.8) + stat_qq_line(colour = "#D55E00") +
    facet_wrap(~group, scales = "free") + labs(title = "Random intercepts QQ", x = "Theoretical", y = "Conditional mode")
  (p1 | p2) / p3 + plot_annotation(title = title) & theme_bw(base_size = 9)
}
for (k in names(res)) {
  p <- diag_plots(res[[k]]$model, paste(k, H[[k]]$label))
  ggsave(file.path("figures/diagnostics", paste0(k, "_diagnostics.png")), p, width = 8, height = 6, dpi = 200)
  ggsave(file.path("figures/diagnostics", paste0(k, "_diagnostics.pdf")), p, width = 8, height = 6)
}
# numeric assumption summaries
assump <- bind_rows(lapply(names(res), function(k) {
  r <- resid(res[[k]]$model)
  data.frame(hypothesis = k, resid_skew = mean((r - mean(r))^3) / sd(r)^3,
             resid_excess_kurtosis = mean((r - mean(r))^4) / sd(r)^4 - 3,
             shapiro_p_subsample = shapiro.test(if (length(r) > 5000) sample(r, 5000) else r)$p.value)
}))
write_csv(assump, "results/assumption_checks.csv")
