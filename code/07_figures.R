# 07: publication figures (PNG 300 dpi and PDF) in figures/.
# Colours: Okabe-Ito blue (#0072B2, control) and vermillion (#D55E00, marked),
# validated for colour vision deficiency; diverging maps use ColorBrewer RdBu
# with a neutral light grey midpoint. Condition identity is always also given
# by labels or position, never colour alone.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork); library(purrr)
})
source("code/functions/montage.R")
dir.create("figures", showWarnings = FALSE)

COL <- c(control = "#0072B2", marked = "#D55E00")
DIV <- c("#2166AC", "#67A9CF", "#D1E5F0", "#F2F2F2", "#FDDBC7", "#EF8A62", "#B2182B")
theme_set(theme_classic(base_size = 9) +
            theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
                  strip.background = element_blank(), strip.text = element_text(face = "bold"),
                  plot.title = element_text(face = "bold", size = 10)))
# Negative numbers are typeset with a true minus sign (U+2212) in tick labels
# and annotations. ggplot2 looks up default position scales by name, so these
# wrappers become the defaults for every continuous axis.
mn <- function(x, d = 3) gsub("-", "\u2212", formatC(x, format = "f", digits = d), fixed = TRUE)
minus_lab <- scales::label_number(style_negative = "minus", drop0trailing = TRUE)
scale_x_continuous <- function(..., labels = minus_lab) ggplot2::scale_x_continuous(..., labels = labels)
scale_y_continuous <- function(..., labels = minus_lab) ggplot2::scale_y_continuous(..., labels = labels)
save_fig <- function(p, name, w, h) {
  ggsave(file.path("figures", paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path("figures", paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
}

ep <- read_csv("results/features_epoch.csv", show_col_types = FALSE)
qc <- read_csv("results/qc_summary.csv", show_col_types = FALSE)
conf <- read_csv("results/confirmatory.csv", show_col_types = FALSE)
rob <- read_csv("results/robustness.csv", show_col_types = FALSE)
chan <- read_csv("results/channel_bandpower.csv", show_col_types = FALSE)
spec <- read_csv("results/region_spectra.csv", show_col_types = FALSE)
locs <- read_csv("results/channel_locations.csv", show_col_types = FALSE)
source("code/functions/models.R")

# ------------------------------------------------------------- 1. QC overview
q <- qc |> mutate(subject = sub("sub-", "", subject),
                  subject = factor(subject, levels = subject[order(available_min / duration_min)]))
p1a <- ggplot(q, aes(y = subject)) +
  geom_col(aes(x = duration_min), fill = "grey85", width = 0.7) +
  geom_col(aes(x = available_min), fill = "#0072B2", width = 0.7) +
  labs(x = "Recording length (min)", y = NULL, title = "A  Recorded data in v1.0.6",
       subtitle = "grey: header length; blue: samples in .fdt")
q2 <- q |> filter(status == "ok")
p1b <- ggplot(q2, aes(y = subject, x = n_bad_channels)) + geom_col(fill = "grey40", width = 0.7) +
  labs(x = "Bad channels (interpolated)", y = NULL, title = "B  Bad channels")
p1c <- q2 |> select(subject, eye = n_ic_eye, muscle = n_ic_muscle, cardiac = n_ic_cardiac) |>
  pivot_longer(-subject, names_to = "type", values_to = "n") |>
  mutate(type = factor(type, levels = c("eye", "muscle", "cardiac"))) |>
  ggplot(aes(y = subject, x = n, fill = type)) + geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = c(eye = "#0072B2", muscle = "#D55E00", cardiac = "#009E73"), name = NULL) +
  labs(x = "ICA components removed (of 70)", y = NULL, title = "C  Removed components") +
  theme(legend.position = "bottom")
p1d <- q2 |> transmute(subject, marked = ifelse(n_marked > 0, n_marked_rejected / n_marked, NA), control = ifelse(n_control > 0, n_control_rejected / n_control, NA)) |>
  pivot_longer(-subject, names_to = "condition", values_to = "frac") |>
  ggplot(aes(y = subject, x = frac, colour = condition, shape = condition)) + geom_point(size = 2) +
  scale_colour_manual(values = COL, name = NULL) + scale_shape_manual(values = c(control = 16, marked = 17), name = NULL) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Epochs rejected after channel repair", y = NULL, title = "D  Rejected epochs", subtitle = "missing point: no epochs") + theme(legend.position = "bottom")
save_fig((p1a | p1b | p1c | p1d) + plot_layout(widths = c(1.4, 0.8, 1, 1)), "fig_qc_overview", 11, 7)

# ------------------------------------------------- 2. QC spectra and ICA metrics
qc_subjects <- c("sub-mit003", "sub-mit039", "sub-mit123")
qs <- map_dfr(qc_subjects, function(s) {
  f <- file.path("results/subjects", paste0(s, ".rds"))
  if (!file.exists(f)) return(NULL)
  o <- readRDS(f)
  tibble(subject = s, f = o$qc_spec$f, unfiltered = o$qc_spec$unfiltered, `filtered, before ICA` = o$qc_spec$before_ica,
         `after ICA and interpolation` = o$qc_spec$after_ica)
}) |> pivot_longer(-c(subject, f), names_to = "stage", values_to = "psd") |> filter(f >= 0.5, f <= 60)
p2a <- ggplot(qs, aes(f, psd, colour = stage, linetype = stage)) + geom_line(linewidth = 0.5) +
  scale_y_log10() + facet_wrap(~subject, nrow = 1) +
  scale_colour_manual(values = c("grey55", "#0072B2", "#D55E00"), name = NULL) +
  scale_linetype_manual(values = c("dotted", "solid", "solid"), name = NULL) +
  labs(x = "Frequency (Hz)", y = expression("Median PSD ("*mu*V^2*"/Hz)"), title = "A  Channel median spectra of the video periods") +
  theme(legend.position = "bottom")
ic <- map_dfr(qc_subjects, function(s) {
  f <- file.path("results/subjects", paste0(s, ".rds"))
  if (!file.exists(f)) return(NULL)
  o <- readRDS(f)
  o$ica_classes |> mutate(subject = s)
}) |> mutate(class = case_when(eye ~ "eye", muscle ~ "muscle", cardiac ~ "cardiac", TRUE ~ "kept"),
             r_eog = pmax(r_veog, r_heog))
p2b <- ggplot(ic, aes(hf_ratio, pmax(focal, outer_frac), colour = class, shape = class)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(size = 1.6, alpha = 0.85) + scale_x_log10() + facet_wrap(~subject, nrow = 1) +
  scale_colour_manual(values = c(kept = "grey60", eye = "#0072B2", muscle = "#D55E00", cardiac = "#009E73"), name = NULL) +
  labs(x = "Power ratio, 20 to 40 Hz over 2 to 20 Hz (log)", y = "max(focality, outer ring fraction)",
       title = "B  ICA component classification") + theme(legend.position = "bottom")
save_fig(p2a / p2b, "fig_qc_subjects", 10, 6.5)

# ------------------------------------------------------------- 3. ratings
vids <- ep |> filter(!rejected, condition == "marked", !is.na(valence)) |> distinct(subject, trial_index, stim, valence, arousal)
allv <- read_csv("results/audit_trials.csv", show_col_types = FALSE) |> filter(!is.na(valence)) |>
  mutate(set = "all rated videos") |>
  bind_rows(vids |> mutate(set = "videos in H2/H3"))
p3a <- allv |> pivot_longer(c(valence, arousal), names_to = "scale") |>
  ggplot(aes(value, fill = set)) + geom_histogram(binwidth = 0.5, boundary = 1, position = "identity", alpha = 0.7, colour = "white", linewidth = 0.2) +
  facet_wrap(~scale) + scale_fill_manual(values = c("grey60", "#0072B2"), name = NULL) +
  labs(x = "Rating (1 to 9)", y = "Videos", title = "A  Rating distributions") + theme(legend.position = "bottom")
p3b <- ggplot(vids, aes(valence, arousal)) + geom_jitter(width = 0.08, height = 0.08, size = 1.4, colour = "#0072B2", alpha = 0.7) +
  labs(x = "Valence", y = "Arousal", title = "B  Valence vs arousal", subtitle = "videos in H2/H3") + coord_equal(xlim = c(1, 9), ylim = c(1, 9))
save_fig(p3a | p3b, "fig_ratings", 9, 4)

# ----------------------------------------------- 4. band power by condition (H1)
h1d <- h1_set(ep, ep$rejected)
sm <- h1d |> group_by(subject, condition) |> summarise(v = mean(ln_alpha_posterior), .groups = "drop")
paired <- sm |> pivot_wider(names_from = condition, values_from = v) |> filter(!is.na(marked), !is.na(control))
cH1 <- conf |> filter(hypothesis == "H1")
p4a <- ggplot(sm |> semi_join(paired, by = "subject"), aes(condition, v)) +
  geom_line(aes(group = subject), colour = "grey75", linewidth = 0.3) +
  geom_point(aes(colour = condition, shape = condition), size = 1.8) +
  scale_colour_manual(values = COL, guide = "none") + scale_shape_manual(values = c(control = 16, marked = 17), guide = "none") +
  scale_x_discrete(limits = c("control", "marked")) +
  labs(x = NULL, y = expression("ln posterior alpha power ("*mu*V^2*")"), title = "A  Posterior alpha, participant means",
       subtitle = sprintf("LMM, marked minus control: %s [%s, %s]", mn(cH1$estimate), mn(cH1$ci_low), mn(cH1$ci_high)))
sp <- spec |> filter(region == "posterior", freq >= 1, freq <= 40) |>
  group_by(condition, freq) |> summarise(m = mean(ln_psd), se = sd(ln_psd) / sqrt(n()), .groups = "drop")
p4b <- ggplot(sp, aes(freq, m, colour = condition, fill = condition)) +
  annotate("rect", xmin = 8, xmax = 13, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_ribbon(aes(ymin = m - se, ymax = m + se), alpha = 0.2, colour = NA) + geom_line(linewidth = 0.6) +
  scale_colour_manual(values = COL, name = NULL) + scale_fill_manual(values = COL, name = NULL) +
  labs(x = "Frequency (Hz)", y = expression("ln PSD ("*mu*V^2*"/Hz)"), title = "B  Posterior spectrum (mean \u00b1 SE over participants)") +
  theme(legend.position = c(0.8, 0.85))
spd <- spec |> filter(freq >= 1, freq <= 40) |> pivot_wider(names_from = condition, values_from = ln_psd) |>
  mutate(d = marked - control) |> filter(!is.na(d)) |>
  group_by(region, freq) |> summarise(m = mean(d), se = sd(d) / sqrt(n()), .groups = "drop")
p4c <- ggplot(spd, aes(freq, m)) + geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = m - 1.96 * se, ymax = m + 1.96 * se), fill = "#D55E00", alpha = 0.2) +
  geom_line(colour = "#D55E00", linewidth = 0.6) + facet_wrap(~region, nrow = 1) +
  labs(x = "Frequency (Hz)", y = "ln PSD, marked minus control", title = "C  Spectral difference by region (mean, 95% CI over participants)")
bp <- read_csv("results/features_long.csv", show_col_types = FALSE) |>
  filter(band %in% c("theta", "alpha", "beta"), region %in% c("frontal", "central", "posterior", "left", "right")) |>
  semi_join(h1d, by = c("subject", "epoch")) |>
  group_by(subject, condition, region, band) |> summarise(v = mean(ln_power), .groups = "drop") |>
  pivot_wider(names_from = condition, values_from = v) |> mutate(d = marked - control) |> filter(!is.na(d)) |>
  mutate(band = factor(band, levels = c("theta", "alpha", "beta")))
p4d <- ggplot(bp, aes(region, d)) + geom_hline(yintercept = 0, colour = "grey50") +
  geom_boxplot(outlier.shape = NA, width = 0.5, colour = "grey40") +
  geom_jitter(width = 0.12, size = 0.8, alpha = 0.6, colour = "#D55E00") + facet_wrap(~band) +
  labs(x = NULL, y = "ln power, marked minus control\n(participant means)", title = "D  Band power difference by region")
save_fig((p4a | p4b) / p4c / p4d + plot_layout(heights = c(1, 0.8, 1)), "fig_bandpower_condition", 10, 10)

# --------------------------------------------------------------- 5. topomaps
topo_xy <- locs |> filter(type == "EEG", !(label %in% OUTER_RING)) |>
  mutate(x = radius * sin(theta * pi / 180), y = radius * cos(theta * pi / 180))
topo_grid <- function(vals, res = 120) {
  d <- topo_xy |> inner_join(vals, by = c("label" = "channel"))
  rmax <- max(sqrt(d$x^2 + d$y^2)) * 1.02
  g <- interp::interp(d$x, d$y, d$value, nx = res, ny = res, linear = TRUE, extrap = FALSE)
  expand.grid(x = g$x, y = g$y) |> mutate(value = as.vector(g$z)) |>
    filter(sqrt(x^2 + y^2) <= rmax, !is.na(value)) |> mutate(rmax = rmax)
}
td <- chan |> pivot_wider(names_from = condition, values_from = c(ln_power, n_epochs)) |>
  mutate(d = ln_power_marked - ln_power_control) |> filter(!is.na(d)) |>
  group_by(channel, band) |> summarise(value = mean(d), t = mean(d) / (sd(d) / sqrt(n())), n = n(), .groups = "drop")
lim <- max(abs(td$value))
grids <- map_dfr(c("theta", "alpha", "beta"), function(b) topo_grid(td |> filter(band == b) |> select(channel, value)) |> mutate(band = b)) |>
  mutate(band = factor(band, levels = c("theta", "alpha", "beta")))
rmax <- max(grids$rmax)
circ <- tibble(a = seq(0, 2 * pi, length.out = 200), x = rmax * cos(a), y = rmax * sin(a))
nose <- tibble(x = c(-0.08, 0, 0.08) * rmax / 0.5, y = c(rmax * 0.99, rmax * 1.12, rmax * 0.99))
post <- topo_xy |> filter(label %in% PREREG_REGIONS$posterior)
p5 <- ggplot(grids, aes(x, y)) + geom_raster(aes(fill = value), interpolate = TRUE) +
  geom_contour(aes(z = value), colour = "grey30", linewidth = 0.15, bins = 8) +
  geom_path(data = circ, linewidth = 0.4) + geom_path(data = nose, linewidth = 0.4) +
  geom_point(data = topo_xy, size = 0.25, colour = "grey20") +
  geom_point(data = post, size = 0.9, shape = 1, colour = "black") +
  scale_fill_gradientn(colours = DIV, limits = c(-lim, lim), name = "ln power\nmarked minus\ncontrol", na.value = "transparent", labels = minus_lab) +
  facet_wrap(~band) + coord_equal() + theme_void(base_size = 9) +
  theme(strip.text = element_text(face = "bold", size = 10), legend.position = "right") +
  labs(title = "Marked minus control, mean over participants (open circles: posterior region)")
save_fig(p5, "fig_topomaps", 9, 3.4)
write_csv(td, "results/topography_differences.csv")

# ------------------------------------------- 6. confirmatory and robustness estimates
rr <- rob |> mutate(analysis = factor(analysis, levels = rev(sort(unique(analysis)))),
                    is_conf = analysis == "R0_confirmatory",
                    hypothesis = factor(hypothesis, levels = c("H1", "H2", "H3a", "H3b"),
                                        labels = c("H1 posterior alpha:\nmarked minus control (ln \u00b5V\u00b2)",
                                                   "H2 FAA ~ valence\n(per SD of valence)",
                                                   "H3a frontal theta ~ arousal\n(ln \u00b5V\u00b2 per SD)",
                                                   "H3b global beta ~ arousal\n(ln \u00b5V\u00b2 per SD)")))
p6 <- ggplot(rr, aes(estimate, analysis)) + geom_vline(xintercept = 0, colour = "grey50") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high, colour = is_conf), width = 0.2, linewidth = 0.5, orientation = "y") +
  geom_point(aes(colour = is_conf, shape = is_conf), size = 2) +
  scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey40"), guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16), guide = "none") +
  facet_wrap(~hypothesis, scales = "free_x", nrow = 1) +
  labs(x = "Estimate with 95% CI (unstandardised)", y = NULL,
       title = "Confirmatory estimates (R0, highlighted) and preregistered robustness analyses")
save_fig(p6, "fig_effects", 11, 3.6)

# -------------------------------------------- 7. rating associations (H2, H3)
rs <- rating_set(ep, ep$rejected)
vm <- rs |> group_by(subject, trial_index, valence, arousal) |>
  summarise(across(c(faa, ln_theta_frontal, ln_beta_global), mean), .groups = "drop") |>
  group_by(subject) |> mutate(across(c(faa, ln_theta_frontal, ln_beta_global), ~ .x - mean(.x), .names = "{.col}_wc")) |> ungroup()
mk <- function(xv, yv, xl, yl, ttl, h) {
  cc <- conf |> filter(hypothesis == h)
  ggplot(vm, aes(.data[[xv]], .data[[yv]])) + geom_hline(yintercept = 0, colour = "grey85") +
    geom_point(size = 1.3, alpha = 0.6, colour = "#0072B2") +
    geom_smooth(method = "lm", formula = y ~ x, colour = "#D55E00", fill = "#D55E00", alpha = 0.15, linewidth = 0.7) +
    labs(x = xl, y = yl, title = ttl, subtitle = sprintf("b = %s [%s, %s], p = %.3f", mn(cc$estimate), mn(cc$ci_low), mn(cc$ci_high), cc$p))
}
p7 <- mk("valence", "faa_wc", "Valence rating", "FAA (centred within participant)", "A  H2", "H2") |
  mk("arousal", "ln_theta_frontal_wc", "Arousal rating", "ln frontal theta\n(centred within participant)", "B  H3a", "H3a") |
  mk("arousal", "ln_beta_global_wc", "Arousal rating", "ln global beta\n(centred within participant)", "C  H3b", "H3b")
save_fig(p7 + plot_annotation(caption = "Points: clip means of marked epochs, centred within participant for display; line: OLS fit for display only."),
         "fig_rating_associations", 11, 3.8)

# ------------------------------------------------------- 8. heart rate (exploratory)
chg <- read_csv("results/click_hr_group.csv", show_col_types = FALSE)
p8 <- ggplot(chg, aes(t, mean)) + geom_hline(yintercept = 0, colour = "grey60") +
  annotate("rect", xmin = -4.5, xmax = -0.5, ymin = -Inf, ymax = Inf, fill = "#D55E00", alpha = 0.08) +
  annotate("rect", xmin = -10, xmax = -6, ymin = -Inf, ymax = Inf, fill = "grey60", alpha = 0.12) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), fill = "#0072B2", alpha = 0.2) +
  geom_line(colour = "#0072B2", linewidth = 0.7) +
  labs(x = "Time relative to click (s)", y = "Heart rate change (bpm)",
       title = "Heart rate around clicks (exploratory)",
       subtitle = "Grey: baseline (\u221210 to \u22126 s); orange: marked window")
save_fig(p8 + theme(plot.title = element_text(size = 8), plot.subtitle = element_text(size = 6.5)), "fig_heart_rate", 3.5, 2.6)

# ---------------------------------------------- 9. diagnostics overview
cat("Figures written to figures/\n")
