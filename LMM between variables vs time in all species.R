rm(list = ls())

## ---- paths ----
setwd("C:/Users/Hp/Desktop/NFU2/Human_demography/")
infile <- "Merge_anthro_variables_NE.csv"
outdir <- "lmm output 3"     # output folder

## ---- packages ----
pkgs <- c("tidyverse","lme4","lmerTest","broom.mixed","ggeffects","emmeans")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# Use asymptotic df in emmeans (more stable with tiny N)
emmeans::emm_options(df = "asymptotic")

## ---- plotting options ----
CI_LEVEL        <- 0.80   # 80% CI (set to 0.95 for 95% if you like)
SHOW_RIBBON     <- TRUE   # draw ribbons when they are reasonable
AUTO_HIDE_RATIO <- 0.90   # hide ribbon if width > 90% of y-range (absurd)

## ---- helpers ----
time_levels <- c("1CE","1700AD","2000BC","4000BC","8000BC")
time_cols <- c("1CE"="#D32F2F","1700AD"="#1976D2","2000BC"="#2E7D32","4000BC"="#7B1FA2","8000BC"="#F39C12")

# scale safely: center-only if SD is tiny
scale_safe <- function(x, tol = 1e-8) {
  s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s < tol) return(x - m)
  (x - m) / s
}

# drop time bins where predictor has no variance (prevents rank deficiency)
drop_const_times <- function(dat, pz){
  dat %>%
    dplyr::group_by(Time) %>%
    dplyr::filter(sd(.data[[pz]], na.rm = TRUE) > 0) %>%
    dplyr::ungroup() %>%
    droplevels()
}

# ensure lower/upper CL exist (fallback using normal quantiles)
ensure_ci <- function(df, level){
  if (!all(c("lower.CL","upper.CL") %in% names(df)) && all(c("estimate","SE") %in% names(df))) {
    z <- qnorm(1 - (1 - level)/2)
    df$lower.CL <- df$estimate - z*df$SE
    df$upper.CL <- df$estimate + z*df$SE
  }
  df
}

# NEW: standardize emtrends() output → columns: Time, estimate, SE, lower.CL, upper.CL
unify_trend_df <- function(x, level = CI_LEVEL) {
  df <- as.data.frame(x)
  
  # Time column
  time_col <- intersect(c("Time","time","TIME"), names(df))
  if (length(time_col) == 0) time_col <- names(df)[1]
  
  # slope (trend) column
  est_col <- grep("\\.trend$|trend$|^emmean$|^estimate$", names(df), value = TRUE)
  if (length(est_col) == 0) {
    num_cols <- names(df)[sapply(df, is.numeric)]
    est_col <- setdiff(num_cols, c("SE","se","Std..Error","std.error","df","lower.CL","upper.CL","asymp.LCL","asymp.UCL"))[1]
  }
  
  # SE column (if present)
  se_col <- intersect(c("SE","se","Std..Error","std.error"), names(df))
  SE <- if (length(se_col)) df[[se_col[1]]] else NA_real_
  
  # CI columns (various names across backends)
  lower_col <- intersect(c("lower.CL","asymp.LCL","LCL"), names(df))
  upper_col <- intersect(c("upper.CL","asymp.UCL","UCL"), names(df))
  
  out <- data.frame(
    Time     = df[[time_col]],
    estimate = df[[est_col]],
    SE       = SE
  )
  
  if (length(lower_col) && length(upper_col)) {
    out$lower.CL <- df[[lower_col[1]]]
    out$upper.CL <- df[[upper_col[1]]]
  } else {
    out <- ensure_ci(out, level)
  }
  out
}

## ---- load & prep ----
frog_raw <- read.csv(infile, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

frog <- frog_raw %>%
  dplyr::rename(
    cropland_Cao_2021 = `cropland Cao et al 2021`,
    cropland_HYDE_v3_5 = `cropland Hyde v.3.5`,
    Years_ago = `Years ago`
  ) %>%
  dplyr::mutate(Time = trimws(Time)) %>%
  dplyr::filter(Time %in% time_levels) %>%
  dplyr::mutate(
    Time  = factor(Time, levels = time_levels, ordered = TRUE),
    logNe = log10(Ne_median)
  )

stopifnot(all(c("ID","Time","Ne_median") %in% names(frog)))

# predictors present in file
predictors <- intersect(
  c("cropland_Cao_2021","cropland_HYDE_v3_5","rf_norice","tot_irri","rurc","urbc","grazing"),
  names(frog)
)

# safe scaling → z_*
for (v in predictors) frog[[paste0("z_", v)]] <- scale_safe(frog[[v]])
zpred <- paste0("z_", predictors)

# if SFS missing, treat all as folded
if (!"SFS" %in% names(frog)) frog$SFS <- "folded"

# outputs
dir.create(outdir, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE)

## ---- core: fit one predictor for one SFS (with fallback) ----
fit_one_predictor <- function(dat, label, pz){
  d <- drop_const_times(dat, pz)
  if (dplyr::n_distinct(d$Time) < 2 || sd(d[[pz]], na.rm = TRUE) == 0) return(NULL)
  
  # 1) Try LMM
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6))
  f_fac <- as.formula(paste0("logNe ~ ", pz, " * Time + (1|ID)"))
  fit <- try(lmer(f_fac, data = d, REML = TRUE, control = ctrl), silent = TRUE)
  model_used <- "lmer"
  
  # 2) Fallback to OLS if singular or failed
  if (inherits(fit, "try-error") || isTRUE(lme4::isSingular(fit, tol = 1e-5))) {
    f_fe <- as.formula(paste0("logNe ~ ", pz, " * Time + ID"))
    fit_try <- try(lm(f_fe, data = d, na.action = na.exclude), silent = TRUE)
    if (inherits(fit_try, "try-error")) {
      fit_try <- lm(as.formula(paste0("logNe ~ ", pz, " * Time")),
                    data = d, na.action = na.exclude)
    }
    fit <- fit_try
    model_used <- "lm"
  }
  
  # 3) Fixed effects table (no CI → avoids non-PD vcov issues)
  fx <- try(broom.mixed::tidy(fit, effects = "fixed", conf.int = FALSE), silent = TRUE)
  if (inherits(fx, "try-error")) {
    sm <- coef(summary(fit))
    fx <- tibble::tibble(term = rownames(sm),
                         estimate = sm[, "Estimate"],
                         std.error = sm[, "Std. Error"],
                         statistic = sm[, if ("t value" %in% colnames(sm)) "t value" else "t"])
  }
  fx$predictor <- sub("^z_","", pz)
  fx$SFS <- label
  fx$model_used <- model_used
  
  # 4) Predicted effect lines by Time with GUARDED CI
  xr <- quantile(d[[pz]], probs = c(0.05, 0.95), na.rm = TRUE, names = FALSE)
  yr <- range(d$logNe, na.rm = TRUE)
  y_pad <- diff(yr) * 0.15
  ylim_use <- c(yr[1] - y_pad, yr[2] + y_pad)
  
  gp <- try(
    ggeffects::ggpredict(
      fit,
      terms  = c(sprintf("%s [%.6f:%.6f]", pz, xr[1], xr[2]), "Time"),
      ci.lvl = CI_LEVEL
    ),
    silent = TRUE
  )
  
  if (!inherits(gp, "try-error")) {
    gdf <- as.data.frame(gp)  # x, predicted, conf.low, conf.high, group (= Time)
    ribbon_too_wide <- any((gdf$conf.high - gdf$conf.low) > (diff(ylim_use) * AUTO_HIDE_RATIO), na.rm = TRUE)
    show_ribbon_now <- isTRUE(SHOW_RIBBON) && !ribbon_too_wide
    
    p <- ggplot2::ggplot(gdf, ggplot2::aes(x = x, y = predicted, color = group))
    if (show_ribbon_now) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = conf.low, ymax = conf.high, fill = group),
        alpha = 0.18, colour = NA
      )
    }
    p <- p +
      ggplot2::geom_line(linewidth = 1.1) +
      ggplot2::scale_color_manual(values = time_cols, drop = FALSE) +
      ggplot2::scale_fill_manual(values = time_cols, drop = FALSE) +
      ggplot2::coord_cartesian(ylim = ylim_use) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste0("Effect of ", sub("^z_","",pz), " × Time (", label, ", ", model_used, ")"),
        x = sub("^z_","",pz), y = "Predicted log10(Ne)", color = "Time", fill = "Time"
      )
    
    ggplot2::ggsave(
      file.path(outdir, "plots", paste0("effect_", sub("^z_","",pz), "_", label, ".png")),
      p, width = 7, height = 5, dpi = 300
    )
  }
  
  # 5) Raw scatter + lm trend, faceted by Time (sanity check)
  base_var <- sub("^z_","", pz)
  rawp <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[base_var]], y = logNe)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, color = "black") +
    ggplot2::facet_wrap(~Time, ncol = 3) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = paste0("Raw: ", base_var, " vs logNe by Time (", label, ")"),
                  x = base_var, y = "log10(Ne)")
  ggplot2::ggsave(file.path(outdir, "plots",
                            paste0("raw_", base_var, "_", label, ".png")),
                  rawp, width = 9, height = 6, dpi = 300)
  
  # 6) Per-time slopes (effect of predictor within each epoch) — unified columns
  slope_tab <- NULL
  slope_try <- try(emmeans::emtrends(fit, ~ Time, var = pz, level = CI_LEVEL), silent = TRUE)
  if (!inherits(slope_try, "try-error")) {
    st <- unify_trend_df(slope_try, level = CI_LEVEL)
    st$predictor <- base_var
    st$SFS <- label
    st$model_used <- model_used
    slope_tab <- st
  }
  
  list(fixed = fx, slopes = slope_tab)
}

## ---- run all predictors per SFS ----
fit_univariate_set <- function(dat, label){
  fixed_all <- list(); slopes_all <- list()
  for (pz in zpred) {
    if (!pz %in% names(dat)) next
    out <- fit_one_predictor(dat, label, pz)
    if (is.null(out)) next
    fixed_all[[pz]]  <- out$fixed
    slopes_all[[pz]] <- out$slopes
  }
  list(
    fixed  = if (length(fixed_all)) dplyr::bind_rows(fixed_all) else NULL,
    slopes = if (length(slopes_all)) dplyr::bind_rows(slopes_all) else NULL
  )
}

spl <- split(frog, frog$SFS)
res <- lapply(names(spl), function(nm) fit_univariate_set(spl[[nm]], nm))
names(res) <- names(spl)

## ---- save tables ----
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
for (nm in names(res)) {
  if (!is.null(res[[nm]]$fixed)) {
    write.csv(res[[nm]]$fixed,
              file.path(outdir, paste0("fixed_effects_", nm, ".csv")),
              row.names = FALSE)
  }
  if (!is.null(res[[nm]]$slopes)) {
    write.csv(res[[nm]]$slopes,
              file.path(outdir, paste0("per_time_slopes_", nm, ".csv")),
              row.names = FALSE)
  }
}

## ---- make forest plots for slopes (per predictor & SFS) ----
make_forest <- function(slopes_df, predictor_label, sfs_label){
  dfp <- subset(slopes_df, predictor == predictor_label)
  if (!nrow(dfp)) return(invisible(NULL))
  dfp$Time <- factor(dfp$Time, levels = time_levels, ordered = TRUE)
  
  g <- ggplot2::ggplot(dfp, ggplot2::aes(x = Time, y = estimate)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower.CL, ymax = upper.CL), width = 0.15) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste0("Slope of Ne vs ", predictor_label, " by Time (", sfs_label, ")"),
      y = "d log10(Ne) / d predictor", x = "Time"
    )
  
  ggplot2::ggsave(file.path(outdir, "plots",
                            paste0("forest_", predictor_label, "_", sfs_label, ".png")),
                  g, width = 7, height = 5, dpi = 300)
}

for (nm in names(res)) {
  if (!is.null(res[[nm]]$slopes)) {
    for (pred in unique(res[[nm]]$slopes$predictor)) {
      make_forest(res[[nm]]$slopes, pred, nm)
    }
  }
}

cat("Done. Results saved to:\n",
    normalizePath(outdir), "\n",
    normalizePath(file.path(outdir, "plots")), "\n")



# Fig 5 companion: species-colored overlap through time
# ============================
library(tidyverse)

# --- paths ---
setwd("C:/Users/Hp/Desktop/NFU2/Human_demography/")
infile <- "Merge_anthro_variables_NE.csv"
outdir <- "lmm output 3/fig5_species"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- config ---
time_levels <- c("1CE","1700AD","2000BC","4000BC","8000BC")
species_order <- c("Dryophytes flaviventris",
                   "Dryophytes japonicus",
                   "Dryophytes immaculatus",
                   "Dryophytes suweonensis")

# match your Ne figure colors
sp_cols <- c(
  "Dryophytes flaviventris" = "#F39C12",  # orange
  "Dryophytes japonicus"    = "#2E7D32",  # green
  "Dryophytes immaculatus"  = "#1976D2",  # blue
  "Dryophytes suweonensis"  = "#D32F2F"   # red
)

# overlap variables to include (only those present will be used)
preds <- c("cropland Cao et al 2021",
           "cropland Hyde v.3.5",
           "rf_norice",
           "tot_irri",
           "rurc",
           "urbc",
           "grazing")

# pretty facet labels
pred_labels <- c(
  "cropland Cao et al 2021" = "Cropland (Cao 2021)",
  "cropland Hyde v.3.5"     = "Cropland (HYDE 3.5)",
  "rf_norice"               = "Rainfed (no rice)",
  "tot_irri"                = "Irrigated (total)",
  "rurc"                    = "Rural built-up",
  "urbc"                    = "Urban built-up",
  "grazing"                 = "Grazing"
)

# --- load & tidy ---
df0 <- read.csv(infile, check.names = FALSE, stringsAsFactors = FALSE)

preds <- preds[preds %in% names(df0)]
stopifnot(length(preds) > 0)

df <- df0 %>%
  mutate(Time = trimws(Time)) %>%
  filter(Time %in% time_levels) %>%
  mutate(Time = factor(Time, levels = time_levels, ordered = TRUE)) %>%
  pivot_longer(all_of(preds), names_to = "variable", values_to = "value") %>%
  # average across SFS within species × time (keeps figure simple)
  group_by(ID, Time, variable) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    ID       = factor(ID, levels = species_order),
    variable = factor(variable, levels = preds, labels = pred_labels[preds])
  ) %>%
  drop_na(ID, value)

# convert to % if data are 0–1
if (max(df$value, na.rm = TRUE) <= 1.05) df <- df %>% mutate(value = 100 * value)

# ============================
# A) Species-colored LINES (facet by variable)
# ============================
p_lines <- ggplot(df, aes(Time, value, color = ID, group = ID)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 2) +
  facet_wrap(~ variable, ncol = 3, labeller = label_wrap_gen(width = 20)) +
  scale_color_manual(values = sp_cols, name = "Species", drop = FALSE) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Species overlap through time (by variable)",
    x = NULL,
    y = "Overlap (%)"
  )

ggsave(file.path(outdir, "fig5_species_lines.png"),
       p_lines, width = 10, height = 6.8, dpi = 300)

# ============================
# B) HEATMAP (facet by variable; species × time)
# ============================
p_heat <- ggplot(df, aes(Time, ID, fill = value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~ variable, ncol = 3, labeller = label_wrap_gen(width = 20)) +
  scale_fill_viridis_c(name = "Overlap (%)", option = "C", direction = 1) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  labs(
    title = "Species × time overlap (by variable)",
    x = NULL, y = NULL
  )

ggsave(file.path(outdir, "fig5_species_heatmap.png"),
       p_heat, width = 10, height = 6.8, dpi = 300)

cat("Saved:\n",
    normalizePath(file.path(outdir, "fig5_species_lines.png")), "\n",
    normalizePath(file.path(outdir, "fig5_species_heatmap.png")), "\n")



# Fig 5 companion — species-coloured overlap vs time (FOLDED ONLY)
library(tidyverse)


setwd("C:/Users/Hp/Desktop/NFU2/Human_demography/")
infile <- "Merge_anthro_variables_NE.csv"
outdir <- "lmm output 3/fig5_species_folded"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- config ---
time_levels <- c("1CE","1700AD","2000BC","4000BC","8000BC")
species_order <- c("Dryophytes flaviventris",
                   "Dryophytes japonicus",
                   "Dryophytes immaculatus",
                   "Dryophytes suweonensis")

# match your Ne palette
sp_cols <- c(
  "Dryophytes flaviventris" = "#F39C12",
  "Dryophytes japonicus"    = "#2E7D32",
  "Dryophytes immaculatus"  = "#1976D2",
  "Dryophytes suweonensis"  = "#D32F2F"
)

# columns to include (only those present are used)
preds <- c("cropland Cao et al 2021","cropland Hyde v.3.5","rf_norice",
           "tot_irri","rurc","urbc","grazing")

nice_lab <- c(
  "cropland Cao et al 2021" = "Cropland (Cao 2021)",
  "cropland Hyde v.3.5"     = "Cropland (HYDE 3.5)",
  "rf_norice"               = "Rainfed (no rice)",
  "tot_irri"                = "Irrigated (total)",
  "rurc"                    = "Rural built-up",
  "urbc"                    = "Urban built-up",
  "grazing"                 = "Grazing"
)

# --- load ---
df0 <- read.csv(infile, check.names = FALSE, stringsAsFactors = FALSE)

# keep only available predictors
preds <- preds[preds %in% names(df0)]
stopifnot(length(preds) > 0)

# default SFS if missing
if (!"SFS" %in% names(df0)) df0$SFS <- "folded"

# ---------- tidy LONG, FOLDED ONLY, median aggregate ----------
df <- df0 |>
  mutate(Time = trimws(Time)) |>
  filter(Time %in% time_levels, SFS == "folded") |>
  mutate(Time = factor(Time, levels = time_levels, ordered = TRUE)) |>
  pivot_longer(all_of(preds), names_to = "variable", values_to = "value") |>
  group_by(ID, Time, variable) |>
  summarise(value = median(value, na.rm = TRUE), .groups = "drop") |>
  mutate(
    ID       = factor(ID, levels = species_order),
    variable = factor(variable, levels = preds, labels = nice_lab[preds])
  ) |>
  drop_na(ID, value)

# ---- auto-fix units: if values look like 0–10000 percentages, divide by 100 ----
maxv <- max(df$value, na.rm = TRUE)
if (is.finite(maxv) && maxv > 100 && maxv <= 10000) {
  df <- df |> mutate(value = value / 100)
}
# convert to % if still in 0–1
if (max(df$value, na.rm = TRUE) <= 1.05) df <- df |> mutate(value = 100 * value)

# =========================
# A) Lines — species colours, facet by variable, FREE y-scale
# =========================
p_lines_free <- ggplot(df, aes(Time, value, color = ID, group = ID)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~ variable, ncol = 3, scales = "free_y",
             labeller = label_wrap_gen(width = 20)) +
  scale_color_manual(values = sp_cols, name = "Species", drop = FALSE) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(title = "Species overlap through time (folded SFS, median)",
       x = NULL, y = "Overlap (%)")

ggsave(file.path(outdir, "fig5_species_lines_freeY_folded.png"),
       p_lines_free, width = 10, height = 6.8, dpi = 300)

# =========================
# B) Lines — NORMALIZED (0–1) within each variable (shape comparison)
# =========================
df_norm <- df |>
  group_by(variable) |>
  mutate(value_norm = (value - min(value, na.rm = TRUE)) /
           (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))) |>
  ungroup()

p_lines_norm <- ggplot(df_norm, aes(Time, value_norm, color = ID, group = ID)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~ variable, ncol = 3, labeller = label_wrap_gen(width = 20)) +
  scale_color_manual(values = sp_cols, name = "Species", drop = FALSE) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(title = "Species overlap through time — normalized (folded SFS, median)",
       x = NULL, y = "Normalized overlap (0–1)")

ggsave(file.path(outdir, "fig5_species_lines_normalized_folded.png"),
       p_lines_norm, width = 10, height = 6.8, dpi = 300)

# =========================
# C) Heatmap — species × time (folded only)
# =========================
p_heat <- ggplot(df, aes(Time, ID, fill = value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~ variable, ncol = 3, labeller = label_wrap_gen(width = 20), scales = "free") +
  scale_fill_viridis_c(name = "Overlap (%)", option = "C", direction = 1) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  labs(title = "Species × time overlap (folded SFS, median)",
       x = NULL, y = NULL)

ggsave(file.path(outdir, "fig5_species_heatmap_folded.png"),
       p_heat, width = 10, height = 6.8, dpi = 300)

cat("Saved:\n",
    normalizePath(file.path(outdir, "fig5_species_lines_freeY_folded.png")), "\n",
    normalizePath(file.path(outdir, "fig5_species_lines_normalized_folded.png")), "\n",
    normalizePath(file.path(outdir, "fig5_species_heatmap_folded.png")), "\n")
