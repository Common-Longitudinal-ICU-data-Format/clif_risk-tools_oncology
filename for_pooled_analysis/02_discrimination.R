# ==============================================================================
# 02_discrimination.R
# AUROC calculations using two complementary approaches:
#   1. Site-level meta-analysis (proper DeLong CIs at each site, then pooled)
#   2. Pooled weighted AUC (from aggregated counts, Hanley-McNeil CIs)
#
# Round-two change: OUTCOME is now a parameter rather than a fixed value.
# The site pipeline emits every AUROC artifact five times, once per outcome
# (composite, nohospice, wardicu, warddeath, hospicedc). Round one read only the
# composite slice, because legacy_view() filters to outcome_key == "composite"
# before anything downstream sees the data. This script now drives off the
# all-outcome staging objects (auroc_enc_all_raw, auroc_h24_all_raw,
# auroc_h12_all_raw, maxscores_ca_all_raw, counts_h24_all_raw, counts_h12_all_raw)
# and computes every outcome at every metric.
#
# legacy_view() itself is unchanged, and every object this script exports under
# its round-one name is still the composite slice, so 05_subgroups.R and
# 06_figures.R run untouched. The per-outcome results are exported alongside
# under *_by_outcome names.
#
# A reproduction self-check (bottom of the approach-1 block) stops the run if
# the composite / main / encounter-level estimates drift from what the
# composite-only path produces, mirroring .verify_composite_reproduction() in
# code/03a_artifacts.R.
#
# Discharges: R02 (ward death as a standalone outcome at the encounter level and
# at 12 and 24 hours), the discrimination half of R27 (composite excluding
# hospice discharge), and R31 (I-squared and tau-squared surfaced in a table
# rather than computed and discarded).
# ==============================================================================

# setup ------------------------------------------------------------------------

library(WeightedROC)
library(metafor)
library(pROC)
library(flextable)
library(officer)
library(here)

## metric specification (data, not code) ---------------------------------------
# Each row names a site-level AUROC source and its matching count source for the
# weighted approach. horizon is NA for the encounter-level row.

METRICS = data.table(
  metric      = c("Encounter max", "12-hour horizon", "24-hour horizon"),
  auroc_obj   = c("auroc_enc_all_raw", "auroc_h12_all_raw", "auroc_h24_all_raw"),
  count_obj   = c("maxscores_ca_all_raw", "counts_h12_all_raw", "counts_h24_all_raw"),
  count_value = c("max_value", "value", "value")
)

METRIC_LEVELS = METRICS$metric

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Calculate weighted AUC from aggregated count data (no expansion needed)
#' Uses WeightedROC package for memory-efficient computation
#' @param df data.table with columns: value (score), outcome (0/1), n (count)
#' @return list with auc, n_total, n_events
calculate_weighted_auc = function(df) {

  if (nrow(df) == 0) return(list(auc = NA_real_, n_total = 0L, n_events = 0L))

  n_total  = sum(df$n, na.rm = TRUE)
  n_events = sum(df$n[df$outcome == 1], na.rm = TRUE)

  if (n_events == 0 || n_events == n_total) {
    return(list(auc = NA_real_, n_total = n_total, n_events = n_events))
  }

  # WeightedROC expects: predictor scores, true labels (0/1), weights
  # We have aggregated data, so we pass value as score, outcome as label, n as weight
  tp_fp = tryCatch(
    WeightedROC(df$value, df$outcome, df$n),
    error = function(e) NULL
  )

  if (is.null(tp_fp)) {
    return(list(auc = NA_real_, n_total = n_total, n_events = n_events))
  }

  auc_val = WeightedAUC(tp_fp)

  list(auc = auc_val, n_total = n_total, n_events = n_events)
}

#' Calculate Hanley-McNeil standard error for AUC
#' Analytical approximation that doesn't require individual observations
#' @param auc the AUC value
#' @param n_pos number of positive cases (events)
#' @param n_neg number of negative cases (non-events)
#' @return standard error of AUC
hanley_mcneil_se = function(auc, n_pos, n_neg) {

  if (is.na(auc) || n_pos <= 0 || n_neg <= 0) return(NA_real_)

  # Hanley & McNeil (1982) approximation
  q1 = auc / (2 - auc)
  q2 = (2 * auc^2) / (1 + auc)

  var_auc = (auc * (1 - auc) + (n_pos - 1) * (q1 - auc^2) + (n_neg - 1) * (q2 - auc^2)) /
    (n_pos * n_neg)

  sqrt(var_auc)
}

#' Meta-analyze site-level AUROCs using random-effects model
#' @param auroc_dt data.table with site, auroc, auroc_se columns
#' @return data.table with pooled estimate, CI, heterogeneity stats
meta_analyze_aurocs = function(auroc_dt) {

  # Filter to valid estimates
  dt = auroc_dt[!is.na(auroc) & !is.na(auroc_se) & auroc_se > 0]

  if (nrow(dt) < 2) {
    if (nrow(dt) == 1) {
      return(data.table(
        k_sites     = 1L,
        auroc       = dt$auroc[1],
        auroc_se    = dt$auroc_se[1],
        ci_lower    = dt$auroc[1] - 1.96 * dt$auroc_se[1],
        ci_upper    = dt$auroc[1] + 1.96 * dt$auroc_se[1],
        tau2        = NA_real_,
        I2          = NA_real_,
        method      = "single_site"
      ))
    }
    return(data.table(
      k_sites = 0L, auroc = NA_real_, auroc_se = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      tau2 = NA_real_, I2 = NA_real_, method = "insufficient_data"
    ))
  }

  # Logit-transform AUROCs before pooling (avoids boundary issues)
  # Delta method: SE_logit = SE_auroc / (auroc * (1 - auroc))
  dt = copy(dt)
  dt[, logit_auroc := log(auroc / (1 - auroc))]
  dt[, logit_se    := auroc_se / (auroc * (1 - auroc))]

  # Random-effects meta-analysis on logit scale
  ma = tryCatch(
    rma(yi = dt$logit_auroc, sei = dt$logit_se, method = "REML"),
    error = function(e) NULL
  )

  if (is.null(ma)) {
    return(data.table(
      k_sites = nrow(dt), auroc = NA_real_, auroc_se = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      tau2 = NA_real_, I2 = NA_real_, method = "meta_failed"
    ))
  }

  # Back-transform from logit scale
  pooled_auroc = plogis(as.numeric(ma$beta))
  ci_lo        = plogis(ma$ci.lb)
  ci_hi        = plogis(ma$ci.ub)
  # Delta method for back-transformed SE
  pooled_se    = ma$se * pooled_auroc * (1 - pooled_auroc)

  data.table(
    k_sites   = nrow(dt),
    auroc     = pooled_auroc,
    auroc_se  = pooled_se,
    ci_lower  = ci_lo,
    ci_upper  = ci_hi,
    tau2      = ma$tau2,
    I2        = ma$I2,
    method    = "random_effects_meta_logit"
  )
}

#' Run meta_analyze_aurocs() over every combination of the key columns.
#' Non-estimable site rows (fewer than MIN_ESTIMABLE events or non-events, per
#' the site-side rule) are dropped before pooling, exactly as legacy_view() did.
#' @param dt          site-level AUROC table (may carry several outcomes)
#' @param metric_lab  label written into the metric column
#' @param key_cols    grouping columns
#' @return one row per combination
meta_over_groups = function(dt, metric_lab,
                            key_cols = c("score_name", "ca_01", "analysis", "outcome_key")) {

  if (is.null(dt) || nrow(dt) == 0) return(data.table())

  d = as.data.table(dt)

  if ("estimable" %in% names(d)) d = d[estimable == TRUE | is.na(estimable)]
  if (nrow(d) == 0) return(data.table())

  combos = unique(d[, ..key_cols])
  out    = vector("list", nrow(combos))

  for (i in seq_len(nrow(combos))) {

    sub = d[combos[i], on = key_cols, nomatch = NULL]
    ma  = meta_analyze_aurocs(sub)

    for (k in key_cols) set(ma, j = k, value = combos[[k]][i])

    ma[, `:=`(
      metric   = metric_lab,
      n_total  = sum(sub$n_obs,    na.rm = TRUE),
      n_events = sum(sub$n_events, na.rm = TRUE)
    )]

    out[[i]] = ma
  }

  rbindlist(out, use.names = TRUE, fill = TRUE)
}

#' Run the pooled weighted AUC over every combination of the key columns.
#' @param dt          count table (maxscores or horizon counts, any outcome)
#' @param metric_lab  label written into the metric column
#' @param value_col   name of the score-value column ("max_value" or "value")
#' @param key_cols    grouping columns
#' @return one row per combination
weighted_over_groups = function(dt, metric_lab, value_col,
                                key_cols = c("score_name", "ca_01", "analysis", "outcome_key")) {

  if (is.null(dt) || nrow(dt) == 0) return(data.table())

  d = as.data.table(dt)

  if (!(value_col %in% names(d))) {
    warning("weighted_over_groups(): '", value_col, "' not found for ", metric_lab)
    return(data.table())
  }

  combos = unique(d[, ..key_cols])
  out    = vector("list", nrow(combos))

  for (i in seq_len(nrow(combos))) {

    sub = d[combos[i], on = key_cols, nomatch = NULL]

    cells  = sub[, .(value = get(value_col), outcome, n)]
    result = calculate_weighted_auc(cells)

    n_pos = result$n_events
    n_neg = result$n_total - result$n_events
    se_hm = hanley_mcneil_se(result$auc, n_pos, n_neg)

    row = data.table(
      metric   = metric_lab,
      auroc    = result$auc,
      auroc_se = se_hm,
      ci_lower = result$auc - 1.96 * se_hm,
      ci_upper = result$auc + 1.96 * se_hm,
      n_total  = result$n_total,
      n_events = result$n_events,
      method   = "weighted_hanley_mcneil"
    )

    for (k in key_cols) set(row, j = k, value = combos[[k]][i])

    out[[i]] = row
  }

  rbindlist(out, use.names = TRUE, fill = TRUE)
}

#' Fetch a staging object by name, returning an empty table when absent.
get_source = function(nm) {
  if (!exists(nm)) return(data.table())
  obj = get(nm)
  if (is.null(obj)) return(data.table())
  as.data.table(obj)
}

# ==============================================================================
# APPROACH 1: SITE-LEVEL META-ANALYSIS
# Uses DeLong CIs computed at each site, then pools via random-effects
# ==============================================================================

message("\n== Approach 1: Site-level meta-analysis (all outcomes) ==")

meta_list = list()

for (m in seq_len(nrow(METRICS))) {

  metric_lab = METRICS$metric[m]
  src        = get_source(METRICS$auroc_obj[m])

  message("\n  -- ", metric_lab, " --")

  if (nrow(src) == 0) {
    message("  WARNING: ", METRICS$auroc_obj[m], " is empty; skipping")
    next
  }

  res = meta_over_groups(src, metric_lab)

  message("    ", nrow(res), " pooled estimates across ",
          uniqueN(res$outcome_key), " outcome(s) and ",
          uniqueN(res$analysis), " variant(s)")

  meta_list[[metric_lab]] = res
}

auroc_meta_all = rbindlist(meta_list, use.names = TRUE, fill = TRUE)

if (nrow(auroc_meta_all) > 0) {
  auroc_meta_all[, approach := "site_meta"]
}

## reproduction self-check -----------------------------------------------------
# The composite / encounter-level estimates produced by the parameterized path
# must equal those produced by pooling the composite-only legacy object. A
# mismatch means the outcome filter, the estimability filter, or the grouping
# changed, and every downstream table would move with it.

message("\n== Reproduction self-check (composite, encounter level) ==")

if (nrow(auroc_meta_all) > 0 && exists("auroc_enc_raw") && nrow(auroc_enc_raw) > 0) {

  ref_check = meta_over_groups(
    auroc_enc_raw, "Encounter max",
    key_cols = c("score_name", "ca_01", "analysis")
  )

  new_check = auroc_meta_all[
    metric == "Encounter max" & outcome_key == "composite",
    .(score_name, ca_01, analysis, auroc, auroc_se)
  ]

  cmp = merge(
    ref_check[, .(score_name, ca_01, analysis, auroc_ref = auroc, se_ref = auroc_se)],
    new_check,
    by = c("score_name", "ca_01", "analysis"),
    all = TRUE
  )

  bad = cmp[
    is.na(auroc) != is.na(auroc_ref) |
      (!is.na(auroc) & abs(auroc - auroc_ref) > 1e-9) |
      (!is.na(auroc_se) & abs(auroc_se - se_ref) > 1e-9)
  ]

  if (nrow(bad) > 0) {
    print(bad)
    stop("Composite reproduction check failed: the parameterized path does not ",
         "reproduce the composite-only estimates.", call. = FALSE)
  }

  message("  Passed: ", nrow(cmp), " composite encounter-level estimates reproduce exactly.")

} else {
  message("  Skipped: no composite encounter-level estimates available to compare.")
}

# ==============================================================================
# APPROACH 2: POOLED WEIGHTED AUC
# Uses aggregated counts with WeightedROC, Hanley-McNeil CIs
# ==============================================================================

message("\n== Approach 2: Pooled weighted AUC (all outcomes) ==")

weighted_list = list()

for (m in seq_len(nrow(METRICS))) {

  metric_lab = METRICS$metric[m]
  src        = get_source(METRICS$count_obj[m])

  message("\n  -- ", metric_lab, " --")

  if (nrow(src) == 0) {
    message("  WARNING: ", METRICS$count_obj[m], " is empty; skipping")
    next
  }

  # The bootstrap variants are per-iteration replicates, not analyses; they
  # belong to eFigure 7 and must never enter the point-estimate tables.
  if ("analysis" %in% names(src)) {
    src = src[!(analysis %in% c("boot", "bootenc"))]
  }

  res = weighted_over_groups(src, metric_lab, METRICS$count_value[m])

  message("    ", nrow(res), " pooled estimates")

  weighted_list[[metric_lab]] = res
}

auroc_weighted_all = rbindlist(weighted_list, use.names = TRUE, fill = TRUE)

if (nrow(auroc_weighted_all) > 0) {
  auroc_weighted_all[, approach := "pooled_weighted"]
}

# ==============================================================================
# COMPARE APPROACHES
# ==============================================================================

message("\n== Comparing approaches ==")

cmp_keys = c("score_name", "ca_01", "analysis", "metric", "outcome_key")

if (nrow(auroc_meta_all) > 0 && nrow(auroc_weighted_all) > 0) {

  comparison = merge(
    auroc_meta_all[, .(score_name, ca_01, analysis, metric, outcome_key,
                       auroc_meta = auroc, se_meta = auroc_se,
                       ci_lo_meta = ci_lower, ci_hi_meta = ci_upper)],
    auroc_weighted_all[, .(score_name, ca_01, analysis, metric, outcome_key,
                           auroc_weighted = auroc, se_weighted = auroc_se,
                           ci_lo_weighted = ci_lower, ci_hi_weighted = ci_upper)],
    by  = cmp_keys,
    all = TRUE
  )

  comparison[, `:=`(
    diff_auroc = auroc_meta - auroc_weighted,
    diff_se    = se_meta - se_weighted
  )]

  message("\n  Comparison of main analysis, composite outcome (meta vs weighted):")
  print(comparison[
    analysis == "main" & outcome_key == "composite",
    .(metric, score_name, ca_01,
      auroc_meta     = round(auroc_meta, 3),
      auroc_weighted = round(auroc_weighted, 3),
      diff           = round(diff_auroc, 4))
  ])

  auroc_comparison = comparison

} else {
  message("  Cannot compare - one or both approaches have no data")
  auroc_comparison = data.table()
}

# ==============================================================================
# PRIMARY RESULTS (use site-level meta if available, else weighted)
# ==============================================================================

message("\n== Assembling primary results ==")

if (nrow(auroc_meta_all) > 0) {
  auroc_primary = copy(auroc_meta_all)
  message("  Using site-level meta-analysis as primary")
} else {
  auroc_primary = copy(auroc_weighted_all)
  message("  Using pooled weighted AUC as primary (no site-level data)")
}

auroc_primary[, `:=`(
  ca_lab       = fifelse(ca_01 == 1L, "Cancer", "Non-cancer"),
  analysis_lab = factor(analysis,    levels = names(ANALYSIS_LABS), labels = ANALYSIS_LABS),
  score_lab    = factor(score_name,  levels = names(SCORE_LABS),    labels = SCORE_LABS),
  outcome_lab  = factor(outcome_key, levels = names(OUTCOME_LABS),  labels = OUTCOME_LABS),
  metric_lab   = factor(metric,      levels = METRIC_LEVELS)
)]

message("  ", nrow(auroc_primary), " estimates across ",
        uniqueN(auroc_primary$outcome_key), " outcomes, ",
        uniqueN(auroc_primary$metric), " metrics, ",
        uniqueN(auroc_primary$analysis), " variants")

# ==============================================================================
# AUROC DIFFERENCES (Cancer vs Non-cancer)
# ==============================================================================

message("\n== Computing AUROC differences ==")

auroc_diff = auroc_primary[, {

  auc_ca    = auroc[ca_01 == 1]
  auc_nonca = auroc[ca_01 == 0]
  se_ca     = auroc_se[ca_01 == 1]
  se_nonca  = auroc_se[ca_01 == 0]

  if (length(auc_ca) == 1 && length(auc_nonca) == 1 &&
      !is.na(auc_ca) && !is.na(auc_nonca)) {

    diff_val = auc_ca - auc_nonca

    # SE of difference (assuming independence between groups)
    se_diff = sqrt(se_ca^2 + se_nonca^2)

    z_val = diff_val / se_diff
    p_val = 2 * pnorm(-abs(z_val))

    list(
      auroc_ca    = auc_ca,
      auroc_nonca = auc_nonca,
      diff_auc    = diff_val,
      se_diff     = se_diff,
      p_value     = p_val
    )
  } else {
    list(
      auroc_ca    = NA_real_,
      auroc_nonca = NA_real_,
      diff_auc    = NA_real_,
      se_diff     = NA_real_,
      p_value     = NA_real_
    )
  }
}, by = .(score_name, analysis, metric, outcome_key)]

auroc_diff[, `:=`(
  sig = fcase(
    is.na(p_value),  "",
    p_value < 0.001, "***",
    p_value < 0.01,  "**",
    p_value < 0.05,  "*",
    default = ""
  ),
  score_lab    = factor(score_name,  levels = names(SCORE_LABS),    labels = SCORE_LABS),
  analysis_lab = factor(analysis,    levels = names(ANALYSIS_LABS), labels = ANALYSIS_LABS),
  outcome_lab  = factor(outcome_key, levels = names(OUTCOME_LABS),  labels = OUTCOME_LABS),
  metric_lab   = factor(metric,      levels = METRIC_LEVELS)
)]

# ==============================================================================
# SUMMARY TABLES
# ==============================================================================

message("\n== Creating summary tables ==")

table_path = function(stem) {
  if (exists("today")) {
    here("output", "tables", paste0(stem, "_", today, ".docx"))
  } else {
    here("output", "tables", paste0(stem, ".docx"))
  }
}

## composite tables (round-one shape, unchanged output) ------------------------

auroc_main_table = auroc_primary[analysis == "main" & outcome_key == "composite", .(
  metric_lab, score_lab, ca_lab,
  auroc_fmt = sprintf("%.3f (%.3f-%.3f)", auroc, ci_lower, ci_upper),
  n         = format_n(n_total),
  n_events  = format_n(n_events)
)]

auroc_main_wide = dcast(
  auroc_main_table,
  metric_lab + score_lab ~ ca_lab,
  value.var = "auroc_fmt"
)

## Difference table ------------------------------------------------------------
##
## Reports the cancer vs non-cancer AUROC gap with a continuous p-value, no
## significance threshold applied. See eMethods for rationale: the random-
## effects meta-analysis SE is dominated by between-site heterogeneity, and a
## dichotomous significance threshold would under-credit the directional
## consistency across the five scores.

auroc_diff_table = auroc_diff[analysis == "main" & outcome_key == "composite", .(
  metric_lab, score_lab,
  diff_fmt = sprintf("%.3f", diff_auc),
  p_fmt    = fifelse(!is.na(p_value) & p_value < 0.001, "<0.001",
                     fifelse(!is.na(p_value), sprintf("%.3f", p_value), "—"))
)]

## Publishable Table: AUROC by cohort with difference and p-value --------------
##
## Main-text table combining per-cohort AUROCs (with 95% CIs) and the cancer
## vs non-cancer difference with p-value. Restricted to encounter-level
## ("Encounter max"); the 24-hour horizon row is shown in eFigure 7 alongside
## the clustered-bootstrap robustness analysis.

auroc_pub_table = merge(
  auroc_main_wide[metric_lab == "Encounter max"],
  auroc_diff_table[metric_lab == "Encounter max"],
  by    = c("metric_lab", "score_lab"),
  all.x = TRUE
)

auroc_pub_table[, metric_lab := NULL]

setcolorder(auroc_pub_table,
            c("score_lab", "Non-cancer", "Cancer", "diff_fmt", "p_fmt"))

setnames(auroc_pub_table,
         old = c("score_lab", "diff_fmt",            "p_fmt"),
         new = c("Score",     "Cancer - Non-cancer", "P value"))

message("  Publishable AUROC table: ", nrow(auroc_pub_table), " rows")

ft_auroc = flextable(auroc_pub_table) |>
  autofit() |>
  align(j = 2:5, align = "center", part = "all") |>
  bold(part = "header")

auroc_path = table_path("table_auroc")
save_as_docx(ft_auroc, path = auroc_path)
message("  Saved: ", auroc_path)

## per-outcome table (R02, R27) ------------------------------------------------
# One block per outcome; within each, one row per score at each metric, with
# both cohorts, the difference, and the p-value. Ward death is the row R02
# promises; the composite excluding hospice discharge is the row R27 turns on.

auroc_outcome_wide = dcast(
  auroc_primary[analysis == "main", .(
    outcome_lab, metric_lab, score_lab, ca_lab,
    auroc_fmt = sprintf("%.3f (%.3f-%.3f)", auroc, ci_lower, ci_upper)
  )],
  outcome_lab + metric_lab + score_lab ~ ca_lab,
  value.var = "auroc_fmt"
)

auroc_outcome_diff = auroc_diff[analysis == "main", .(
  outcome_lab, metric_lab, score_lab,
  diff_fmt = sprintf("%.3f", diff_auc),
  p_fmt    = fifelse(!is.na(p_value) & p_value < 0.001, "<0.001",
                     fifelse(!is.na(p_value), sprintf("%.3f", p_value), "—")),
  n_events = format_n(0L)
)]

auroc_outcome_events = auroc_primary[
  analysis == "main" & ca_01 == 1L,
  .(outcome_lab, metric_lab, score_lab, ca_events = format_n(n_events))
]

auroc_by_outcome = merge(
  auroc_outcome_wide,
  auroc_outcome_diff[, .(outcome_lab, metric_lab, score_lab, diff_fmt, p_fmt)],
  by    = c("outcome_lab", "metric_lab", "score_lab"),
  all.x = TRUE
)

auroc_by_outcome = merge(
  auroc_by_outcome, auroc_outcome_events,
  by    = c("outcome_lab", "metric_lab", "score_lab"),
  all.x = TRUE
)

setorder(auroc_by_outcome, outcome_lab, metric_lab, score_lab)

setcolorder(auroc_by_outcome,
            c("outcome_lab", "metric_lab", "score_lab",
              "Non-cancer", "Cancer", "diff_fmt", "p_fmt", "ca_events"))

setnames(auroc_by_outcome,
         old = c("outcome_lab", "metric_lab", "score_lab",
                 "diff_fmt", "p_fmt", "ca_events"),
         new = c("Outcome", "Metric", "Score",
                 "Cancer - Non-cancer", "P value", "Cancer events"))

ft_by_outcome = flextable(auroc_by_outcome) |>
  autofit() |>
  bold(part = "header")

by_outcome_path = table_path("table_auroc_by_outcome")
save_as_docx(ft_by_outcome, path = by_outcome_path)
message("  Saved: ", by_outcome_path, " (", nrow(auroc_by_outcome), " rows)")

## heterogeneity table (R31) ---------------------------------------------------
# I-squared and tau-squared were already computed inside meta_analyze_aurocs()
# and then discarded. R31 promises both as summary measures of between-site
# heterogeneity, so they are surfaced here for the composite outcome at every
# metric. tau-squared is on the logit scale, which the footnote must state.

het_table = auroc_primary[analysis == "main" & outcome_key == "composite", .(
  metric_lab, score_lab, ca_lab,
  k_sites,
  i2_fmt   = fifelse(is.na(I2),   "—", sprintf("%.1f%%", I2)),
  tau2_fmt = fifelse(is.na(tau2), "—", sprintf("%.4f", tau2))
)]

het_wide = dcast(
  het_table,
  metric_lab + score_lab ~ ca_lab,
  value.var = c("k_sites", "i2_fmt", "tau2_fmt")
)

setorder(het_wide, metric_lab, score_lab)

ft_het = flextable(het_wide) |>
  autofit() |>
  bold(part = "header")

het_path = table_path("table_heterogeneity")
save_as_docx(ft_het, path = het_path)
message("  Saved: ", het_path)

# ==============================================================================
# EXPORTS
# ==============================================================================

message("\n== Discrimination analysis complete ==")

if (nrow(auroc_meta_all) > 0) {
  message("  Site-level meta-analysis estimates: ", nrow(auroc_meta_all))
}
if (nrow(auroc_weighted_all) > 0) {
  message("  Pooled weighted estimates: ", nrow(auroc_weighted_all))
}
message("  Primary results (all outcomes): ", nrow(auroc_primary))

# Add proper encounter Ns from VARIANT_N (not n_total which may be observations)
# This ensures figures show correct encounter counts
if (exists("VARIANT_N") && exists("COHORT_N")) {

  message("\n== Adding canonical encounter Ns ==")

  auroc_primary[analysis == "main" & ca_01 == 0, n_encounters := COHORT_N["0"]]
  auroc_primary[analysis == "main" & ca_01 == 1, n_encounters := COHORT_N["1"]]

  for (v in names(VARIANT_N)) {
    if (v != "main") {
      auroc_primary[analysis == v, n_encounters := VARIANT_N[v]]
    }
  }

  message("  Added n_encounters column to auroc_primary")
}

## composite slice under the round-one names ----------------------------------
# 05_subgroups.R and 06_figures.R consume these and assume one outcome per row.
# Handing them the full per-outcome table would silently multiply every row, so
# the legacy names stay composite-only and the per-outcome results travel
# separately.

auroc_primary_composite = auroc_primary[outcome_key == "composite"]
auroc_diff_composite    = auroc_diff[outcome_key == "composite"]

auroc_results_final    = auroc_primary_composite
auroc_meta_final       = auroc_meta_all[outcome_key == "composite"]
auroc_weighted_final   = auroc_weighted_all[outcome_key == "composite"]
auroc_comparison_final = auroc_comparison
auroc_diff_final       = auroc_diff_composite
auroc_main_table_final = auroc_main_wide
auroc_diff_table_final = auroc_diff_table
auroc_pub_table_final  = auroc_pub_table

## per-outcome objects ---------------------------------------------------------

auroc_primary_by_outcome = auroc_primary
auroc_diff_by_outcome    = auroc_diff
auroc_by_outcome_final   = auroc_by_outcome
het_table_final          = het_wide

# Also export site-level data for forest plots
if (exists("auroc_enc_raw")) auroc_site_enc = auroc_enc_raw
if (exists("auroc_h24_raw")) auroc_site_h24 = auroc_h24_raw
if (exists("auroc_h12_raw")) auroc_site_h12 = auroc_h12_raw

gc()
