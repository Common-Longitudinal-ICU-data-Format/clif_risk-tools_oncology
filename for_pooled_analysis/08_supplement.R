# ==============================================================================
# 08_supplement.R
# Pooled analyses for artifact families that 00_load.R reads and that no other
# script consumed. Each block discharges a specific reviewer response.
#
#   S1  first-ward-observation score distributions          R06
#   S2  eTable 1: health systems, hospitals, types          R08, R18
#   S3  admission diagnoses by chapter and code stem        R09
#   S4  number needed to evaluate, lead time, crossing      R25
#   S5  monitoring intensity and component missingness      R30
#   S6  between-site descriptive heterogeneity              R31
#
# Every block guards on its source object, reports what it found, and writes a
# Word table to output/tables. A missing family produces a message and an empty
# object rather than an error, so a partial site collection still runs.
# ==============================================================================

# setup ------------------------------------------------------------------------

library(flextable)
library(officer)
library(here)

message("\n== 08_supplement.R ==")

sup_path = function(stem) {
  if (exists("today")) {
    here("output", "tables", paste0(stem, "_", today, ".docx"))
  } else {
    here("output", "tables", paste0(stem, ".docx"))
  }
}

save_supplement = function(dt, stem, label, footnote = NULL) {

  if (is.null(dt) || nrow(dt) == 0) {
    message("  ", label, ": nothing to write")
    return(invisible(NULL))
  }

  ft = flextable(as.data.frame(dt)) |>
    bold(part = "header")

  if (!is.null(footnote)) {
    ft = add_footer_lines(ft, values = footnote)
  }

  ft = autofit(ft)

  p = sup_path(stem)
  save_as_docx(ft, path = p)
  message("  Saved: ", p, " (", nrow(dt), " rows)")
  invisible(p)
}

## NEWS positivity footnote ----------------------------------------------------
# Every table built from the sesp / crossclass artifacts applies the full NEWS
# rule (aggregate >= 5 OR any single parameter scoring 3). The threshold-
# equivalence tables in 04_meta.R sweep the aggregate score alone, so NEWS
# operating characteristics at threshold 5 differ between the two families.
# The same sentence is attached to both so the reader sees the reason.

NEWS_RULE_FOOTNOTE = paste0(
  "NEWS positivity applies the full published rule: aggregate score \u22655 or ",
  "any single parameter scoring 3. The threshold-equivalence tables sweep the ",
  "aggregate score alone and therefore report different NEWS operating ",
  "characteristics at threshold 5."
)

have = function(nm) exists(nm) && nrow(get(nm)) > 0

## weighted summary helpers ----------------------------------------------------
# The site artifacts are count tables, so every distributional summary has to be
# computed from (value, n) pairs rather than from observations.

weighted_mean_from_counts = function(value, n) {
  sum(as.numeric(value) * n) / sum(n)
}

weighted_sd_from_counts = function(value, n) {
  m  = weighted_mean_from_counts(value, n)
  nn = sum(n)
  if (nn <= 1) return(NA_real_)
  sqrt(sum(n * (as.numeric(value) - m)^2) / (nn - 1))
}

weighted_quantile_from_counts = function(value, n, probs) {
  o   = order(value)
  v   = as.numeric(value)[o]
  w   = as.numeric(n)[o]
  cum = cumsum(w) / sum(w)
  vapply(probs, function(p) v[which(cum >= p)[1]], numeric(1))
}

fmt_pct = function(x, digits = 1) {
  fifelse(is.na(x), "—", sprintf(paste0("%.", digits, "f%%"), 100 * x))
}

#' Drop rows from a site that still ships the round-one schema.
#' load_from_catalog() rbinds with fill = TRUE, so a site whose artifacts predate
#' the round-two rewrite contributes rows with NA in every new column rather than
#' failing. Those rows would surface as a phantom "<= NA h" column and a blank
#' crossing definition. Detect, name the sites, and drop.
drop_round_one = function(dt, required, label) {

  d = as.data.table(dt)

  missing_cols = setdiff(required, names(d))
  if (length(missing_cols) > 0L) {
    message("  ", label, ": missing column(s) ",
            paste(missing_cols, collapse = ", "), "; skipping")
    return(data.table())
  }

  stale = Reduce(`|`, lapply(required, function(cc) is.na(d[[cc]])))

  if (any(stale)) {
    stale_sites = sort(unique(d$site[stale]))
    message("  ", label, ": dropped ", sum(stale), " round-one row(s) from ",
            paste(stale_sites, collapse = ", "),
            " — these site(s) must re-run 03b_leadtime.R")
  }

  d[!stale]
}

# ==============================================================================
# S1. FIRST-WARD-OBSERVATION SCORE DISTRIBUTIONS  (R06)
# ==============================================================================
# Reviewer 1 asked for a subanalysis by emergency-department triage category.
# CLIF carries no triage acuity, so the response offers the score at the first
# ward observation as a measure of presenting physiology at the moment the
# monitoring period begins. firstscore_raw holds (score_name, ca_01, value, n).

message("\n-- S1: first-ward-observation score distributions (R06) --")

firstscore_summary = data.table()

if (have("firstscore_raw")) {

  fs = as.data.table(firstscore_raw)[!is.na(value)]

  fs_pooled = fs[, .(n = sum(n)), by = .(score_name, ca_01, value)]

  firstscore_summary = fs_pooled[, {
    q = weighted_quantile_from_counts(value, n, c(0.25, 0.50, 0.75))
    .(
      n_encounters = sum(n),
      mean_score   = weighted_mean_from_counts(value, n),
      sd_score     = weighted_sd_from_counts(value, n),
      median_score = q[2],
      q25          = q[1],
      q75          = q[3]
    )
  }, by = .(score_name, ca_01)]

  firstscore_summary[, `:=`(
    score_lab = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS),
    ca_lab    = fifelse(ca_01 == 1L, "Cancer", "Non-cancer")
  )]

  setorder(firstscore_summary, score_lab, -ca_01)

  firstscore_table = firstscore_summary[, .(
    Score      = score_lab,
    Cohort     = ca_lab,
    N          = format_n(n_encounters),
    `Mean (SD)`     = sprintf("%.2f (%.2f)", mean_score, sd_score),
    `Median (IQR)`  = sprintf("%.0f (%.0f-%.0f)", median_score, q25, q75)
  )]

  save_supplement(firstscore_table, "table_s_firstscore", "First-ward score")

} else {
  message("  firstscore_raw unavailable; skipping")
  firstscore_table = data.table()
}

# ==============================================================================
# S2. eTABLE 1: HEALTH SYSTEMS, HOSPITALS, AND HOSPITAL TYPES  (R08, R18)
# ==============================================================================
# eTable 1 columns, in order: Health system, Study period, Hospitals (n),
# hospital-type composition, Hospitalizations (n), Adult inpatient beds, Adult
# ICU beds, IRB. Three sources feed it.
#
# 1. site_metadata.yaml (repository root), hand-maintained. One key per site
#    token with:
#      health_system         display name
#      hospital_ids          list of literal ADT hospital_id values, lowercase
#      n_hospitals           OR a bare count, for a site whose ids are not yet
#                            known to the coordinating center
#      date_range, adult_inpatient_beds, adult_icu_beds, irb
#    hospital_ids mirrors the site's config `hospital_ids` key, which
#    01_cohort.R enforces with a stop. When a list is present it defines the
#    hospital count; when only n_hospitals is present, that number is reported;
#    when neither is present the observed count is reported.
#
# 2. leadtime_median_raw, unit == "hospital": one unit_id per hospital that
#    contributed at least one primary-outcome event. This is the observed set,
#    checked against the declared list or count. Undeclared observed ids are
#    excluded and named on the console; declared-but-unobserved ids are named.
#
# 3. hospital_types_raw: counts by ADT hospital_type (academic / community /
#    LTACH / unspecified). Round-two 01_cohort.R restricts this tally to the
#    declared list after cohort pruning, so for new uploads its total equals
#    the declared count. Uploads predating that change (hopkins) can carry a
#    larger total; the console names the mismatch and the table reports the
#    declared count.
#
# Hospitalizations (n) is the final flow row per site (flow_data_raw), i.e. the
# analytic ED-admit cohort, matching the submitted eTable 1 definition.

message("\n-- S2: eTable 1 (R08, R18) --")

hospital_type_table = data.table()
total_hospitals     = NA_integer_

## 1. site metadata --------------------------------------------------------------

meta_file = here("site_metadata.yaml")
site_meta = data.table()
declared  = data.table()

scalar_or_na = function(x) {
  if (is.null(x) || length(x) == 0L) NA_character_ else as.character(x[[1]])
}

if (file.exists(meta_file)) {

  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("site_metadata.yaml is present but the yaml package is not installed.",
         call. = FALSE)
  }

  meta_list = yaml::read_yaml(meta_file)

  if (!is.list(meta_list) || is.null(names(meta_list)) || any(!nzchar(names(meta_list)))) {
    stop("site_metadata.yaml must be a mapping keyed by site token.", call. = FALSE)
  }

  site_meta = rbindlist(lapply(names(meta_list), function(s) {
    m   = meta_list[[s]]
    ids = m$hospital_ids
    data.table(
      site                 = tolower(trimws(s)),
      health_system        = scalar_or_na(m$health_system),
      n_declared_list      = if (is.null(ids)) NA_integer_ else length(unlist(ids)),
      n_declared_count     = if (is.null(m$n_hospitals)) NA_integer_ else as.integer(m$n_hospitals),
      date_range           = scalar_or_na(m$date_range),
      adult_inpatient_beds = scalar_or_na(m$adult_inpatient_beds),
      adult_icu_beds       = scalar_or_na(m$adult_icu_beds),
      irb                  = scalar_or_na(m$irb)
    )
  }))

  declared = rbindlist(lapply(names(meta_list), function(s) {
    ids = meta_list[[s]]$hospital_ids
    if (is.null(ids) || length(ids) == 0L) return(NULL)
    data.table(
      site        = tolower(trimws(s)),
      hospital_id = tolower(trimws(as.character(unlist(ids))))
    )
  }))

  if (nrow(declared) > 0) {
    declared = declared[nzchar(hospital_id)]
    dup = declared[, .N, by = .(site, hospital_id)][N > 1]
    if (nrow(dup) > 0) {
      stop("site_metadata.yaml lists duplicate hospital_ids for: ",
           paste(unique(dup$site), collapse = ", "), call. = FALSE)
    }
  }

  both = site_meta[!is.na(n_declared_list) & !is.na(n_declared_count) &
                     n_declared_list != n_declared_count]
  if (nrow(both) > 0) {
    stop("site_metadata.yaml gives both hospital_ids and a different n_hospitals for: ",
         paste(both$site, collapse = ", "), call. = FALSE)
  }

  message("  Read site_metadata.yaml: ", nrow(site_meta), " site(s); id lists for ",
          length(unique(declared$site)), ", counts only for ",
          site_meta[is.na(n_declared_list) & !is.na(n_declared_count), .N])

} else {
  message("  site_metadata.yaml not found at ", meta_file,
          "; health-system names, declared hospitals, bed counts, IRB numbers, ",
          "and date ranges left for manual entry")
}

## 2. observed hospitals and the reconciled count -------------------------------

hospital_counts = data.table()

if (have("leadtime_median_raw")) {

  lm_units = drop_round_one(leadtime_median_raw,
                            c("crossing_def", "unit", "unit_id", "n_events", "median_h"),
                            "leadtime_median (hospital count)")

  lm_units = lm_units[unit == "hospital"]
  lm_units[, unit_id := tolower(trimws(as.character(unit_id)))]

  observed = unique(lm_units[, .(site, hospital_id = unit_id)])

  sites_seen = sort(unique(c(observed$site, site_meta$site)))
  rows = vector("list", length(sites_seen))

  for (i in seq_along(sites_seen)) {

    s    = sites_seen[i]
    obs  = sort(observed[site == s, hospital_id])
    dec  = sort(declared[site == s, hospital_id])
    n_ct = if (nrow(site_meta) > 0 && s %in% site_meta$site) {
      site_meta[site == s, n_declared_count]
    } else NA_integer_

    if (length(obs) == 0L) {
      message("  ", s, ": in site_metadata.yaml but no artifacts loaded; skipped")
      next
    }

    if (length(dec) > 0L) {

      undeclared = setdiff(obs, dec)
      unobserved = setdiff(dec, obs)
      if (length(undeclared) > 0L) {
        message("  ", s, ": observed but not declared, EXCLUDED from count: ",
                paste(undeclared, collapse = ", "))
      }
      if (length(unobserved) > 0L) {
        message("  ", s, ": declared but no primary-outcome events in lead-time artifact: ",
                paste(unobserved, collapse = ", "))
      }
      n_rep = length(dec)
      src   = "declared list"
      message("  ", s, ": ", n_rep, " hospital(s) from declared list [",
              paste(dec, collapse = ", "), "]")

    } else if (!is.na(n_ct)) {

      if (length(obs) != n_ct) {
        message("  ", s, ": declared n_hospitals = ", n_ct, " but ", length(obs),
                " observed [", paste(obs, collapse = ", "),
                "]; reporting the declared count. Obtain the id list from the site.")
      }
      n_rep = n_ct
      src   = "declared count"
      message("  ", s, ": ", n_rep, " hospital(s) from declared count (observed ",
              length(obs), ")")

    } else {

      n_rep = length(obs)
      src   = "observed"
      message("  ", s, ": ", n_rep, " observed hospital(s), nothing declared [",
              paste(obs, collapse = ", "), "]")
    }

    rows[[i]] = data.table(
      site            = s,
      hospitals_total = as.integer(n_rep),
      hospitals_obs   = length(obs),
      count_source    = src
    )
  }

  hospital_counts = rbindlist(rows)

  total_hospitals = sum(hospital_counts$hospitals_total)
  message("  Pooled hospital count: ", total_hospitals,
          " across ", nrow(hospital_counts), " health system(s)")

} else {
  message("  leadtime_median_raw unavailable; hospital count not computed")
}

## 3. hospital-type composition from ADT ---------------------------------------

ht_wide   = data.table()
type_cols = character(0)

if (have("hospital_types_raw")) {

  ht = as.data.table(hospital_types_raw)[, .(n_hospitals = sum(n_hospitals)),
                                         by = .(site, hospital_type)]

  ht_wide = dcast(ht, site ~ hospital_type, value.var = "n_hospitals", fill = 0L)

  type_cols = setdiff(names(ht_wide), "site")
  ht_wide[, hospitals_adt := rowSums(.SD), .SDcols = type_cols]

  message("  ADT hospital_types total: ", sum(ht_wide$hospitals_adt), "; by type: ",
          paste(type_cols, "=", colSums(ht_wide[, ..type_cols]), collapse = ", "))

} else {
  message("  hospital_types_raw unavailable; type composition not available")
}

## 4. hospitalizations (final analytic cohort per site) --------------------------

site_hosp_n = data.table()

if (exists("flow_data_raw") && nrow(flow_data_raw) > 0 && exists("FLOW_FINAL_STEP")) {
  site_hosp_n = as.data.table(flow_data_raw)[step %like% FLOW_FINAL_STEP, .(
    site,
    hospitalizations = n_remaining_ca + n_remaining_no
  )]
} else {
  message("  flow_data_raw / FLOW_FINAL_STEP unavailable; hospitalizations column blank")
}

## 5. assemble -----------------------------------------------------------------

if (nrow(hospital_counts) > 0 || nrow(ht_wide) > 0) {

  et1 = if (nrow(hospital_counts) > 0) copy(hospital_counts) else data.table(site = ht_wide$site)

  if (nrow(ht_wide) > 0)     et1 = merge(et1, ht_wide,     by = "site", all.x = TRUE)
  if (nrow(site_hosp_n) > 0) et1 = merge(et1, site_hosp_n, by = "site", all.x = TRUE)
  if (nrow(site_meta) > 0)   et1 = merge(et1, site_meta,   by = "site", all.x = TRUE)

  for (mc in c("hospitals_total", "hospitals_adt", "hospitalizations")) {
    if (!mc %in% names(et1)) et1[, (mc) := NA_integer_]
  }
  for (mc in c("health_system", "date_range", "adult_inpatient_beds", "adult_icu_beds", "irb")) {
    if (!mc %in% names(et1)) et1[, (mc) := NA_character_]
  }

  ## QC: ADT type total versus reported count --------------------------------------
  mismatch = et1[!is.na(hospitals_adt) & !is.na(hospitals_total) &
                   hospitals_adt != hospitals_total]
  if (nrow(mismatch) > 0) {
    warning("hospital_types total differs from the reported hospital count at: ",
            paste0(mismatch$site, " (types ", mismatch$hospitals_adt, " vs reported ",
                   mismatch$hospitals_total, ")", collapse = "; "),
            ". Expected only for uploads predating the round-two 01_cohort.R tally; ",
            "otherwise the site's declared list or ADT needs correction.",
            call. = FALSE)
  }

  et1[, health_system := fifelse(is.na(health_system), toupper(site), health_system)]
  setorder(et1, health_system)

  keep_cols = c("health_system", "date_range", "hospitals_total", type_cols,
                "hospitalizations", "adult_inpatient_beds", "adult_icu_beds", "irb")
  keep_cols = intersect(keep_cols, names(et1))

  hospital_type_table = et1[, ..keep_cols]

  if ("hospitalizations" %in% names(hospital_type_table)) {
    hospital_type_table[, hospitalizations :=
                          fifelse(is.na(hospitalizations), "",
                                  format_n(hospitalizations))]
  }
  for (cc in intersect(c("adult_inpatient_beds", "adult_icu_beds"), names(hospital_type_table))) {
    v = suppressWarnings(as.numeric(hospital_type_table[[cc]]))
    hospital_type_table[, (cc) := fifelse(is.na(v), "", format_n(v))]
  }

  pretty = c(
    health_system        = "Health system",
    date_range           = "Study period",
    hospitals_total      = "Hospitals, n",
    academic             = "Academic",
    community            = "Community",
    LTACH                = "Long-term acute care",
    unspecified          = "Type unspecified",
    hospitalizations     = "Hospitalizations, n",
    adult_inpatient_beds = "Adult inpatient beds, n",
    adult_icu_beds       = "Adult ICU beds, n",
    irb                  = "IRB"
  )
  old_names = names(pretty)[names(pretty) %in% names(hospital_type_table)]
  setnames(hospital_type_table, old = old_names, new = unname(pretty[old_names]))

  save_supplement(
    hospital_type_table, "table_s_hospital_types", "eTable 1",
    footnote = paste0(
      "Hospitals, n counts hospitals contributing encounters to the analytic cohort. ",
      "Hospital type is the majority-vote classification of each hospital in the ",
      "site's admission-discharge-transfer table. Hospitalizations, n is the final ",
      "analytic cohort (emergency-department admissions to a ward) at each health system."
    )
  )

} else {
  message("  Neither hospital counts nor hospital_types available; skipping eTable 1")
}

# ==============================================================================
# S3. ADMISSION DIAGNOSES  (R09)
# ==============================================================================
# Reviewer 2 asked why patients were hospitalized and whether the reasons differ
# by cancer status. CLIF carries no harmonized chief complaint, so the response
# offers primary hospital diagnoses by ICD-10-CM chapter and by the most
# frequent three-character stems. Site tallies below six encounters are already
# suppressed at export (02_scores.R), so pooled percentages are computed over
# the retained rows and the denominator is stated.

message("\n-- S3: admission diagnoses (R09) --")

adm_chapter_table = data.table()
adm_stem_table    = data.table()

if (have("adm_dx_chapter_raw")) {

  ch = as.data.table(adm_dx_chapter_raw)[, .(n = sum(n)), by = .(ca_01, chapter)]
  ch[, pct := n / sum(n), by = ca_01]

  # ICD-10-CM chapter labels keyed on the leading letter of the code. Chapter D
  # straddles neoplasms (D00-D49) and blood/immune disease (D50-D89); the label
  # names both. Letters not in the map print as the bare letter.
  ICD_CHAPTER_LABS = c(
    A = "A: Infectious and parasitic diseases",
    B = "B: Infectious and parasitic diseases",
    C = "C: Malignant neoplasms",
    D = "D: In situ, benign, or uncertain neoplasms; blood and immune disorders",
    E = "E: Endocrine, nutritional, and metabolic diseases",
    F = "F: Mental and behavioral disorders",
    G = "G: Diseases of the nervous system",
    H = "H: Diseases of the eye and ear",
    I = "I: Diseases of the circulatory system",
    J = "J: Diseases of the respiratory system",
    K = "K: Diseases of the digestive system",
    L = "L: Diseases of the skin",
    M = "M: Musculoskeletal and connective tissue diseases",
    N = "N: Diseases of the genitourinary system",
    O = "O: Pregnancy, childbirth, and the puerperium",
    P = "P: Perinatal conditions",
    Q = "Q: Congenital malformations",
    R = "R: Symptoms, signs, and abnormal findings",
    S = "S: Injury",
    T = "T: Injury, poisoning, and complications of care",
    U = "U: Codes for special purposes (COVID-19)",
    V = "V: External causes of morbidity",
    W = "W: External causes of morbidity",
    X = "X: External causes of morbidity",
    Y = "Y: External causes of morbidity",
    Z = "Z: Factors influencing health status"
  )

  ch[, chapter_lab := fcoalesce(unname(ICD_CHAPTER_LABS[as.character(chapter)]),
                                as.character(chapter))]

  adm_chapter_table = dcast(
    ch[, .(chapter, chapter_lab, ca_01,
           cell = paste0(format_n(n), " (", fmt_pct(pct), ")"))],
    chapter + chapter_lab ~ ca_01,
    value.var = "cell",
    fill = "—"
  )

  setorder(adm_chapter_table, chapter)
  adm_chapter_table[, chapter := NULL]
  setnames(adm_chapter_table, "chapter_lab", "chapter")

  setnames(adm_chapter_table, c("0", "1"), c("No Cancer", "Cancer"), skip_absent = TRUE)
  setnames(adm_chapter_table, "chapter", "ICD-10-CM chapter")

  save_supplement(adm_chapter_table, "table_s_admission_dx_chapter", "Admission dx chapter")

} else {
  message("  adm_dx_chapter_raw unavailable; skipping chapter table")
}

if (have("adm_dx_stem_raw")) {

  st = as.data.table(adm_dx_stem_raw)[, .(n = sum(n)), by = .(ca_01, code_stem)]
  st[, pct := n / sum(n), by = ca_01]

  # top 15 stems within each cohort, unioned so both columns are comparable
  top_stems = unique(c(
    st[ca_01 == 1][order(-n)][1:min(15L, .N)]$code_stem,
    st[ca_01 == 0][order(-n)][1:min(15L, .N)]$code_stem
  ))
  top_stems = top_stems[!is.na(top_stems)]

  adm_stem_table = dcast(
    st[code_stem %in% top_stems,
       .(code_stem, ca_01, cell = paste0(format_n(n), " (", fmt_pct(pct), ")"))],
    code_stem ~ ca_01,
    value.var = "cell",
    fill = "—"
  )

  setnames(adm_stem_table, c("0", "1"), c("No Cancer", "Cancer"), skip_absent = TRUE)
  setnames(adm_stem_table, "code_stem", "ICD-10-CM code")

  save_supplement(adm_stem_table, "table_s_admission_dx_stem", "Admission dx stem")

} else {
  message("  adm_dx_stem_raw unavailable; skipping stem table")
}

# ==============================================================================
# S4. NUMBER NEEDED TO EVALUATE, LEAD TIME, AND CROSSING  (R25)
# ==============================================================================
# Reviewer 4 asked for number needed to evaluate and for time from threshold
# crossing to deterioration. Number needed to evaluate is 1/PPV at each score's
# standard threshold and derives from the sesp cells already pooled.
#
# Lead time arrives in two artifacts. leadtime holds cumulative counts at fixed
# thresholds, which pool by summation and give the encounter-weighted "within x
# hours" proportions directly. leadtime_median holds exact medians computed at
# each site for the health system and for each hospital_id; those are NOT
# pooled, and the range across units is the reported statistic. Both carry
# crossing_def: `first` is the conventional first-alert-of-the-encounter
# definition, and `final_onset` is the onset of the last uninterrupted positive
# run, which is the warning time a clinician means. crossclass gives the 2x2 of
# ever-crossed against ever-had-the-event.
#
# No mean is reported. The distribution is bounded below at zero with a heavy
# right tail, and the round-one mean of 182 hours with an SD of 255 described it
# badly. The sufficient statistics are still exported, so a mean can be
# recovered if a reviewer asks.

message("\n-- S4: NNE, lead time, threshold crossing (R25) --")

nne_table = data.table()

if (have("sesp_raw")) {

  sesp_pool = as.data.table(sesp_raw)[, .(
    n_total = sum(n_total),
    tp      = sum(tp),
    fp      = sum(fp),
    tn      = sum(tn),
    fn      = sum(fn),
    n_pos   = sum(n_pos)
  ), by = .(score_name, ca_01)]

  sesp_pool[, `:=`(
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    ppv         = tp / (tp + fp),
    positivity  = n_pos / n_total
  )]

  # Number needed to evaluate: alerts that must be worked up per true event.
  sesp_pool[, nne := fifelse(ppv > 0, 1 / ppv, NA_real_)]

  sesp_pool[, `:=`(
    score_lab = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS),
    ca_lab    = fifelse(ca_01 == 1L, "Cancer", "Non-cancer")
  )]

  setorder(sesp_pool, score_lab, -ca_01)

  nne_table = sesp_pool[, .(
    Score           = score_lab,
    Cohort          = ca_lab,
    `Alert rate`    = fmt_pct(positivity),
    Sensitivity     = fmt_pct(sensitivity),
    PPV             = fmt_pct(ppv),
    `NNE (1/PPV)`   = fifelse(is.na(nne), "—", sprintf("%.1f", nne))
  )]

  save_supplement(nne_table, "table_s_nne", "Number needed to evaluate",
                  footnote = NEWS_RULE_FOOTNOTE)

} else {
  message("  sesp_raw unavailable; skipping NNE table")
}

leadtime_table        = data.table()
leadtime_median_table = data.table()

## S4a. cumulative lead-time distribution (poolable) ---------------------------
# leadtime_raw holds n_at_or_below and n_total at fixed thresholds, so pooling is
# summation with no distributional assumption. "Within x" comes straight from the
# ratio; "beyond x" is its complement. Reported alongside the alert rate, because
# a long warning time from a score that flags most encounters is a different
# claim from the same warning time at low positivity.

if (have("leadtime_raw")) {

  lt = drop_round_one(leadtime_raw,
                      c("crossing_def", "threshold_h", "n_at_or_below", "n_total"),
                      "leadtime")

  if (nrow(lt) == 0L) {
    message("  leadtime: no round-two rows; skipping distribution table")
  }

  lt = lt[, .(
    n_at_or_below = sum(n_at_or_below),
    n_total       = sum(n_total)
  ), by = .(score_name, ca_01, outcome_key, crossing_def, threshold_h)]

  lt[, pct_within := n_at_or_below / n_total]

  lt[, `:=`(
    score_lab   = factor(score_name,  levels = names(SCORE_LABS),   labels = SCORE_LABS),
    ca_lab      = fifelse(ca_01 == 1L, "Cancer", "Non-cancer"),
    outcome_lab = factor(outcome_key, levels = names(OUTCOME_LABS), labels = OUTCOME_LABS),
    def_lab     = fifelse(crossing_def == "first",
                          "First crossing", "Final positive run")
  )]

  setorder(lt, outcome_lab, def_lab, score_lab, -ca_01, threshold_h)

  leadtime_table = dcast(
    lt[, .(outcome_lab, def_lab, score_lab, ca_lab,
           th = paste0("\u2264 ", threshold_h, " h"),
           cell = fmt_pct(pct_within))],
    outcome_lab + def_lab + score_lab + ca_lab ~ factor(th, levels = unique(th)),
    value.var = "cell",
    fill = "—"
  )

  setnames(leadtime_table,
           c("outcome_lab", "def_lab", "score_lab", "ca_lab"),
           c("Outcome", "Crossing", "Score", "Cohort"))

  save_supplement(leadtime_table, "table_s_leadtime", "Lead-time distribution")

} else {
  message("  leadtime_raw unavailable; skipping lead-time distribution table")
}

## S4b. per-unit medians, reported as a range ----------------------------------
# Each median is exact, computed on line-level data at its own site. Sites are
# not pooled; the reported statistic is the range across units. Hospital-level
# units are filtered by LEADTIME_EVENT_FLOOR first, because the range is the
# statistic most sensitive to the smallest contributors and some health systems
# report nine or more hospitals.

LEADTIME_EVENT_FLOOR = 20L

if (have("leadtime_median_raw")) {

  lm = drop_round_one(leadtime_median_raw,
                      c("crossing_def", "unit", "unit_id", "n_events", "median_h"),
                      "leadtime_median")

  n_dropped = lm[unit == "hospital" & n_events < LEADTIME_EVENT_FLOOR, .N]
  if (n_dropped > 0L) {
    message("  Dropped ", n_dropped, " hospital row(s) below ",
            LEADTIME_EVENT_FLOOR, " events")
  }

  lm = lm[unit == "health_system" | n_events >= LEADTIME_EVENT_FLOOR]

  lm_range = lm[, .(
    k_units  = .N,
    n_events = sum(n_events),
    md_lo    = min(median_h),
    md_md    = median(median_h),
    md_hi    = max(median_h)
  ), by = .(score_name, ca_01, outcome_key, crossing_def, unit)]

  lm_range[, `:=`(
    score_lab   = factor(score_name,  levels = names(SCORE_LABS),   labels = SCORE_LABS),
    ca_lab      = fifelse(ca_01 == 1L, "Cancer", "Non-cancer"),
    outcome_lab = factor(outcome_key, levels = names(OUTCOME_LABS), labels = OUTCOME_LABS),
    def_lab     = fifelse(crossing_def == "first",
                          "First crossing", "Final positive run"),
    unit_lab    = fifelse(unit == "health_system", "Health system", "Hospital")
  )]

  setorder(lm_range, outcome_lab, def_lab, unit_lab, score_lab, -ca_01)

  leadtime_median_table = lm_range[, .(
    Outcome                     = outcome_lab,
    Crossing                    = def_lab,
    Unit                        = unit_lab,
    Score                       = score_lab,
    Cohort                      = ca_lab,
    Units                       = k_units,
    Events                      = format_n(n_events),
    `Median of unit medians, h` = sprintf("%.1f", md_md),
    `Range across units, h`     = sprintf("%.1f-%.1f", md_lo, md_hi)
  )]

  save_supplement(leadtime_median_table, "table_s_leadtime_median", "Lead-time medians")

} else {
  message("  leadtime_median_raw unavailable; skipping lead-time median table")
}

crossclass_table = data.table()

if (have("crossclass_raw")) {

  cc = as.data.table(crossclass_raw)[, .(n = sum(n)),
                                     by = .(score_name, ca_01, outcome_key, crossed, event)]

  cc_wide = cc[, .(
    n_total       = sum(n),
    crossed_event = sum(n[crossed == 1L & event == 1L]),
    crossed_no    = sum(n[crossed == 1L & event == 0L]),
    nocross_event = sum(n[crossed == 0L & event == 1L]),
    nocross_no    = sum(n[crossed == 0L & event == 0L])
  ), by = .(score_name, ca_01, outcome_key)]

  cc_wide[, `:=`(
    pct_crossed        = (crossed_event + crossed_no) / n_total,
    pct_event_crossed  = crossed_event / (crossed_event + nocross_event),
    pct_missed         = nocross_event / (crossed_event + nocross_event)
  )]

  cc_wide[, `:=`(
    score_lab   = factor(score_name,  levels = names(SCORE_LABS),   labels = SCORE_LABS),
    ca_lab      = fifelse(ca_01 == 1L, "Cancer", "Non-cancer"),
    outcome_lab = factor(outcome_key, levels = names(OUTCOME_LABS), labels = OUTCOME_LABS)
  )]

  setorder(cc_wide, outcome_lab, score_lab, -ca_01)

  crossclass_table = cc_wide[, .(
    Outcome                      = outcome_lab,
    Score                        = score_lab,
    Cohort                       = ca_lab,
    `Ever crossed threshold`     = fmt_pct(pct_crossed),
    `Events preceded by crossing` = fmt_pct(pct_event_crossed),
    `Events never flagged`       = fmt_pct(pct_missed)
  )]

  save_supplement(crossclass_table, "table_s_crossclass", "Threshold crossing",
                  footnote = NEWS_RULE_FOOTNOTE)

} else {
  message("  crossclass_raw unavailable; skipping crossing table")
}

# ==============================================================================
# S5. MONITORING INTENSITY AND COMPONENT MISSINGNESS  (R30)
# ==============================================================================
# Reviewer 4 observed that carry-forward plus treat-missing-as-normal interacts
# with monitoring intensity, which may differ by cancer status. The monitoring
# artifact reports, per encounter and measure, the number of distinct
# measurement times within the at-risk ward window and the implied rate per 24
# hours; missing_vlab reports the share of encounters with no value at all.

message("\n-- S5: monitoring intensity and missingness (R30) --")

monitoring_table = data.table()

if (have("monitoring_raw")) {

  mon = as.data.table(monitoring_raw)[, .(
    n_enc         = sum(n_enc),
    sum_obs       = sum(sum_obs),
    sumsq_obs     = sum(sumsq_obs),
    sum_ward_h    = sum(sum_ward_hours),
    sum_rate      = sum(sum_rate_per24h),
    sumsq_rate    = sum(sumsq_rate_per24h)
  ), by = .(ca_01, measure)]

  mon[, `:=`(
    mean_obs   = sum_obs  / n_enc,
    sd_obs     = calculate_sd_from_sums(sum_obs,  sumsq_obs,  n_enc),
    mean_rate  = sum_rate / n_enc,
    sd_rate    = calculate_sd_from_sums(sum_rate, sumsq_rate, n_enc)
  )]

  mon[, ca_lab := fifelse(ca_01 == 1L, "Cancer", "Non-cancer")]
  setorder(mon, measure, -ca_01)

  monitoring_table = mon[, .(
    Measure                       = measure,
    Cohort                        = ca_lab,
    Encounters                    = format_n(n_enc),
    `Measurements, mean (SD)`     = sprintf("%.1f (%.1f)", mean_obs,  sd_obs),
    `Per 24 h, mean (SD)`         = sprintf("%.2f (%.2f)", mean_rate, sd_rate)
  )]

  save_supplement(monitoring_table, "table_s_monitoring", "Monitoring intensity")

} else {
  message("  monitoring_raw unavailable; skipping monitoring table")
}

monitoring_bins_table = data.table()

if (have("monitoring_bins_raw")) {

  mb = as.data.table(monitoring_bins_raw)[, .(n_enc = sum(n_enc)),
                                          by = .(ca_01, measure, rate_bin)]
  mb[, pct := n_enc / sum(n_enc), by = .(ca_01, measure)]

  monitoring_bins_table = dcast(
    mb[, .(measure, rate_bin, ca_01, cell = fmt_pct(pct))],
    measure + rate_bin ~ ca_01,
    value.var = "cell",
    fill = "—"
  )

  setnames(monitoring_bins_table, c("0", "1"), c("No Cancer", "Cancer"), skip_absent = TRUE)
  setnames(monitoring_bins_table, c("measure", "rate_bin"),
           c("Measure", "Measurements per 24 h"))

  save_supplement(monitoring_bins_table, "table_s_monitoring_bins", "Monitoring bins")

} else {
  message("  monitoring_bins_raw unavailable; skipping monitoring-bin table")
}

missingness_table = data.table()

if (have("missing_vlab_raw")) {

  mv = as.data.table(missing_vlab_raw)[, .(
    n_total   = sum(n_total),
    n_missing = sum(n_missing)
  ), by = .(variable, ca_01)]

  mv[, pct_missing := n_missing / n_total]
  mv[, ca_lab := fifelse(ca_01 == 1L, "Cancer", "Non-cancer")]

  missingness_table = dcast(
    mv[, .(variable, ca_lab,
           cell = paste0(format_n(n_missing), " (", fmt_pct(pct_missing), ")"))],
    variable ~ ca_lab,
    value.var = "cell",
    fill = "—"
  )

  setnames(missingness_table, "variable", "Component")

  save_supplement(missingness_table, "table_s_missingness", "Component missingness")

} else {
  message("  missing_vlab_raw unavailable; skipping missingness table")
}

# ==============================================================================
# S6. BETWEEN-SITE DESCRIPTIVE HETEROGENEITY  (R31)
# ==============================================================================
# R31 promises, alongside I-squared and tau-squared (produced in
# 02_discrimination.R), the range across sites of the primary outcome rate and
# of standard-threshold alert positivity. Both come from the per-site sesp
# cells: the outcome rate is invariant across scores within a site, so one score
# supplies it.

message("\n-- S6: between-site descriptive heterogeneity (R31) --")

site_range_table = data.table()

if (have("sesp_raw")) {

  ss = as.data.table(sesp_raw)

  ## outcome rate by site (one score as reference) ------------------------------

  ref_score = ss$score_name[1]

  rate_by_site = ss[score_name == ref_score, .(
    outcome_rate = sum(n_outcome) / sum(n_total)
  ), by = .(site, ca_01)]

  rate_range = rate_by_site[, .(
    metric = "Primary outcome rate",
    score  = "—",
    k      = .N,
    lo     = min(outcome_rate),
    md     = median(outcome_rate),
    hi     = max(outcome_rate)
  ), by = ca_01]

  ## alert positivity by site and score ----------------------------------------

  pos_by_site = ss[, .(positivity = sum(n_pos) / sum(n_total)),
                   by = .(site, ca_01, score_name)]

  pos_range = pos_by_site[, .(
    metric = "Alert positivity at standard threshold",
    k      = .N,
    lo     = min(positivity),
    md     = median(positivity),
    hi     = max(positivity)
  ), by = .(ca_01, score_name)]

  setnames(pos_range, "score_name", "score")

  site_range = rbindlist(list(rate_range, pos_range), use.names = TRUE, fill = TRUE)

  site_range[, `:=`(
    ca_lab    = fifelse(ca_01 == 1L, "Cancer", "Non-cancer"),
    score_lab = fifelse(score == "—", "—",
                        as.character(factor(score,
                                            levels = names(SCORE_LABS),
                                            labels = SCORE_LABS)))
  )]

  setorder(site_range, metric, score_lab, -ca_01)

  site_range_table = site_range[, .(
    Measure               = metric,
    Score                 = score_lab,
    Cohort                = ca_lab,
    Sites                 = k,
    `Median (range)`      = paste0(fmt_pct(md), " (", fmt_pct(lo), "-", fmt_pct(hi), ")")
  )]

  save_supplement(site_range_table, "table_s_site_range", "Between-site ranges",
                  footnote = NEWS_RULE_FOOTNOTE)

} else {
  message("  sesp_raw unavailable; skipping between-site range table")
}

# ==============================================================================
# EXPORTS
# ==============================================================================

message("\n== 08_supplement.R complete ==")

firstscore_table_final      = firstscore_table
hospital_type_table_final   = hospital_type_table
adm_chapter_table_final     = adm_chapter_table
adm_stem_table_final        = adm_stem_table
nne_table_final             = nne_table
leadtime_table_final        = leadtime_table
leadtime_median_table_final = leadtime_median_table
crossclass_table_final      = crossclass_table
monitoring_table_final      = monitoring_table
monitoring_bins_table_final = monitoring_bins_table
missingness_table_final     = missingness_table
site_range_table_final      = site_range_table

gc()
