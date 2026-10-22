# Paper 2, Part A (public data only): normative muscle ageing, reliability of change, muscle-ageing acceleration (MAA)
# and pre-diagnostic divergence of lung-cancer cases (resected set p2 + non-resected set p2b) from matched cancer-free
# controls. Design: docs/paper2/P2_00 (public-data version).
#
# Inputs : results/p2/p2_paraspinal.csv (pipeline/p2_metrics.py), data/nlst_serial/{p2,p2b}_{series,cases,controls}.csv
#          (P2_SETS=p2 restricts to the resected set), data/nlst/nlst_prsn_public.csv (can_scr: mode of detection)
# Outputs: results/p2/partA_*.csv, results/p2/partA_summary.txt, results/p2/p2_maa.csv (one row per participant)
# The primary exposure is chosen here by the rule pre-specified in docs/paper2/P2_01 section 4.2, using controls only;
# no outcome data exist in this project when this runs.
# Usage: Rscript stats/p2_partA.R [out_dir]   (out_dir defaults to results/p2)
#        env P2_CTRL=_riskset P2_QC=<review csv> P2_BOOT=0 Rscript stats/p2_partA.R results/p2/sens_riskset
#        (docs/paper2/P2_04 sensitivity analyses; defaults reproduce the original analysis)
.libPaths(c("env/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(lme4); library(splines); library(survival) })
set.seed(20261009)
args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args)) args[1] else "results/p2"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, "partA_summary.txt"), split = TRUE)
cat("Paper 2 Part A run", format(Sys.time()), "\n\n")

# ---------------- data ----------------
met <- fread(Sys.getenv("P2_MET", "results/p2/p2_paraspinal.csv"))
if (!"error" %in% names(met)) met[, error := NA_character_]
met <- unique(met[order(!is.na(error))], by = "case_id")   # one row per scan (error-free row preferred)
SETS <- strsplit(Sys.getenv("P2_SETS", "p2,p2b"), ",")[[1]]
rd <- function(stem, cc) rbindlist(lapply(SETS, function(s) {
  f <- sprintf("data/nlst_serial/%s_%s.csv", s, stem)
  if (file.exists(f)) fread(f, colClasses = list(character = cc))[, set_name := s]
}), fill = TRUE)
ser <- rd("series", "pid")
cas <- rd("cases", c("pid", "rounds"))
# P2_CTRL=_riskset reads the risk-set control sets (pipeline/p2_resample_riskset.py) and the series of new controls
CTRL <- Sys.getenv("P2_CTRL", "")
ctl <- rd(paste0("controls", CTRL), c("pid", "matched_case"))
if (CTRL != "") ser <- unique(rbind(ser, rd(paste0("series", CTRL, "_new"), "pid"), fill = TRUE), by = "case_id")
cas[, resected := set_name == "p2"]
prs <- fread("data/nlst/nlst_prsn_public.csv", select = c("pid", "can_scr"), colClasses = list(character = "pid"))
cas <- merge(cas, prs, by = "pid", all.x = TRUE)
cat("sets:", paste(SETS, collapse = ", "), sprintf("| cases %d (resected %d) | controls %d\n", nrow(cas), sum(cas$resected), nrow(ctl)))

ppl <- rbind(
  cas[, .(pid, group = "case", set = pid, age0 = age, sex = gender, smoke = cigsmok, index_days = candx_days,
          scr_days0, scr_days1, scr_days2, stage7 = de_stag_7thed, morph = lc_morph, lesionsize)],
  ctl[, .(pid, group = "control", set = matched_case, age0 = age, sex = gender, smoke = cigsmok,
          index_days = index_days, scr_days0, scr_days1, scr_days2, stage7 = NA, morph = NA, lesionsize = NA)])
ppl[, `:=`(sex = factor(sex, 1:2, c("male", "female")), smoke = factor(smoke, 0:1, c("former", "current")))]

d <- merge(ser[, .(case_id, pid, round, kernel, vendor, kernel_consistent)], met, by = "case_id")
d <- merge(d, ppl, by = "pid")
d[, scan_days := fifelse(round == 0, scr_days0, fifelse(round == 1, scr_days1, scr_days2))]
n_raw <- d[, .N]
# sensitivity option: the independent segmentation model (TotalSegmentator `total` autochthon) replaces paraspinal
SEG <- Sys.getenv("P2_SEG", "psm")
if (SEG == "aut") for (k in c("available", "edge", "vol_cm3", "hu", "fat_pct")) d[[paste0("psm_", k)]] <- d[[paste0("aut_", k)]]
# visual QC review (results/p2/p2_qc_review.csv): scans marked qc_pass = 0 are excluded, with the reason recorded there
qcf <- Sys.getenv("P2_QC", "results/p2/p2_qc_review.csv")
if (file.exists(qcf)) {
  qc <- fread(qcf)
  d <- d[!case_id %in% qc[qc_pass == 0, case_id]]
  cat(sprintf("visual QC: %d scans excluded (%s)\n", sum(qc$qc_pass == 0), qcf))
}
d <- d[is.na(error) & psm_available == 1 & psm_edge == 0 & !is.na(scan_days) & !is.na(psm_hu)]
cat(sprintf("scans measured %d; passing QC (paraspinal complete over T6-T11, not at image edge) %d\n", n_raw, nrow(d)))

d[, psm_lvol := 100 * log(psm_vol_cm3)]          # slope = approx % per year
d[, noise := log(liver_sd_raw)]
d[is.na(noise), noise := median(d$noise, na.rm = TRUE)]
d[, noise_c := noise - median(noise)]
d[, age0c := age0 - 62]
d[, vendor := factor(vendor)]
d[, vendor0 := vendor[which.min(scan_days)], by = pid]   # person-level vendor (first scan) for person-level models
d[, group := factor(group, c("control", "case"))]
# time since the participant's first measured scan (years) and time relative to the case's diagnosis day
d[, t := (scan_days - min(scan_days)) / 365.25, by = pid]
d[, tau := (scan_days - index_days) / 365.25]
d <- d[group == "control" | scan_days <= index_days]
FILTER <- Sys.getenv("P2_FILTER", "none")
if (FILTER == "kernel_consistent") d <- d[kernel_consistent == 1]
if (FILTER == "gap90") d <- d[group == "control" | tau <= -90 / 365.25]   # drop case scans within 90 days of diagnosis
if (FILTER == "three_rounds") d <- d[pid %in% d[, .N, by = pid][N == 3, pid]]
cat(sprintf("sensitivity settings: segmentation = %s, filter = %s\n", SEG, FILTER))
keep <- d[, .N, by = pid][N >= 2, pid]
d <- d[pid %in% keep]
cat(sprintf("participants with >=2 QC-passing scans: cases %d, controls %d\n\n",
            uniqueN(d[group == "case", pid]), uniqueN(d[group == "control", pid])))
fwrite(data.table(item = c("scans_measured", "scans_after_qc_and_filters", "cases_analysed", "controls_analysed",
                           "case_scans", "control_scans", "cases_3_scans", "controls_3_scans"),
                  n = c(n_raw, nrow(d), uniqueN(d[group == "case", pid]), uniqueN(d[group == "control", pid]),
                        d[group == "case", .N], d[group == "control", .N],
                        sum(d[group == "case", .N, by = pid]$N == 3), sum(d[group == "control", .N, by = pid]$N == 3))),
       file.path(OUT, "partA_counts.csv"))

# predict from a control-fitted lm; a vendor never seen among those controls gets the most common vendor's effect
safe_predict <- function(fm, nd) {
  nd <- copy(as.data.table(nd))
  lv <- fm$xlevels$vendor0
  if (!is.null(lv)) nd[!(as.character(vendor0) %in% lv), vendor0 := names(which.max(table(fm$model$vendor0)))]
  predict(fm, newdata = nd)
}

MEAS <- data.table(var = c("psm_hu", "psm_lvol", "psm_fat_pct", "T10_trab_hu", "lung_perc15_hu"),
                   label = c("Paraspinal lean-muscle radiodensity (HU)", "Paraspinal volume (100 x log cm3)",
                             "Paraspinal intramuscular fat (%)", "T10 trabecular bone (HU)",
                             "Lung density, 15th percentile (HU; emphysema proxy, positive control)"),
                   ageing_sign = c(-1, -1, 1, -1, -1), candidate = c(TRUE, TRUE, TRUE, FALSE, FALSE))
MEAS <- MEAS[var %in% names(d)]

# ---------------- A1 + A2: normative model and reliability in controls ----------------
ctrl <- d[group == "control"]
rel_fun <- function(dd, v) {
  dd <- dd[!is.na(get(v))]
  f <- lmer(as.formula(paste(v, "~ t + age0c + sex + smoke + vendor + noise_c + (1 + t | pid)")), data = dd,
            control = lmerControl(optimizer = "bobyqa", check.conv.singular = "ignore"))
  vc <- as.data.frame(VarCorr(f))
  tau0 <- vc$vcov[which(vc$grp == "pid" & vc$var1 == "(Intercept)" & is.na(vc$var2))]
  tau1 <- vc$vcov[which(vc$grp == "pid" & vc$var1 == "t" & is.na(vc$var2))]
  s2 <- sigma(f)^2
  sxx <- dd[, .(sxx = sum((t - mean(t))^2)), by = pid]$sxx
  c(slope_true_sd = sqrt(tau1), level_true_sd = sqrt(tau0), resid_sd = sqrt(s2), level_icc = tau0 / (tau0 + s2),
    slope_reliability = mean(tau1 / (tau1 + s2 / sxx)), mdc95_2y = 1.96 * sqrt(2) * sqrt(s2))
}
norm_rows <- list(); rel_rows <- list()
for (i in seq_len(nrow(MEAS))) {
  v <- MEAS$var[i]
  dd <- ctrl[!is.na(get(v))]
  f <- lmer(as.formula(paste(v, "~ t * (age0c + sex + smoke) + vendor + noise_c + (1 + t | pid)")), data = dd,
            control = lmerControl(optimizer = "bobyqa", check.conv.singular = "ignore"))
  co <- summary(f)$coefficients
  ci <- confint(f, method = "Wald")
  keep_terms <- grep("^t$|^t:|noise_c", rownames(co), value = TRUE)
  norm_rows[[v]] <- data.table(measure = v, term = keep_terms, estimate = co[keep_terms, 1],
                               lo = ci[keep_terms, 1], hi = ci[keep_terms, 2], se = co[keep_terms, 2])
  # reliability: point estimate, participant bootstrap, and subgroups
  pt <- rel_fun(dd, v)
  ids <- unique(dd$pid)
  boot <- replicate(max(1, as.integer(Sys.getenv("P2_BOOT", "200"))), if (Sys.getenv("P2_BOOT", "200") == "0") NA_real_ else {
    s <- sample(ids, replace = TRUE)
    bd <- dd[data.table(pid = s, new = seq_along(s)), on = "pid", allow.cartesian = TRUE][, pid := paste(pid, new)]
    tryCatch(rel_fun(bd, v)["slope_reliability"], error = function(e) NA_real_)
  })
  rel_rows[[v]] <- cbind(data.table(measure = v, subset = "all controls", n = length(ids)), as.data.table(as.list(pt)),
                         data.table(rel_lo = quantile(boot, .025, na.rm = TRUE), rel_hi = quantile(boot, .975, na.rm = TRUE)))
  for (sub in list(list("kernel consistent", dd[kernel_consistent == 1]),
                   list("kernel switched", dd[kernel_consistent == 0]),
                   list("GE", dd[vendor == "GE"]), list("Siemens", dd[vendor == "SE"]))) {
    if (uniqueN(sub[[2]]$pid) < 40) next
    r <- tryCatch(rel_fun(sub[[2]], v), error = function(e) NULL)
    if (!is.null(r)) rel_rows[[paste(v, sub[[1]])]] <- cbind(data.table(measure = v, subset = sub[[1]],
                                                                        n = uniqueN(sub[[2]]$pid)), as.data.table(as.list(r)))
  }
}
norm <- rbindlist(norm_rows); rel <- rbindlist(rel_rows, fill = TRUE)
fwrite(norm, file.path(OUT, "partA_normative_model.csv")); fwrite(rel, file.path(OUT, "partA_reliability.csv"))
cat("A1 normative annual change (controls; reference = male, former smoker, age 62, median noise):\n")
print(norm[term == "t" | grepl("^t:", term), .(measure, term, est = round(estimate, 3),
                                               ci = sprintf("%.3f to %.3f", lo, hi))])
cat("\nA2 reliability of within-person change (controls):\n")
print(rel[, .(measure, subset, n, slope_rel = round(slope_reliability, 2),
              rel_ci = ifelse(is.na(rel_lo), "", sprintf("%.2f-%.2f", rel_lo, rel_hi)),
              level_icc = round(level_icc, 2), resid_sd = round(resid_sd, 2), mdc95_2y = round(mdc95_2y, 2))])

# cross-model agreement of change: paraspinal (abdominal_muscles task) vs autochthon (total task)
d[, aut_lvol := 100 * log(aut_vol_cm3)]
slope_of <- function(y, t) if (sum(!is.na(y)) >= 2) coef(lm(y ~ t))[2] else NA_real_
xm <- d[, .(psm_hu = slope_of(psm_hu, t), aut_hu = slope_of(aut_hu, t),
            psm_lvol = slope_of(psm_lvol, t), aut_lvol = slope_of(aut_lvol, t)), by = pid]
xm <- merge(xm, d[, .(n_scans = .N), by = pid], by = "pid")
xm[, psm_fat_pct := d[, .(v = slope_of(psm_fat_pct, t)), by = pid][match(xm$pid, pid), v]]
xm[, aut_fat_pct := d[, .(v = slope_of(aut_fat_pct, t)), by = pid][match(xm$pid, pid), v]]
fwrite(xm, file.path(OUT, "partA_crossmodel_slopes.csv"))
cat(sprintf("\nCross-model agreement of individual slopes (Spearman): radiodensity %.2f, volume %.2f\n",
            cor(xm$psm_hu, xm$aut_hu, use = "c", method = "s"), cor(xm$psm_lvol, xm$aut_lvol, use = "c", method = "s")))

# ---------------- pre-specified primary-exposure rule (P2_01 section 4.2) ----------------
r_all <- rel[subset == "all controls" & measure %in% MEAS[candidate == TRUE, var]]
if (r_all[measure == "psm_hu", slope_reliability] >= 0.40) {
  primary <- "psm_hu"; why <- "psm_hu slope reliability >= 0.40"
} else if (any(r_all$slope_reliability >= 0.40)) {
  primary <- r_all[which.max(slope_reliability), measure]; why <- "psm_hu < 0.40; most reliable candidate >= 0.40"
} else {
  primary <- NA_character_; why <- "no candidate reaches 0.40: individual change not reliably measurable (fallback)"
}
cat(sprintf("\nPRIMARY EXPOSURE RULE: %s -> %s\n", why, ifelse(is.na(primary), "level + slope joint model", primary)))
writeLines(c(paste("primary", primary), paste("reason", why), paste("run", format(Sys.time()))),
           file.path(OUT, "partA_primary_exposure.txt"))

# ---------------- MAA: observed slope vs normative slope for the same round pattern ----------------
# Noise-adjusted values (the normative noise coefficient), then OLS slope per participant.
maa_rows <- list()
for (i in seq_len(nrow(MEAS))) {
  v <- MEAS$var[i]; sg <- MEAS$ageing_sign[i]
  gam_noise <- norm[measure == v & term == "noise_c", estimate]
  dv <- d[!is.na(get(v))][, y := get(v) - gam_noise * noise_c]
  pp <- dv[, .(pattern = paste(sort(round), collapse = ""), n_scans = .N, span = max(t) - min(t),
               slope = if (.N >= 2) coef(lm(y ~ t))[2] else NA_real_, last = y[which.max(t)]),
           by = .(pid, group, set, age0c, sex, smoke, vendor0)]
  # normative slope for each pattern: controls re-measured on the same subset of rounds
  cpat <- rbindlist(lapply(c("01", "02", "12", "012"), function(p) {
    rr <- as.integer(strsplit(p, "")[[1]])
    dv[group == "control" & round %in% rr][, if (.N == length(rr)) .(slope = coef(lm(y ~ t))[2]),
                                           by = .(pid, age0c, sex, smoke, vendor0)][, pattern := p]
  }))
  for (p in unique(pp$pattern)) {
    cp <- cpat[pattern == p]
    if (nrow(cp) < 50) next
    fm <- lm(slope ~ age0c + sex + smoke + vendor0, data = cp)
    rsd <- sigma(fm)
    idx <- pp$pattern == p
    pp[idx, exp_slope := safe_predict(fm, pp[idx])]
    pp[idx, maa_z := sg * (slope - exp_slope) / rsd]
    pp[idx, norm_sd := rsd]
  }
  # "excess years of normative ageing per calendar year", relative to the mean control slope for the pattern
  mean_ctrl <- cpat[, .(mslope = mean(slope)), by = pattern]
  pp <- merge(pp, mean_ctrl, by = "pattern", all.x = TRUE)
  pp[, excess_years_per_year := (slope - exp_slope) / abs(mslope)]
  maa_rows[[v]] <- pp[, .(pid, measure = v, pattern, n_scans, span, slope, exp_slope, maa_z, excess_years_per_year,
                          last_value = last)]
}
maa <- rbindlist(maa_rows)
fwrite(maa, file.path(OUT, "p2_maa.csv"))

# analysis-level long data (one row per QC-passing scan) for figures and secondary analyses
fwrite(d[, c("case_id", "pid", "group", "set", "round", "t", "tau", "scan_days", "index_days", "age0", "sex", "smoke",
             "vendor", "vendor0", "kernel", "kernel_consistent", "noise_c", MEAS$var, "aut_hu", "aut_fat_pct",
             "aut_lvol"), with = FALSE], file.path(OUT, "partA_long.csv"))

# ---------------- aligned level difference: case scan minus its matched controls' mean in the same round ----------------
al <- rbindlist(lapply(MEAS$var, function(v) {
  gam <- norm[measure == v & term == "noise_c", estimate]
  a <- d[!is.na(get(v)), .(pid, group, set, round, tau, y = get(v) - gam * noise_c)]
  cm <- a[group == "control", .(cmean = mean(y), nctl = .N), by = .(set, round)]
  a[group == "case"][cm, on = .(set, round), nomatch = 0][, .(pid, round, tau, diff = y - cmean, nctl)][, measure := v]
}))
al[, tau_bin := cut(tau, c(-Inf, -4, -3, -2, -1, 0.001), right = FALSE,
                    labels = c("<= -4", "-4 to -3", "-3 to -2", "-2 to -1", "-1 to 0"))]
als <- al[, .(n = .N, mean_tau = mean(tau), mean_diff = mean(diff), lo = mean(diff) - 1.96 * sd(diff) / sqrt(.N),
              hi = mean(diff) + 1.96 * sd(diff) / sqrt(.N)), by = .(measure, tau_bin)][order(measure, tau_bin)]
fwrite(al, file.path(OUT, "partA_aligned_scans.csv")); fwrite(als, file.path(OUT, "partA_aligned_bins.csv"))

# ---------------- A3: pre-diagnostic divergence of cases from matched controls ----------------
cat("\nA3 case vs matched control (cases' scans before diagnosis; tau = years relative to the case's diagnosis day)\n")
a3_rows <- list()
mz <- merge(maa, unique(d[, .(pid, group, set, age0c, sex, smoke)]), by = "pid")
for (v in MEAS$var) {
  dv <- d[!is.na(get(v))]
  f <- lmer(as.formula(paste(v, "~ group * tau + age0c + sex + smoke + vendor + noise_c + (1 + tau | pid) + (1 | set)")),
            data = dv, control = lmerControl(optimizer = "bobyqa", check.conv.singular = "ignore"))
  co <- summary(f)$coefficients; ci <- confint(f, method = "Wald")
  for (tm in c("groupcase", "groupcase:tau")) {
    z <- co[tm, 1] / co[tm, 2]
    a3_rows[[paste(v, tm)]] <- data.table(measure = v, analysis = "mixed model", term = tm, estimate = co[tm, 1],
                                          lo = ci[tm, 1], hi = ci[tm, 2], p = 2 * pnorm(-abs(z)))
  }
  # same model additionally adjusted for lung density (emphysema proxy) at each scan
  if (v != "lung_perc15_hu" && "lung_perc15_hu" %in% names(dv)) {
    dl <- dv[!is.na(lung_perc15_hu)][, lung_c := (lung_perc15_hu - median(lung_perc15_hu)) / 10]
    fl <- lmer(as.formula(paste(v, "~ group * tau + age0c + sex + smoke + vendor + noise_c + lung_c + (1 + tau | pid) + (1 | set)")),
               data = dl, control = lmerControl(optimizer = "bobyqa", check.conv.singular = "ignore"))
    co2 <- summary(fl)$coefficients; ci2 <- confint(fl, method = "Wald")
    for (tm in c("groupcase", "groupcase:tau")) {
      a3_rows[[paste(v, tm, "lungadj")]] <- data.table(measure = v, analysis = "mixed model, adjusted for lung density",
                                                       term = tm, estimate = co2[tm, 1], lo = ci2[tm, 1], hi = ci2[tm, 2],
                                                       p = 2 * pnorm(-abs(co2[tm, 1] / co2[tm, 2])))
    }
  }
  # MAA of cases vs their matched controls: conditional logistic regression on matched sets
  m <- mz[measure == v & !is.na(maa_z)]
  cl <- clogit(I(group == "case") ~ maa_z + strata(set), data = m)
  s <- summary(cl)$conf.int
  a3_rows[[paste(v, "clogit")]] <- data.table(measure = v, analysis = "conditional logistic", term = "MAA per SD",
                                              estimate = s[1, 1], lo = s[1, 3], hi = s[1, 4],
                                              p = summary(cl)$coefficients[1, 5])
}
# interval-level: annualised change per interval, cases' intervals classified by time from interval end to diagnosis
iv <- rbindlist(lapply(MEAS$var, function(v) {
  gam_noise <- norm[measure == v & term == "noise_c", estimate]
  dv <- d[!is.na(get(v))][order(pid, t)][, y := get(v) - gam_noise * noise_c]
  dv[, .(r0 = head(round, -1), r1 = tail(round, -1), rate = diff(y) / diff(t), end_tau = tail(tau, -1)),
     by = .(pid, group, set, age0c, sex, smoke, vendor0)][, measure := v]
}))
iv[, sign := MEAS$ageing_sign[match(measure, MEAS$var)]]
# standardise each interval rate against controls' rates over the same pair of rounds; the control reference is built
# for every pair (0-1, 1-2 and 0-2), because some cases were imaged at T0 and T2 only
cref <- rbindlist(lapply(MEAS$var, function(v) {
  gam_noise <- norm[measure == v & term == "noise_c", estimate]
  dv <- d[group == "control" & !is.na(get(v))][, y := get(v) - gam_noise * noise_c]
  rbindlist(lapply(list(c(0, 1), c(1, 2), c(0, 2)), function(pr) {
    w <- dcast(dv[round %in% pr], pid + age0c + sex + smoke + vendor0 ~ round, value.var = c("y", "t"))
    setnames(w, c(paste0("y_", pr), paste0("t_", pr)), c("ya", "yb", "ta", "tb"))
    w[!is.na(ya) & !is.na(yb), .(pid, age0c, sex, smoke, vendor0, rate = (yb - ya) / (tb - ta))][, `:=`(r0 = pr[1], r1 = pr[2])]
  }))[, measure := v]
}))
iv[, z := {
  m_ <- .BY$measure; a_ <- .BY$r0; b_ <- .BY$r1
  cr <- cref[measure == m_ & r0 == a_ & r1 == b_]
  fm <- lm(rate ~ age0c + sex + smoke + vendor0, data = cr)
  sign * (rate - safe_predict(fm, .SD)) / sigma(fm)
}, by = .(measure, r0, r1)]
iv[, window := fifelse(group == "control", "control",
                fifelse(end_tau > -1, "case: ends <1 y before dx",
                        fifelse(end_tau > -2, "case: ends 1-2 y before dx", "case: ends >=2 y before dx")))]
ivs <- iv[, .(n = .N, mean_z = mean(z), lo = mean(z) - 1.96 * sd(z) / sqrt(.N), hi = mean(z) + 1.96 * sd(z) / sqrt(.N)),
          by = .(measure, window)][order(measure, window)]
fwrite(ivs, file.path(OUT, "partA_interval_windows.csv"))
fwrite(iv, file.path(OUT, "partA_intervals.csv"))      # interval-level z (for figures and subgroup analyses)
# case-vs-control difference in interval z by window (matched set and participant as random effects)
for (v in MEAS$var) {
  ivv <- iv[measure == v]
  ivv[, window := relevel(factor(window), ref = "control")]
  f <- lmer(z ~ window + (1 | set) + (1 | pid), data = ivv,
            control = lmerControl(check.conv.singular = "ignore"))
  co <- summary(f)$coefficients
  for (tm in grep("^window", rownames(co), value = TRUE)) {
    a3_rows[[paste(v, tm)]] <- data.table(measure = v, analysis = "interval z: case window vs controls", term = tm,
                                          estimate = co[tm, 1], lo = co[tm, 1] - 1.96 * co[tm, 2],
                                          hi = co[tm, 1] + 1.96 * co[tm, 2], p = 2 * pnorm(-abs(co[tm, 1] / co[tm, 2])))
  }
}
a3 <- rbindlist(a3_rows)
fwrite(a3, file.path(OUT, "partA_divergence.csv"))
print(a3[, .(measure, analysis, term, est = round(estimate, 3), ci = sprintf("%.3f to %.3f", lo, hi), p = signif(p, 2))])
cat("\nInterval z (ageing direction positive; controls centred at 0 by construction):\n")
print(ivs[, .(measure, window, n, mean_z = round(mean_z, 2), ci = sprintf("%.2f to %.2f", lo, hi))])

# dose-response among cases: stage, tumour size, histology, mode of detection, resection set
cs <- merge(mz[group == "case"], cas[, .(pid, stage7 = de_stag_7thed, morph = lc_morph, lesionsize, can_scr, resected)],
            by = "pid")
cs[, stage4 := fcase(stage7 %in% c(110, 120), "I", stage7 %in% c(210, 220), "II", stage7 %in% c(310, 320), "III",
                     stage7 == 400, "IV", default = NA_character_)]
cs[, histo := fcase(morph %in% c(8041, 8042, 8043, 8044, 8045), "small cell",
                    morph %in% c(8070, 8071, 8072, 8073, 8074), "squamous",
                    morph %in% c(8140, 8250, 8251, 8252, 8253, 8254, 8255, 8260, 8310, 8323, 8480, 8481, 8490, 8550,
                                 8551, 8560), "adenocarcinoma", default = "other")]
cs[, detection := fcase(can_scr == 1, "screen-detected", can_scr %in% c(2, 3, 4), "interval/post-screening",
                        default = NA_character_)]
cs[, size_mm := suppressWarnings(as.numeric(lesionsize))]
cs[size_mm > 200 | size_mm <= 0, size_mm := NA]
summ <- function(x, by) x[!is.na(get(by)), .(n = sum(!is.na(maa_z)), mean_maa = mean(maa_z, na.rm = TRUE),
                                             se = sd(maa_z, na.rm = TRUE) / sqrt(sum(!is.na(maa_z)))),
                          by = c("measure", by)][, `:=`(factor = by)][]
strat <- rbindlist(lapply(c("stage4", "histo", "detection", "resected"), function(b) {
  o <- summ(cs, b); setnames(o, b, "stratum"); o[, stratum := as.character(stratum)] }))
fwrite(strat, file.path(OUT, "partA_case_strata.csv"))
cat("\nCase MAA by stage, histology, detection mode and resection set (controls = 0 by construction):\n")
print(strat[order(measure, factor, stratum), .(measure, factor, stratum, n, mean_maa = round(mean_maa, 2), se = round(se, 2))])
tr <- rbindlist(lapply(MEAS$var, function(v) {
  x <- cs[measure == v & !is.na(maa_z)]
  a <- summary(lm(maa_z ~ as.integer(factor(stage4, c("I", "II", "III", "IV"))) + age0c + sex + smoke, data = x))$coefficients
  b <- summary(lm(maa_z ~ I(size_mm / 10) + age0c + sex + smoke, data = x))$coefficients
  data.table(measure = v, test = c("stage trend (per stage group)", "tumour size (per 10 mm)"),
             estimate = c(a[2, 1], b[2, 1]), se = c(a[2, 2], b[2, 2]), p = c(a[2, 4], b[2, 4]))
}))
fwrite(tr, file.path(OUT, "partA_dose_response.csv"))
cat("\nDose-response among cases (MAA z):\n"); print(tr[, .(measure, test, est = round(estimate, 3), se = round(se, 3), p = signif(p, 2))])
sink()
