# Paper 2: sensitivity analyses of docs/paper2/P2_04_review_sensitivity_prespec.md (internal review, 2026-10-15).
# Refuses to run if that plan changed after its hash was recorded (results/p2/rev_prespec_sha256.txt).
# Needs the Part A runs results/p2 (original), results/p2/sens_riskset and results/p2/sens_novisualqc.
#   Rscript stats/p2_review_sens.R
# Outputs: results/p2/rev_sens_estimates.csv, rev_sens_summary.txt
.libPaths(c("env/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(survival) })
PLAN <- "docs/paper2/P2_04_review_sensitivity_prespec.md"
rec <- strsplit(readLines("results/p2/rev_prespec_sha256.txt")[1], " +")[[1]][1]
now <- strsplit(system(paste("sha256sum", PLAN), intern = TRUE), " +")[[1]][1]
if (rec != now) stop("analysis plan changed after its hash was recorded")
sink("results/p2/rev_sens_summary.txt", split = TRUE)
cat("Paper 2 review sensitivity analyses, run", format(Sys.time()), "\nplan", PLAN, "sha256", now, "\n\n")

met <- fread("results/p2/p2_paraspinal.csv", select = c("case_id", "lung_vol_l", "lung_laa950_pct"))
met <- unique(met, by = "case_id")
cases <- rbindlist(lapply(c("p2", "p2b"), function(s)
  fread(sprintf("data/nlst_serial/%s_cases.csv", s), colClasses = list(character = c("pid", "rounds")))))
three <- cases[nchar(rounds) == 3, pid]

COMP <- c(T10_trab_hu = 1, psm_lvol = 1, lung_perc15_hu = 1, psm_fat_pct = -1, psm_hu = 1)
safe_predict <- function(fm, nd) {
  nd <- copy(nd)
  lv <- fm$xlevels$vendor0
  nd[!(as.character(vendor0) %in% lv), vendor0 := names(which.max(table(fm$model$vendor0)))]
  predict(fm, newdata = nd)
}
# baseline scores as in stats/p2_phenotype.R; lung = "perc15" (as published), "voladj" (log lung volume added to the
# reference model) or "laa950" (log(%LAA-950 + 0.1), reversed; noise_c enters the reference model because the normative
# model has no noise coefficient for it)
baseline <- function(dir, lung = "perc15") {
  long <- fread(file.path(dir, "partA_long.csv"), colClasses = list(character = c("pid", "set")))
  norm <- fread(file.path(dir, "partA_normative_model.csv"))
  long[, age0c := age0 - 62]
  long[, `:=`(sex = factor(sex), smoke = factor(smoke), vendor0 = factor(vendor0))]
  long <- merge(long, met, by = "case_id", all.x = TRUE)
  gam <- function(v) norm[measure == v & term == "noise_c", estimate]
  for (v in names(COMP)) long[, (paste0("y_", v)) := get(v) - gam(v) * noise_c]
  long[, lvol := log(lung_vol_l)]
  if (lung == "laa950") long[, y_lung_perc15_hu := -log(lung_laa950_pct + 0.1)]
  s <- long[order(pid, t)][, .SD[1], by = pid]
  for (v in names(COMP)) {
    y <- paste0("y_", v)
    rhs <- "age0c + sex + smoke + vendor0"
    if (v == "lung_perc15_hu" && lung == "voladj") rhs <- paste(rhs, "+ lvol")
    if (v == "lung_perc15_hu" && lung == "laa950") rhs <- paste(rhs, "+ noise_c")
    ok <- !is.na(s[[y]]) & (if (v == "lung_perc15_hu" && lung == "voladj") !is.na(s$lvol) else TRUE)
    fm <- lm(as.formula(paste(y, "~", rhs)), data = s[group == "control" & ok])
    s[, (paste0("z_", v)) := NA_real_]
    s[ok, (paste0("z_", v)) := COMP[[v]] * (get(y) - safe_predict(fm, s[ok])) / sigma(fm)]
  }
  std <- function(x) x / sd(x[s$group == "control"], na.rm = TRUE)
  s[, S := std(rowMeans(cbind(z_T10_trab_hu, z_psm_lvol, z_lung_perc15_hu)))]
  s[, S_msk := std(rowMeans(cbind(z_T10_trab_hu, z_psm_lvol)))]
  s[, case := as.integer(group == "case")]
  d <- s[!is.na(S)]
  d <- d[set %in% d[case == 1, set] & set %in% d[case == 0, set]]
  d[, `:=`(negS = -S, negS_msk = -S_msk, neg_zL = -z_lung_perc15_hu)]
  d
}

est <- list()
fit <- function(d, f, term, analysis, label) {
  m <- clogit(as.formula(f), data = d)
  ci <- summary(m)$conf.int; co <- summary(m)$coefficients
  est[[length(est) + 1]] <<- data.table(analysis = analysis, model = label, or = ci[term, "exp(coef)"],
                                       lo = ci[term, "lower .95"], hi = ci[term, "upper .95"],
                                       p = co[term, "Pr(>|z|)"], n_cases = sum(d$case), n_controls = sum(1 - d$case))
}
score_models <- function(d, analysis, components = TRUE) {
  fit(d, "case ~ negS + strata(set)", "negS", analysis, "Primary score, per 1 SD lower")
  if (components) for (v in c("T10_trab_hu", "psm_lvol", "lung_perc15_hu", "psm_fat_pct", "psm_hu")) {
    dv <- copy(d)[, negz := -get(paste0("z_", v))]
    fit(dv[!is.na(negz)], "case ~ negz + strata(set)", "negz", analysis, paste0("Component: ", v))
  }
  fit(d, "case ~ negS_msk + neg_zL + strata(set)", "negS_msk", analysis, "Musculoskeletal score adjusted for lung z")
}
conc <- function(d, var, B = 1000, seed = 20261015) {
  pr <- merge(d[case == 1, .(set, xc = get(var))], d[case == 0, .(set, xk = get(var))], by = "set",
              allow.cartesian = TRUE)
  pr[, w := fifelse(xc < xk, 1, fifelse(xc == xk, 0.5, 0))]
  ps <- pr[, .(w = sum(w), n = .N), by = set]
  set.seed(seed)
  idx <- replicate(B, sample(nrow(ps), replace = TRUE), simplify = FALSE)
  list(ps = ps, est = sum(ps$w) / sum(ps$n), idx = idx)
}
add_conc <- function(d, analysis, var, label) {
  k <- conc(d, var)
  bs <- sapply(k$idx, function(i) sum(k$ps$w[i]) / sum(k$ps$n[i]))
  est[[length(est) + 1]] <<- data.table(analysis = analysis, model = label, or = k$est, lo = quantile(bs, .025),
                                       hi = quantile(bs, .975), p = NA_real_, n_cases = sum(d$case),
                                       n_controls = sum(1 - d$case))
  k
}

# 1-2. risk-set controls; no visual-QC exclusions
for (run in list(c("results/p2/sens_riskset", "1. Risk-set controls"),
                 c("results/p2/sens_novisualqc", "2. No visual-QC exclusions"))) {
  d <- baseline(run[1])
  score_models(d, run[2])
  add_conc(d, run[2], "S", "Concordance: score")
  dv <- fread(file.path(run[1], "partA_divergence.csv"))[measure == "psm_fat_pct" & grepl("interval z", analysis)]
  for (i in seq_len(nrow(dv))) est[[length(est) + 1]] <- data.table(
    analysis = run[2], model = paste("Primary hypothesis, intramuscular fat interval z,", dv$term[i]),
    or = dv$estimate[i], lo = dv$lo[i], hi = dv$hi[i], p = dv$p[i], n_cases = NA_integer_, n_controls = NA_integer_)
}
# 3. sets whose case had three pre-diagnostic rounds (original data)
d0 <- baseline("results/p2")
fit(d0[set %in% three], "case ~ negS + strata(set)", "negS", "3. Case with three pre-diagnostic rounds",
    "Primary score, per 1 SD lower")
# 4. discrimination
a4 <- "4. Discrimination (original data)"
kS <- add_conc(d0, a4, "S", "Concordance: score")
kL <- add_conc(d0, a4, "z_lung_perc15_hu", "Concordance: lung density z alone")
add_conc(d0, a4, "S_msk", "Concordance: musculoskeletal score alone")
dif <- sapply(seq_along(kS$idx), function(b) {
  i <- kS$idx[[b]]
  sum(kS$ps$w[i]) / sum(kS$ps$n[i]) - sum(kL$ps$w[i]) / sum(kL$ps$n[i])
})
stopifnot(identical(kS$ps$set, kL$ps$set))
est[[length(est) + 1]] <- data.table(analysis = a4, model = "Concordance difference: score minus lung alone",
                                     or = kS$est - kL$est, lo = quantile(dif, .025), hi = quantile(dif, .975),
                                     p = NA_real_, n_cases = sum(d0$case), n_controls = sum(1 - d0$case))
# 5-6. lung density adjusted for lung volume; %LAA-950
for (run in list(c("voladj", "5. Lung density adjusted for lung volume"), c("laa950", "6. %LAA below -950 HU"))) {
  d <- baseline("results/p2", run[1])
  fit(d, "case ~ neg_zL + strata(set)", "neg_zL", run[2], "Lung z alone, per 1 SD lower")
  score_models(d, run[2], components = FALSE)
}
# 7. E-values (risk-ratio approximation)
ev <- function(rr) ifelse(rr <= 1, 1, rr + sqrt(rr * (rr - 1)))
E <- rbindlist(est, fill = TRUE)
prim <- fread("results/p2/pheno_estimates.csv")
for (lab in c("Primary: composite score (bone + muscle + lung), per 1 SD lower",
              "Musculoskeletal score, adjusted for lung density z")) {
  r <- prim[analysis == lab]
  E <- rbind(E, data.table(analysis = "7. E-value", model = paste0(lab, " (OR ", sprintf("%.2f", r$or), ")"),
                           or = ev(r$or), lo = ev(r$lo), hi = NA_real_, p = NA_real_), fill = TRUE)
}
fwrite(E, "results/p2/rev_sens_estimates.csv")
print(E[, .(analysis, model, est = sprintf("%.3f (%.3f-%.3f)", or, lo, hi), p = signif(p, 2), n_cases, n_controls)],
      nrows = 200)
cat("\n7. E-value rows: 'or' = E-value of the point estimate, 'lo' = E-value of the lower confidence limit\n")
sink()
