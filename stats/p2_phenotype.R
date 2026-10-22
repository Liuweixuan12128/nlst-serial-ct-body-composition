# Paper 2: CT host-phenotype score (bone + muscle + lung reserve) and later lung cancer diagnosis.
# Implements docs/paper2/P2_03_phenotype_prespec.md exactly; refuses to run if that plan changed after its hash was
# recorded (results/p2/pheno_prespec_sha256.txt). Run once, from the project root, after stats/p2_partA.R:
#   Rscript stats/p2_phenotype.R
# Outputs: results/p2/pheno_scores.csv, pheno_estimates.csv, pheno_quartiles.csv, pheno_summary.txt
.libPaths(c("env/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(survival) })
OUT <- "results/p2"
PLAN <- "docs/paper2/P2_03_phenotype_prespec.md"
rec <- strsplit(readLines(file.path(OUT, "pheno_prespec_sha256.txt"))[1], " +")[[1]][1]
now <- strsplit(system(paste("sha256sum", PLAN), intern = TRUE), " +")[[1]][1]
if (rec != now) stop("analysis plan changed after its hash was recorded")
sink(file.path(OUT, "pheno_summary.txt"), split = TRUE)
cat("Paper 2 CT host-phenotype score, run", format(Sys.time()), "\nplan", PLAN, "sha256", now, "\n\n")

long <- fread(file.path(OUT, "partA_long.csv"), colClasses = list(character = c("pid", "set")))
norm <- fread(file.path(OUT, "partA_normative_model.csv"))
long[, age0c := age0 - 62]
long[, `:=`(sex = factor(sex), smoke = factor(smoke), vendor0 = factor(vendor0))]
gam <- function(v) norm[measure == v & term == "noise_c", estimate]
# orientation: +1 = higher value means more reserve
COMP <- c(T10_trab_hu = 1, psm_lvol = 1, lung_perc15_hu = 1, psm_fat_pct = -1, psm_hu = 1)
for (v in names(COMP)) long[, (paste0("y_", v)) := get(v) - gam(v) * noise_c]

# case-level descriptors (set = case pid)
cas <- rbindlist(lapply(c("p2", "p2b"), function(s)
  fread(sprintf("data/nlst_serial/%s_cases.csv", s), colClasses = list(character = c("pid", "rounds")))[, series := s]),
  fill = TRUE)
prs <- fread("data/nlst/nlst_prsn_public.csv", select = c("pid", "can_scr"), colClasses = list(character = "pid"))
cas <- merge(cas, prs, by = "pid", all.x = TRUE)
grp <- cas[, .(set = pid, series = fifelse(series == "p2", "Resected", "Not resected"),
               stage2 = fcase(de_stag_7thed %in% c(110, 120, 210, 220), "Stage I\u2013II",
                              de_stag_7thed %in% c(310, 320, 400), "Stage III\u2013IV", default = NA_character_),
               histo = fcase(lc_morph %in% c(8041, 8042, 8043, 8044, 8045), "Small cell",
                             lc_morph %in% c(8070, 8071, 8072, 8073, 8074), "Squamous",
                             lc_morph %in% c(8140, 8250, 8251, 8252, 8253, 8254, 8255, 8260, 8310, 8323, 8480, 8481,
                                             8490, 8550, 8551, 8560), "Adenocarcinoma", default = "Other"),
               detection = fcase(can_scr == 1, "Screen-detected", can_scr %in% c(2, 3, 4), "Interval/post-screening",
                                 default = NA_character_),
               sex_s = fifelse(gender == 1, "Men", "Women"),
               smoke_s = fifelse(cigsmok == 1, "Current smoker", "Former smoker"))]

safe_predict <- function(fm, nd) {
  nd <- copy(nd)
  lv <- fm$xlevels$vendor0
  nd[!(as.character(vendor0) %in% lv), vendor0 := names(which.max(table(fm$model$vendor0)))]
  predict(fm, newdata = nd)
}

# z-scores against controls' scans of the same kind (baseline or last), then the scores
build <- function(which = c("baseline", "last")) {
  which <- match.arg(which)
  s <- long[order(pid, t)][, if (which == "baseline") .SD[1] else .SD[.N], by = pid]
  for (v in names(COMP)) {
    y <- paste0("y_", v)
    ok <- !is.na(s[[y]])
    fm <- lm(as.formula(paste(y, "~ age0c + sex + smoke + vendor0")), data = s[group == "control" & ok])
    s[, (paste0("z_", v)) := NA_real_]
    s[ok, (paste0("z_", v)) := COMP[[v]] * (get(y) - safe_predict(fm, s[ok])) / sigma(fm)]
  }
  std <- function(x) x / sd(x[s$group == "control"], na.rm = TRUE)
  s[, S := std(rowMeans(cbind(z_T10_trab_hu, z_psm_lvol, z_lung_perc15_hu)))]           # NA unless all three
  s[, S_msk := std(rowMeans(cbind(z_T10_trab_hu, z_psm_lvol)))]
  s[, S5 := std(rowMeans(cbind(z_T10_trab_hu, z_psm_lvol, z_lung_perc15_hu, z_psm_fat_pct, z_psm_hu)))]
  s[, case := as.integer(group == "case")]
  s <- merge(s, grp, by = "set", all.x = TRUE)
  lead <- s[case == 1, .(set, lead = -tau)]
  s <- merge(s, lead, by = "set", all.x = TRUE)
  s[, lead3 := cut(lead, c(-Inf, 2, 4, Inf), right = FALSE, labels = c("<2 y", "2\u2013<4 y", "\u22654 y"))]
  s
}

est <- list()
fit <- function(d, f, label, term, population) {
  m <- clogit(as.formula(f), data = d)
  co <- summary(m)$coefficients; ci <- summary(m)$conf.int
  est[[length(est) + 1]] <<- data.table(analysis = label, population = population, term = term,
                                       or = ci[term, "exp(coef)"], lo = ci[term, "lower .95"],
                                       hi = ci[term, "upper .95"], p = co[term, "Pr(>|z|)"],
                                       n_cases = sum(d$case), n_controls = sum(1 - d$case),
                                       n_sets = uniqueN(d[case == 1, set]))
  invisible(m)
}
lr_int <- function(d, fac, label) {
  d <- d[!is.na(get(fac))]
  d <- d[set %in% d[case == 1, set] & set %in% d[case == 0, set]]
  d[, fac_ := factor(get(fac))]
  m0 <- clogit(case ~ negS + strata(set), data = d)
  m1 <- clogit(case ~ negS + negS:fac_ + strata(set), data = d)
  lr <- 2 * (m1$loglik[2] - m0$loglik[2])
  dfree <- sum(!is.na(coef(m1))) - 1
  p_lr <- pchisq(lr, dfree, lower.tail = FALSE)
  for (lv in levels(d$fac_)) {
    fit(d[fac_ == lv], "case ~ negS + strata(set)", paste0(label, ": ", lv), "negS", "baseline")
    est[[length(est)]] <<- copy(est[[length(est)]])[, p_interaction := p_lr]
  }
  p_lr
}

b <- build("baseline")
fwrite(b[, .(pid, set, group, series, tau, lead, z_T10_trab_hu, z_psm_lvol, z_lung_perc15_hu, z_psm_fat_pct, z_psm_hu,
             S, S_msk, S5)], file.path(OUT, "pheno_scores.csv"))
d <- b[!is.na(S)]
d <- d[set %in% d[case == 1, set] & set %in% d[case == 0, set]]
d[, `:=`(negS = -S, negS_msk = -S_msk, negS5 = -S5, neg_zL = -z_lung_perc15_hu)]
cat(sprintf("baseline scans with all three components, complete matched sets: %d cases, %d controls, %d sets\n",
            sum(d$case), sum(1 - d$case), uniqueN(d[case == 1, set])))
cat(sprintf("case baseline scan before diagnosis: median %.2f y (IQR %.2f-%.2f)\n\n",
            median(d[case == 1, lead]), quantile(d[case == 1, lead], .25), quantile(d[case == 1, lead], .75)))

# primary
fit(d, "case ~ negS + strata(set)", "Primary: composite score (bone + muscle + lung), per 1 SD lower", "negS", "baseline")
# 1. quartiles (control cut-points), Q4 = highest reserve = reference
qc <- quantile(d[case == 0, S], c(.25, .5, .75))
d[, Q := findInterval(S, qc) + 1L]           # Q1 = lowest reserve ... Q4 = highest
d[, Qf := relevel(factor(Q, 1:4, paste0("Q", 1:4)), ref = "Q4")]
mq <- clogit(case ~ Qf + strata(set), data = d)
ciq <- summary(mq)$conf.int; coq <- summary(mq)$coefficients
d[, Qrev := 5L - Q]                          # 1 = highest reserve ... 4 = lowest, so OR > 1 = risk rises
mt <- clogit(case ~ Qrev + strata(set), data = d)
p_trend <- summary(mt)$coefficients["Qrev", "Pr(>|z|)"]
qt <- data.table(quartile = paste0("Q", 1:4),
                 score_range = c(sprintf("< %.2f", qc[1]), sprintf("%.2f to %.2f", qc[1], qc[2]),
                                 sprintf("%.2f to %.2f", qc[2], qc[3]), sprintf("\u2265 %.2f", qc[3])),
                 n_cases = sapply(1:4, function(k) d[Q == k & case == 1, .N]),
                 n_controls = sapply(1:4, function(k) d[Q == k & case == 0, .N]),
                 or = c(ciq[paste0("QfQ", 1:3), "exp(coef)"], 1), lo = c(ciq[paste0("QfQ", 1:3), "lower .95"], NA),
                 hi = c(ciq[paste0("QfQ", 1:3), "upper .95"], NA), p = c(coq[paste0("QfQ", 1:3), "Pr(>|z|)"], NA),
                 or_per_quartile = exp(coef(mt)[["Qrev"]]), p_trend = p_trend)
fwrite(qt, file.path(OUT, "pheno_quartiles.csv"))
# 2. components
for (v in c("T10_trab_hu", "psm_lvol", "lung_perc15_hu", "psm_fat_pct", "psm_hu")) {
  dv <- copy(d)[, negz := -get(paste0("z_", v))]
  fit(dv[!is.na(negz)], "case ~ negz + strata(set)", paste0("Component: ", v, ", per 1 SD lower reserve"), "negz",
      "baseline")
}
# 3. musculoskeletal score, alone and adjusted for lung density
fit(d, "case ~ negS_msk + strata(set)", "Musculoskeletal score (bone + muscle), per 1 SD lower", "negS_msk", "baseline")
fit(d, "case ~ negS_msk + neg_zL + strata(set)", "Musculoskeletal score, adjusted for lung density z", "negS_msk",
    "baseline")
est[[length(est)]] <- copy(est[[length(est)]])[, note := "same model; lung density term reported in next row"]
fit(d, "case ~ negS_msk + neg_zL + strata(set)", "Lung density z, adjusted for musculoskeletal score", "neg_zL",
    "baseline")
# 4. five-measure score
d5 <- d[!is.na(S5)]
fit(d5, "case ~ negS5 + strata(set)", "Five-measure score (adds intramuscular fat, lean radiodensity), per 1 SD lower",
    "negS5", "baseline")
# 5-8. replication, lead time, subgroups
p_int <- c(series = lr_int(d, "series", "Series"), lead = lr_int(d, "lead3", "Lead time"),
           histo = lr_int(d, "histo", "Histology"), stage = lr_int(d, "stage2", "Stage"),
           detection = lr_int(d, "detection", "Detection"), sex = lr_int(d, "sex_s", "Sex"),
           smoke = lr_int(d, "smoke_s", "Smoking"))
# 7. latest pre-diagnostic scan
l <- build("last")[!is.na(S)]
l <- l[set %in% l[case == 1, set] & set %in% l[case == 0, set]][, negS := -S]
fit(l, "case ~ negS + strata(set)", "Composite score at latest pre-diagnostic scan, per 1 SD lower", "negS", "last")
cat(sprintf("latest scan: case median %.2f y before diagnosis\n", median(l[case == 1, -tau])))
# 9. matched-set concordance (case has lower score than its control)
pairs <- merge(d[case == 1, .(set, Sc = S)], d[case == 0, .(set, Sk = S)], by = "set", allow.cartesian = TRUE)
pairs[, w := fifelse(Sc < Sk, 1, fifelse(Sc == Sk, 0.5, 0))]
ps <- pairs[, .(w = sum(w), n = .N), by = set]
cidx <- sum(ps$w) / sum(ps$n)
set.seed(20261011)
bs <- replicate(1000, { k <- sample(nrow(ps), replace = TRUE); sum(ps$w[k]) / sum(ps$n[k]) })
est[[length(est) + 1]] <- data.table(analysis = "Matched-set concordance (case lower than its control)",
                                     population = "baseline", term = "C", or = cidx, lo = quantile(bs, .025),
                                     hi = quantile(bs, .975), p = NA_real_, n_cases = sum(d$case),
                                     n_controls = sum(1 - d$case), n_sets = uniqueN(d$set))

E <- rbindlist(est, fill = TRUE)
fwrite(E, file.path(OUT, "pheno_estimates.csv"))
cat("\nAll estimates (OR per 1 SD lower reserve unless stated):\n")
print(E[, .(analysis, OR = sprintf("%.2f (%.2f-%.2f)", or, lo, hi), p = signif(p, 2), p_int = signif(p_interaction, 2),
            n_cases)], nrows = 200)
cat("\nQuartiles (control cut-points; Q4 = highest reserve = reference):\n"); print(qt)
cat("\nInteraction P:\n"); print(signif(p_int, 2))
sink()
