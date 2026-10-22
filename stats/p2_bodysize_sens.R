# Paper 2: body-size sensitivity analysis of docs/paper2/P2_05_bodysize_sensitivity_prespec.md (2026-10-18).
# Refuses to run if that plan changed after its hash was recorded (results/p2/bodysize_prespec_sha256.txt).
#   Rscript stats/p2_bodysize_sens.R
# Outputs: results/p2/bodysize_sens_estimates.csv, bodysize_sens_summary.txt
.libPaths(c("env/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(survival) })
PLAN <- "docs/paper2/P2_05_bodysize_sensitivity_prespec.md"
rec <- strsplit(readLines("results/p2/bodysize_prespec_sha256.txt")[1], " +")[[1]][1]
now <- strsplit(system(paste("sha256sum", PLAN), intern = TRUE), " +")[[1]][1]
if (rec != now) stop("analysis plan changed after its hash was recorded")
sink("results/p2/bodysize_sens_summary.txt", split = TRUE)
cat("Paper 2 body-size sensitivity analysis, run", format(Sys.time()), "\nplan", PLAN, "sha256", now, "\n\n")

met <- fread("results/p2/p2_paraspinal.csv", select = c("case_id", "psm_nslices", "slice_mm", "lung_laa950_pct"))
met <- unique(met, by = "case_id")[, slab_mm := psm_nslices * slice_mm]

COMP <- c(T10_trab_hu = 1, psm_lvol = 1, lung_perc15_hu = 1)
safe_predict <- function(fm, nd) {
  nd <- copy(nd)
  lv <- fm$xlevels$vendor0
  nd[!(as.character(vendor0) %in% lv), vendor0 := names(which.max(table(fm$model$vendor0)))]
  predict(fm, newdata = nd)
}
# baseline z-scores as in stats/p2_phenotype.R and stats/p2_review_sens.R (original control sets and QC)
long <- fread("results/p2/partA_long.csv", colClasses = list(character = c("pid", "set")))
norm <- fread("results/p2/partA_normative_model.csv")
long[, age0c := age0 - 62]
long[, `:=`(sex = factor(sex), smoke = factor(smoke), vendor0 = factor(vendor0))]
long <- merge(long, met, by = "case_id", all.x = TRUE)
gam <- function(v) norm[measure == v & term == "noise_c", estimate]
for (v in names(COMP)) long[, (paste0("y_", v)) := get(v) - gam(v) * noise_c]
long[, `:=`(y_laa = -log(lung_laa950_pct + 0.1), y_slab = -log(slab_mm), lslab = log(slab_mm))]
s <- long[order(pid, t)][, .SD[1], by = pid]
zfit <- function(y, extra = "") {
  rhs <- paste("age0c + sex + smoke + vendor0", extra)
  ok <- !is.na(s[[y]]) & !is.na(s$lslab)
  fm <- lm(as.formula(paste(y, "~", rhs)), data = s[group == "control" & ok])
  z <- rep(NA_real_, nrow(s))
  z[ok] <- (s[[y]][ok] - safe_predict(fm, s[ok])) / sigma(fm)
  z
}
s[, `:=`(z_bone = zfit("y_T10_trab_hu"), z_vol = zfit("y_psm_lvol"), z_vol_bs = zfit("y_psm_lvol", "+ lslab"),
         z_lung = zfit("y_lung_perc15_hu"), z_laa = zfit("y_laa", "+ noise_c"), z_slab = zfit("y_slab"))]
std <- function(x) x / sd(x[s$group == "control"], na.rm = TRUE)
s[, `:=`(S = std(rowMeans(cbind(z_bone, z_vol, z_lung))), S_bs = std(rowMeans(cbind(z_bone, z_vol_bs, z_lung))),
         M = std(rowMeans(cbind(z_bone, z_vol))), M_bs = std(rowMeans(cbind(z_bone, z_vol_bs))))]
s[, case := as.integer(group == "case")]
d <- s[!is.na(S_bs) & !is.na(z_laa)]
d <- d[set %in% d[case == 1, set] & set %in% d[case == 0, set]]

est <- list()
fit <- function(f, term, analysis, label) {
  m <- clogit(as.formula(f), data = d)
  ci <- summary(m)$conf.int; co <- summary(m)$coefficients
  est[[length(est) + 1]] <<- data.table(analysis = analysis, model = label, or = ci[term, "exp(coef)"],
                                       lo = ci[term, "lower .95"], hi = ci[term, "upper .95"],
                                       p = co[term, "Pr(>|z|)"], n_cases = sum(d$case), n_controls = sum(1 - d$case))
}
d[, `:=`(nz_slab = -z_slab, nz_vol = -z_vol, nz_vol_bs = -z_vol_bs, nS = -S, nS_bs = -S_bs, nM = -M, nM_bs = -M_bs,
         nz_lung = -z_lung, nz_laa = -z_laa)]
A1 <- "8. Body size: slab length"; A2 <- "8. Body size: muscle volume adjusted for slab length"
fit("case ~ nz_slab + strata(set)", "nz_slab", A1, "T6-T11 slab length, per 1 SD shorter")
fit("case ~ nz_vol + strata(set)", "nz_vol", A1, "Muscle volume (as published), same sets")
fit("case ~ nS + strata(set)", "nS", A1, "Primary score (as published), same sets")
fit("case ~ nz_vol_bs + strata(set)", "nz_vol_bs", A2, "Muscle volume, per 1 SD lower")
fit("case ~ nS_bs + strata(set)", "nS_bs", A2, "Primary score, per 1 SD lower")
fit("case ~ nM_bs + strata(set)", "nM_bs", A2, "Bone + muscle score alone")
fit("case ~ nM_bs + nz_lung + strata(set)", "nM_bs", A2, "Bone + muscle score adjusted for lung density z")
fit("case ~ nM_bs + nz_laa + strata(set)", "nM_bs", A2, "Bone + muscle score adjusted for %LAA-950 z")
E <- rbindlist(est)
fwrite(E, "results/p2/bodysize_sens_estimates.csv")
cat("first-scan slab length (mm), mean by group and sex\n")
print(dcast(s, sex ~ group, value.var = "slab_mm", fun.aggregate = function(x) round(mean(x, na.rm = TRUE), 1)))
cat("\ncorrelation of z_vol with z_vol_bs in controls:", round(cor(d[case == 0, z_vol], d[case == 0, z_vol_bs]), 3), "\n\n")
print(E[, .(analysis, model, est = sprintf("%.3f (%.3f-%.3f)", or, lo, hi), p = signif(p, 2), n_cases, n_controls)])
sink()
