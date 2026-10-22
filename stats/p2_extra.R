# Paper 2 secondary analyses for the figures (run after stats/p2_partA.R; reads its exported analysis data).
#  1. subgroup analysis of excess ageing (primary measure): case vs matched controls within case-defined subgroups,
#     all windows pooled and the <1-year window, with an interaction likelihood-ratio test
#  2. level offsets (case minus matched controls) within stage, histology and detection subgroups, every measure
#  3. normative centile charts (quantile regression on age, by sex; controls' first scan) and each case's centile
# Outputs: results/p2/extra_subgroups.csv, extra_offsets_by_subgroup.csv, extra_centile_curves.csv,
#          extra_case_centiles.csv
.libPaths(c("env/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(lme4); library(quantreg); library(splines) })
OUT <- "results/p2"
long <- fread(file.path(OUT, "partA_long.csv"), colClasses = list(character = c("pid", "set")))
iv <- fread(file.path(OUT, "partA_intervals.csv"), colClasses = list(character = c("pid", "set")))
prim <- strsplit(readLines(file.path(OUT, "partA_primary_exposure.txt"))[1], " ")[[1]][2]
cas <- rbindlist(lapply(c("p2", "p2b"), function(s)
  fread(sprintf("data/nlst_serial/%s_cases.csv", s), colClasses = list(character = c("pid", "rounds")))[, resected := s == "p2"]),
  fill = TRUE)
prs <- fread("data/nlst/nlst_prsn_public.csv", select = c("pid", "can_scr"), colClasses = list(character = "pid"))
cas <- merge(cas, prs, by = "pid", all.x = TRUE)
grp <- cas[, .(set = pid,
               age_dx = fifelse(age + candx_days / 365.25 >= 65, "\u226565 y at diagnosis", "<65 y at diagnosis"),
               sex_s = fifelse(gender == 1, "Men", "Women"),
               smoke_s = fifelse(cigsmok == 1, "Current smoker", "Former smoker"),
               stage = fcase(de_stag_7thed %in% c(110, 120), "I", de_stag_7thed %in% c(210, 220), "II",
                             de_stag_7thed %in% c(310, 320), "III", de_stag_7thed == 400, "IV", default = NA_character_),
               histo = fcase(lc_morph %in% c(8041, 8042, 8043, 8044, 8045), "Small cell",
                             lc_morph %in% c(8070, 8071, 8072, 8073, 8074), "Squamous",
                             lc_morph %in% c(8140, 8250, 8251, 8252, 8253, 8254, 8255, 8260, 8310, 8323, 8480, 8481,
                                             8490, 8550, 8551, 8560), "Adenocarcinoma", default = "Other"),
               detection = fcase(can_scr == 1, "Screen-detected", can_scr %in% c(2, 3, 4), "Interval/post-screening",
                                 default = NA_character_),
               resected = fifelse(resected, "Resected", "Not resected"))]
grp[, stage2 := fcase(stage %in% c("I", "II"), "Stage I\u2013II", stage %in% c("III", "IV"), "Stage III\u2013IV",
                      default = NA_character_)]
vend <- unique(long[group == "case", .(set = pid, vendor = fcase(vendor0 == "GE", "GE", vendor0 == "SE", "Siemens",
                                                              default = "Philips/Toshiba"))])
grp <- merge(grp, vend, by = "set", all.x = TRUE)
ctl <- lmerControl(check.conv.singular = "ignore")

# ---------------- 1. subgroups
FACT <- c(age_dx = "Age", sex_s = "Sex", smoke_s = "Smoking", stage2 = "Stage", histo = "Histology",
          detection = "Detection", resected = "Resection", vendor = "Scanner vendor")
x <- merge(iv[measure == prim], grp, by = "set")
x[, group := factor(group, c("control", "case"))]
sub_rows <- list()
for (f in names(FACT)) {
  xf <- x[!is.na(get(f))]
  xf[, sub := factor(get(f))]
  for (win in c("all", "<1 y")) {
    xw <- if (win == "all") xf else xf[group == "control" | end_tau > -1]
    for (lv in levels(xf$sub)) {
      m <- lmer(z ~ group + (1 | set) + (1 | pid), data = xw[sub == lv], control = ctl)
      co <- summary(m)$coefficients["groupcase", ]
      sub_rows[[paste(f, win, lv)]] <- data.table(factor = FACT[[f]], subgroup = lv, window = win,
                                                 n_cases = uniqueN(xw[sub == lv & group == "case", pid]),
                                                 estimate = co[1], lo = co[1] - 1.96 * co[2], hi = co[1] + 1.96 * co[2])
    }
    m1 <- lmer(z ~ group * sub + (1 | set) + (1 | pid), data = xw, REML = FALSE, control = ctl)
    m0 <- lmer(z ~ group + sub + (1 | set) + (1 | pid), data = xw, REML = FALSE, control = ctl)
    pint <- anova(m0, m1)[2, "Pr(>Chisq)"]
    for (lv in levels(xf$sub)) sub_rows[[paste(f, win, lv)]][, p_interaction := pint]
  }
}
sg <- rbindlist(sub_rows)
fwrite(sg, file.path(OUT, "extra_subgroups.csv"))
cat("subgroups (", prim, "): interaction P range", sprintf("%.2f-%.2f", min(sg$p_interaction), max(sg$p_interaction)), "\n")

# ---------------- 2. level offsets within tumour subgroups
long[, age0c := age0 - 62]
long[, group := factor(group, c("control", "case"))]
MEAS <- c("psm_fat_pct", "psm_hu", "psm_lvol", "T10_trab_hu", "lung_perc15_hu")
off_rows <- list()
lg <- merge(long, grp[, .(set, stage, histo, detection)], by = "set")
for (f in c("stage", "histo", "detection")) {
  for (lv in sort(unique(na.omit(lg[[f]])))) {
    for (v in MEAS) {
      dv <- lg[get(f) == lv & !is.na(get(v))]
      m <- lmer(as.formula(paste(v, "~ group * tau + age0c + sex + smoke + vendor + noise_c + (1 + tau | pid) + (1 | set)")),
                data = dv, control = lmerControl(optimizer = "bobyqa", check.conv.singular = "ignore"))
      co <- summary(m)$coefficients["groupcase", ]
      off_rows[[paste(f, lv, v)]] <- data.table(factor = f, subgroup = lv, measure = v,
                                               n_cases = uniqueN(dv[group == "case", pid]),
                                               estimate = co[1], lo = co[1] - 1.96 * co[2], hi = co[1] + 1.96 * co[2])
    }
  }
}
fwrite(rbindlist(off_rows), file.path(OUT, "extra_offsets_by_subgroup.csv"))
cat("offsets by subgroup written\n")

# ---------------- 3. centile charts (controls' first scan) and case centiles
first <- long[order(pid, t)][, .SD[1], by = pid]
first[, age_scan := age0 + t]   # t = 0 on the first scan; age at entry is the age at T0
first[, vol_cm3 := exp(psm_lvol / 100)]
CM <- c("psm_fat_pct", "psm_hu", "vol_cm3", "T10_trab_hu")
taus <- c(0.05, 0.25, 0.5, 0.75, 0.95)
grid <- seq(55, 75, by = 0.5)
fine <- seq(0.01, 0.99, by = 0.01)
curves <- list(); cc <- list()
for (v in CM) for (sx in c("male", "female")) {
  ref <- first[group == "control" & sex == sx & !is.na(get(v))]
  fit <- rq(as.formula(paste(v, "~ bs(age_scan, df = 3, Boundary.knots = c(54, 76))")), tau = fine, data = ref)
  pg <- predict(fit, newdata = data.frame(age_scan = grid))
  for (tq in taus) curves[[paste(v, sx, tq)]] <- data.table(measure = v, sex = sx, tau = tq, age = grid,
                                                             value = pg[, which.min(abs(fine - tq))])
  for (g in c("case", "control")) {
    q <- first[group == g & sex == sx & !is.na(get(v))]
    pq <- predict(fit, newdata = data.frame(age_scan = q$age_scan))
    pq <- t(apply(pq, 1, sort))                       # monotone in tau (removes quantile crossing)
    cen <- vapply(seq_len(nrow(q)), function(i) mean(pq[i, ] <= q[[v]][i]), 0)
    cc[[paste(v, sx, g)]] <- data.table(pid = q$pid, group = g, sex = sx, measure = v, centile = cen)
  }
}
fwrite(rbindlist(curves), file.path(OUT, "extra_centile_curves.csv"))
cen <- rbindlist(cc)
fwrite(cen, file.path(OUT, "extra_case_centiles.csv"))
print(cen[, .(n = .N, median_centile = median(centile), below_25th = mean(centile < 0.25)), by = .(measure, group)])
