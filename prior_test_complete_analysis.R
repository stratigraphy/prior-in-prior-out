# Prior in, prior out: the complete analysis of the 56 published estimates of the precession constant k ####
# This one script makes every run, number, table and figure of the paper, in this order:
#   Part 1. The prior test: TimeOptBMCMC (R package astrochron) ten times for each of the 56 estimates, 560 runs: nine
#           runs with the prior mean of k at -2 to +2 SD from the TimeOptB (Farhat et al., 2022) value in steps of
#           0.5 SD, and one run with that value and a five times wider SD (the wide run).
#   Part 2. The readings: the 560 runs combined into the posterior of every run and the reading of every estimate,
#           with the geological floor.
#   Part 3. Table 1 and Table S1 as printed.
#   Part 4. The curve of k against age, with modBchron() of the R package WaverideR.
#   Part 5. The authors' own methods (TimeOptMCMC, or TimeOptBMCMC for the Hakatai Shale) with their own prior of k,
#           and with that prior shifted by -2 and +2 SD: 66 runs for the 22 estimates not published with TimeOptB,
#           and Table S2.
#   Part 6. The synthetic test: 24 synthetic records with a known k, four TimeOptBMCMC runs each (96 runs), and
#           Figure S2.
#   Part 7. Figures 1 to 3, S1, S3 and S4.
# Folders, relative to the project folder (the folder of prior_test.Rproj and this script; open prior_test.Rproj in
# RStudio, then source this script or run it part by part):
#   data/                   inputs (see below); Part 6 writes the 24 synthetic records to data/synthetic/
#   output/figures/         all figures (PNG, 300 dpi)
#   output/tables/          all tables as printed (CSV): Table 1, Table S1, Table S2 and the synthetic test
#   output/saved_results/   all results: the saved runs (runs_prior_test/, runs_authors/ and runs_synthetic/), the
#                           readings, the curves of k against age and the Bchron models (RDS and CSV)
# Files needed in data/:
#   prior_test_runs.csv  one line per run of Part 1: the record, its study and published method, its age, the
#                        published value, the TimeOptB prior of k at its age, the data file, the sedimentation-rate
#                        range, the TimeOptBMCMC options, the random seed, the prior mean and SD of k of this run and
#                        the name of the output file.
#   authors_runs.csv     one line per run of Part 5: the record, the method, the data file, the sedimentation-rate
#                        range, the TimeOptBMCMC options, the mean frequencies g1-g5 of TimeOptMCMC, the authors' prior
#                        of k, the shift, the prior mean and SD of k of this run, the seed, the number of samples and
#                        the name of the output file.
#   <record>.txt         the series of each record (depth in m, value), as prepared from the original study.
#   AstroGeo22.txt       the Farhat et al. (2022) solution.
# How to run:
#   All runs are made on one core, one run after another, and every run is saved in output/saved_results/ as soon as
#   it is finished. A run that is already saved is not run again, so the script can be stopped and
#   started again at any time, and each record (each run_record() or run_authors_record() line) can be run on its
#   own, also on different computers: copy the saved .rds files into the same folders afterwards. Parts 2 to 4 and 7
#   need only the saved runs of Part 1, so with these runs in output/saved_results/runs_prior_test/ they can be run
#   without Parts 5 and 6.
#   One TimeOptBMCMC run takes 15 minutes to several hours (median about 50 minutes); the 560 runs of Part 1 take
#   about four weeks on one core. A TimeOptMCMC run of Part 5 (100,000 samples) takes several hours to several days.
#   Settings: 50,000 samples per TimeOptBMCMC run, burn-in 20% (Part 5: the first quarter of the chain).

# Packages and settings ####
library(astrochron)
library(WaverideR)
# folders (relative to the project folder)
data_dir <- "data"
if (!dir.exists(data_dir)) {
  stop("Folder data/ not found. Open prior_test.Rproj, or set the working directory to the folder that holds ",
       "this script, and run it again.")
}
figures_dir <- file.path("output", "figures")
tables_dir <- file.path("output", "tables")
results_dir <- file.path("output", "saved_results")
runs_dir <- file.path(results_dir, "runs_prior_test")
authors_dir <- file.path(results_dir, "runs_authors")
synthetic_dir <- file.path(results_dir, "runs_synthetic")
for (d in c(figures_dir, tables_dir, results_dir, runs_dir, authors_dir, synthetic_dir,
            file.path(data_dir, "synthetic"))) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}
# run settings
runs <- read.csv(file.path(data_dir, "prior_test_runs.csv"), stringsAsFactors = FALSE)
runs$shift <- as.character(runs$shift)
authors_runs <- read.csv(file.path(data_dir, "authors_runs.csv"), stringsAsFactors = FALSE)

# Saving the figures ####
# Every figure is saved twice in output/figures/, at 300 dpi: as PNG and as TIFF with LZW compression (TIFF is one of
# the formats AGU asks for with a revised manuscript).
open_figure <- function(name, device, width, height) {
  if (device == "png") {
    png(file.path(figures_dir, paste0(name, ".png")), width = width, height = height, units = "in", res = 300)
  } else {
    tiff(file.path(figures_dir, paste0(name, ".tif")), width = width, height = height, units = "in", res = 300,
         compression = "lzw")
  }
}
n_samples <- 50000
k_now <- 50.475838
pm <- "\u00b1"
# synthetic records: known k, age, sedimentation rate (cm/kyr) and sample step (m), as Sopelana L (Wu et al., 2024)
k_true <- 52.15
synthetic_age_Ma <- 70
sed_rate <- 3
step_m <- 0.0731

# The TimeOptB prior of k and how a run changes it ####
# Where the prior comes from:
#   timeOptBMCMC() has no argument for the prior of k. Inside timeOptBMCMC(), the priors of all astronomical
#   parameters (the secular frequencies g1-g5, s1-s6 and the precession constant k) are taken from a function inside
#   the astrochron package, .getpriorpdfastropara(tGa, mvopt). It returns the TimeOptB prior (Farhat et al., 2022) at
#   the age tGa (in billion years) as a list; the elements used here are:
#     mu      the mean of each parameter;
#     sigma   the standard deviation of each parameter;
#     snmean  the mean used by timeOptBMCMC() for the resolution check (mean +/- 2 sigma); for k it equals mu.
#   The position of k in these vectors depends on mvopt (which s-terms are included) and is given by the astrochron
#   function .getmvindex("k", mvopt).
# How a run gets its own prior of k:
#   Before a run, use_k_prior() builds a new version of .getpriorpdfastropara() that first calls the original
#   function and then overwrites only the three values of k (mu, snmean and sigma) with the prior mean and SD of k
#   of this run. All other priors stay exactly as in astrochron. assignInNamespace() then puts this new version inside
#   the astrochron package, under the same name, so that timeOptBMCMC() uses it. Nothing is changed on disk: the
#   change lasts only for this R session.
#   After the run, use_original_prior() puts the original function back with assignInNamespace(), so that every run
#   starts from the unchanged astrochron prior. If a run stops with an error, the changed prior stays in place until
#   the next run sets its own prior, until use_original_prior() is run, or until R is restarted.

# A copy of the original astrochron function, kept to build the new prior and to put it back after each run
original_prior <- astrochron:::.getpriorpdfastropara

# Put a prior of k with mean k_mean and standard deviation k_sd inside astrochron
use_k_prior <- function(k_mean, k_sd, mvopt) {
  # position of k in the vectors of the prior
  ik <- astrochron:::.getmvindex("k", mvopt)
  # the new version: the original TimeOptB prior, with only the mean and SD of k replaced
  new_prior <- function(tGa, mvopt) {
    q <- original_prior(tGa, mvopt)
    q$mu[ik] <- k_mean
    q$snmean[ik] <- k_mean
    q$sigma[ik] <- k_sd
    q
  }
  # replace .getpriorpdfastropara inside the astrochron package by the new version
  assignInNamespace(".getpriorpdfastropara", new_prior, ns = "astrochron")
}

# Put the original TimeOptB prior back inside astrochron
use_original_prior <- function() {
  assignInNamespace(".getpriorpdfastropara", original_prior, ns = "astrochron")
}

# One TimeOptBMCMC run (one line of prior_test_runs.csv) ####
# The random seed of each run is given in prior_test_runs.csv. It is 1 for all runs except Wumishan at +2 SD, where
# the run with seed 1 stopped with an error inside timeOptBMCMC() and the run was made with seed 2.
run_bmcmc <- function(i) {
  run <- runs[i, ]
  file <- file.path(runs_dir, run$output_file)
  # a complete run that is already saved is not run again
  if (file.exists(file) && length(readRDS(file)$k) == n_samples) {
    cat(run$output_file, "is already saved\n")
    return(invisible(file))
  }
  dat <- read.table(run$data_file, header = TRUE)
  cat(sprintf("%s, shift %s: prior of k %.2f +/- %.2f arcsec/yr, seed %d\n", run$record, run$shift, run$prior_mean,
              run$prior_sd, run$seed))
  t0 <- Sys.time()
  # 1. the prior of k of this run
  use_k_prior(run$prior_mean, run$prior_sd, run$mvopt)
  # 2. the run
  set.seed(run$seed)
  result <- timeOptBMCMC(dat, age = run$age_Ma, env = TRUE, mvopt = run$mvopt, sedmin = run$sedmin,
                         sedmax = run$sedmax, ftol = run$ftol, roll = 10^7, nsamples = n_samples, detrend = TRUE,
                         genplot = FALSE, check = run$check, verbose = FALSE)
  # 3. the original prior back
  use_original_prior()
  # 4. save the samples of k and their log-likelihood, with the settings of the run
  summary <- data.frame(record = run$record, shift = run$shift, prior_mean = run$prior_mean,
                        prior_sd = run$prior_sd, n_samples = n_samples, seed = run$seed,
                        minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(list(summary = summary, k = result$mcmcres$mv[, "k"], loglik = result$mcmcres$mv[, "logpdf"]), file)
  cat(sprintf("  saved %s after %.1f minutes\n", run$output_file, summary$minutes))
  invisible(file)
}

# The ten runs of one record ####
run_record <- function(record) {
  lines <- which(runs$record == record)
  cat("\n", record, " (", runs$study[lines[1]], ", ", runs$age_Ma[lines[1]], " Ma): ", length(lines), " runs\n",
      sep = "")
  for (i in lines) {
    run_bmcmc(i)
  }
}

# Runs with the authors' own methods (Table S2) ####
# SD of the prior of g1-g5, and the step sizes of the sampler (as in the TimeOptMCMC example of astrochron)
g_sd <- c(0.125, 0.015, 0.15, 0.15, 0.00002)
step_sizes <- c(0.2, 0.25, 0.35, 0.35, 0.8, 0.85, 0.6, 0.35, 0.9, 0.35, 0.9) / 20

# One run with the authors' method (one line of authors_runs.csv) ####
# TimeOptMCMC checks before it starts that the shortest astronomical period stays resolved at k = kAve - 2 kSd. When
# kAve - 2 kSd is zero or negative (Dales Gorge at -2 SD: 45 - 2 x 25 < 0), this check cannot be made and
# timeOptMCMC() stops; such a run is not made, and the pull in Table S2 is computed from the runs at 0 and +2 SD.
run_authors <- function(i) {
  run <- authors_runs[i, ]
  file <- file.path(authors_dir, run$output_file)
  # a complete run that is already saved is not run again
  if (file.exists(file) && length(readRDS(file)$k) == run$n_samples) {
    cat(run$output_file, "is already saved\n")
    return(invisible(file))
  }
  if (run$method == "TimeOptMCMC" && run$prior_mean - 2 * run$prior_sd <= 0) {
    cat(run$output_file, "is not run: the prior mean minus 2 SD of k is not positive\n")
    return(invisible(NA))
  }
  dat <- read.table(run$data_file, header = TRUE)
  cat(sprintf("%s (%s), shift %s: prior of k %.2f +/- %.2f arcsec/yr\n", run$record, run$method, run$shift,
              run$prior_mean, run$prior_sd))
  t0 <- Sys.time()
  set.seed(run$seed)
  if (run$method == "TimeOptMCMC") {
    result <- timeOptMCMC(dat, iopt = 1, sedmin = run$sedmin, sedmax = run$sedmax,
                          sedstart = (run$sedmin + run$sedmax) / 2,
                          gAve = c(run$g1, run$g2, run$g3, run$g4, run$g5), gSd = g_sd, gstart = c(-1, -1, -1, -1, -1),
                          kAve = run$prior_mean, kSd = run$prior_sd, kstart = -1, rhomin = 0, rhomax = 0.9999,
                          rhostart = -1, sigmamin = NULL, sigmamax = NULL, sigmastart = 1, nsamples = run$n_samples,
                          epsilon = step_sizes, ran = TRUE, burnin = -1, genplot = FALSE, verbose = FALSE)
    k <- result[, "k"]
  } else {
    use_k_prior(run$prior_mean, run$prior_sd, run$mvopt)
    result <- timeOptBMCMC(dat, age = run$age_Ma, env = TRUE, mvopt = run$mvopt, sedmin = run$sedmin,
                           sedmax = run$sedmax, ftol = run$ftol, roll = 10^7, nsamples = run$n_samples,
                           detrend = TRUE, genplot = FALSE, check = TRUE, verbose = FALSE)
    use_original_prior()
    k <- result$mcmcres$mv[, "k"]
  }
  summary <- data.frame(record = run$record, method = run$method, shift = run$shift, prior_mean = run$prior_mean,
                        prior_sd = run$prior_sd, n_samples = run$n_samples, seed = run$seed,
                        minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(list(summary = summary, k = k), file)
  cat(sprintf("  saved %s after %.1f minutes\n", run$output_file, summary$minutes))
  invisible(file)
}

# The three runs of one record with the authors' method ####
run_authors_record <- function(record) {
  lines <- which(authors_runs$record == record)
  cat("\n", record, " (", authors_runs$study[lines[1]], ", ", authors_runs$method[lines[1]], "): ",
      length(lines), " runs\n", sep = "")
  for (i in lines) {
    run_authors(i)
  }
}

# The TimeOptB prior at 70 Ma: mean frequencies g1-g5 and the SD of the prior of k ####
prior_70 <- astrochron:::.getpriorpdfastropara(synthetic_age_Ma / 1000, 0)
g <- prior_70$mu[1:5]
ik <- astrochron:::.getmvindex("k", 0)
k_prior_sd <- prior_70$sigma[ik]

# The 24 synthetic records of the synthetic test (Text S2) ####
make_record <- function(case, snr, seed) {
  length_m <- 25
  if (case == "short") {
    length_m <- 14
  }
  # time in years, one value per sample
  t_yr <- seq(0, length_m / (sed_rate / 100) * 1000, by = step_m / (sed_rate / 100) * 1000)
  set.seed(seed)
  # climatic precession: the four main k + g terms (amplitudes of La2004), random phases
  amp_p <- c(0.0189, 0.0163, 0.0130, 0.0088)
  g_p <- c(g[5], g[2], g[4], g[3])
  precession <- rep(0, length(t_yr))
  for (i in 1:4) {
    precession <- precession + amp_p[i] * sin(2 * pi * (k_true + g_p[i]) * t_yr / 1296000 + runif(1, 0, 2 * pi))
  }
  # obliquity: k + s3, k + s4 and k + s6, random phases
  s_terms <- c(-26.33, -18.85, -17.75)
  amp_o <- c(0.3, 1, 0.6)
  obliquity <- rep(0, length(t_yr))
  for (i in 1:3) {
    obliquity <- obliquity + amp_o[i] * sin(2 * pi * (k_true + s_terms[i]) * t_yr / 1296000 + runif(1, 0, 2 * pi))
  }
  # eccentricity as the envelope of precession
  envelope <- hilbert(cbind(t_yr, precession), genplot = FALSE, verbose = FALSE)[, 2]
  signal <- as.numeric(scale(envelope)) + as.numeric(scale(precession)) + 0.5 * as.numeric(scale(obliquity))
  signal <- as.numeric(scale(signal))
  # AR(1) noise with the same standard deviation as the signal, divided by the signal-to-noise ratio
  noise <- as.numeric(scale(arima.sim(list(ar = 0.6), n = length(signal))))
  value <- signal + noise / snr
  data.frame(depth = t_yr / 1000 * (sed_rate / 100), value = value)
}
synthetic_records <- NULL
for (case in c("long", "short")) {
  for (snr in c(0.3, 0.5, 1, 2)) {
    for (seed in 1:3) {
      name <- paste0(case, "_snr", snr, "_s", seed)
      file <- file.path(data_dir, "synthetic", paste0(name, ".txt"))
      if (!file.exists(file)) {
        rec <- make_record(case, snr, seed)
        writeLines(c("depth value", sprintf("%.17g %.17g", rec$depth, rec$value)), file)
      }
      synthetic_records <- rbind(synthetic_records, data.frame(record = name, case = case, snr = snr, seed = seed,
                                                               data_file = file))
    }
  }
}

# The four runs of each synthetic record ####
synthetic_runs <- NULL
for (i in seq_len(nrow(synthetic_records))) {
  for (shift in c("-2", "0", "2", "w")) {
    prior_mean <- k_true
    prior_sd <- k_prior_sd
    if (shift == "w") {
      prior_sd <- 5 * k_prior_sd
    } else {
      prior_mean <- k_true + as.numeric(shift) * k_prior_sd
    }
    synthetic_runs <- rbind(synthetic_runs, data.frame(synthetic_records[i, ], shift = shift, prior_mean = prior_mean,
                                                       prior_sd = prior_sd,
                                                       output_file = paste0(synthetic_records$record[i], "_", shift,
                                                                            ".rds")))
  }
}
write.csv(synthetic_runs, file.path(results_dir, "synthetic_runs.csv"), row.names = FALSE)

# One TimeOptBMCMC run of a synthetic record ####
run_synthetic <- function(i) {
  run <- synthetic_runs[i, ]
  file <- file.path(synthetic_dir, run$output_file)
  # a complete run that is already saved is not run again
  if (file.exists(file) && length(readRDS(file)$k) == n_samples) {
    cat(run$output_file, "is already saved\n")
    return(invisible(file))
  }
  dat <- read.table(run$data_file, header = TRUE)
  cat(sprintf("%s, shift %s: prior of k %.2f +/- %.2f arcsec/yr\n", run$record, run$shift, run$prior_mean,
              run$prior_sd))
  t0 <- Sys.time()
  use_k_prior(run$prior_mean, run$prior_sd, 0)
  set.seed(1)
  result <- timeOptBMCMC(dat, age = synthetic_age_Ma, env = TRUE, mvopt = 0, sedmin = 2, sedmax = 4, ftol = 0.005,
                         roll = 10^7, nsamples = n_samples, detrend = TRUE, genplot = FALSE, check = TRUE,
                         verbose = FALSE)
  use_original_prior()
  summary <- data.frame(record = run$record, shift = run$shift, prior_mean = run$prior_mean,
                        prior_sd = run$prior_sd, n_samples = n_samples, seed = 1,
                        minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(list(summary = summary, k = result$mcmcres$mv[, "k"]), file)
  cat(sprintf("  saved %s after %.1f minutes\n", run$output_file, summary$minutes))
  invisible(file)
}

# Option: all 721 runs at once on a multicore computer (switched off) ####
# Instead of one core, the runs of Parts 1, 5 and 6 can be spread over several cores of one computer. This is switched
# off: to use it, remove the # at the start of the lines below (not of this explanation). Because a saved run is not
# run again, the run_record(), run_authors_record() and run_synthetic() lines of Parts 1, 5 and 6 then only report
# that every run is already saved, and the rest of the script continues as usual.
# WARNING: every core runs its own chain at full load. The computer can become slow or freeze while other large
# programs are open, and a computer with little memory can run out of it (each core needs a few hundred MB). Keep
# n_cores well below the number of cores of the computer and close other large programs.
#
# library(parallel)
# library(foreach)
# library(doParallel)
# n_cores <- 4
# warning("All 721 runs are started on ", n_cores, " cores. The computer may become slow or freeze while other ",
#         "large programs are open, and it can run out of memory. Close other programs before you continue.")
# cl <- makeCluster(n_cores, type = "PSOCK")
# registerDoParallel(cl)
# clusterExport(cl, ls(globalenv()))
# jobs <- rbind(data.frame(kind = "prior test", i = seq_len(nrow(runs))),
#               data.frame(kind = "authors", i = seq_len(nrow(authors_runs))),
#               data.frame(kind = "synthetic", i = seq_len(nrow(synthetic_runs))))
# done <- foreach(j = seq_len(nrow(jobs)), .packages = "astrochron") %dopar% {
#   if (jobs$kind[j] == "prior test") {
#     run_bmcmc(jobs$i[j])
#   } else if (jobs$kind[j] == "authors") {
#     run_authors(jobs$i[j])
#   } else {
#     run_synthetic(jobs$i[j])
#   }
# }
# stopCluster(cl)

# Part 1: the 560 runs of the prior test, one record after another, on one core ####
# 1 Dales Gorge Member (Zhou et al., 2022), 2465 Ma ####
run_record("DalesGorgeZhou")

# 2 Dagushi (Fan et al., 2026), 1787 Ma ####
run_record("Dagushi")

# 3 Chuanlinggou (Shi et al., 2023), 1640 Ma ####
run_record("ChuanlinggouShi")

# 4 Chuanlinggou (Zhou et al., 2024), 1634 Ma ####
run_record("Chuanlinggou")

# 5 Wumishan (Zhou et al., 2024), 1480 Ma ####
run_record("Wumishan")

# 6 Xiamaling (Meyers & Malinverno, 2018), 1400 Ma ####
run_record("Xiamaling")

# 7 Hakatai Shale (Lantink et al., 2026), 1250 Ma ####
run_record("Hakatai")

# 8 Yemahe (Zhou et al., 2024), 1215 Ma ####
run_record("Yemahe")

# 9 Datangpo (Zhou et al., 2022), 655 Ma ####
run_record("Datangpo")

# 10 Datangpo (Wu et al., 2024), 650 Ma ####
run_record("DatangpoWu")

# 11 Fiq (Wu et al., 2024), 640 Ma ####
run_record("Fiq")

# 12 Muraykhah (Wu et al., 2024), 600 Ma ####
run_record("Muraykhah")

# 13 Wuhe (Wu et al., 2024), 570 Ma ####
run_record("Wuhe")

# 14 Zhengjiatang (Huang et al., 2024), 570 Ma ####
run_record("Zhengjiatang")

# 15 Qiongzhusi (Zhang et al., 2022), 526 Ma ####
run_record("Qiongzhusi")

# 16 Luoyixi (Wu et al., 2024), 500 Ma ####
run_record("Luoyixi")

# 17 Alum (Wu et al., 2024), 492 Ma ####
run_record("Alum")

# 18 Liangjiashan (Huang et al., 2024), 470 Ma ####
# the wide run has check = FALSE in prior_test_runs.csv: the resolution check leaves no sedimentation rate
run_record("Liangjiashan")

# 19 Liangjiashan (Wu et al., 2024), 470 Ma ####
run_record("LiangjiashanWu")

# 20 Pingliang (Zhou et al., 2022), 455 Ma ####
run_record("Pingliang")

# 21 Vaureal (Wu et al., 2024), 450 Ma ####
run_record("Vaureal")

# 22 Anticosti (Huang et al., 2024), 448 Ma ####
run_record("Anticosti")

# 23 Wanhe S1 (Wu et al., 2024), 445 Ma ####
run_record("WanheS1")

# 24 Wanhe S2 (Wu et al., 2024), 445 Ma ####
run_record("WanheS2")

# 25 Paslek (Wu et al., 2024), 435 Ma ####
run_record("Paslek")

# 26 Pozar-CS (Huang et al., 2024), 410 Ma ####
run_record("PozarCS")

# 27 Pod Barrandovem (Wu et al., 2024), 405 Ma ####
run_record("PodBarrandovem")

# 28 H32 (Huang et al., 2024), 375 Ma ####
run_record("H32")

# 29 Sinsin (Wu et al., 2024), 374 Ma ####
run_record("Sinsin")

# 30 Lali (Wu et al., 2024), 360 Ma ####
run_record("Lali")

# 31 Luokun S1 (Wu et al., 2024), 325 Ma ####
run_record("LuokunS1")

# 32 Luokun S3 (Wu et al., 2024), 310 Ma ####
run_record("LuokunS3")

# 33 Naqing (Wu et al., 2024), 308 Ma ####
run_record("Naqing")

# 34 Fengcheng CAL (Wu et al., 2024), 300 Ma ####
run_record("FengchengCAL")

# 35 Fengcheng GR (Wu et al., 2024), 300 Ma ####
run_record("FengchengGR")

# 36 Estelle Mine (Wu et al., 2024), 297 Ma ####
run_record("EstelleMine")

# 37 Lucaogou (Huang et al., 2024), 290 Ma ####
run_record("Lucaogou")

# 38 Wujiaping (Zhou et al., 2022), 259 Ma ####
run_record("Wujiaping")

# 39 Baikouquan (Wu et al., 2024), 250 Ma ####
run_record("Baikouquan")

# 40 Guandao (Huang et al., 2024), 245 Ma ####
run_record("Guandao")

# 41 Chang (Wu et al., 2024), 240 Ma ####
run_record("Chang")

# 42 Newark (Wu et al., 2024), 220 Ma ####
run_record("Newark")

# 43 St Audrie's Bay A (Wu et al., 2024), 210 Ma ####
run_record("StAudrieA")

# 44 St Audrie's Bay L (Wu et al., 2024), 210 Ma ####
run_record("StAudrieL")

# 45 Tanba (Wu et al., 2024), 203 Ma ####
run_record("Tanba")

# 46 Mochras (Wu et al., 2024), 190 Ma ####
run_record("Mochras")

# 47 Torre de Busi (Wu et al., 2024), 160 Ma ####
run_record("TorreDeBusi")

# 48 Shahezi (Wu et al., 2024), 135 Ma ####
run_record("Shahezi")

# 49 Serre Chaitieu (Wu et al., 2024), 119 Ma ####
run_record("SerreChaitieu")

# 50 Songke (Wu et al., 2024), 82 Ma ####
run_record("Songke")

# 51 Sopelana L (Wu et al., 2024), 70 Ma ####
run_record("SopelanaL")

# 52 Sopelana MS (Wu et al., 2024), 70 Ma ####
run_record("SopelanaMS")

# 53 Walvis Ridge (Meyers & Malinverno, 2018), 55 Ma ####
run_record("WalvisRidge")

# 54 Walvis Ridge (Wu et al., 2024), 55 Ma ####
run_record("WalvisWu")

# 55 Newfoundland Ridge (De Vleeschouwer et al., 2023), 41 Ma ####
run_record("NewfoundlandRidge")

# 56 Shahejie (Ge et al., 2026), 39.23 Ma ####
run_record("Shahejie")

# Part 2: the readings: combine the 560 runs into the posterior of every run and the reading of every estimate ####
# Nothing is reweighted or interpolated: every number is a statistic of the samples of one run.
records <- runs[!duplicated(runs$record), c("record", "study", "published_method", "age_Ma", "published_mean",
                                            "published_sd", "farhat_mean", "farhat_sd")]
nine_shifts <- c("-2", "-1.5", "-1", "-0.5", "0", "0.5", "1", "1.5", "2")

# Samples of one run after the burn-in (20%), and whether the two halves of the chain agree ####
read_run <- function(record, shift) {
  file <- file.path(runs_dir, paste0(record, "_bmcmc_", shift, ".rds"))
  if (!file.exists(file)) {
    return(NULL)
  }
  z <- readRDS(file)
  n <- length(z$k)
  k <- z$k[(floor(0.2 * n) + 1):n]
  h <- floor(length(k) / 2)
  mixed <- abs(mean(k[1:h]) - mean(k[(h + 1):length(k)])) <= max(sd(k), 0.05)
  list(k = k, prior_mean = z$summary$prior_mean, prior_sd = z$summary$prior_sd, mixed = mixed)
}

# Posterior of every run of one record inside the plausible range (lowest to highest plausible k) ####
record_stats <- function(record, lowest, highest) {
  out <- NULL
  for (s in c(nine_shifts, "w")) {
    z <- read_run(record, s)
    if (is.null(z)) {
      next
    }
    keep <- z$k >= lowest & z$k <= highest
    kk <- z$k[keep]
    k_mean <- NA
    k_sd <- NA
    k_median <- NA
    if (length(kk) > 1) {
      k_mean <- mean(kk)
      k_sd <- sd(kk)
      k_median <- median(kk)
    }
    out <- rbind(out, data.frame(record = record, shift = s, prior_mean = z$prior_mean, prior_sd = z$prior_sd,
                                 n = length(kk), share_removed = 1 - mean(keep), mean = k_mean, sd = k_sd,
                                 median = k_median, mixed = z$mixed))
  }
  out
}

# Reading of one record from its nine runs ####
# A chain stays in the solution where it settles. The runs are therefore grouped into solutions: sorted by posterior
# median, two neighbouring runs belong to the same solution when their medians differ by less than 3 posterior SD of
# either run. The main solution is the one most runs land in.
#   X  more than half of the samples of the centred run, or all samples of the wide run, are dropped;
#   P  the run medians rise with the prior mean (Spearman 0.9 or more) across several solutions or with a pull over
#      all runs of 0.5 or more, or the runs of the main solution have a pull of 0.9 or more, or of 0.7 or more with a
#      shrinkage of 0.65 or more;
#   D  at least seven of the nine runs land in the main solution; strong when all nine do, the pull is below 0.5 and
#      the shrinkage below 0.65, weak otherwise;
#   none  no solution holds seven or more runs.
# pull: least-squares slope of the posterior mean against the prior mean (runs of the main solution); shrinkage:
# median posterior SD of these runs / prior SD.
read_record <- function(nine, prior_sd, wide) {
  nine <- nine[!is.na(nine$mean), ]
  if (nrow(nine) < 3) {
    return(list(pull = NA, shrink = NA, n_sol = NA, main_share = NA, main_k = NA, main_sd = NA, unmixed = NA,
                reading = "incomplete"))
  }
  nine <- nine[order(nine$prior_mean), ]
  o <- order(nine$median)
  cl <- integer(nrow(nine))
  cl[o[1]] <- 1
  for (j in seq_along(o)[-1]) {
    a1 <- o[j - 1]
    a2 <- o[j]
    gap <- nine$median[a2] - nine$median[a1]
    if (gap > 3 * max(nine$sd[a1], nine$sd[a2])) {
      cl[a2] <- cl[a1] + 1
    } else {
      cl[a2] <- cl[a1]
    }
  }
  n_sol <- max(cl)
  runs_per_solution <- table(cl)
  main <- as.integer(names(runs_per_solution)[which.max(runs_per_solution)])
  m <- nine[cl == main, ]
  main_share <- nrow(m) / nrow(nine)
  pull <- NA
  if (nrow(m) >= 3) {
    pull <- unname(coef(lm(mean ~ prior_mean, data = m))[2])
  }
  pull_all <- unname(coef(lm(mean ~ prior_mean, data = nine))[2])
  shrink <- median(m$sd) / prior_sd
  rho <- suppressWarnings(cor(nine$prior_mean, nine$median, method = "spearman"))
  centre <- nine[nine$shift == "0", ]
  no_wide <- !is.null(wide) && wide$share_removed > 0.99
  if ((nrow(centre) == 1 && centre$share_removed > 0.5) || no_wide) {
    reading <- "X"
  } else if (!is.na(rho) && rho >= 0.9 && (n_sol >= 2 || pull_all >= 0.5)) {
    reading <- "P"
  } else if (nrow(m) >= 7 && !is.na(pull) && (pull >= 0.9 || (pull >= 0.7 && shrink >= 0.65))) {
    reading <- "P"
  } else if (nrow(m) >= 7) {
    if (main_share == 1 && pull < 0.5 && shrink < 0.65) {
      reading <- "D strong"
    } else {
      reading <- "D weak"
    }
  } else {
    reading <- "none"
  }
  main_sd <- median(m$sd)
  if (nrow(m) > 1) {
    main_sd <- sqrt(median(m$sd)^2 + sd(m$mean)^2)
  }
  if (is.na(pull)) {
    pull <- pull_all
  }
  list(pull = pull, shrink = shrink, n_sol = n_sol, main_share = main_share, main_k = mean(m$mean),
       main_sd = main_sd, unmixed = sum(!nine$mixed), reading = reading)
}

# Reading of all 56 records with a given lower limit of k for each record ####
read_all <- function(lowest_k) {
  table1 <- NULL
  table_runs <- NULL
  for (i in seq_len(nrow(records))) {
    r <- records$record[i]
    highest <- records$farhat_mean[i] + 5 * records$farhat_sd[i]
    lowest <- max(k_now, lowest_k[r])
    rs <- record_stats(r, lowest, highest)
    table_runs <- rbind(table_runs, rs)
    nine <- rs[rs$shift != "w", ]
    wide <- rs[rs$shift == "w", ]
    if (nrow(wide) == 0) {
      wide <- NULL
    }
    rd <- read_record(nine, records$farhat_sd[i], wide)
    wide_mean <- NA
    wide_sd <- NA
    if (!is.null(wide)) {
      wide_mean <- wide$mean
      wide_sd <- wide$sd
    }
    table1 <- rbind(table1, data.frame(record = r, study = records$study[i],
                                       published_method = records$published_method[i], age_Ma = records$age_Ma[i],
                                       published_mean = records$published_mean[i],
                                       published_sd = records$published_sd[i],
                                       farhat_mean = records$farhat_mean[i], farhat_sd = records$farhat_sd[i],
                                       floor = lowest_k[r], n_runs = sum(!is.na(nine$mean)), wide_mean = wide_mean,
                                       wide_sd = wide_sd, pull = rd$pull, shrink = rd$shrink, n_sol = rd$n_sol,
                                       main_share = rd$main_share, main_k = rd$main_k, main_sd = rd$main_sd,
                                       unmixed = rd$unmixed, reading = rd$reading))
  }
  list(table1 = table1, runs = table_runs)
}

# Step 1: readings inside the plausible window only, to find the anchors of the geological floor ####
# Anchors are estimates read as set by the data (strong) whose main solution agrees with the published value within
# 2 SD; the SD of an anchor includes a method term of 2% of k. Naqing is not used: it lies above Fengcheng GR, which is
# only 8 Myr younger, and the tidal history does not allow such a rise.
no_floor <- rep(k_now, nrow(records))
names(no_floor) <- records$record
step1 <- read_all(no_floor)$table1
step1$sd_m <- sqrt(step1$main_sd^2 + (0.02 * step1$main_k)^2)
step1$agree <- abs(step1$main_k - step1$published_mean) <= 2 * sqrt(step1$sd_m^2 + step1$published_sd^2)
anchors <- step1[step1$reading == "D strong" & step1$agree & step1$record != "Naqing", ]
cat("\nanchors of the geological floor:", paste(anchors$record, collapse = ", "), "\n")

# Step 2: the geological floor of every record and the final readings ####
# k only increases back in time, so a record cannot have a lower k than a younger anchor: the floor at an age is the
# highest value of (anchor k - 3 SD) over all younger anchors.
floor_k <- no_floor
for (i in seq_len(nrow(records))) {
  younger <- anchors[anchors$age_Ma < records$age_Ma[i], ]
  if (nrow(younger) > 0) {
    floor_k[records$record[i]] <- max(k_now, max(younger$main_k - 3 * younger$sd_m))
  }
}
final <- read_all(floor_k)
table1 <- final$table1
table_runs <- final$runs

# Step 3: a solution at the model value (M) ####
# For a record read P or X: the runs, ordered by prior mean, are cut into segments where the median jumps by more
# than 1.5 times the prior step and 1.5 times the smallest posterior SD; segments whose medians lie within 3
# posterior SD of each other are merged. A segment counts when it holds at least three runs, its median lies within
# 1.5 prior SD of the Farhat et al. (2022) value, and its median rises by less than 0.5 per unit of prior mean; of
# several, the one closest to the model value is used. The data then hold a solution at the model value, but nothing
# in the inversion picks it.
table1$model_k <- NA
table1$model_sd <- NA
table1$model_runs <- NA
for (i in which(table1$reading %in% c("P", "X"))) {
  r <- table1$record[i]
  q <- table_runs[table_runs$record == r & table_runs$shift != "w" & !is.na(table_runs$mean), ]
  if (nrow(q) < 3) {
    next
  }
  q <- q[order(q$prior_mean), ]
  seg <- integer(nrow(q))
  seg[1] <- 1
  for (j in seq_len(nrow(q))[-1]) {
    dm <- q$median[j] - q$median[j - 1]
    dp <- q$prior_mean[j] - q$prior_mean[j - 1]
    if (abs(dm) > 1.5 * dp && abs(dm) > 1.5 * min(q$sd[j], q$sd[j - 1])) {
      seg[j] <- seg[j - 1] + 1
    } else {
      seg[j] <- seg[j - 1]
    }
  }
  seg_median <- numeric(max(seg))
  seg_sd <- numeric(max(seg))
  for (g in seq_along(seg_median)) {
    seg_median[g] <- median(q$median[seg == g])
    seg_sd[g] <- median(q$sd[seg == g])
  }
  grp <- seq_along(seg_median)
  for (g1 in seq_along(seg_median)) {
    for (g2 in seq_along(seg_median)) {
      if (g2 > g1 && abs(seg_median[g1] - seg_median[g2]) < 3 * min(seg_sd[g1], seg_sd[g2])) {
        grp[grp == grp[g2]] <- grp[g1]
      }
    }
  }
  q$sol <- grp[seg]
  if (length(unique(q$sol)) < 2) {
    next
  }
  m0 <- table1$farhat_mean[i]
  s0 <- table1$farhat_sd[i]
  best <- Inf
  for (g in unique(q$sol)) {
    m <- q[q$sol == g, ]
    if (nrow(m) < 3 || abs(median(m$median) - m0) > 1.5 * s0) {
      next
    }
    slope <- unname(coef(lm(median ~ prior_mean, data = m))[2])
    if (slope < 0.5 && abs(median(m$median) - m0) < best) {
      best <- abs(median(m$median) - m0)
      table1$reading[i] <- "M"
      table1$model_k[i] <- mean(m$median)
      table1$model_sd[i] <- sqrt(median(m$sd)^2 + sd(m$median)^2)
      table1$model_runs[i] <- nrow(m)
      table1$main_k[i] <- table1$model_k[i]
      table1$main_sd[i] <- table1$model_sd[i]
      table1$main_share[i] <- nrow(m) / nrow(q)
    }
  }
}

# Step 4: readings set by a stated observation (O) ####
# Hakatai Shale (read P by the rules above): one solution lies at the model value, and Lantink et al. (2026) chose it
# for a stated observation (the obliquity-scale peak at about 87 cm fits s6 + k). Its k is that of the runs in the
# solution that holds the published value.
i <- which(table1$record == "Hakatai")
q <- table_runs[table_runs$record == "Hakatai" & table_runs$shift != "w" & !is.na(table_runs$mean), ]
j0 <- which.min(abs(q$median - table1$published_mean[i]))
near <- abs(q$median - q$median[j0]) < 3 * q$sd[j0]
table1$reading[i] <- "O"
table1$main_k[i] <- mean(q$mean[near])
table1$main_sd[i] <- median(q$sd[near])
if (sum(near) > 1) {
  table1$main_sd[i] <- sqrt(median(q$sd[near])^2 + sd(q$mean[near])^2)
}
table1$main_share[i] <- sum(near) / nrow(q)
# Dales Gorge Member (read M by the rules above): the Joffre Member directly above it gives k = 108.6 +/- 8.5 arcsec/yr
# (2 SD; Lantink et al., 2022). With the floor rule of the anchors (mean - 3 SD, SD including the 2% method term),
# Joffre sets a lower limit of 94.3 arcsec/yr; the k of the estimate is that of the runs above this limit.
joffre_mean <- 108.6
joffre_sd <- 8.5 / 2
joffre_low <- joffre_mean - 3 * sqrt(joffre_sd^2 + (0.02 * joffre_mean)^2)
i <- which(table1$record == "DalesGorgeZhou")
q <- table_runs[table_runs$record == "DalesGorgeZhou" & table_runs$shift != "w" & !is.na(table_runs$mean), ]
above <- q$median >= joffre_low
table1$reading[i] <- "O"
table1$main_k[i] <- mean(q$mean[above])
table1$main_sd[i] <- median(q$sd[above])
if (sum(above) > 1) {
  table1$main_sd[i] <- sqrt(median(q$sd[above])^2 + sd(q$mean[above])^2)
}
table1$main_share[i] <- sum(above) / nrow(q)

# Results ####
write.csv(table1, file.path(results_dir, "readings.csv"), row.names = FALSE)
write.csv(table_runs, file.path(results_dir, "posterior_of_every_run.csv"), row.names = FALSE)
saveRDS(list(readings = table1, posterior_of_every_run = table_runs, anchors = anchors, floor = floor_k),
        file.path(results_dir, "readings.rds"))
print(table1[, c("record", "age_Ma", "published_mean", "floor", "main_k", "main_sd", "wide_mean", "wide_sd", "pull",
                 "shrink", "reading")], digits = 3, row.names = FALSE)
print(table(table1$reading))

# Part 3: Table 1 and Table S1 as printed ####
# Names of the records as printed ####
printed_name <- c(PozarCS = "Po\u017e\u00e1r-CS", DalesGorgeZhou = "Dales Gorge", ChuanlinggouShi = "Chuanlinggou",
                  DatangpoWu = "Datangpo", LiangjiashanWu = "Liangjiashan", WanheS1 = "Wanhe S1",
                  WanheS2 = "Wanhe S2", Paslek = "Pas\u0142ek", PodBarrandovem = "Pod Barrandovem",
                  LuokunS1 = "Luokun S1", LuokunS3 = "Luokun S3", FengchengCAL = "Fengcheng CAL",
                  FengchengGR = "Fengcheng GR", EstelleMine = "Estelle Mine", StAudrieA = "St Audrie's Bay A",
                  StAudrieL = "St Audrie's Bay L", TorreDeBusi = "Torre de Busi", SerreChaitieu = "Serre Chaitieu",
                  SopelanaL = "Sopelana L", SopelanaMS = "Sopelana MS", WalvisRidge = "Walvis Ridge",
                  WalvisWu = "Walvis Ridge", NewfoundlandRidge = "Newfoundland Ridge", Vaureal = "Vaur\u00e9al")
name_of <- function(record) {
  if (record %in% names(printed_name)) {
    return(printed_name[[record]])
  }
  record
}

# Readings as printed ####
printed_reading <- c("D strong" = "\u2713 D", "D weak" = "\u2713 D (w)", "O" = "\u2713 O", "M" = "~ M",
                     "P" = "\u2717 P", "X" = "\u2717 X", "none" = "-")

# Mean +/- SD as printed, "-" when there is no value ####
mean_sd <- function(m, s, digits) {
  if (is.na(m)) {
    return("-")
  }
  paste(formatC(m, format = "f", digits = digits), pm, formatC(s, format = "f", digits = digits))
}

# Table 1 as printed ####
t1 <- NULL
for (i in seq_len(nrow(table1))) {
  r <- table1[i, ]
  runs_in_main <- round(r$main_share * r$n_runs)
  t1 <- rbind(t1, data.frame(Record = name_of(r$record), Study = r$study, Age_Ma = round(r$age_Ma),
                             Published_k = mean_sd(r$published_mean, r$published_sd, 2),
                             Wide_prior = mean_sd(r$wide_mean, r$wide_sd, 1),
                             Main_solution = paste0(mean_sd(r$main_k, r$main_sd, 1), " (", runs_in_main, "/",
                                                    r$n_runs, ")"),
                             Pull = formatC(r$pull, format = "f", digits = 2),
                             Shrink = formatC(r$shrink, format = "f", digits = 2),
                             Reading = printed_reading[[r$reading]]))
}
write.csv(t1, file.path(tables_dir, "Table1.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Table S1 as printed ####
shifts <- c("-2", "-1.5", "-1", "-0.5", "0", "0.5", "1", "1.5", "2", "w")
s1 <- NULL
for (i in seq_len(nrow(table1))) {
  r <- table1[i, ]
  line <- data.frame(Record = name_of(r$record), Floor = formatC(r$floor, format = "f", digits = 1))
  for (s in shifts) {
    q <- table_runs[table_runs$record == r$record & table_runs$shift == s, ]
    value <- "-"
    if (nrow(q) == 1) {
      value <- mean_sd(q$mean, q$sd, 1)
    }
    line[[paste0("run_", s)]] <- value
  }
  line$Reading <- printed_reading[[r$reading]]
  s1 <- rbind(s1, line)
}
write.csv(s1, file.path(tables_dir, "Table_S1.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(t1, row.names = FALSE)

# Part 4: the curve of k against age (Figures S3 and S4) ####
# The curve is built from the estimates read as set by the data (D) or by an observation (O) in Part 2, with the
# monotone Bayesian age-depth model of Bchron (Haslett & Parnell, 2008; modified by Trayler et al., 2020), as function
# modBchron() of the R package WaverideR (0.6.0).
# The curve is made twice: from all these estimates (Figures S3 and S4) and from the strong ones only (Figure S3).
# Why an age-depth model:
#   k only increases back in time, just as age only increases with depth in a sedimentary section. modBchron() is
#   therefore used with "age" = k (arcsec/yr) and "position" = 5000 - age (Ma), so that the youngest estimate lies
#   highest in the "section" and has the lowest k. The model gives a curve that can only rise with age.
# How an estimate enters:
#   k = k of the main solution of the runs; SD = the SD of that solution combined with a method term of 2% of k
#   (analyses of the same record by different studies differ by 2-5%). A weak reading (D weak, and the two O
#   readings) enters with twice that SD. Each estimate also gets an outlier part: a uniform distribution over its
#   plausible range (present-day k to twice the Farhat et al. (2022) value at its age), which holds 10% of the weight
#   of a strong estimate and 25% of a weak one, so that one estimate cannot force the curve.
#   Fixed points: the present-day k (La2004/AstroGeo22) and the present rate from lunar laser ranging (Williams &
#   Boggs, 2016), entered at 10 Ma. The formation of the Moon (4425 Ma) gives k of 2000 arcsec/yr or more; it enters
#   only as the upper bound of the model.
# Output in output/saved_results/: kcurve_all.csv and kcurve_strong.csv (age, median and 68/95% bands of k, every
# 10 Myr), and the Bchron models kcurve_all.rds and kcurve_strong.rds.
# Time: about 7 minutes per curve.
n_mc <- 20000
n_copy <- 20
position_top <- 5000

# The Farhat et al. (2022) solution (AstroGeo22), for the plausible range of each estimate ####
astrogeo <- read.table(file.path(data_dir, "AstroGeo22.txt"), skip = 1)
astrogeo_age_Ma <- 1000 * astrogeo[, 1]
astrogeo_k <- astrogeo[, 11]

# The inputs of the curve ####
used <- table1[table1$reading %in% c("D strong", "D weak", "O"), ]
inputs <- data.frame(id = used$record, age_Ma = used$age_Ma, age_unc_Ma = 0.01 * used$age_Ma, k = used$main_k,
                     k_sd = sqrt(used$main_sd^2 + (0.02 * used$main_k)^2), p_out = 0.1, fixed = FALSE,
                     strong = used$reading == "D strong")
for (i in seq_len(nrow(inputs))) {
  if (!inputs$strong[i]) {
    inputs$k_sd[i] <- 2 * inputs$k_sd[i]
    inputs$p_out[i] <- 0.25
  }
}
present <- data.frame(id = "present", age_Ma = 0, age_unc_Ma = 0, k = 50.467727, k_sd = 0.001, p_out = 0,
                      fixed = TRUE, strong = TRUE)
llr <- data.frame(id = "LLR", age_Ma = 10, age_unc_Ma = 0, k = 50.702297, k_sd = 0.03, p_out = 0, fixed = TRUE,
                  strong = TRUE)
inputs <- rbind(present, llr, inputs)
inputs <- inputs[order(inputs$age_Ma), ]
k_moon <- 21000 - 19000

# One curve from a set of inputs ####
# Every estimate is entered as n_copy "dates" at its position: copies of its Gaussian distribution in proportion to
# its weight (1 - outlier share) and copies of the uniform outlier distribution in proportion to the outlier share.
make_curve <- function(inp, output_name) {
  ages <- NULL
  sds <- NULL
  types <- NULL
  positions <- NULL
  thicknesses <- NULL
  ids <- NULL
  for (i in seq_len(nrow(inp))) {
    n_before <- length(ages)
    n_gauss <- max(1, round((1 - inp$p_out[i]) * n_copy))
    ages <- c(ages, rep(inp$k[i], n_gauss))
    sds <- c(sds, rep(inp$k_sd[i], n_gauss))
    types <- c(types, rep("G", n_gauss))
    if (!inp$fixed[i] && inp$p_out[i] > 0) {
      n_uniform <- max(1, round(inp$p_out[i] * n_copy))
      k_high <- 2 * approx(astrogeo_age_Ma, astrogeo_k, inp$age_Ma[i], rule = 2)$y
      ages <- c(ages, rep((k_now + k_high) / 2, n_uniform))
      sds <- c(sds, rep((k_high - k_now) / 2, n_uniform))
      types <- c(types, rep("U", n_uniform))
    }
    n_new <- length(ages) - n_before
    positions <- c(positions, rep(position_top - inp$age_Ma[i], n_new))
    thicknesses <- c(thicknesses, rep(inp$age_unc_Ma[i], n_new))
    ids <- c(ids, rep(inp$id[i], n_new))
  }
  o <- order(positions)
  set.seed(1)
  fit <- modBchron(ages = ages[o], ageSds = sds[o], positions = positions[o], positionThicknesses = thicknesses[o],
                   ids = ids[o], distTypes = types[o], MC = n_mc, burn = round(n_mc / 4),
                   predictPositions = position_top - seq(4540, 0, by = -10), truncateUp = 50,
                   truncateDown = k_moon, engine = "R")
  # samples of the curve after the burn-in: one row per sample, one column per predicted age
  m <- fit$model
  if (nrow(m) == length(fit$predictPositions)) {
    m <- t(m)
  }
  m <- m[(round(n_mc / 4) + 1):nrow(m), , drop = FALSE]
  bands <- NULL
  for (j in seq_len(ncol(m))) {
    bands <- rbind(bands, quantile(m[, j], c(0.025, 0.16, 0.5, 0.84, 0.975), na.rm = TRUE))
  }
  curve <- data.frame(age_Ma = position_top - fit$predictPositions, q025 = bands[, 1], q16 = bands[, 2],
                      median = bands[, 3], q84 = bands[, 4], q975 = bands[, 5])
  curve <- curve[order(curve$age_Ma), ]
  write.csv(curve, paste0(output_name, ".csv"), row.names = FALSE)
  saveRDS(fit, paste0(output_name, ".rds"))
  curve
}

# The curve from all estimates set by the data or by an observation ####
curve_all <- make_curve(inputs, file.path(results_dir, "kcurve_all"))
print(curve_all[curve_all$age_Ma %in% c(0, 350, 510, 1000, 1640, 2470), ], row.names = FALSE)

# The curve from the strong estimates only ####
curve_strong <- make_curve(inputs[inputs$strong, ], file.path(results_dir, "kcurve_strong"))
print(curve_strong[curve_strong$age_Ma %in% c(0, 350, 510, 1000), ], row.names = FALSE)

# Part 5: the 66 runs with the authors' own methods, one record after another, on one core ####
# 1 Dales Gorge Member (Zhou et al., 2022) ####
run_authors_record("DalesGorgeZhou")

# 2 Dagushi (Fan et al., 2026) ####
run_authors_record("Dagushi")

# 3 Chuanlinggou (Shi et al., 2023) ####
run_authors_record("ChuanlinggouShi")

# 4 Chuanlinggou (Zhou et al., 2024) ####
run_authors_record("Chuanlinggou")

# 5 Wumishan (Zhou et al., 2024) ####
run_authors_record("Wumishan")

# 6 Xiamaling (Meyers & Malinverno, 2018) ####
run_authors_record("Xiamaling")

# 7 Hakatai Shale (Lantink et al., 2026; TimeOptBMCMC) ####
run_authors_record("Hakatai")

# 8 Yemahe (Zhou et al., 2024) ####
run_authors_record("Yemahe")

# 9 Datangpo (Zhou et al., 2022) ####
run_authors_record("Datangpo")

# 10 Zhengjiatang (Huang et al., 2024) ####
run_authors_record("Zhengjiatang")

# 11 Qiongzhusi (Zhang et al., 2022) ####
run_authors_record("Qiongzhusi")

# 12 Liangjiashan (Huang et al., 2024) ####
run_authors_record("Liangjiashan")

# 13 Pingliang (Zhou et al., 2022) ####
run_authors_record("Pingliang")

# 14 Anticosti (Huang et al., 2024) ####
run_authors_record("Anticosti")

# 15 Pozar-CS (Huang et al., 2024) ####
run_authors_record("PozarCS")

# 16 H32 (Huang et al., 2024) ####
run_authors_record("H32")

# 17 Lucaogou (Huang et al., 2024) ####
run_authors_record("Lucaogou")

# 18 Wujiaping (Zhou et al., 2022) ####
run_authors_record("Wujiaping")

# 19 Guandao (Huang et al., 2024) ####
run_authors_record("Guandao")

# 20 Walvis Ridge (Meyers & Malinverno, 2018) ####
run_authors_record("WalvisRidge")

# 21 Newfoundland Ridge (De Vleeschouwer et al., 2023) ####
run_authors_record("NewfoundlandRidge")

# 22 Shahejie (Ge et al., 2026) ####
run_authors_record("Shahejie")

# Table S2 ####
# Authors' method: posterior mean +/- SD of the run at the centred prior (burn-in: the first quarter of the chain);
# pull = (mean of the highest run - mean of the lowest run) / (prior mean of the highest run - prior mean of the lowest
# run), from the runs present. *: the means of the two halves of the chain differ by more than one posterior SD
# (at least 0.05 arcsec/yr) in at least one run.
# TimeOptBMCMC: the same numbers from the runs of Part 1 at -2, 0 and +2 SD (after the plausible window and the
# geological floor; Part 2); if one of these runs has no sample left, the nearest run present is used.
runs_bmcmc <- table_runs
s2 <- NULL
for (record in unique(authors_runs$record)) {
  lines <- authors_runs[authors_runs$record == record, ]
  # authors' method
  present <- NULL
  for (j in seq_len(nrow(lines))) {
    file <- file.path(authors_dir, lines$output_file[j])
    if (!file.exists(file)) {
      next
    }
    k <- readRDS(file)$k
    k <- k[(length(k) %/% 4 + 1):length(k)]
    h <- length(k) %/% 2
    mixed <- abs(mean(k[1:h]) - mean(k[(h + 1):length(k)])) <= max(sd(k), 0.05)
    present <- rbind(present, data.frame(shift = lines$shift[j], prior_mean = lines$prior_mean[j], mean = mean(k),
                                         sd = sd(k), mixed = mixed))
  }
  authors_post <- "-"
  authors_pull <- "-"
  if (!is.null(present) && nrow(present) >= 2) {
    present <- present[order(present$shift), ]
    star <- ""
    if (!all(present$mixed)) {
      star <- "*"
    }
    centre <- present[present$shift == 0, ]
    if (nrow(centre) == 1) {
      authors_post <- paste0(sprintf("%.1f", centre$mean), " ", pm, " ", sprintf("%.1f", centre$sd), star)
    }
    n <- nrow(present)
    authors_pull <- sprintf("%.2f", (present$mean[n] - present$mean[1]) /
                              (present$prior_mean[n] - present$prior_mean[1]))
  }
  # TimeOptBMCMC runs of Part 1
  q <- runs_bmcmc[runs_bmcmc$record == record & runs_bmcmc$shift != "w" & !is.na(runs_bmcmc$mean), ]
  q$s <- as.numeric(q$shift)
  centre <- q[q$s == 0, ]
  lowest <- q[which.min(q$s), ]
  highest <- q[which.max(q$s), ]
  star <- ""
  if (!all(c(lowest$mixed, centre$mixed, highest$mixed))) {
    star <- "*"
  }
  bmcmc_post <- "-"
  if (nrow(centre) == 1) {
    bmcmc_post <- paste0(sprintf("%.1f", centre$mean), " ", pm, " ", sprintf("%.1f", centre$sd), star)
  }
  bmcmc_pull <- sprintf("%.2f", (highest$mean - lowest$mean) / (highest$prior_mean - lowest$prior_mean))
  t1 <- table1[table1$record == record, ]
  s2 <- rbind(s2, data.frame(record = record, age_Ma = round(t1$age_Ma), method = lines$method[1],
                             published = paste(sprintf("%.2f", t1$published_mean), pm,
                                               sprintf("%.2f", t1$published_sd)),
                             authors_post = authors_post, authors_pull = authors_pull, bmcmc_post = bmcmc_post,
                             bmcmc_pull = bmcmc_pull))
}
write.csv(s2, file.path(tables_dir, "Table_S2.csv"), row.names = FALSE, fileEncoding = "UTF-8")
authors_pull <- suppressWarnings(as.numeric(s2$authors_pull))
cat("authors' methods: pull above 0.7 in", sum(authors_pull > 0.7, na.rm = TRUE), "of", sum(!is.na(authors_pull)),
    "records; median pull", round(median(authors_pull, na.rm = TRUE), 2), "\n")
print(s2, row.names = FALSE)

# Part 6: the synthetic test: the 96 runs, one after another, on one core ####
for (i in seq_len(nrow(synthetic_runs))) {
  run_synthetic(i)
}

# Summary of the synthetic test and Figure S2 ####
# Posterior of each run after the burn-in (20%); samples below the present-day k are kept, so that a record whose
# data point below the present-day k shows it.
posterior_of <- function(record, shift) {
  file <- file.path(synthetic_dir, paste0(record, "_", shift, ".rds"))
  if (!file.exists(file)) {
    return(c(NA, NA))
  }
  k <- readRDS(file)$k
  k <- k[(floor(0.2 * length(k)) + 1):length(k)]
  c(mean(k), sd(k))
}
summary_table <- NULL
for (i in seq_len(nrow(synthetic_records))) {
  r <- synthetic_records$record[i]
  wide <- posterior_of(r, "w")
  centred <- posterior_of(r, "0")
  low <- posterior_of(r, "-2")
  high <- posterior_of(r, "2")
  pull <- (high[1] - low[1]) / (4 * k_prior_sd)
  summary_table <- rbind(summary_table, data.frame(synthetic_records[i, c("record", "case", "snr", "seed")],
                                                   wide_mean = wide[1], wide_sd = wide[2],
                                                   below_present = wide[1] < k_now, centred_mean = centred[1],
                                                   centred_sd = centred[2], pull = pull))
}
write.csv(summary_table, file.path(tables_dir, "Table_synthetic_test.csv"), row.names = FALSE)
print(summary_table, digits = 3, row.names = FALSE)

for (device in c("png", "tiff")) {
  open_figure("Figure_S2", device, 10, 3.6)
  par(mfrow = c(1, 3), mar = c(4, 4, 2.2, 0.5), mgp = c(2.3, 0.7, 0))
  snr_values <- sort(unique(summary_table$snr))
  x <- match(summary_table$snr, snr_values)
  x[summary_table$case == "long"] <- x[summary_table$case == "long"] - 0.12
  x[summary_table$case == "short"] <- x[summary_table$case == "short"] + 0.12
  symbol <- rep(16, nrow(summary_table))
  symbol[summary_table$case == "short"] <- 17
  colour <- rep("black", nrow(summary_table))
  colour[summary_table$below_present] <- "#C0392B"
  # (a) the run with the wide prior
  plot(x, summary_table$wide_mean, pch = symbol, col = colour, xaxt = "n", xlab = "Signal-to-noise ratio (SD)",
       ylab = "k, run with the wide prior (arcsec/yr)", main = "(a) Where do the data alone put k?",
       ylim = range(c(summary_table$wide_mean, k_now - 1, k_true + 1), na.rm = TRUE))
  axis(1, at = seq_along(snr_values), labels = snr_values)
  abline(h = k_true, lty = 2)
  abline(h = k_now, lty = 4, col = "grey40")
  legend("topleft", bty = "n", cex = 0.8, pch = c(16, 17, 16, NA, NA), lty = c(NA, NA, NA, 2, 4),
         col = c("black", "black", "#C0392B", "black", "grey40"),
         legend = c("25 m record", "14 m record", "Below present-day k", "True k", "Present-day k"))
  # (b) the pull
  plot(x, summary_table$pull, pch = symbol, xaxt = "n", ylim = c(0, 1.4), xlab = "Signal-to-noise ratio (SD)",
       ylab = "Pull (posterior shift / prior shift)", main = "(b) Does the posterior follow a shifted prior?")
  axis(1, at = seq_along(snr_values), labels = snr_values)
  abline(h = c(0, 1), lty = 3, col = "grey50")
  # (c) the run with the prior centred on the true value
  plot(x, summary_table$centred_mean, pch = symbol, xaxt = "n", xlab = "Signal-to-noise ratio (SD)",
       ylab = "k, prior centred on the true value (arcsec/yr)", main = "(c) Prior centred on the true value",
       ylim = range(c(summary_table$centred_mean, k_true - 1, k_true + 1), na.rm = TRUE))
  axis(1, at = seq_along(snr_values), labels = snr_values)
  abline(h = k_true, lty = 2)
  dev.off()
}

# Part 7: Figures 1 to 3, S1, S3 and S4 ####
# Nothing is reweighted or interpolated between runs: every run is drawn from its own samples.
# Names of the records in the panel titles ####
panel_name <- c(PozarCS = "Po\u017e\u00e1r-CS", DalesGorgeZhou = "Dales Gorge", ChuanlinggouShi = "Chuanlinggou (Shi)",
                Chuanlinggou = "Chuanlinggou (Zhou)", PodBarrandovem = "Pod Barrandovem",
                StAudrieA = "St Audrie's Bay A", StAudrieL = "St Audrie's Bay L", WanheS1 = "Wanhe S1",
                WanheS2 = "Wanhe S2", LuokunS1 = "Luokun S1", LuokunS3 = "Luokun S3", FengchengCAL = "Fengcheng CAL",
                FengchengGR = "Fengcheng GR", EstelleMine = "Estelle Mine", TorreDeBusi = "Torre de Busi",
                SerreChaitieu = "Serre Chaitieu", SopelanaL = "Sopelana L", SopelanaMS = "Sopelana MS",
                WalvisRidge = "Walvis Ridge (M&M)", WalvisWu = "Walvis Ridge (Wu)",
                NewfoundlandRidge = "Newfoundland Ridge", DatangpoWu = "Datangpo (Wu)", Datangpo = "Datangpo (Zhou)",
                LiangjiashanWu = "Liangjiashan (Wu)", Liangjiashan = "Liangjiashan (Huang)", Paslek = "Pas\u0142ek",
                Vaureal = "Vaur\u00e9al")
panel_name_of <- function(record) {
  if (record %in% names(panel_name)) {
    return(panel_name[[record]])
  }
  record
}

# Badges of the readings ####
badge_text <- c("D strong" = "\u2713 D", "D weak" = "\u2713 D (w)", "O" = "\u2713 O", "M" = "~ M", "P" = "\u2717 P",
                "X" = "\u2717 X", "none" = "-")
badge_colour <- c("D strong" = "#239B56", "D weak" = "#239B56", "O" = "#E67E22", "M" = "#F5B041", "P" = "#C0392B",
                  "X" = "#C0392B", "none" = "#B7950B")
draw_badge <- function(reading, cex_b) {
  label <- badge_text[[reading]]
  usr <- par("usr")
  bx <- usr[1] + 0.03 * diff(usr[1:2])
  by <- usr[4] - 0.03 * diff(usr[3:4])
  bw <- strwidth(label, cex = cex_b) + 0.04 * diff(usr[1:2])
  bh <- 1.6 * strheight(label, cex = cex_b)
  rect(bx, by - bh, bx + bw, by, col = badge_colour[[reading]], border = NA)
  text(bx + 0.02 * diff(usr[1:2]), by - bh / 2, label, adj = c(0, 0.5), cex = cex_b, col = "white", font = 2)
}

# Samples of k of one run after the burn-in (20%) ####
load_run <- function(record, shift) {
  file <- file.path(runs_dir, paste0(record, "_bmcmc_", shift, ".rds"))
  if (!file.exists(file)) {
    return(NULL)
  }
  z <- readRDS(file)
  if (!is.numeric(z$k) || length(z$k) < 1000) {
    return(NULL)
  }
  n <- length(z$k)
  z$k[(floor(0.2 * n) + 1):n]
}

# Axis range of the panel of one record (same scale on both axes) ####
panel_range <- function(q, runs) {
  m0 <- q$farhat_mean
  s0 <- q$farhat_sd
  limits <- NULL
  for (k in runs) {
    if (length(k) > 50) {
      limits <- c(limits, quantile(k, c(q$probs_low, q$probs_high)))
    }
  }
  rng <- range(c(m0 - 2.3 * s0, m0 + 2.3 * s0, limits, q$published_mean))
  c(max(rng[1], max(k_now, m0 - 5 * s0) - 0.5), min(rng[2], m0 + 5 * s0))
}

# Samples of the nine runs and of the wide run of one record inside its plausible range ####
record_samples <- function(q) {
  highest <- q$farhat_mean + 5 * q$farhat_sd
  lowest <- max(k_now, q$floor)
  runs <- list()
  for (s in nine_shifts) {
    k <- load_run(q$record, s)
    if (!is.null(k)) {
      runs[[s]] <- k[k >= lowest & k <= highest]
    }
  }
  wide <- numeric(0)
  k <- load_run(q$record, "w")
  if (!is.null(k)) {
    wide <- k[k >= lowest & k <= highest]
  }
  list(runs = runs, wide = wide)
}

# The Farhat et al. (2022) solution (AstroGeo22) and the TimeOptB prior of k ####
astrogeo <- read.table(file.path(data_dir, "AstroGeo22.txt"), skip = 1)
names(astrogeo) <- c("T", "a", "amin", "amax", "lod", "lodmin", "lodmax", "obl", "oblmin", "oblmax", "pre", "premin",
                     "premax")
astrogeo$age <- 1000 * astrogeo$T
astrogeo$lo <- pmin(astrogeo$premin, astrogeo$premax)
astrogeo$hi <- pmax(astrogeo$premin, astrogeo$premax)
astrogeo_2700 <- astrogeo[astrogeo$age <= 2700, ]
prior_age <- seq(0, 2600, by = 10)
prior_mean <- rep(NA, length(prior_age))
prior_sd <- rep(NA, length(prior_age))
for (i in seq_along(prior_age)) {
  p <- astrochron:::.getpriorpdfastropara(max(prior_age[i], 1) / 1000, 2)
  prior_mean[i] <- p$mu[6]
  prior_sd[i] <- p$sigma[6]
}
prior_mean[1] <- k_now

# Broken time axis of Figures 1 and S3: 0-700 Ma over 60% of the width, 700-2600 Ma over 40% ####
t_break <- 700
tx <- function(a) {
  ifelse(a <= t_break, 0.6 * a / t_break, 0.6 + 0.4 * (a - t_break) / (2600 - t_break))
}
axis_time <- function(break_marks) {
  at1 <- seq(0, 700, by = 100)
  at2 <- seq(1000, 2600, by = 400)
  axis(1, at = tx(c(at1, at2)), labels = c(at1, at2))
  if (break_marks) {
    for (dx in c(-0.006, 0.006)) {
      segments(tx(t_break) + dx - 0.004, par("usr")[3] - 0.02 * diff(par("usr")[3:4]), tx(t_break) + dx + 0.004,
               par("usr")[3] + 0.02 * diff(par("usr")[3:4]), xpd = NA, lwd = 1.2)
    }
  }
  abline(v = tx(t_break), lty = 3, col = "grey60")
}

# Geological time scale above Figures 1 and S3 (eras and periods, ICS colours from WaverideR) ####
gts <- GTS_info[GTS_info$type %in% c("Era", "Period"), ]
gts <- gts[gts$top < 2600, ]
gts$bot <- pmin(gts$bot, 2600)
gts_bar <- function(ylo, yhi) {
  h <- (yhi - ylo) / 2
  for (lev in c("Period", "Era")) {
    g <- gts[gts$type == lev, ]
    y0 <- ylo
    if (lev == "Era") {
      y0 <- ylo + h
    }
    for (i in seq_len(nrow(g))) {
      x0 <- tx(g$top[i])
      x1 <- tx(g$bot[i])
      rect(x0, y0, x1, y0 + h, col = rgb(g$R[i], g$G[i], g$B[i], maxColorValue = 255), border = "grey30",
           lwd = 0.4, xpd = NA)
      label <- sub("_Period", "", g$name[i])
      if (strwidth(label, cex = 0.55) > 0.95 * (x1 - x0)) {
        label <- substr(label, 1, 2)
      }
      if (strwidth(label, cex = 0.55) < 0.95 * (x1 - x0)) {
        text((x0 + x1) / 2, y0 + h / 2, label, cex = 0.55, xpd = NA)
      }
    }
  }
}

# Background of Figures 1 and S3: TimeOptB prior, Farhat et al. (2022) solution and present-day k ####
time_background <- function() {
  axis(2, at = c(50, 55, 60, 70, 80, 90, 100, 110, 120), las = 1)
  polygon(tx(c(prior_age, rev(prior_age))), c(prior_mean - 2 * prior_sd, rev(prior_mean + 2 * prior_sd)),
          col = "#E5E7E9", border = NA)
  lines(tx(prior_age), prior_mean, lty = 2, col = "grey45")
  polygon(tx(c(astrogeo_2700$age, rev(astrogeo_2700$age))), c(astrogeo_2700$lo, rev(astrogeo_2700$hi)),
          col = "#2C3E5040", border = NA)
  lines(tx(astrogeo_2700$age), astrogeo_2700$pre, lwd = 1.4, col = "black")
}

# Figure 1: the published estimates ####
method_pch <- c("TimeOptMCMC" = 17, "TimeOptB" = 16, "TimeOptBMCMC" = 15)
method_colour <- c("TimeOptMCMC" = "#7D3C98", "TimeOptB" = "#2E86C1", "TimeOptBMCMC" = "#B9770E")
for (device in c("png", "tiff")) {
  open_figure("Figure1", device, 10, 5.2)
  layout(matrix(1:2, ncol = 1), heights = c(1, 0.16))
  par(mar = c(3.8, 4.2, 2.6, 1), mgp = c(2.4, 0.7, 0))
  plot(NA, xlim = c(0, 1), ylim = c(48, 125), log = "y", xaxs = "i", xaxt = "n", yaxt = "n", xlab = "Age (Ma)",
       ylab = "k (arcsec/yr, log scale)")
  axis_time(TRUE)
  time_background()
  abline(h = k_now, lty = 4, col = "grey40")
  colours <- method_colour[table1$published_method]
  arrows(tx(table1$age_Ma), table1$published_mean - table1$published_sd, tx(table1$age_Ma),
         table1$published_mean + table1$published_sd, angle = 90, code = 3, length = 0.012, col = colours)
  points(tx(table1$age_Ma), table1$published_mean, pch = method_pch[table1$published_method], col = colours, cex = 1)
  box()
  u <- par("usr")
  gts_bar(10^u[4] * 1.005, 10^u[4] * 1.10)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", ncol = 3, bty = "n", cex = 0.8,
         pch = c(16, 17, 15, NA, NA, NA, NA), lty = c(NA, NA, NA, 1, 2, 4, NA),
         col = c(method_colour[c("TimeOptB", "TimeOptMCMC", "TimeOptBMCMC")], "black", "grey45", "grey40", NA),
         lwd = c(NA, NA, NA, 1.6, 1, 1, NA),
         legend = c("TimeOptB", "TimeOptMCMC", "TimeOptBMCMC", "Farhat et al. (2022a) solution and envelope",
                    "TimeOptB prior and 2 SD band", "Present-day k", ""))
  dev.off()
}

# Figures 2 and 3: one column per TimeOptBMCMC run ####
# x: the prior mean of k of the nine runs (each run fills a column 0.5 SD wide, centred on its prior mean); y: k.
# The colour of a cell is the density of the run's samples relative to the highest density of that run: darkest
# within 2 times of the maximum, then 10, 100, 10^3 and 10^4 times less. Points: posterior mean of each run.
blues <- c("#DEEBF7", "#C6DBEF", "#9ECAE1", "#4292C6", "#08519C")
levels_log10 <- c(-4, -3, -2, -1, log10(0.5))
matrix_panel <- function(i, label, cex_t = 0.5, cex_b = 0.5) {
  q <- table1[i, ]
  m0 <- q$farhat_mean
  s0 <- q$farhat_sd
  samples <- record_samples(q)
  q$probs_low <- 0.002
  q$probs_high <- 0.998
  rng <- panel_range(q, samples$runs)
  plot(NA, xlim = rng, ylim = rng, xlab = "", ylab = "", las = 1, xaxs = "i", yaxs = "i")
  ny <- 300
  yg <- seq(rng[1], rng[2], length.out = ny + 1)
  for (s in names(samples$runs)) {
    k <- samples$runs[[s]]
    if (length(k) < 50) {
      next
    }
    x0 <- m0 + (as.numeric(s) - 0.25) * s0
    x1 <- m0 + (as.numeric(s) + 0.25) * s0
    d <- density(k, n = 1024, from = rng[1], to = rng[2], bw = max(bw.nrd0(k), diff(rng) / 400))
    dy <- approx(d$x, d$y, xout = (yg[-1] + yg[-(ny + 1)]) / 2)$y
    relative <- log10(pmax(dy, 1e-300) / max(dy))
    class <- findInterval(relative, levels_log10)
    for (j in which(class > 0)) {
      rect(x0, yg[j], x1, yg[j + 1], col = blues[class[j]], border = blues[class[j]], lwd = 0.3)
    }
    points((x0 + x1) / 2, mean(k), pch = 16, cex = 0.3, col = "black")
  }
  if (q$floor > rng[1]) {
    rect(rng[1] - 100, rng[1] - 100, rng[2] + 100, q$floor, col = "grey80", border = NA)
  }
  abline(0, 1, lty = 3, col = "grey20")
  abline(h = m0, v = m0, lty = 5, col = "#1E8449")
  points(rng[1], q$published_mean, pch = 15, col = "#C0392B", cex = 0.9, xpd = NA)
  wide <- samples$wide
  if (length(wide) > 50) {
    segments(rng[2], mean(wide) - sd(wide), rng[2], mean(wide) + sd(wide), col = "grey35", lwd = 2.5, xpd = NA,
             lend = 1)
    points(rng[2], mean(wide), pch = 23, bg = "white", col = "grey25", cex = 0.6, xpd = NA)
  }
  box()
  mtext(paste0(label, "  ", panel_name_of(q$record), ", ", round(q$age_Ma), " Ma"), side = 3, line = 0.25, cex = cex_t,
        font = 2)
  draw_badge(q$reading, cex_b)
}
matrix_legend <- function(cex) {
  legend("center", ncol = 3, bty = "n", cex = cex,
         pch = c(22, 22, 22, 22, 22, 16, 15, 23, NA, NA, 22, 22, 22, 22, 22),
         pt.bg = c(rev(blues), NA, NA, "white", NA, NA, "grey80", "#239B56", "#E67E22", "#F5B041", "#C0392B"),
         pt.cex = c(1.6, 1.6, 1.6, 1.6, 1.6, 0.6, 0.9, 0.8, NA, NA, 1.6, 1.6, 1.6, 1.6, 1.6),
         lty = c(NA, NA, NA, NA, NA, NA, NA, NA, 3, 5, NA, NA, NA, NA, NA),
         col = c(rep("grey40", 5), "black", "#C0392B", "grey25", "grey20", "#1E8449", "grey80", "#239B56", "#E67E22",
                 "#F5B041", "#C0392B"),
         legend = c("Density of one run: within 2x of its maximum", "10x less", "100x less", "1,000x less",
                    "10,000x less", "Posterior mean of the run", "Published k (left edge)",
                    paste("Run with a wide prior, mean", pm, "1 SD (right edge)"), "Posterior = prior mean (1:1 line)",
                    "Farhat et al. (2022a) value at that age", "Below the geological floor",
                    "\u2713 D: set by the data (w: weak)", "\u2713 O: a stated observation picks the solution",
                    "~ M: a solution at the model value (an observation is needed)",
                    "\u2717 P: prior in, prior out; X: no plausible answer"))
}
n_records <- nrow(table1)
n_columns <- 7
n_rows <- ceiling(n_records / n_columns)

# Figure 2: four examples ####
examples <- c("Chang", "WalvisRidge", "TorreDeBusi", "Hakatai")
for (device in c("png", "tiff")) {
  open_figure("Figure2", device, 10, 3.9)
  layout(rbind(1:4, rep(5, 4)), heights = c(2.7, 1.2))
  par(mar = c(2.6, 3.2, 1.8, 0.9), mgp = c(1.6, 0.4, 0), tcl = -0.25, cex.axis = 0.75, pty = "s")
  for (j in seq_along(examples)) {
    matrix_panel(which(table1$record == examples[j]), paste0("(", letters[j], ")"), cex_t = 0.75, cex_b = 0.7)
    mtext("prior mean of k", side = 1, line = 1.5, cex = 0.6)
    if (j == 1) {
      mtext("posterior of k (arcsec/yr)", side = 2, line = 2.1, cex = 0.6)
    }
  }
  par(mar = c(0, 0, 0, 0), pty = "m")
  plot.new()
  matrix_legend(0.72)
  dev.off()
}

# Figure 3: all 56 estimates ####
for (device in c("png", "tiff")) {
  open_figure("Figure3", device, 10, 1.5 * n_rows + 1.1)
  layout(rbind(matrix(1:(n_columns * n_rows), nrow = n_rows, byrow = TRUE), rep(n_columns * n_rows + 1, n_columns)),
         heights = c(rep(1.5, n_rows), 1.1))
  par(mar = c(1.3, 1.9, 1.5, 0.6), mgp = c(1.2, 0.3, 0), tcl = -0.2, cex.axis = 0.58, pty = "s")
  for (i in seq_len(n_records)) {
    matrix_panel(i, paste0(LETTERS[(i - 1) %/% n_columns + 1], (i - 1) %% n_columns + 1))
  }
  for (j in seq_len(n_columns * n_rows - n_records)) {
    plot.new()
  }
  par(mar = c(0, 0, 0, 0), pty = "m")
  plot.new()
  matrix_legend(0.6)
  mtext(paste("x: prior mean of k (arcsec/yr), nine TimeOptBMCMC runs at -2 to +2 SD (one column each); ",
              "y: posterior of k (arcsec/yr)"), side = 3, line = -0.9, cex = 0.6)
  dev.off()
}

# Figure S1: the runs of Figure 3 as violins ####
violin_panel <- function(i, label, cex_t = 0.5, cex_b = 0.5) {
  q <- table1[i, ]
  m0 <- q$farhat_mean
  s0 <- q$farhat_sd
  samples <- record_samples(q)
  q$probs_low <- 0.005
  q$probs_high <- 0.995
  rng <- panel_range(q, samples$runs)
  plot(NA, xlim = rng, ylim = rng, xlab = "", ylab = "", las = 1)
  if (q$floor > rng[1]) {
    rect(rng[1] - 100, rng[1] - 100, rng[2] + 100, q$floor, col = "grey85", border = NA)
  }
  wide <- samples$wide
  if (length(wide) > 50) {
    rect(rng[1] - 100, mean(wide) - sd(wide), rng[2] + 100, mean(wide) + sd(wide), col = "#D5D8DC80", border = NA)
    abline(h = mean(wide), col = "grey50", lwd = 0.8)
  }
  abline(0, 1, lty = 3, col = "grey30")
  abline(h = m0, v = m0, lty = 5, col = "#1E8449")
  half_width <- 0.22 * s0
  for (s in names(samples$runs)) {
    k <- samples$runs[[s]]
    if (length(k) < 50) {
      next
    }
    x <- m0 + as.numeric(s) * s0
    d <- density(k, n = 256, from = min(k), to = max(k))
    w <- half_width * d$y / max(d$y)
    polygon(c(x - w, rev(x + w)), c(d$x, rev(d$x)), col = "#5DADE2", border = "#1F618D", lwd = 0.3)
    segments(x, mean(k) - sd(k), x, mean(k) + sd(k), lwd = 0.7)
    points(x, mean(k), pch = 16, cex = 0.35)
  }
  points(rng[1], q$published_mean, pch = 15, col = "#C0392B", cex = 0.9, xpd = NA)
  box()
  mtext(paste0(label, "  ", panel_name_of(q$record), ", ", round(q$age_Ma), " Ma"), side = 3, line = 0.25, cex = cex_t,
        font = 2)
  draw_badge(q$reading, cex_b)
}
for (device in c("png", "tiff")) {
  open_figure("Figure_S1", device, 10, 1.5 * n_rows + 0.9)
  layout(rbind(matrix(1:(n_columns * n_rows), nrow = n_rows, byrow = TRUE), rep(n_columns * n_rows + 1, n_columns)),
         heights = c(rep(1.5, n_rows), 0.9))
  par(mar = c(1.3, 1.9, 1.5, 0.4), mgp = c(1.2, 0.3, 0), tcl = -0.2, cex.axis = 0.58, pty = "s")
  for (i in seq_len(n_records)) {
    violin_panel(i, paste0(LETTERS[(i - 1) %/% n_columns + 1], (i - 1) %% n_columns + 1))
  }
  for (j in seq_len(n_columns * n_rows - n_records)) {
    plot.new()
  }
  par(mar = c(0, 0, 0, 0), pty = "m")
  plot.new()
  legend("center", ncol = 3, bty = "n", cex = 0.62,
         pch = c(22, 16, 22, 15, NA, NA, 22),
         pt.bg = c("#5DADE2", NA, "#D5D8DC", NA, NA, NA, "grey85"), pt.cex = c(1.6, 0.6, 1.6, 0.9, NA, NA, 1.6),
         lty = c(NA, 1, NA, NA, 3, 5, NA),
         col = c("#1F618D", "black", "grey50", "#C0392B", "grey30", "#1E8449", "grey85"),
         legend = c("Posterior of one TimeOptBMCMC run, at its prior mean (width: density)",
                    paste("Posterior mean", pm, "1 SD of that run"),
                    paste("Run with a wide prior (5 SD): mean", pm, "1 SD"),
                    "Published k (left edge)", "Posterior = prior mean (1:1 line)",
                    "Farhat et al. (2022a) value at that age", "Below the geological floor (samples dropped)"))
  mtext(paste("x: prior mean of k (arcsec/yr), nine TimeOptBMCMC runs at -2 to +2 SD;  y: posterior of k (arcsec/yr);",
              "same scale on both axes"), side = 3, line = -0.9, cex = 0.6)
  dev.off()
}

# Figure S3: the estimates set by the data or by an observation, and the curve of k against age ####
used <- table1[table1$reading %in% c("D strong", "D weak", "O"), ]
weak <- used$reading != "D strong"
used_pch <- method_pch[used$published_method]
for (i in which(weak)) {
  used_pch[i] <- c("16" = 1, "17" = 2, "15" = 0)[[as.character(used_pch[i])]]
}
colour_data <- "#239B56"
used_colour <- rep(colour_data, nrow(used))
used_colour[used$reading == "O"] <- "#E67E22"
oldest <- max(used$age_Ma)
for (device in c("png", "tiff")) {
  open_figure("Figure_S3", device, 10, 6.4)
  layout(matrix(c(1, 2), ncol = 1), heights = c(1, 0.26))
  par(mar = c(3.8, 4.2, 3.2, 1), mgp = c(2.3, 0.7, 0))
  plot(NA, xlim = c(0, 1), ylim = c(48, 125), log = "y", xaxs = "i", xaxt = "n", yaxt = "n", xlab = "Age (Ma)",
       ylab = "k (arcsec/yr, log scale)")
  axis_time(FALSE)
  time_background()
  c4 <- curve_all[curve_all$age_Ma <= 2600, ]
  polygon(tx(c(c4$age_Ma, rev(c4$age_Ma))), c(c4$q025, rev(c4$q975)), col = "#D3540030", border = NA)
  lines(tx(c4$age_Ma[c4$age_Ma <= oldest]), c4$median[c4$age_Ma <= oldest], lwd = 1.6, col = "#D35400")
  lines(tx(c4$age_Ma[c4$age_Ma >= oldest]), c4$median[c4$age_Ma >= oldest], lwd = 1.4, lty = 2, col = "#D35400")
  cs <- curve_strong[curve_strong$age_Ma <= max(used$age_Ma[!weak]), ]
  lines(tx(cs$age_Ma), cs$q025, lty = 3, lwd = 0.8, col = "#6C3483")
  lines(tx(cs$age_Ma), cs$q975, lty = 3, lwd = 0.8, col = "#6C3483")
  lines(tx(cs$age_Ma), cs$median, lwd = 1.3, lty = 5, col = "#6C3483")
  abline(h = k_now, lty = 4, col = "grey40")
  segments(tx(used$age_Ma), used$published_mean, tx(used$age_Ma), used$main_k, col = "grey70", lwd = 0.8)
  points(tx(used$age_Ma), used$published_mean, pch = 4, col = "grey55", cex = 0.6)
  arrows(tx(used$age_Ma), used$main_k - used$main_sd, tx(used$age_Ma), used$main_k + used$main_sd, angle = 90,
         code = 3, length = 0.012, col = used_colour)
  points(tx(used$age_Ma), used$main_k, pch = used_pch, col = used_colour, bg = "white", cex = 1)
  # solutions at the model value (M): shown, but not used for the curve
  model <- table1[table1$reading == "M", ]
  arrows(tx(model$age_Ma), model$main_k - model$main_sd, tx(model$age_Ma), model$main_k + model$main_sd, angle = 90,
         code = 3, length = 0.012, col = "#F5B041")
  points(tx(model$age_Ma), model$main_k, pch = 5, col = "#D68910", bg = "white", cex = 0.9)
  box()
  u <- par("usr")
  gts_bar(10^u[4] * 1.005, 10^u[4] * 1.10)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", ncol = 2, bty = "n", cex = 0.72,
         pch = c(16, 17, 0, 5, 1, 4, NA, NA, NA, NA),
         col = c(colour_data, colour_data, "#E67E22", "#D68910", colour_data, "grey55", "black", "#D35400", "#6C3483",
                 "grey45"),
         lty = c(NA, NA, NA, NA, NA, NA, 1, 1, 5, 2), lwd = c(NA, NA, NA, NA, NA, NA, 1.4, 1.6, 1.3, 1),
         legend = c("Set by the data: published with TimeOptB", "Set by the data: published with TimeOptMCMC",
                    "Orange: set by an observation (Hakatai, Dales Gorge; open: weak)",
                    "Solution at the model value, not picked by the data (M; not in the curve)", "Open: weak reading",
                    "Published value (joined to the value of the runs)", "Farhat et al. (2022a), with envelope",
                    "This study: k(t) from these estimates (median, 95%; weak ones penalised)",
                    "Same, strong readings only (dotted: 95%)", "TimeOptB prior, 2 SD"))
  dev.off()
}

# Figure S4: the curve from today to the formation of the Earth ####
for (device in c("png", "tiff")) {
  open_figure("Figure_S4", device, 10, 5.2)
  par(mar = c(4, 4.5, 1, 1), mgp = c(2.5, 0.7, 0))
  plot(NA, xlim = c(0, 4540), ylim = c(48, 3000), log = "y", xaxs = "i", xlab = "Age (Ma)",
       ylab = "k (arcsec/yr, log scale)", las = 1, yaxt = "n")
  axis(2, at = c(50, 70, 100, 200, 500, 1000, 2000), las = 1)
  polygon(c(astrogeo$age, rev(astrogeo$age)), c(astrogeo$lo, rev(astrogeo$hi)), col = "#2C3E5040", border = NA)
  lines(astrogeo$age, astrogeo$pre, lwd = 1.4)
  polygon(c(curve_all$age_Ma, rev(curve_all$age_Ma)), c(curve_all$q025, rev(curve_all$q975)), col = "#D3540030",
          border = NA)
  lines(curve_all$age_Ma[curve_all$age_Ma <= oldest], curve_all$median[curve_all$age_Ma <= oldest], lwd = 1.6,
        col = "#D35400")
  lines(curve_all$age_Ma[curve_all$age_Ma >= oldest], curve_all$median[curve_all$age_Ma >= oldest], lwd = 1.4, lty = 2,
        col = "#D35400")
  arrows(used$age_Ma, used$main_k - used$main_sd, used$age_Ma, used$main_k + used$main_sd, angle = 90, code = 3,
         length = 0.02, col = used_colour)
  points(used$age_Ma, used$main_k, pch = 16, cex = 0.7, col = used_colour)
  abline(h = k_now, lty = 4, col = "grey40")
  abline(v = 4425, lty = 3, col = "grey40")
  text(4400, 1500, "Moon formation (4425 Ma), k of 2000 or more:\nused only as the upper bound of the curve", adj = 1,
       cex = 0.7)
  legend("topleft", bty = "n", cex = 0.75, lty = c(1, 1, 2, NA), pch = c(NA, NA, NA, 16),
         lwd = c(1.4, 1.6, 1.4, NA), col = c("black", "#D35400", "#D35400", colour_data),
         legend = c("Farhat et al. (2022a, 2022b), AstroGeo22, with envelope", "This study: monotone k(t), median and 95%",
                    paste0("Beyond the oldest estimate (", oldest, " Ma): extrapolation"),
                    paste("Estimates set by the data (", pm, " 1 SD)", sep = "")))
  dev.off()
}
