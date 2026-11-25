# Prior in, prior out

[![DOI](https://zenodo.org/badge/1470449640.svg)](https://doi.org/10.5281/zenodo.24508760)

R code and inputs of the re-analysis of 56 published cyclostratigraphic estimates of Earth's precession constant *k*:

> Arts, M. (submitted). Prior in, prior out: a re-analysis of 56 cyclostratigraphic estimates of Earth's precession
> constant. *Paleoceanography and Paleoclimatology* (Research Letter).

## How to run

1. Install R with the packages **astrochron** (1.5 or later, for `timeOptBMCMC()`) and **WaverideR** (0.6.0 or later,
   for `modBchron()`).
2. Open `prior_test.Rproj` in RStudio. This makes this folder the working folder.
3. Open one of the two scripts and source it, or run it part by part:
   - `prior_test_from_saved_runs.R` gives every table and figure from the saved runs. Only the curve of *k* against
     age (Bchron, about 7 minutes per curve) is run again. Put the saved runs from the Zenodo record of the runs in
     `output/saved_results/` first (see below).
   - `prior_test_complete_analysis.R` makes all 721 runs (560 of the prior test, 65 with the authors' methods and 96 of
     the synthetic test) and every number, table and figure. The runs are made one after another on one core; the 560
     runs of the prior test take about four weeks. Every run is saved as soon as it is finished, and a complete saved
     run is not run again, so the script can be stopped and started again. An option to make the runs on several cores
     is included but switched off.

The seven parts of the complete script:

| Part | What it does |
|---|---|
| 1 | The prior test: 560 TimeOptBMCMC runs, ten for each of the 56 estimates |
| 2 | The readings (set by the data, prior in prior out, no plausible answer, solution at the model value, set by an observation) |
| 3 | Table 1 and Table S1 |
| 4 | The curve of *k* against age (Bchron, `modBchron()` of WaverideR) |
| 5 | The authors' own methods with their own priors: 65 runs (Dales Gorge at −2 SD cannot start), and Table S2 |
| 6 | The synthetic test: 24 synthetic records, 96 runs, and Figure S2 |
| 7 | Figures 1 to 3, S1, S3 and S4 |

## Folders

| Folder | Contents |
|---|---|
| `data/` | Inputs: `prior_test_runs.csv` (settings of the 560 runs), `authors_runs.csv` (settings of the 66 runs with the authors' methods, of which 65 can start), the 56 prepared series (`<record>.txt`: depth in m, value) and `AstroGeo22.txt` (Farhat et al., 2022b). The complete script writes the synthetic records to `data/synthetic/`. |
| `output/figures/` | All figures, PNG and TIFF (300 dpi) |
| `output/tables/` | All tables as printed (CSV) |
| `output/saved_results/` | All results: the saved runs (`runs_prior_test/`, `runs_authors/`, `runs_synthetic/`), the readings, the curves of *k* against age and the Bchron models (RDS and CSV) |

The saved runs (about 300 MB) are not in this repository; they are archived in a separate Zenodo record
(https://doi.org/10.5281/zenodo.24509195). Unzip that record into `output/saved_results/`.

## Data

The series are those of the cited studies, prepared with the data window, interpolation and sedimentation-rate range
of each original study; please cite the original studies when using them. The Požár-CS series is that of Da Silva et
al. (2016).

## Licence

GPL-3.0 (see `LICENSE`).
