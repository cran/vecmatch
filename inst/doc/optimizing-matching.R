## ----include = FALSE----------------------------------------------------------
library(vecmatch)
data(cancer) # or data("cancer", package = "vecmatch")
formula_cancer <- formula(status ~ age * sex)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.asp = 0.8,
  echo = TRUE,
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  cache = TRUE # so that estimate_gps doesn’t re-run on every knit
)

## ----include=FALSE------------------------------------------------------------
library(vecmatch)

formula_cancer <- formula(status ~ age * sex)

## -----------------------------------------------------------------------------
opt_args <- make_opt_args(
  data            = cancer,
  formula         = formula_cancer,
  reference       = c("control", "adenoma", "crc_benign", "crc_malignant"),
  gps_method      = c("m1", "m7", "m8"),
  matching_method = c("fullopt", "nnm"),
  caliper         = seq(0.01, 5, 0.01),
  cluster         = 1:3,
  ratio           = 1:3,
  min_controls    = 1:3,
  max_controls    = 1:3
)

opt_args

## ----warning=FALSE, message=FALSE---------------------------------------------
library(future)
library(doFuture)

## 1. Register future as the foreach backend
doFuture::registerDoFuture()

## 2. Choose parallel strategy
old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
future::plan(
  future::multisession,
  workers = 4
)

## 3. Seeding before calling optimize_gps()
set.seed(167894)
seed_before <- .Random.seed

opt_results <- optimize_gps(
  data     = cancer,
  formula  = formula_cancer,
  opt_args = opt_args,
  n_iter   = 6000
)

summary(opt_results)

## -----------------------------------------------------------------------------
select_results <- select_opt(opt_results,
  perc_matched = c("adenoma", "crc_malignant"),
  smd_variables = "age",
  smd_type = "max"
)

summary(select_results)

## ----estimate_gps-------------------------------------------------------------
# estimating gps
final_match <- run_selected_matching(select_results,
  cancer,
  formula_cancer,
  smd_group = "0.05-0.10"
)

summary(final_match)

## -----------------------------------------------------------------------------
balqual(final_match,
  formula = formula_cancer,
  type = "smd",
  statistic = "max",
  round = 3,
  cutoffs = 0.2
)

all.equal(seed_before, .Random.seed)

