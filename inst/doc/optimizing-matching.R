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
  reference       = c("control", "adenoma", "crc_beningn", "crc_malignant"),
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
set.seed(167894)
seed_before <- .Random.seed

opt_results <- optimize_gps(
  data     = cancer,
  formula  = formula_cancer,
  opt_args = opt_args,
  n_iter   = 1500
)

## -----------------------------------------------------------------------------
select_results <- select_opt(opt_results,
  perc_matched = c("adenoma", "crc_malignant"),
  smd_variables = "age",
  smd_type = "max"
)

## -----------------------------------------------------------------------------
param_df <- attr(select_results, "param_df")

subset(param_df, smd_group == "0.10-0.15")

## ----estimate_gps-------------------------------------------------------------
# estimating gps
gps_mat <- estimate_gps(formula_cancer,
  cancer,
  method = "multinom",
  link = "generalized_logit",
  reference = "adenoma"
)

## -----------------------------------------------------------------------------
# csr with refitting
csr_mat <- csregion(gps_mat)

# matching
matched_df <- match_gps(csr_mat,
  method = "nnm",
  caliper = 4.29,
  reference = "adenoma",
  ratio = 2,
  order = "asc",
  replace = TRUE,
  ties = TRUE,
  kmeans_cluster = 1
)

## -----------------------------------------------------------------------------
balqual(matched_df,
  formula = formula_cancer,
  type = "smd",
  statistic = "max",
  round = 3,
  cutoffs = 0.2
)

all.equal(seed_before, .Random.seed)

