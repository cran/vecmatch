## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.asp = 0.8,
  echo = TRUE
)

## ----message=FALSE------------------------------------------------------------
library(vecmatch)
library(ggplot2)

raincloud(cancer,
          age,
          status,
          significance = "t_test",
          sig_label_color = TRUE,
          sig_label_size = 3,
          limits = c(10, 120)) +
    scale_y_continuous(breaks = seq(10, 100, 10))

## ----message=FALSE------------------------------------------------------------
mosaic(cancer,
       status,
       sex,
       group_counts = TRUE,
       significance = TRUE)

## -----------------------------------------------------------------------------
formula_cancer <- status ~ age * sex
gps_matrix <- estimate_gps(formula_cancer,
                           cancer,
                           method = 'multinom',
                           reference = 'control')
head(gps_matrix, 7)

## -----------------------------------------------------------------------------
csr_matrix <- csregion(gps_matrix)

## -----------------------------------------------------------------------------
dim(gps_matrix)
dim(csr_matrix)

## -----------------------------------------------------------------------------
set.seed(164373)
matched_cancer <- match_gps(csr_matrix,
                            caliper = 0.21,
                            kmeans_cluster = 2,
                            reference = 'control',
                            replace = TRUE,
                            method = 'fullopt',
                            order = 'desc')

## -----------------------------------------------------------------------------
balqual(matched_cancer,
        formula_cancer,
        type = 'smd',
        statistic = 'max',
        round = 4)

## -----------------------------------------------------------------------------
  matched_cancer$dataset <- 'matched'
  unmatched_cancer <- cancer
  unmatched_cancer$dataset <- 'unmatched'
  data_full <- rbind(matched_cancer, unmatched_cancer)

## ----message=FALSE, fig.asp=1.5-----------------------------------------------
raincloud(data_full,
          age,
          status,
          dataset,
          significance = "t_test",
          sig_label_color = TRUE,
          sig_label_size = 3,
          limits = c(10, 120)) +
    scale_y_continuous(breaks = seq(10, 100, 10))

## ----message=FALSE, fig.asp=1.5-----------------------------------------------
mosaic(data_full,
       status,
       sex,
       dataset,
       group_counts = TRUE,
       significance = TRUE)

