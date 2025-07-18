#' @title Filter the data based on common support region
#'
#' @description The `csregion()` function estimates the boundaries of the
#'   rectangular common support region, as defined by Lopez and Gutman (2017),
#'   and filters the matrix of generalized propensity scores based on these
#'   boundaries. The function returns a matrix of observations whose generalized
#'   propensity scores lie within the treatment group-specific boundaries.
#'
#' @param gps_matrix An object of classes `gps` and `data.frame` (e.g., created
#'   by the `estimate_gps()` function). The first column corresponds to the
#'   treatment or grouping variable, while the other columns represent the
#'   treatment assignment probabilities calculated separately for each
#'   hypotetical treatment group. The number of columns should therefore be
#'   equal to the number of unique levels of the treatment variable plus one
#'   (for the treatment variable itself). The number of rows should correspond
#'   to the number of subjects for which generalized propensity scores were
#'   estimated.
#'
#' @param borders A character string specifying how to handle observations at
#'   the edges of the Common Support Region (CSR). Acceptable values are
#'   `"include"` and `"exclude"`. If `"include"` is selected (default),
#'   observations with Generalized Propensity Scores (GPS) exactly equal to the
#'   CSR boundaries are retained for further analysis. This corresponds to a
#'   non-strict inequality: \code{lower_bound <= GPS <= upper_bound}. If
#'   `"exclude"` is selected, observations lying exactly on the CSR boundaries
#'   are removed. This corresponds to a strict inequality: \code{lower_bound <
#'   GPS < upper_bound}. Using `"exclude"` will typically result in a slightly
#'   smaller matched sample size compared to `"include"`, but may be preferred
#'   for more conservative matching.
#'
#' @param refit Logical. If `TRUE` (default), the model used to estimate the GPS
#'   is refitted after excluding samples outside the common support region,
#'   using the same formula and method as in the original `estimate_gps()` call.
#'   If `FALSE`, the model is not refitted, but still only samples within the
#'   CSR are retained. Refitting is recommended, as suggested by Lopez and
#'   Gutman (2017).
#'
#' @return A numeric matrix similar to the one returned by `estimate_gps()`,
#'   but with the number of rows reduced to exclude those observations that do
#'   not fit within the common support region (CSR) boundaries. The returned
#'   object also possesses additional attributes that summarize the calculation
#'   process of the CSR boundaries:
#'  * `filter_matrix` - A logical matrix with the same dimensions as the
#'  gps-part of `gps_matrix`, indicating which treatment assignment
#'  probabilities fall within the CSR boundaries,
#'  * `filter_vector` - A vector indicating whether each observation was kept
#'  (`TRUE`) or removed (`FALSE`), essentially a row-wise
#'  sum of `filter_matrix`,
#'  * `csr_summary` - A summary of the CSR calculation process, including
#'  details of the boundaries and the number of observations filtered.
#'  * `csr_data` - The original dataset used for the estimation of generalized
#'  propensity scores (`original_data` attribute of the `gps` object) filtered
#'  by the `filter_vector`
#'
#' @examples
#' # We could estimate simples generalized propensity scores for the `iris`
#' # dataset
#' gps <- estimate_gps(Species ~ Sepal.Length, data = iris)
#'
#' # And then define the common support region boundaries using `csregion()`
#' gps_csr <- csregion(gps)
#'
#' # The additional information of the CSR-calculation process are
#' # accessible through the attributes described in the `*Value*` section
#' attr(gps_csr, "filter_matrix")
#' attr(gps_csr, "csr_summary")
#' attr(gps_csr, "csr_data")
#'
#' @export
csregion <- function(gps_matrix,
                     borders = "include",
                     refit = TRUE) {
  csr_data <- attr(gps_matrix, "original_data")

  .chk_cond(
    "gps" %nin% class(gps_matrix),
    "The `gps_matrix` argument must be of class `gps`."
  )

  # check borders arg
  chk::chk_character(borders)
  chk::chk_length(borders, length = 1)
  .chk_cond(
    borders %nin% c("include", "exclude"),
    'The `borders` argument can only take one of the
            following values: "include", "exclude".'
  )

  # check refit arg
  chk::chk_flag(refit)

  ## Calculating the csr_low
  csr_low <- apply(
    stats::aggregate(. ~ treatment,
      data = gps_matrix,
      FUN = function(x) min(x)
    )[, -1],
    2,
    max
  )

  ## Calculating the csr_high
  csr_high <- apply(
    stats::aggregate(. ~ treatment,
      data = gps_matrix,
      FUN = function(x) max(x)
    )[, -1],
    2,
    min
  )

  ## filter out the unvalid observations
  filter_matrix <- mapply(function(df_col, vec_low, vec_high) {
    switch(borders,
      include = vec_low <= df_col & df_col <= vec_high,
      exclude = vec_low < df_col & df_col < vec_high
    )
  }, gps_matrix[, 2:ncol(gps_matrix)], csr_low, csr_high)

  ## summarizing the logical matrix to a subset vector
  filter_vector <- apply(filter_matrix, 1, all)

  ## defining the number of negatives
  n_negative <- sum(!filter_vector)
  n_negative_matrix <- colSums(!filter_matrix)

  ## refitting the gps_matrix
  if (refit) {
    # 1.) Filter out observations from original data
    # 2.) Change function call
    # 3.) Evaluate function call
    # 4.) Overwrite gps_matrix

    # Limiting original data only to the csr
    csr_filtered <- csr_data[filter_vector, ]

    # Saving and changing function call for estimate_gps
    estimate_call <- attr(gps_matrix, "function_call")
    estimate_call$data <- csr_filtered

    # Evaluating the changed call
    gps_matrix <- tryCatch(
      eval(estimate_call),
      error = function(e) {
        warning("Refitting caused an error and will be ignored.
                Setting refit = FALSE.")
        refit <<- FALSE # Assign to parent env (if refit is defined there)
        csr_filtered # Return filtered data
      }
    )
  } else {
    ## subsetting the gps_matrix
    gps_matrix <- subset(gps_matrix, filter_vector)
  }



  ## drop unused levels of treatment variable
  gps_matrix[, "treatment"] <- droplevels(gps_matrix[, "treatment"])

  ## detect low number of observations in groups and print a warning
  if (any(table(gps_matrix[, "treatment"]) < 20)) {
    chk::wrn(strwrap("Some groups have fewer than 20 observations, which may
    impact the performance of the matching process. Consider using
    `replace = TRUE`in `match_gps()` to address this.",
      prefix = " ", initial = ""
    ))
  }

  ## assembling the print results list
  csr_summary <- data.frame(
    treatment = colnames(gps_matrix)[2:ncol(gps_matrix)],
    csr_low = csr_low,
    csr_high = csr_high,
    n_negative_matrix = n_negative_matrix
  )

  res <- list(
    data_csr = csr_summary,
    excluded = n_negative
  )

  ## Setting a new class for the results
  csres <- structure(res, class = "csres")

  ## print custom output
  show_csres(csres)

  ## adding attributes to csr output
  attr(gps_matrix, "filter_matrix") <- filter_matrix
  attr(gps_matrix, "filter_vector") <- filter_vector
  attr(gps_matrix, "csr_summary") <- csr_summary
  attr(gps_matrix, "csr_data") <- csr_data[filter_vector, ]

  # Assign attributes and class
  class(gps_matrix) <- c("data.frame", "gps", "csr")

  # return the gps_matrix
  return(invisible(gps_matrix))
}

# Function to show the contents of an object of class 'csres'
show_csres <- function(object) {
  # Function to print a data frame in a table-like format
  print_table <- function(df, colnames_df) {
    # Print column headers
    cat(paste(sprintf("%-15s", colnames_df), collapse = " | "), "\n")
    cat(paste(rep("-", 17 * ncol(df)), collapse = ""), "\n")

    # Print each row of the data frame
    apply(df, 1, function(row) {
      cat(paste(sprintf("%-15s", row), collapse = " | "), "\n")
    })
  }

  # Print header and content
  cat("\n")
  cat("Rectangular CSR Borders Evaluation", "\n")
  cat("==================================\n\n")

  print_table(object$data_csr,
    colnames = c(
      "Treatment", "Lower CSR limit",
      "Upper CSR limit", "Number excluded"
    )
  )

  cat("\n")
  cat("===================================================\n")
  cat("The total number of excluded observations is:\t", object$excluded, "\n")
  cat(strwrap("Note: You can view the summary of the
              CSR calculation using the `attr()` function.",
    prefix = " ", initial = ""
  ))
}
