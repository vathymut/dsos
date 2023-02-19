#' @noRd
#' @keywords  internal
rudirichlet <- function(n) {
  # Function adapted (inspired) from bayesboot package
  weights <- stats::rexp(n, rate = 1)
  weights <- weights / sum(weights)
  return(weights * n)
}

#' @noRd
#' @keywords  internal
draw_bb <- function(os_train, os_test) {
  # _bb suffix stands for bayesian bootstrap
  n_train <- length(os_train)
  n_test <- length(os_test)
  w_train <- rudirichlet(n_train)
  w_test <- rudirichlet(n_test)
  weight <- c(w_train, w_test)
  stat_bb <- wauc_from_os(os_train, os_test, weight = weight)
  return(stat_bb)
}

#' @noRd
#' @keywords  internal
repeat_fn <- function(os_train, os_test, fn, n_pt = 4e3) {
  wauc_dist <- future.apply::future_replicate(
    n_pt,
    fn(os_train, os_test)
  )
  return(wauc_dist)
}

#' @noRd
#' @keywords  internal
wauc_bb <- function(os_train, os_test, n_pt = 4e3) {
  wauc_dist <- repeat_fn(
    os_train,
    os_test,
    fn = draw_bb,
    n_pt = n_pt
  )
  return(wauc_dist)
}

#' @noRd
#' @keywords  internal
draw_bb_and_perm <- function(os_train, os_test) {
  n_train <- length(os_train)
  n_test <- length(os_test)
  w_train <- rudirichlet(n_train)
  w_test <- rudirichlet(n_test)
  weight <- c(w_train, w_test)
  wauc_bb <- wauc_from_os(os_train, os_test, weight = weight)
  shuffled <- shuffle_os(c(os_train, os_test), n_test)
  wauc_perm <- wauc_from_os(shuffled$train, shuffled$test, weight = NULL)
  return(c(permuted = wauc_perm, posterior = wauc_bb))
}

#' @noRd
#' @keywords  internal
wauc_samples <- function(os_train, os_test, n_pt = 4e3) {
  wauc_dist <- repeat_fn(
    os_train,
    os_test,
    fn = draw_bb_and_perm,
    n_pt = n_pt
  )
  return(list(permuted = wauc_dist[1, ], posterior = wauc_dist[2, ]))
}

#' @title
#' Bayesian Test from Outlier Scores
#'
#' @param os_train Outlier scores in training (reference) set.
#' @param os_test Outlier scores in test set.
#' @param n_pt The number of permutations.
#' @param threshold Threshold for adverse shift. Defaults to 1 / 12,
#' the asymptotic value of the test statistic when the two samples are drawn
#' from the same distribution.
#'
#' @inherit pt_from_os description
#' @inheritSection at_from_os Notes
#'
#' @return
#' A named list of class \code{outlier.bayes} containing:
#' \itemize{
#'    \item \code{posterior}: Posterior distribution of WAUC test statistic
#'    \item \code{threshold}: WAUC threshold for adverse shift
#'    \item \code{adverse_probability}: probability of adverse shift
#'    \item \code{bayes_factor}: Bayes factor
#'    \item \code{outlier_scores}: outlier scores from training and test set
#' }
#'
#' @references Kamulete, V. M. (2023).
#' \emph{Are you OK? A Bayesian test for adverse shift}.
#' Manuscript in preparation.
#'
#' @references Johnson, V. E. (2005).
#' \emph{Bayes factors based on test statistics}.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 67(5), 689-701.
#'
#' @references Gu, J., Ghosal, S., & Roy, A. (2008).
#' \emph{Bayesian bootstrap estimation of ROC curve}.
#' Statistics in medicine, 27(26), 5407-5420.
#'
#' @details
#' The posterior distribution of the test statistic is based on \code{n_pt}
#' (boostrap) permutations. The method uses the Bayesian bootstrap as a
#' resampling procedure as in Gu et al (2008). Johnson (2005) shows to
#' leverage (turn) a test statistic into a Bayes factor. The test statistic
#' is the weighted AUC (WAUC).
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' os_train <- rnorm(n = 100)
#' os_test <- rnorm(n = 100)
#' bayes_test <- bf_from_os(os_train, os_test)
#' bayes_test
#' # To run in parallel on local cluster, uncomment the next two lines.
#' # library(future)
#' # future::plan(future::multisession)
#' parallel_test <- bf_from_os(os_train, os_test)
#' parallel_test
#' }
#'
#' @family bayesian-test
#'
#' @export
bf_from_os <- function(os_train,
                       os_test,
                       n_pt = 4e3,
                       threshold = 1 / 12) {
  posterior <- wauc_bb(os_train, os_test, n_pt = n_pt)
  adverse_prob <- 1 - stats::ecdf(posterior)(threshold)
  bayes_factor <- adverse_prob / (1 - adverse_prob)
  result <- list(
    posterior = posterior,
    threshold = threshold,
    adverse_probability = adverse_prob,
    bayes_factor = bayes_factor,
    outlier_scores = list(train = os_train, test = os_test)
  )
  class(result) <- "outlier.bayes"
  return(result)
}

#' @title
#' Bayesian and Frequentist Test from Outlier Scores
#'
#' @inherit pt_from_os description
#' @inheritSection at_from_os Notes
#' @inheritParams bf_from_os
#'
#' @return
#' A list of factors (BF) for 3 different test specifications:
#' \itemize{
#'    \item \code{frequentist}: Frequentist BF.
#'    \item \code{bayes_noperm}: Bayestion BF test with asymptotic threshold.
#'    \item \code{bayes_perm}: Bayestion BF with exchangeable threshold.
#' }
#'
#' @details
#' This compares the Bayesian to the frequentist approach for convenience.
#' The Bayesian test mimics `bf_from_os()` and the frequentist one,
#' `pt_from_os()`. The Bayesian test computes Bayes factors based on the
#' asymptotic (defaults to 1/12) and the exchangeable threshold. The latter
#' calculates the threshold as the median weighted AUC (WAUC) after \code{n_pt}
#' permutations assuming outlier scores are exchangeable. This is recommended
#' for small samples. The frequentist test converts the one-sided (one-tailed)
#' p-value to the Bayes factor - see \code{as_bf} function.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' set.seed(12345)
#' os_train <- rnorm(n = 100)
#' os_test <- rnorm(n = 100)
#' bayes_test <- bf_compare(os_train, os_test)
#' bayes_test
#' # To run in parallel on local cluster, uncomment the next two lines.
#' # library(future)
#' # future::plan(future::multisession)
#' parallel_test <- bf_compare(os_train, os_test)
#' parallel_test
#' }
#'
#' @family bayesian-test
#'
#' @seealso
#' [bf_from_os()] for bayes factor, the Bayesian test.
#' [pt_from_os()] for p-value, the frequentist test.
#'
#' @export
bf_compare <- function(os_train,
                       os_test,
                       threshold = 1 / 12,
                       n_pt = 4e3) {
  draws <- wauc_samples(os_train, os_test, n_pt = n_pt)
  # Get p-value
  permuted <- draws$permuted
  test_stat <- wauc_from_os(os_train, os_test)
  p_value <- 1 - stats::ecdf(permuted)(test_stat)
  freq_bf <- as_bf(p_value)
  # Get bayes factor from asymptotic threshold
  posterior <- draws$posterior
  cdf_fn <- stats::ecdf(posterior)
  asym_prob <- 1 - cdf_fn(threshold)
  asym_bf <- asym_prob / (1 - asym_prob)
  # Get bayes factor from exchangeable threshold
  pt_threshold <- stats::quantile(permuted, probs = 0.5)
  pt_prob <- 1 - cdf_fn(pt_threshold)
  pt_bf <- pt_prob / (1 - pt_prob)
  bf <- list(bayes_perm = pt_bf, bayes_noperm = asym_bf)
  bf[["frequentist"]] <- freq_bf
  return(bf)
}

#' @title
#' Convert P-value to Bayes Factor
#'
#' @param pvalue P-value.
#'
#' @return Bayes Factor (scalar value).
#'
#' @references Marsman, M., & Wagenmakers, E. J. (2017).
#' \emph{Three insights from a Bayesian interpretation of the one-sided P value}.
#' Educational and Psychological Measurement, 77(3), 529-539.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' bf_from_pvalue <- as_bf(pvalue = 0.5)
#' bf_from_pvalue
#' }
#'
#' @family bayesian-test
#'
#' @seealso
#' [as_pvalue()] to convert Bayes factor to p-value.
#'
#' @export
as_bf <- function(pvalue) {
  bf <- exp(stats::qlogis(pvalue))
  inv_bf <- 1. / bf
  return(inv_bf)
}

#' @title
#' Convert Bayes Factor to P-value
#'
#' @param bf Bayes factor.
#'
#' @return p-value (scalar value).
#'
#' @references Marsman, M., & Wagenmakers, E. J. (2017).
#' \emph{Three insights from a Bayesian interpretation of the one-sided P value}.
#' Educational and Psychological Measurement, 77(3), 529-539.
#'
#' @examples
#' \donttest{
#' library(dsos)
#' pvalue_from_bf <- as_pvalue(bf = 1)
#' pvalue_from_bf
#' }
#'
#' @family bayesian-test
#'
#' @seealso
#' [as_bf()] to convert p-value to Bayes factor.
#'
#' @export
as_pvalue <- function(bf) {
  inv_bf <- 1. / bf
  pvalue <- stats::plogis(log(inv_bf))
  return(pvalue)
}
