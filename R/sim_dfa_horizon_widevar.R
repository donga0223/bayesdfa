#' Simulate from a DFA
#'
#' @param num_trends The number of trends.
#' @param num_years The number of years.
#' @param num_ts The number of timeseries.
#' @param num_horizons The number of horizons.
#' @param phi AR(1) coefficients.
#' @param loadings_matrix A loadings matrix. The number of rows should match the
#'   number of timeseries and the number of columns should match the number of
#'   trends. Note that this loadings matrix will be internally manipulated by
#'   setting some elements to 0 and constraining some elements to 1 so that the
#'   model can be fitted. See [fit_dfa()]. See the outfit element `Z` in
#'   the returned list is to see the manipulated loadings matrix. If not
#'   specified, a random matrix `~ N(0, 1)` is used.
#' @param sigma A vector of standard deviations on the observation error. Should
#'   be of the same length as the number of trends. If not specified, random
#'   numbers are used `rlnorm(1, meanlog = log(0.2), 0.1)`.
#' @param sigmaH A covariance H by H matrix of horizon. Diagonal part increasing by h and correlated by H.
#' @param rho Assume sigmaH has AR(1) correlation matrix.
#' @param varIndx Indices of unique observation variances. Defaults to `c(1, 1,
#'   1, 1)`. Unique observation error variances would be specified as `c(1, 2, 3,
#'   4)` in the case of 4 time series.
#' @param extreme_value Value added to the random walk in the extreme time step.
#'   Defaults to not included.
#' @param n.outlier
#' @param extreme_loc Location of single extreme event in the process. The same
#'   for all processes, and defaults to `round(n_t/2)` where `n_t` is the time
#'   series length
#' @param nu_fixed Nu is the degrees of freedom parameter for the
#'   t-distribution, defaults to 100, which is effectively normal.
#' @param user_supplied_deviations An optional matrix of deviations for the trend
#'   random walks. Columns are for trends and rows are for each time step.
#' @export
#' @return A list with the following elements: `y_sim` is the simulated data,
#'   pred is the true underlying data without observation error added, `x` is
#'   the underlying trends, `Z` is the manipulated loadings matrix that is fed
#'   to the model.
#' @importFrom stats rlnorm rnorm rt
#' @examples
#' x <- sim_dfa(num_trends = 2)
#' names(x)
#' matplot(t(x$y_sim), type = "l")
#' matplot(t(x$x), type = "l")
#'
#' set.seed(42)
#' x <- sim_dfa(extreme_value = -4, extreme_loc = 10)
#' matplot(t(x$x), type = "l")
#' abline(v = 10)
#' matplot(t(x$pred), type = "l")
#' abline(v = 10)
#'
#' set.seed(42)
#' x <- sim_dfa()
#' matplot(t(x$x), type = "l")
#' abline(v = 10)
#' matplot(t(x$pred), type = "l")
#' abline(v = 10)
#' @export

sim_dfa_horizon_widevar <- function (num_trends = 1, num_years = 20, num_ts = 4, num_horizons = 2, phi = c(0.8, 0.7)
                                     ,zeta = 0.5
                                     , loadings_matrix = matrix(nrow = num_ts, ncol = num_trends, rnorm(num_ts * num_trends, 0, 1))
                                     , sigma = rlnorm(1, meanlog = log(10), 1)
                                     , sigmaH = diag(num_horizons), rho = 0.8
                                     , varIndx = rep(1, num_ts), n.outlier = NULL
                                     , extreme_loc = NULL, nu_fixed = 100, user_supplied_deviations = NULL) 
{
  y_ignore <- matrix(rnorm(num_ts * num_years), nrow = num_ts, 
                     ncol = num_years)
  d <- fit_dfa(y_ignore, num_trends = num_trends, estimation = "none", 
               scale = "center", varIndx = varIndx, nu_fixed = nu_fixed, 
               trend_model = "rw")
  Z <- loadings_matrix
  y <- vector(mode = "numeric", length = d$sampling_args$data$N)
  for (k in seq_len(d$sampling_args$data$K)) {
    Z[k, k] <- abs(Z[k, k])
  }
  for (k in seq_len(d$sampling_args$data$K)) {
    for (p in seq_len(d$sampling_args$data$P)) {
      if (p < k) 
        Z[p, k] <- 0
    }
  }
  x <- matrix(nrow = d$sampling_args$data$K, ncol = d$sampling_args$data$N)
  for (k in seq_len(d$sampling_args$data$K)) {
    if (!is.null(user_supplied_deviations)) {
      devs <- user_supplied_deviations[, k]
    }
    else {
      devs <- rt(d$sampling_args$data$N, df = d$sampling_args$data$nu_fixed)
    }
    x[k, 1] <- rnorm(1, 0, 1)
    if (is.null(n.outlier)) { #extreme_value
      for (t in seq(2, d$sampling_args$data$N)) {
        x[k, t] <- phi[k]*x[k, t - 1] + devs[t]
      }
    }
    else {
      if(is.null(extreme_loc)){
        extreme_loc <- sample(c(5:(num_years-5)), n.outlier)
        extreme_loc <- sort(extreme_loc)
      }else{
        extreme_loc <- extreme_loc
      }
      
      for(t in 2:(extreme_loc[1]-1)){
        x[k, t] <- phi[k]*x[k, t - 1] + devs[t]
      }
      if(n.outlier > 1){
        for(kk in 1:(n.outlier-1)){
          x[k, extreme_loc[kk]] <- phi[k]*x[k, (extreme_loc[kk] - 1)] *runif(1, 3, 7) + devs[t]
          for(t in seq(extreme_loc[kk]+1,extreme_loc[(kk+1)]-1)){
            x[k, t] <- phi[k]*x[k, t - 1] + devs[t]
          }
        }
      }
      x[k, extreme_loc[n.outlier]] <- phi[k]*x[k, (extreme_loc[n.outlier] - 1)] *runif(1, 3, 7) + devs[t]
      for(t in seq(extreme_loc[n.outlier]+1, d$sampling_args$data$N)){
        x[k, t] <- phi[k]*x[k, t - 1] + devs[t]
      }
    }
  }
  
  if(num_horizons == 1){
    R <- 1
  }else{
    R <- matrix(NA, nrow = num_horizons, ncol = num_horizons)
    diag(R) <- 1
    for(i in 1:(num_horizons-1)){
      for(j in ((i+1):num_horizons)){
        R[i,j] <- R[j,i] <- exp(-0.4+0.025*max(i,j)-0.125*abs(i-j))
      }
    }
  }
  
  
  sigmah <- 1+zeta*sqrt(0:(num_horizons-1))
  sigmaH <- sigma*diag(sigmah)%*%R%*%diag(sigmah)
  
  SigmaHL <- bdiag(replicate(d$sampling_args$data$P, sigmaH, simplify=FALSE))
  
  
  pred <- Z %*% x
  A <- diag(d$sampling_args$data$P)%x%matrix(1, nrow = num_horizons, ncol = 1)
  A.pred <- A%*%pred
  y_sim <- matrix(NA, nrow = dim(A.pred)[1], ncol = dim(A.pred)[2])
  for(i in 1:(d$sampling_args$data$N)){
    y_sim[,i] <- mvrnorm(1, A.pred[,i], SigmaHL)
  }
  
  
  list(y_sim = y_sim, pred = pred, x = x, Z = Z, phi = phi, sigma = sigma, sigmaH = sigmaH, SigmaHL = SigmaHL, extreme_loc = extreme_loc)
}



