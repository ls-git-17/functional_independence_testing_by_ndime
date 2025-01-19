# load needed packages
library(fda)


# basis expansion in Fourier and B-spline bases ---------------------------

# arguments:
#  x - a list of m\times n data.frames of data, whose each column is a discretized version 
#      of a function and rows correspond to design time points. The ith element of this list 
#      contains the data of ith feature, i = 1, ..., d.
#  Bk - a vector of numbers of basis functions in basis expansion of the data
#  int - a vector of two elements representing the interval [a,b]
#  time - a vector of design time points
# values:
#  matrix of coefficients of size sum(Bk)\times n

basis_expansion_fourier <- function(x, Bk, int, time) {
  nn <- ncol(x[[1]])
  pp <- length(x)
  if (pp != length(Bk)) {
    stop("length of Bk should be the same as length of x")
  }
  
  data_fd_F <- NULL
  for (ii_x_fd in seq_len(pp)) {
    fbasis1 <- create.fourier.basis(rangeval = int, Bk[ii_x_fd])
    data_fd_F <- rbind(data_fd_F, smooth.basis(time, x[[ii_x_fd]], fbasis1)$fd$coefs)
  }
  
  return(data_fd_F)
}

basis_expansion_bspline <- function(x, Bk, int, time) {
  nn <- ncol(x[[1]])
  pp <- length(x)
  if (pp != length(Bk)) {
    stop("length of Bk should be the same as length of x")
  }
  
  data_fd_B <- NULL
  for (ii_x_fd in seq_len(pp)) {
    fbasis2 <- create.bspline.basis(rangeval = int, Bk[ii_x_fd])
    data_fd_B <- rbind(data_fd_B, smooth.basis(time, x[[ii_x_fd]], fbasis2)$fd$coefs)
  }
  
  return(data_fd_B)
}


# NDIME and mNDIME tests --------------------------------------------------


# basic functions ---------------------------------------------------------

# K-matrix for estimator of NDIME
# arguments:
#  x - n\times d matrix of n observations of d-dimensional random vector
# values:
#  a K matrix in the formula for estimator of NDIME
KK <- function(x) {
  odl <- dist(x)
  sigma <- median(odl)
  
  l_return <- array(0, c(nrow(x), nrow(x)))
  const_1 <- 2 * (sigma^2) / (1 + sigma^2)
  const_2 <- 2 * (sigma^4) / (1 + sigma^2) + 2 * const_1
  
  norm_2 <- as.matrix(odl^2)
  il_sk <- x %*% t(x)
  
  l_return <- exp(-(norm_2 + const_1 * il_sk) / const_2)
  
  return(l_return)
}

# NDIME test statistic
# arguments:
#  x - n\times d_1 matrix of n observations of d_1-dimensional random vector
#  y - n\times d_2 matrix of n observations of d_2-dimensional random vector
#  HH - a centering matrix H
# values:
#  value of NDIME test statistic
NDIME <- function(x, y, HH) {
  return(sum(diag(KK(x) %*% HH %*% KK(y) %*% HH)) / nrow(x)^2)
}

# marginal test statistic based on test statistic implemented in func argument
# arguments:
#  x - n\times d_1 matrix of n observations of d_1-dimensional random vector
#  y - n\times d_2 matrix of n observations of d_2-dimensional random vector, y has to be different than x
#  func - a function for calculating the value of a test statistic, e.g., NDIME()
#  ... - additional arguments to func function
# values:
#  value of marginal test statistic
marginal_test_statistic <- function(x, y, func, ...) {
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)
  if (n != nrow(y)) {
    stop("different number of observations in x and y")
  }
  
  temp <- 0
  for (i in seq_len(p)) {
    for (j in seq_len(q)) {
      temp <- temp + func(matrix(x[, i], ncol = 1), matrix(y[, j], ncol = 1), ...)
    }
  }
  
  return(temp)
}


# functions for pairwise independence testing -----------------------------


# permutation test based on test statistic given in func argument
# arguments:
#  x - n\times d_1 matrix of n observations of d_1-dimensional random vector
#  y - n\times d_2 matrix of n observations of d_2-dimensional random vector
#  func - a function for calculating the value of a test statistic, e.g., NDIME()
#  n_perm - a number of permutation samples used
# values:
#  p-value of a test
perm_test_ndime <- function(x, y, func, n_perm = 1000) {
  n <- nrow(x)
  HH <- diag(1, n) - (array(1, n) %*% t(array(1, n))) / n
  test_statistic <- func(x, y, HH)
  
  results <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    y_perm <- y[sample(seq_len(n), n), ]
    results[i] <- func(x, y_perm, HH)
  }
  
  p_value <- mean(results > test_statistic)
  
  return(p_value)
}

# permutation test for marginal version of the given test statistic based on 
# test statistic given in func argument
# arguments:
#  x - n\times d_1 matrix of n observations of d_1-dimensional random vector
#  y - n\times d_2 matrix of n observations of d_2-dimensional random vector
#  func - a function for calculating the value of a test statistic, e.g., NDIME()
#  n_perm - a number of permutation samples used
# values:
#  p-value of a test
perm_test_ndime_m <- function(x, y, func, n_perm = 1000)
{
  n <- nrow(x)
  HH <- diag(1, n) - (array(1, n) %*% t(array(1, n))) / n
  test_statistic <- marginal_test_statistic(x, y, func, HH)
  
  results <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    y_perm <- y[sample(seq_len(n), n), ]
    results[i] <- marginal_test_statistic(x, y_perm, func, HH)
  }
  
  p_value <- mean(results > test_statistic)
  
  return(p_value)
}


# functions for mutual independence testing -------------------------------


# R and S methods
# arguments:
#  x - a list of n\times m matrices of vector data. The ith element of this list 
#      contains the data of ith group, i = 1,..., L.
#  func - function for test statistic, e.g., NDIME()
# values:
#  R or S

r_method_ndime <- function(x, func) {
  kk <- length(x)
  
  temp_fun <- matrix(0, nrow = kk - 1, ncol = 2)
  for (c in seq_len(kk - 1)) {
    x2 <- x[[c+1]]
    if (c < kk - 1) {
      for (i_c in (c + 2):kk) {
        x2 <- cbind(x2, x[[i_c]])
      }
    }
    temp_fun[c, ] <- func(x[[c]], x2)
  }
  
  return(colMeans(temp_fun^2))
}

s_method_ndime <- function(x, func) {
  kk <- length(x)
  
  temp_fun <- matrix(0, nrow = kk, ncol = 2)
  for (c in seq_len(kk)) {
    temp_list <- x[-c]
    x2 <- temp_list[[1]]
    for (i_c in 2:(kk - 1)) {
      x2 <- cbind(x2, temp_list[[i_c]])
    }
    temp_fun[c, ] <- func(x[[c]], x2)
  }
  
  return(colMeans(temp_fun^2))
}

# test statistics for perm_test_r_s_ndime() function given below
# arguments:
#  x - n\times d_1 matrix of n observations of d_1-dimensional random vector
#  y - n\times d_2 matrix of n observations of d_2-dimensional random vector
# values:
#  a vector of two test statistics NDIME and mNDIME
test_statistics_ndime <- function(x, y) {
  n <- nrow(x)
  HH <- diag(1, n) - (array(1, n) %*% t(array(1, n))) / n
  test_stat_ndime <- c(NDIME(x, y, HH), marginal_test_statistic(x, y, NDIME, HH))
  
  return(test_stat_ndime)
}

# permutation test for more than two groups
# arguments:
#  x - a list of n\times m matrices of vector data. The ith element of this list 
#      contains the data of ith group, i = 1,..., L.
#  meth - a function for R or S method (one of the functions r_method_ndime and 
#         s_method_ndime)
#  func - function for test statistics, e.g., test_statistics_ndime()
#  n_perm - a number of permutation samples used
# values:
#  p-values of tests
perm_test_r_s_ndime <- function(x, meth, func, n_perm = 1000) {
  kk <- length(x)
  nn <- nrow(x[[kk]])
  T_0 <- meth(x, func)
  
  T_perm <- matrix(0, nrow = n_perm, ncol = 2)
  x_perm <- x
  for (i in seq_len(n_perm)) {
    for (j in 2:kk) {
      x_perm[[j]] <- x[[j]][sample(nn), ]
    }
    T_perm[i, ] <- meth(x_perm, func)
  }
  
  return(colMeans(T_perm > matrix(T_0, nrow = n_perm, ncol = 2, byrow = TRUE)))
}
