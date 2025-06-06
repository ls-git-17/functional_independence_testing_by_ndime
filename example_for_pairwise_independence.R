# load functions for independence tests
source("functional_ndime_mdime.R")

# load needed packages
library(mvtnorm)
library(Matrix)

# simulated data generation as in Case 1, Setting 2, Scenario 1 in Section 3 of the paper

# hyperparameters for simulated data
# number of functional vectors
kk <- 2
# dimensions of functional data (d1, d2)
pp <- c(3, 3)
# number of observations
nn <- 20
# time points
ss <- 50
time <- seq(0, 1, length.out = ss)
# numbers of basis functions in basis expansion (K1, ..., Kp1, L1, ..., Lp2)
klm <- vector("list", kk)
for (i_klm in seq_len(kk)) klm[[i_klm]] <- rep(5, pp[i_klm])
p <- sum(klm[[1]])
# basis values
fourier_basis <- create.fourier.basis(rangeval = c(min(time), max(time)), 5)
basis_values <- eval.basis(time, fourier_basis)

# generating simulated data
sigma_mat <- diag(p)
x <- rmvnorm(nn, rep(0, p), sigma_mat)
y <- x^2
alpha_beta_gamma <- cbind(x, y)
x_list <- vector("list", pp[1])
y_list <- vector("list", pp[2])
for (i_x_list in seq_len(pp[1])) {
  x_list[[i_x_list]] <- matrix(0, nrow = ss, ncol = nn)
}
for (i_y_list in seq_len(pp[2])) {
  y_list[[i_y_list]] <- matrix(0, nrow = ss, ncol = nn)
}
for (i_nn in seq_len(nn)) {
  data_i_nn <- matrix(0, nrow = sum(pp), ncol = ss)
  for (i_ss in seq_len(ss)) {
    Phi <- vector("list", sum(pp))
    for (i_Phi in seq_len(length(Phi))) {
      Phi[[i_Phi]] <- basis_values[i_ss, ]
    }
    Phi_mat <- t(as.matrix(bdiag(Phi)))
    data_i_nn[, i_ss] <- as.numeric(Phi_mat %*% alpha_beta_gamma[i_nn, ])
  }
  for (i_x_list in seq_len(pp[1])) {
    rx <- max(data_i_nn[i_x_list, ]) - min(data_i_nn[i_x_list, ])
    x_list[[i_x_list]][, i_nn] <- data_i_nn[i_x_list, ] + rnorm(ss, 0, 0.025 * rx)
  }
  for (i_y_list in seq_len(pp[2])) {
    ry <- max(data_i_nn[i_y_list + pp[1], ]) - min(data_i_nn[i_y_list + pp[1], ])
    y_list[[i_y_list]][, i_nn] <- data_i_nn[i_y_list + pp[1], ] + rnorm(ss, 0, 0.025 * ry)
  }
}

# plot for trajectories of simulated data
par(mfrow = c(2, 3), mar = c(4, 4, 2, 0.1))
for (i in 1:3) {
  matplot(x_list[[i]], type = "l", col = 1, lty = 1, xlab = "Time", ylab = "Value", main = paste("Variable X1", i, sep = "_"))
}
for (i in 1:3) {
  matplot(y_list[[i]], type = "l", col = 1, lty = 1, xlab = "Time", ylab = "Value", main = paste("Variable X2", i, sep = "_"))
}

# basis expansion of simulated data in B-spline basis
x_basis <- basis_expansion_bspline(x_list, Bk = klm[[1]], int = range(time), time = time)
y_basis <- basis_expansion_bspline(y_list, Bk = klm[[2]], int = range(time), time = time)
x <- t(x_basis)
y <- t(y_basis)

# functional NDIME permutation test
perm_test_ndime(x, y, NDIME, n_perm = 1000)
# functional mNDIME permutation test
perm_test_ndime_m(x, y, NDIME, n_perm = 1000)
