rm(list=ls())
library(mvtnorm);library(ks);library(plot3D)
# ----------------------------------------------------------- Functions
norm2.mix = function(x, y, mu1, sig1, mu2, sig2, p){
  t = cbind(x, y)
  tmp = p*dmvnorm(t, mu1, sig1) + (1-p)*dmvnorm(t, mu2, sig2)
  return (tmp)
}
trapezoidal_integ_2d = function(fxy, x_width=0.05, y_width=0.05){
  n = dim(fxy)[1]
  fxy_sum = sum(fxy[1:(n-1), 1:(n-1)] + fxy[2:n, 1:(n-1)]+ fxy[1:(n-1), 2:n] + fxy[2:n, 2:n])
  return (1/4 * x_width * y_width * fxy_sum)
}
# ----------------------------------------------------------- Simulation Settings
s_num = 500
x1_grid = seq(-3,3,0.05)
x2_grid = seq(-3,3,0.05)
x_grids = expand.grid(x1_grid, x2_grid)
n_vec = c(200, 400, 800)
S = matrix(c(1, 1, 1, 2),2,2)
m1 = c(-1, -1)
m2 = c(1, 1)
s1 = diag((2/3)^2, 2)
s2 = (2/3)^2*S
p = 1/3
# ----------------------------------------------------------- For Results
n_grid = length(x1_grid)
fx_true = outer(x1_grid, x2_grid, norm2.mix, mu1=m1, sig1=s1, mu2=m2, sig2=s2, p=p)
fx_trues = cbind(fx_true, fx_true, fx_true)
idx_list = list(1:n_grid, (n_grid+1):(n_grid*2), (n_grid*2+1):(n_grid*3))
x_idx1 = which((x_grids[, 1] == -1) + (x_grids[, 2] == 1) == 2)
x_idx2 = which((x_grids[, 1] == 0) + (x_grids[, 2] == 0) == 2)
x_idx3 = which((x_grids[, 1] == 1) + (x_grids[, 2] == 1) == 2)
bw_methods = c('h*I', 'h*diag(1,2)', 'h*S')
# ----------------------------------------------------------- Simulation 
rst_bias = rst_var = rst_f = rst_est = NULL
for(n_pos in 1:length(n_vec)){
  n = n_vec[n_pos]
  cat('n=', n, '\n')
  h= n^(-1/6)
  bw_list = list(diag(h, 2), h*diag(c(1,2)), h*S)
  est_500 = NULL
  fx_hat_500 = f_hat_sq_500 = matrix(0, n_grid*n_grid, 3)
  for(s_pos in 1:s_num){
    if(s_pos %% 10 == 0){print(s_pos)}
    I = rbinom(n, 1, p)
    X = I*rmvnorm(n, m1, s1) + (1-I)*rmvnorm(n, m2, s2)
    fx_hat_mat = est_vec = NULL
    for(bw_pos in 1:3){
      d = kde(X, bw_list[[bw_pos]], eval.points = x_grids)$estimate
      fx_hat_mat = cbind(fx_hat_mat, d)
      est_vec = c(est_vec, c(d[x_idx1], d[x_idx2], d[x_idx3]))
    }
    fx_hat_500 = fx_hat_500 + fx_hat_mat
    f_hat_sq_500 = f_hat_sq_500 + fx_hat_mat^2
    est_500[[s_pos]] = est_vec
  }
  e_f_hat = matrix(fx_hat_500 / s_num, n_grid, n_grid*3)
  e_f_hat_sq = matrix(f_hat_sq_500 / s_num, n_grid, n_grid*3)
  rst_f[[n_pos]] = e_f_hat
  rst_bias[[n_pos]] = (e_f_hat - fx_trues) ^ 2
  rst_var[[n_pos]] =  e_f_hat_sq - (e_f_hat ^ 2)
  rst_est[[n_pos]] = do.call(rbind, est_500)
}
# ----------------------------------------------------------------- Result 1,2,3.
result_bias = result_var = matrix(0, length(n_vec), length(bw_methods), 
                                  dimnames = list(paste("n=", n_vec), bw_methods))
par(mfrow=c(1, 2))
for(i in 1:length(n_vec)){
  persp3D(x1_grid, x2_grid, fx_true, main="fx_true")
  for(j in 1:length(bw_methods)){
    persp3D(x1_grid, x2_grid, rst_f[[i]][, idx_list[[j]]], main=paste('n=', n_vec[i], bw_methods[j]))
    result_bias[i, j] = trapezoidal_integ_2d(rst_bias[[i]][, idx_list[[j]]])
    result_var[i, j] = trapezoidal_integ_2d(rst_var[[i]][, idx_list[[j]]])
  }
}
result_bias;result_var
# ----------------------------------------------------------------- Result 4.
est_name = as.vector(sapply(bw_methods, paste, c('x=(-1, 1)', 'x=(0, 0)', 'x=(1, 1)')))
x_mat = matrix(c(-1, 0, 1, 1, 0, 1), 3, 2)
par(mfrow=c(3, 3))
for(x in 1:dim(x_mat)[1]){
  fx = norm2.mix(x_mat[x, 1], x_mat[x, 2], mu1=m1, sig1=s1, mu2=m2, sig2=s2, p=p)
  for(n in 1:length(n_vec)){
    for(i in c(0, 3, 6)){
      idx = x+i
      max = ifelse(max(rst_est[[n]][, idx])<fx, fx, max(rst_est[[n]][, idx]))
      min = ifelse(min(rst_est[[n]][, idx])>fx, fx, min(rst_est[[n]][, idx]))
      hist(rst_est[[n]][, idx], main = paste('n=', n_vec[n], est_name[idx]), xlab='', xlim=c(max, min))
      abline(v=fx, col='red')
    }
  }
}
