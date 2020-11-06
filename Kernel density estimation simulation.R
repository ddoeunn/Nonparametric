rm(list=ls())
norm_mix = function(x, mu1=-1, mu2=1, sig1=2/3, sig2=2/3, p=1/3){
  return (p*dnorm(x,mu1,sig1)+(1-p)*dnorm(x,mu2,sig2)) 
}
norm_deri_2 = function(x, mu1=-1, mu2=1, sig1=2/3, sig2=2/3, p=1/3){
  tmp1 = p/sig1^3 * exp(-(x-mu1)^2/(2*sig1^2))*(1-(x-mu1)^2/(sig1^2))
  tmp2 = (1-p)/sig2^3 * exp(-(x-mu2)^2/(2*sig2^2))*(1-(x-mu2)^2/(sig2^2))
  return(-(tmp1+tmp2)/sqrt(2*pi))
}
trapezoidal_integ = function(y, width=0.05){
  n = length(y)
  return (sum((y[1:(n-1)] + y[2:n]) * width / 2))
}
# ------------------------------------------------------------------------- Simulation Setting
x_grid = seq(-3,3,0.05) # grid points
x_vec = c(-1, 0, 1)
bandwidths = list(0, 'nrd0', 'ucv')
n_vec = c(200, 400, 800)
s_num = 500
m1 = -1; m2 = 1; sd1 = 2/3; sd2 = 2/3; p = 1/3
n_grid = length(x_grid) 
fx_true = norm_mix(x_grid) # true density f(x)
n_bw = length(bandwidths)
integ_f_deri_2 = integrate(function(x){norm_deri_2(x)^2}, -Inf, Inf)$value 
h_opts = (1/(2*sqrt(pi)) / integ_f_deri_2) ^ (1/5) * n_vec^(-1/5)
# ------------------------------------------------------------------------- Simulation
set.seed(0)
par(mfrow=c(1, 2))
rst_bias = rst_var = rst_f = rst_est = NULL
for(n_pos in 1:length(n_vec)){
  n = n_vec[n_pos]
  bandwidths[[1]] = h_opts[n_pos]
  cat('n=', n, '\n')
  fx_hat_list = f_hat_sq_list = est_list = bw_list = NULL
  for(s_pos in 1:s_num){
    if(s_pos %% 50 == 0){print(s_pos)}
    I = rbinom(n, 1, p)
    data = I * rnorm(n, m1, sd1) + (1-I) * rnorm(n, m2, sd2) # Scenario 2
    fx_hat_mat = est_mat = bw_mat = NULL
    for(bw_pos in 1:n_bw){
      d = density(data, bw=bandwidths[[bw_pos]], kernel='gaussian', n=n_grid, from=-3, to=3)
      fx_hat = d$y
      est = fx_hat[c(which(d$x == -1), which(d$x == 0), which(d$x == 1))]
      fx_hat_mat = cbind(fx_hat_mat, fx_hat)
      est_mat = c(est_mat, est)
      bw_mat = cbind(bw_mat, d$bw)
    }
    fx_hat_list[[s_pos]] = fx_hat_mat
    f_hat_sq_list[[s_pos]] = fx_hat_mat ^ 2
    est_list[[s_pos]] = est_mat
    bw_list[[s_pos]] = bw_mat
  }
  e_f_hat = Reduce('+', fx_hat_list) / s_num
  e_f_hat_sq = Reduce('+', f_hat_sq_list) / s_num
  rst_bias[[n_pos]] = (e_f_hat - fx_true) ^ 2
  rst_var[[n_pos]] =  e_f_hat_sq - (e_f_hat ^ 2)
  rst_f[[n_pos]] = e_f_hat
  rst_est[[n_pos]] = do.call(rbind, est_list)
  # ------------------------------------------------------------------------- Result 5. 
  rst_bw = do.call(rbind, bw_list)
  hist(rst_bw[, 2], main = paste('n=', n,'bw selector=',  bandwidths[[2]]), xlab='', xlim=c(h_opts[n_pos], max(rst_bw[, 2])))
  abline(v=h_opts[n_pos], col='red')
  hist(rst_bw[, 3], main = paste('n=', n,'bw selector=',  bandwidths[[3]]), xlab='')
  abline(v=h_opts[n_pos], col='red')
}
# ------------------------------------------------------------------------- Result 1, 2 
bw_methods = c('h_opt', 'nrd0', 'ucv')
result_bias = t(sapply(rst_bias, function(x){apply(x, 2, trapezoidal_integ)}))
result_var = t(sapply(rst_var, function(x){apply(x, 2, trapezoidal_integ)}))
rownames(result_bias) = rownames(result_var) = paste('n=', n_vec, sep='')
colnames(result_bias) = colnames(result_var) = bw_methods
result_bias; result_var
# ------------------------------------------------------------------------- Result 3. 
par(mfrow=c(1, 3))
for(i in 1:length(n_vec)){
  for(j in 1:n_bw){
    title = paste('n=', n_vec[i],  bw_methods[j])
    plot(x_grid, rst_f[[i]][, j], type='l', main=title, ylab='', ylim = c(0, .4))
    lines(x_grid, fx_true, col=2)
  }
}
# ------------------------------------------------------------------------- Result 4.
est_name = as.vector(sapply(bw_methods, paste, paste('(x=', x_vec, ')', sep=''))) 
fx_vec = norm_mix(x_vec)
for(x in 1:length(x_vec)){
  for(n in 1:length(n_vec)){
    for(i in c(0, 3, 6)){
      idx = x+i
      hist(rst_est[[n]][, idx], main = paste('n=', n_vec[n], est_name[idx]), xlab='')
      abline(v=fx_vec[x], col='red')
    }
  }
}
