rm(list=ls())
norm_mix = function(x, mu1=-1, mu2=1, sig1=2/3, sig2=2/3, p=1/3){
  return (p*dnorm(x,mu1,sig1)+(1-p)*dnorm(x,mu2,sig2)) # pdf
}

get_breaks = function(data, nclass){
  labs = levels(cut(data, nclass, dig.lab=10))
  lower = as.numeric(sub("\\((.+),.*", "\\1", labs))[1]
  upper = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", labs))               
  return (c(lower, upper))
}

den_hist = function(x, object){
  br = object$breaks
  idx = sum(br < x)
  ret = ifelse(idx == 0 | idx == length(br), 0, object$density[idx])
  return (ret)
}

trapezoidal_integ = function(y){
  n = dim(y)[1]
  ret = colSums((y[1:(n-1), ] + y[2:n, ]) * 0.05 / 2)
  return (ret)
}
# ------------------------------------------------------------------------- Simulation Setting
s_num = 500
x_grid = seq(-3,3,0.05) # grid points
n_grid = length(x_grid) 
x_vec = c(-1, 0, 1)
n_vec = c(200, 400, 800)
bin_width = c('Sturges', 'Scott', 'FD')
n_bin_method = length(bin_width)
m1 = -1
m2 = 1
sd1 = sd2 =2/3
p = 1/3
fx_true = norm_mix(x_grid) # true density f(x)
est_name = as.vector(sapply(bin_width, paste, paste('(x=', x_vec, ')', sep=''))) # for problem 4.
# ------------------------------------------------------------------------- Simulation
set.seed(2020)
result_bias = result_var = result_f = result_m = result_est =  NULL
for(n_pos in 1:length(n_vec)){
  f_hat_list = f_hat_sq_list = est_list = mse_list =NULL 
  n = n_vec[n_pos]
  cat('n=', n, '\n')
  for(s_pos in 1:s_num){
    if(s_pos %% 50 == 0){print(s_pos)}
    I = rbinom(n, 1, p)
    data = I * rnorm(n, m1, sd1) + (1-I) * rnorm(n, m2, sd2)
    nclass = c(nclass.Sturges(data), nclass.scott(data), nclass.FD(data))
    f_hat = est = NULL
    for(i in 1:n_bin_method){
      breaks = get_breaks(data, nclass[i])
      hist = hist(data, breaks = breaks, plot=F)
      est = c(est, sapply(x_vec, den_hist, hist))
      f_hat = cbind(f_hat, sapply(x_grid, den_hist, hist))
    } 
    mse_list[[s_pos]] = (f_hat - fx_true) ^ 2
    f_hat_list[[s_pos]] = f_hat
    est_list[[s_pos]] = est
    f_hat_sq_list[[s_pos]] = f_hat ^ 2
  }
  e_f_hat = Reduce('+', f_hat_list) / s_num
  e_f_hat_sq = Reduce('+', f_hat_sq_list) / s_num
  result_m[[n_pos]] = Reduce('+', mse_list) / s_num 
  result_f[[n_pos]] = e_f_hat
  result_bias[[n_pos]] = (e_f_hat - fx_true) ^ 2
  result_var[[n_pos]] =  e_f_hat_sq - (e_f_hat ^ 2)
  result_est[[n_pos]] = do.call(rbind, est_list)
}
# ------------------------------------------------------------------------- Result 1, 2
bias_result = t(sapply(result_bias, trapezoidal_integ))
var_result = t(sapply(result_var, trapezoidal_integ))
mse_result = t(sapply(result_m, trapezoidal_integ))
colnames(mse_result) = colnames(bias_result) = colnames(var_result) = bin_width
rownames(mse_result) = rownames(bias_result) = rownames(var_result) = paste('n=', n_vec)
bias_result # 1. 
var_result # 2.
mse_result; round(bias_result + var_result, 10) == round(mse_result, 10) # check (1+2)
# ------------------------------------------------------------------------- Result 3
par(mfrow=c(3, 3))
for(i in 1:length(n_vec)){ # 3.
  for(j in 1:n_bin_method){
    title = paste('n=', n_vec[i],  bin_width[j])
    plot(x_grid, result_f[[i]][, j], type='l', main=title, ylab='', ylim = c(0, .4))
    lines(x_grid, fx_true, col=2)
  }
}
# ------------------------------------------------------------------------- Result 4
fx_vec = norm_mix(x_vec)
for(x in 1:length(x_vec)){
  for(n in 1:length(n_vec)){
    for(i in c(0, 3, 6)){
      idx = x+i
      hist(result_est[[n]][, idx], main = paste('n=', n_vec[n], est_name[idx]), xlab='')
      abline(v=fx_vec[x], col='red')
    }
  }
}

