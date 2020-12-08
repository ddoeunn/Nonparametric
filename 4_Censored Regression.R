rm(list=ls())
get_h = function(i, k, X){
  if((i+k) > length(X)){
    h = X[i] - X[i-k]
  }else if((i-k) <= 0){
    h = X[i+k] - X[i]
  }else{
    h = (X[i+k] - X[i-k])/2
  }
  return (h)
}

transform_y = function(X, Z, D, k){
  trans_y = c()
  for(i in 1:length(X)){
    h = get_h(i, k, X)
    sum_idx = Z[i] < Z
    if(sum(D[sum_idx]) == 0){
      trans_y[i] = Z[i]
    }else{
      dn = dnorm((X[i] - X[sum_idx]) / (h+1e-20)) * D[sum_idx]
      up = Z[sum_idx] * dn
      trans_y[i] = D[i]*Z[i] + (1-D[i]) * sum(up) / (sum(dn)+1e-20)
    }
  }
  return (trans_y)
}

get_m_x_hat = function(x, X_data, Y_star, k){
  l_idx = which.min(abs(X_data - x))
  K_h = dnorm((x - X_data) / get_h(l_idx, k, X_data))
  s_n2 = sum(K_h * (x - X_data)**2)
  s_n1 = sum(K_h * (x - X_data))
  w_i = K_h * (s_n2 - (x - X_data) * s_n1)
  m_x_hat = sum(w_i * Y_star) / (sum(w_i)+1e-20)
  return (m_x_hat)
}

leave_one_out = function(X, Y_star, D, k){
  cv_k = 0
  for(i in 1:length(Y_star)){
    m_x_i = get_m_x_hat(X[i], X[-i], Y_star[-i], k)
    cv_k = cv_k + sum((Y_star[i] - m_x_i)^2)
  }
  return (cv_k)
}

get_opt_k = function(X, Z, D){
  k_max = as.integer((length(X) - 1) / 2)
  cv_k_vec = rep(0, k_max)
  for(k in 1:k_max){
    cat('k = ', k, '\n')
    Y_star = transform_y(X, Z, D, k)
    cv_k_vec[k] = leave_one_out(X, Y_star, D, k)
  }
  opt_k = which.min(cv_k_vec)
  main = paste('Cross-validation Curve (optimal k = ', opt_k, ')', sep='')
  plot(cv_k_vec, pch=20, xlab='k', ylab="CV(k)", main=main)
  print(cv_k_vec)
  return (opt_k)
}

fit_mx = function(x, X_data, Y_star, k){
  m_X = sapply(x, get_m_x_hat, X_data, Y_star, k)
  return (m_X)
}

# ------------------------------------------------------- Simulation Data
set.seed(2020)
datasize = 200
X = sort(runif(datasize))
Y = 4.5 - 64*X^2*(1-X)^2 - 16*(X -0.5)^2 + 0.25*rnorm(datasize)
cx = I((0 <= X) & (X <= 0.5))*3*(1.25-abs(4*X-1)) + I(X > 0.5)*3*(1.25-abs(4*X-3))
C = rexp(cx)
Z = mapply(min, Y, C)
D = (Y <= C)

k_opt = get_opt_k(X, Z, D)                # get optimal k
trans_y = transform_y(X, Z, D, k_opt)     # transform data
m_X = fit_mx(X, X, trans_y, k_opt)        # prediction
m_X2 = fit_mx(X, X, Z, k_opt)

# ------------------------------------------------------- Visualization
plot(X, trans_y, pch=as.integer(!D)+1, col=as.integer(!D)+3,
     cex=.7, ylab='y*', xlab='x', main = 'Transformed by Local Average',
     ylim=c(min(Y) - 0.5, max(Y) + 1))
lines(X, 4.5 - 64*X^2*(1-X)^2 - 16*(X -0.5)^2)
lines(X, m_X, col=6)
lines(X, m_X2, col=7)

plot(X, Y, ylim=c(min(Y) - 0.5, max(Y) + 1), cex=0.5, xlab='x', ylab='y',
     main='Unobserved Simulated Data')
plot(X, Z, pch=as.integer(!D)+1, col=as.integer(!D)+3, cex=.5, 
     ylab='z', xlab='x',
     ylim=c(min(Y) - 0.5, max(Y) + 1), main='Observed Simulated Data')
legend('topright', legend=c('non-censored', 'censored'), 
       col=c(3, 4), pch=c(1, 2), bty = "n")

# ------------------------------------------------------- Real Data
library(survival)
stanford2 = survival::stanford2
data = stanford2[!is.na(stanford2$t5),]
x_vec = data$age
z_vec = log1p(sqrt(data$time))
d_vec = data$status

opt_k = get_opt_k(x_vec, z_vec, d_vec)            # get optimal k
y_star = transform_y(x_vec, z_vec, d_vec, opt_k)  # transform data
m_X = fit_mx(x_vec, x_vec, y_star, opt_k)         # prediction
m_X_2 = fit_mx(x_vec, x_vec, z_vec, opt_k)

# ------------------------------------------------------- Visualization
plot(x_vec, z_vec, pch=as.integer(!d_vec)+1, col=as.integer(!d_vec)+3,
     cex=.5, ylab='log(Time)', xlab='x', ylim=c(-1, 6),
     main = 'Stanford Heart Transplant Data')

plot(x_vec, y_star, pch=as.integer(!d_vec)+1, col=as.integer(!d_vec)+3,
     cex=.5, ylab='y*', xlab='x', ylim=c(-1, 6),
     main = 'Transformed by Local Average')
slope = lm(z_vec[x_vec >= 48]~x_vec[x_vec >= 48])$coef[2]
y_true =  mean(z_vec[x_vec < 48]) + slope*(x_vec-48)*(x_vec>48)
lines(x_vec, y_true)
lines(x_vec, m_X, col=6)
lines(x_vec, m_X_2, col=7)




