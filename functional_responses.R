# writes functional responses used in tree hole model

# density-dependent mortality
dd_mu = function(mu_L, K_L, x1, x2) mu_L*x1*(1 + (x1 + x2)/K_L)

# density-dependent cell division
dd_div = function(delta, K_R, R) delta*R*(1-R/K_R)

# temperature-dependent functions
briere = function(x, c, T0, Tm){
  if((x < T0) | (x > Tm))
    0.0
  else
    c*x*(x-T0)*sqrt(Tm-x)
}

quadratic = function(x, c, T0, Tm, offset = 0.0){
  if((x < T0) | (x > Tm))
    offset
  else
    c*(x-T0)*(x-Tm)
}

inverted_quadratic = function(x, c, T0, Tm, offset = 0.33){
  if((x < T0) | (x > Tm))
    offset
  else
    1.0/(c*(x-T0)*(x-Tm))
}

truncated_linear = function(x, m, z, trunc.temp, trunc.val){
  if(x < trunc.temp)
    trunc.val
  else if (x > z/m)
    0.0
  else 
    -m*x+z
}

# trait thermal response functions
# using parameter values fitted to Cx. tarsalis

# larval mortality (mu_L)
mu_L = function(temp){
  -log(quadratic(temp, -2.12e-03, 5.9, 43.1, offset = 0.01))/30
}

x = c(0:45)
mu_plot = c()
for (i in 1:length(x)) mu_plot[i] = mu_L(x[i])
plot(x, 1/mu_plot, type = "l")

# larval development rate (d_L)
d_L = function(temp){
  briere(temp, 4.12e-05, 4.3, 39.9)
}

x = c(0:45)
dL_plot = c()
for (i in 1:length(x)) dL_plot[i] = d_L(x[i])
plot(x, dL_plot, type = "l")

# trophont mortality rate (mu_R) and theront mortality rate (mu_H)
# this one is totally made up, for now
mu_R = function(temp){
  inverted_quadratic(temp, -2.42e-01, 0, 32, offset = 5)*10
}

mu_H = function(temp){
  inverted_quadratic(temp, -2.42e-01, 0, 32, offset = 5)*10
}

x = c(0:45)
mr_plot = c()
for (i in 1:length(x)) mr_plot[i] = mu_R(x[i])
plot(x, 1/mr_plot, type = "l")

# cell division rate (delta)
# this one is totally made up too
delta = function(temp){
  briere(temp, 1e-03, 0, 32)
}

x = c(0:45)
delta_plot = c()
for (i in 1:length(x)) delta_plot[i] = delta(x[i])
plot(x, delta_plot, type = "l")

# consumption rate of trophonts by larvae (f_LR)
# totally made up function
f_LR = function(temp){
  briere(temp, 2e-06, 10, 32)
}

x = c(0:45)
fLR_plot = c()
for (i in 1:length(x)) fLR_plot[i] = f_LR(x[i])
plot(x, fLR_plot, type = "l")
