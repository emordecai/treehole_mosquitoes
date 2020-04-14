# Simulate mosquito-ciliate model

# Load packages
library(deSolve)
library(RColorBrewer)

# Fixed (temperature-independent) parameters
fixed.pars = c(
  'K_L' = 500,
  'K_R' = 1000,
  'beta' = 0.001,
  'bc' = 0.5,
  'n' = 0.25,
  'm' = 50,
  'alpha_L' = 1/100
)

# Load model
source("mosquito-ciliate_model_T.R")

# Initial conditions
inits = c(1000, 0, 0, 0, 100, 0)

# Simulation time and temperature sequence
times = seq(1,365, by = 0.1)
sinfun = function(times, trange = 10, tmean = 20) trange*sin(2*pi*(times + 365*3/4)/365) + tmean
temp = sinfun(times, 5, 20)
# plot(times, temp, type = "l")

# Parameters
parameters = c(
  d_L,
  mu_L,
  mu_R,
  mu_H,
  delta,
  f_LR
)

# Run simulation
state = inits
names(state) = c("Lu", "Li","Pu", "Pi", "R", "H")
out = ode(y = state, times = times, func = ode_larvae_T, parms = parameters, method = "rk4")
out2 = as.data.frame(out)
head(out2)
# pdf(file = "output/sim_output.pdf", width = 6, height = 6)
cols = brewer.pal(6, "Paired")
matplot(out2[,1], out2[,2:7], type = 'l', xlab = 'day', ylab = 'number',
        lty = c(rep(1,4), rep(2,2)), col = cols, lwd = 2)
legend('topright', legend = names(state), lty = c(rep(1,4), rep(2,2)), 
       col = cols, bty = 'n', lwd = 2)
# dev.off()

### Calculate some outputs of interest
# Final prevalence at each life stage
lp = with(tail(out2,1),
     c(
       'pupal prevalence' = Pi/(Pi+Pu),
       'theront prevalence' = H/(H+R),
       'total pupae' = Pi+Pu
       ))

print(lp)
