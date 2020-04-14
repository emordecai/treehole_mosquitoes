# Simulate mosquito-ciliate model

# Load packages
library(deSolve)
library(RColorBrewer)

# Load model
source("mosquito-ciliate_model.R")

# Initial conditions
inits = c(1000, 0, 0, 0, 100, 0)

# Parameters
parameters = c(
  #'q' = 1,
  #'p_EL' = 0.9,
  #'c_L' = 0.1,
  'd_L' = 1/100,
  'beta' = 0.001,
  'bc' = 0.5,
  #'d_P' = 1/5,
  'mu_L' = 1/100,
  #'mu_P' = 1/100,
  #'mu_A' = 1/20,
  'mu_R' = 1/4,
  'mu_H' = 1/2,
  'delta' = 1.5,
  'n' = 0.25,
  'm' = 50,
  'alpha_L' = 1/100,
  'K_R' = 1000,
  'K_L' = 500,
  'f_LR' = 0.001
  #'a_h' = 0.05
)

# Run simulation
times = seq(1,365, by = 0.1)
state = inits
names(state) = c("Lu", "Li","Pu", "Pi", "R", "H")
out = ode(y = state, times = times, func = ode_larvae, parms = parameters, method = "rk4")
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
