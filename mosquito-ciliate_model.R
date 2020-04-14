# Model of mosquito - ciliate interactions
# within a season and within a tree hole
# as a function of temperature

source('functional_responses.R')

### Mosquito stages
# L_u: uninfected larvae (all stages, for now)
# L_i: infected larvae (all stages, for now)
# P_u: uninfected pupae
# P_i: infected pupae

### Ciliate stages
# R: trophonts (free-living non-parasitic ciliates)
# H: theronts (free-living parasitic ciliates)
# C: cysts (dessication-resistant resting stages)

### Parameters or functions
# mu_L, mu_R, mu_H: temperature-dependent background mortality rates for each life stage
# alpha_L: additional parasite-induced mortality of infected larvae
# K_R: carrying capacity of trophonts within a tree hole
# K_L: carrying capacity of mosquito larvae within a tree hole
# d_L: temperature-dependent larval development rate
# beta: contact rate between larvae and theronts
# delta: temperature-dependent cell division rates for free-living ciliates
# n: fraction of daughter cells that are theronts
# m: number of trophonts produced per killed infected larva
# f_LR: linear consumption rate of larval mosquitoes on trophonts, temperature-dependent



# a model that tracks larval and pupal dynamics
ode_larvae = function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    dLu = -dd_mu(mu_L, K_L, Lu, Li) - beta*H*Lu - d_L*Lu  
    dLi = -dd_mu(mu_L, K_L, Li, Lu) - alpha_L*Li + beta*H*Lu - d_L*Li
    dR = dd_div(delta, K_R, R) - n*R + (alpha_L + mu_L)*Li*m - mu_R*R - (Lu + Li)*f_LR*R
    dH = n*R - mu_H*H - beta*H*(Lu+Li)
    dPu = d_L*Lu
    dPi = d_L*Li
    list(c(dLu, dLi, dPu, dPi, dR, dH))
  })
}



### Key assumptions embedded in these equations:
# Transmission only occurs in larval stage
# Larvae and pupae do not clear infection
# Trophont cell division is logistic
# Infection is intensity-independent
# Parasite-killed larvae produce trophonts
# Larvae are an external (between-season) input
# Larval mortality rate is density-dependent
# Carrying capacity and consumption rate are independent of trophont density
# Fraction of daughter cells that are theronts is constant and 
  # independent of the larval density
# Larval mosquitoes filter-feed on trophonts with a linear functional response
# Larval development does not depend on trophont consumption because most
  # larval food comes from other sources (biofilm)
