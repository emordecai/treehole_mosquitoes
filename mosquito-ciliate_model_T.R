# Model of mosquito - ciliate interactions
# within a season and within a tree hole
# as a function of temperature

### Mosquito stages
# L_u: uninfected larvae (all stages, for now)
# L_i: infected larvae (all stages, for now)
# P_u: uninfected pupae
# P_i: infected pupae

### Ciliate stages
# R: trophonts (free-living non-parasitic ciliates)
# H: theronts (free-living parasitic ciliates)

### Parameters or functions
# mu_L, mu_R, mu_H: temperature-dependent background mortality rates for each life stage
# alpha_L: additional parasite-induced mortality of infected larvae
# K_R: carrying capacity of trophonts within a tree hole
# K_L: inverse strength of density-dependence for mosquito larvae within a tree hole
# d_L: temperature-dependent larval development rate
# beta: per-contact transmission probability for larvae that consume theronts
# delta: temperature-dependent cell division rates for free-living ciliates
# n: fraction of daughter cells that are theronts
# m: number of trophonts produced per killed infected larva
# f_LR: linear consumption rate of larval mosquitoes on trophonts, temperature-dependent


# a model that tracks larval and pupal dynamics
ode_larvae_T = function(t, state, parameters){
  with(as.list(c(state, fixed.pars, parameters)), {
    dLu = -dd_mu(mu_L(temp[t]), K_L, Lu, Li) - beta*H*Lu - d_L(temp[t])*Lu 
    dLi = -dd_mu(mu_L(temp[t]), K_L, Li, Lu) - alpha_L*Li + beta*H*Lu - d_L(temp[t])*Li
    dR = dd_div(delta(temp[t]), K_R, R) - n*R + (alpha_L + mu_L(temp[t]))*Li*m - mu_R(temp[t])*R - (Lu + Li)*f_LR(temp[t])*R
    dH = n*R - mu_H(temp[t])*H - beta*H*(Lu+Li)
    dPu = d_L(temp[t])*Lu
    dPi = d_L(temp[t])*Li
    list(c(dLu, dLi, dPu, dPi, dR, dH))
  })
}


### Key assumptions embedded in these equations:
# Transmission only occurs in larval stage
# Larvae, pupae, and adults do not clear infection
# Trophont cell division is logistic
# Infection is intensity-independent
# Parasite-killed larvae produce trophonts
# Larvae are an external (between-season) input
# Infected adults are castrated
# Larval mortality rate is density-dependent
# Carrying capacity and consumption rate are independent of trophont density
# Fraction of daughter cells that are theronts is constant and 
  # independent of the larval density
# Larval mosquitoes filter-feed on trophonts with a linear functional response
# Larval development does not depend on trophont consumption because most
  # larval food comes from other sources (biofilm)
