
# -------- Unbiased p-value calculations

# TODO : cut dT from the parameter list
unbiased_p_val <- function(t_obs, pT, intT, E_T = NA){
  # t_obs : Observed value of statistic T
  # pT : CDF of T
  # dT : PDF of T
  # intT : t * f_T(t), needed to calulcate some integrals
  # E_T : Expected value of T
  
  if(is.na(E_T)){
    # Calculate E_T
    E_T = integrate(f = intT, lower = -Inf, upper = Inf)$value
  }
  
  # Ensuring that t_obs is single-valued
  if(length(t_obs) > 1){
    # TODO : Why this handling ?
    return(sapply(t_obs, find_t_dual_and_calc_tail_prob, pT = pT, intT = intT, E_T = E_T))
  }
  else{
    return(find_t_dual_and_calc_tail_prob(t_obs, pT, intT, E_T))
  }
}

find_t_dual_and_calc_tail_prob <- function(t_obs, pT, intT, E_T){
  # Checking if t_obs is left or right critical value
  if(t_obs < E_T){
    # Finding the 0 of the dual equation, t_dual
    t_dual = uniroot(dual_equation_high, t_obs = t_obs, pT = pT, intT = intT, E_T = E_T, interval = c(E_T, E_T + 1e-2), extendInt = "upX")$root

    # calculating tails
    p_val = 1 + pT(t_obs) - pT(t_dual)
    return(p_val)
  }
  else{
    # Finding the 0 of the dual equation, t_dual
    t_dual = uniroot(dual_equation_low, t_obs = t_obs, pT = pT, intT = intT, E_T = E_T, interval = c(E_T -1e-2, E_T), extendInt = "upX")$root
    # Calculating tails
    p_val = 1 + pT(t_dual) - pT(t_obs)
    return(p_val)
  }
}

dual_equation_high <- function(t_prop, t_obs, pT, dT, intT, E_T){
  # Calculate integral(t * f_T) from t_obs to t_prop
  E_T_short = integrate(f = intT, lower = t_obs, upper = t_prop)$value
  E_T_short - (pT(t_prop) - pT(t_obs)) * E_T # TODO: Can save pT(t_obs) from one iteration to the next
}

dual_equation_low<- function(t_prop, t_obs, pT, dT, intT, E_T){
  # Calculate integral(t * f_T) from t_prop to t_obs
  E_T_short = integrate(f = intT, lower = t_prop, upper = t_obs)$value
  E_T_short - (pT(t_obs) - pT(t_prop)) * E_T # TODO: Can save pT(t_obs) from one iteration to the next
}


# ------- Tail p-value

tail_p_val<-function(t_obs, pT){
  pT_obs = pT(t_obs)
  if(pT_obs < .5){
    return(2 * pT_obs)
  }
  else{
    return(2*(1- pT_obs))
  }
}

# ------- Point p-value

point_p_val<- function(t_obs, dT, pT, mode_T = NA){
  # Finding the mode if not provided
  if(is.na(mode_T)){
    mode_T = optimize(dT, interval = c(-100, 100), maximum = TRUE)$maximum
  }
  # Are we below or above the mode
  dT_obs = dT(t_obs)
  
  dT_match <- function(t) {dT(t) - dT_obs}
  
  if(t_obs < mode_T){
    # calculate the point which has the same density at the other side
    t_dual = uniroot(dT_match, lower = mode_T, upper = mode_T + 1e-2, extendInt = "downX")$root
    # calcualte the probability of being above/below
    return(1 - pT(t_dual) + pT(t_obs))
  }
  else{
    # calculate the point which has the same density at the other side
    t_dual = uniroot(dT_match, lower = mode_T - 1e-2, upper = mode_T, extendInt = "upX")$root
    # calcualte the probability of being above/below
    return(1 - pT(t_obs) + pT(t_dual))
  }
}

# ----- Calculating power

# 1. Find the critical values of t such that the p-values are below the significance level
find_crit_values <- function(method, significance_value,  ..., mid_point){
  
  modified_func <- function(t){method(t, ...) - significance_value}
  # TODO : Update so we get automatic guess for mid-point
  t_1 <- uniroot(modified_func, interval = c(-25, mid_point))$root
  t_2 <- uniroot(modified_func, interval = c(mid_point, 25))$root
  return(c(t_1, t_2))
}

# 2. Calculate the CDF of pT being above or below certain points with alpha_true
power_from_crit_values_gamma <- function(t1, t2, alpha_true){
  return(1 - plgamma(t2, shape = alpha_true) + plgamma(t1, shape = alpha_true))
}
