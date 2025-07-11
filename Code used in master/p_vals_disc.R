
# ---- Unbiased p-value

unbiased_p_val <- function(t_obs, dT, pT, intT, E_T = NA, bndry_min = -Inf, bndry_max = Inf){
  # t_obs : Observed value of statistic T
  # pT : CDF of T
  # intT : t * f_T(t), needed to calulcate some expected values
  # E_T : Expected value of T
  
  if(is.na(E_T)){
    # Calculate E_T
    E_T = sum(sapply(seq(bndry_min, bndry_max), intT))
  }
  
  # Assuming that t_obs is single-valued
  
  return(find_t_dual_and_calc_tail_prob(t_obs, dT = dT, pT = pT, intT = intT, E_T = E_T, bndry_min = bndry_min, bndry_max = bndry_max))
  
}

find_t_dual_and_calc_tail_prob <- function(t_obs, dT, pT, intT, E_T, bndry_min, bndry_max){
  # Checking if t_obs is left or right critical value
  if(t_obs < E_T){
    # Finding the 0 of the dual equation, t_dual
    # TODO : readable if we are at a upcrossing or downcrossing
    dual_vals = bisection_discrete_sign(dual_equation_high, lower = ceiling(E_T), upper = bndry_max, t_obs = t_obs, dT = dT, pT = pT, intT = intT, E_T = E_T)
    t_dual = dual_vals[1]
    gamma_dual = dual_vals[2]
    
    # calculating tails
    paste0("The dual values :",dual_vals)
    # 1 - pT(t_dual) = P(T > t_dual)
    p_val = 1 - pT(t_dual) + pT(t_obs) + dT(t_dual) * gamma_dual
    return(p_val)
  }
  else{
    # Finding the 0 of the dual equation, t_dual
    dual_vals = bisection_discrete_sign(dual_equation_low, lower = bndry_min, upper = floor(E_T), t_obs, dT = dT, pT = pT, intT = intT, E_T = E_T)
    t_dual = dual_vals[1]
    gamma_dual = dual_vals[2]
    
    # Calculating tails
    paste0("The dual values :",dual_vals)
    # 1 - pT(t_obs-1) = P(T >= t_obs)
    # pT(t_dual -1) = P(T < t_dual)
    p_val = 1  - pT(t_obs - 1) + pT(t_dual -1) + dT(t_dual) * gamma_dual
    return(p_val)
  }
}

bisection_discrete_sign <- function(f, lower, upper, ...){
  endCond = FALSE
  
  lower = round(lower)
  upper = round(upper)
  print("yo debugger")
  print(c(lower, upper))
  fl = f(lower, ...)
  fu = f(upper, ...)
  
  if(fl * fu > 0){
    print("Edge case")
    # want to check if we want to add the whole tail
    
    # left side
    fl_minus = f(lower -1, ...)
    if(fl_minus * fl <= 0){
      if(fl >= 0){
        gamma = calc_gamma(fl, lower, ...)
        dual_vals = c(lower, gamma)
        return(dual_vals)
      }
      else{
        gamma = calc_gamma(fl_minus, lower-1, ...)
        dual_vals = c(lower-1, gamma)
        return(dual_vals)
      }
    }
    
    # right side
    fu_plus = f(upper+1, ...)
    if(fu * fu_plus <= 0){
      if(fu >= 0){
        gamma = calc_gamma(fu, upper, ...)
        dual_vals = c(upper, gamma)
        return(dual_vals)
      }
      else{
        gamma = calc_gamma(fu, upper+1, ...)
        dual_vals = c(upper+1, gamma)
        return(dual_vals)
      }
    }
    
    # Other handling
    stop("Error : The endpoints are not of opposite sign")
  }
  while(!endCond){
    mid = round((lower + upper)/ 2)
    fm = f(mid, ...)
    if(fm * fl > 0){
      lower = mid
      fl = fm
    }
    else{
      upper = mid
      fu = fm
    }
    print(c(lower, mid, upper))
    if(upper - lower == 1){
      endCond = TRUE
    }
  }
  if(fl <= 0){
    gamma = calc_gamma(fl, lower, ...)
    dual_vals = c(lower, gamma)
    return(dual_vals)
  }
  else{
    gamma = calc_gamma(fu, upper, ...)
    dual_vals = c(upper, gamma)
    return(dual_vals)
  }
  # add [1] if low value is wanted and [2] if high value is wanted
}

calc_gamma <- function(diff, crit_val, t_obs, dT, pT, intT, E_T, tol = 1e-8){
  dCrit = dT(crit_val)
  
  # diff = (crit_val - E_T)*dCrit*(1-gamma)
  gamma = 1 - abs(diff)/(abs(crit_val - E_T)*dCrit)
  
  print("calc_gamma start")
  print("Critical point prob :")
  print(dCrit)
  print("Critcal value distance :")
  print(crit_val - E_T)
  print("Difference to be made up :")
  print(diff)
  print("Calculated gamma :")
  print(gamma)
  print("calc_gamma end")
  
  # rare occurence handling
  ## the case where dCrit = 0
  if(dCrit == 0){
    gamma = 0
    return(gamma)
  }
  ## the case where the critical value is the mean value
  if(abs(crit_val - E_T) < tol){
    # TODO: Ad hoc handling
    print(diff)
    if(diff == 0){
      gamma = 0
      return(gamma)
    }
    else{
      gamma = 1
      return(gamma)
    }
  }
  
  # normal handling
  if(gamma > 1){
    print("Error incoming")
    print(dCrit)
    print(crit_val - E_T)
    print(diff)
    print(gamma)
    stop("Error : gamma too large")
  }
  else{
    return(gamma)
  }
}

dual_equation_high <- function(t_prop, t_obs, dT, pT, intT, E_T){
  # Calculating E[T*phi] - E[T]*E[phi], whole values C1 and C2, C1 = T_obs
  if(t_prop == t_obs +1){
    # handling for when you have nothign in between (otherwise the seq function will mess with things)
    return(0)
  }
  
  # Calculate the sum(t * P(t)) from t_obs to t_prop
  E_T_short = sum(sapply(seq(t_obs+1, t_prop-1), intT))
  
  # P(x <= T <= y) = P(T <= y) - P(T <= x-1) = pT(y) - pT(x-1)
  E_T_short - (pT(t_prop-1) - pT(t_obs)) * E_T # TODO: Can save pT(t_obs) from one iteration to the next
}

dual_equation_low<- function(t_prop, t_obs, dT, pT, intT, E_T){
  # Calculating E[T*phi] - E[T]*E[phi], whole values C1 and C2, C2 = T_obs
  if(t_prop == t_obs -1){
    # handling for when you have nothign in between (otherwise the seq function will mess with things)
    return(0)
  }
  
  # Calculate the sum(t * P(t)) from t_prop+1 to t_obs (since we include the whole observed val)
  E_T_short = sum(sapply(seq(t_prop+1, t_obs-1), intT))
  
  # P(x <= T <= y) = P(T <= y) - P(T <= x-1) = pT(y) - pT(x-1)
  -(E_T_short - (pT(t_obs-1) - pT(t_prop)) * E_T) # TODO: Can save pT(t_obs) from one iteration to the next
}

# ----- Tail p-value

# Assuming point mass function
tail_p_val <- function(t_obs, dT, pT, bndry_min, bndry_max){
  left_tail = pT(t_obs)
  right_tail = 1 - left_tail + dT(t_obs)
  
  if(left_tail < right_tail){
    right_tail_dual = 0 # Starting proposal
    t_prop =  bndry_max 
    tail_add = 0 # Initializing so it doesnt have to be created in loop
    right_tail_prop = right_tail_dual
    
    # Stop the loop when you add something that makes the right_tail proposal 
    # larger than the left tail
    while(right_tail_prop <= left_tail){
      # accept the suggestion if it is lower than left tail
      right_tail_dual = right_tail_prop 
      
      # Check the next step and stop loop if larger than the left tail
      tail_add = dT(t_prop)
      right_tail_prop = right_tail_prop + tail_add
      t_prop = t_prop - 1
    }
    t_dual = t_prop + 1
    return(left_tail + right_tail_dual)
  }
  else{
    left_tail_dual = 0 # Starting proposal
    t_prop =  bndry_min 
    left_tail_prop = left_tail_dual
    
    # Stop the loop when you add something that makes the left_tail-proposal 
    # larger than the right tail
    while(left_tail_prop <= right_tail){
      # accept the suggestion if it is lower than right tail
      left_tail_dual = left_tail_prop 
      
      # Check the next step and stop loop if larger than the left tail
      tail_add = dT(t_prop)
      left_tail_prop = left_tail_prop + tail_add
      t_prop = t_prop + 1
    }
    t_dual = t_prop - 1
    return(right_tail + left_tail_dual)
  }
}

# ----- Point p-value

# TODO : Effectivice this computation
point_p_val <- function(t_obs, dT, bndry_min, bndry_max){
  
  S = seq(bndry_min, bndry_max)
  p_val = 0
  
  dT_obs = dT(t_obs)
  for(y in S){
    dTy = dT(y)
    if(dTy <= dT_obs){
      p_val = p_val + dTy
    }
  }
  p_val
}