



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