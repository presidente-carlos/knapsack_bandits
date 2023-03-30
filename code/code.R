#--------------------------------------#
# Title: BwW
# Author: Carlos Gonzalez
# Date: March 2023
#--------------------------------------#

# Contact:
# Please email: carlos.gonzalezperez@economics.ox.ac.uk shall you encounter
# any mistake or have any doubt regarding the code

## Description:
#··············
# This code creates the BwW function which implements the BwW Algorithm in Gonzalez 2023
# Additionally, it generates auxiliary functions for DGP, distance calculation, etc.
# It analyzes and produce relevant graphs for the uniform-linear degenerate case


# 0. Libraries
#-------------
library(dplyr)           # Data manipulation
library(ptinpoly)        # Point \in 3D polytope
library(partitions)      # Compute probability simplex space
library(lpSolve)         # LP solving
library(ggplot2)         # Graphing
library(future)          # Parallel computing
library(future.apply)    # User friendly parallel computing
library(tictoc)          # Time management

# 1. Data preparation (set S and prob space)
#-------------------------------------------

# 1.1. Set S
#-----------

# We first create set S taking as input budget constraints
# The idea is that, given the nature of the problem each vertex is the solution
# of a linear problem for some square subset of the A matrix
# where A includes budget constraints and feasibility constraints (resources >=0)

create_set2 = function(budget, K, d, eps){
  
  A = matrix(c(1, 0, 0,                                # Average rewards
               0, 1, 0,                                # Average wage costs
               0, 0, budget[1] / budget[2],            # Average space consumption
               1, 0, 0,                                # Feasibility constraints
               0, 1, 0,
               0, 0, 1), ncol = d) 
  
  # Resources averages are B = budget[1] constrained
  # d = 3 (as rewards are considered a resource)
  # rewards are unconstrained (i.e. B_r = 1)
  b_vector = c(1, rep(budget[1] / K * (1-eps), d-1), rep(0, d))

  # Solve the iteration of LP for each square A-subset combination
  comb_vector = as.vector(combn(1:nrow(A), 3, simplify = FALSE))
  vertices = list()
  t = 1
    for (j in 1:length(comb_vector)){
      interim_A = A[comb_vector[[j]],]
      if (det(interim_A)!=0){ # If A is not singular and a unique solution exists
        interim_b = b_vector[comb_vector[[j]]]
        vertices[[t]] = solve(interim_A, interim_b) # Solve the LP
      t = t+1
      }
    }
  # Clean and return
  vertices_clean = do.call(rbind, vertices) |> as.matrix()
  colnames(vertices_clean) =  c("x", "y", "z")
  vertices_clean
}

# Load faces data
# It happens to match the indexes of our vertices (cube)
data("faces")


# 1.2. Probability simplex space 
#--------------------------------

# We discretize the probability simplex for time considerations

simplex_partition = function(m, prob_gap){
  
  # Find "all" partitions of P
  n_perm = 1 / prob_gap
  C = t(restrictedparts(m, n_perm, include.zero = T))
  perm_unique = function(vec){
    unique(permutations(vec, k = length(vec)))
  }
  all_P = apply(C, 1, perm_unique)
  do.call(rbind, all_P) / m
}

# 2. Confidence Bound functions
#------------------------------

# Radian Operator
rad = function(nu, N, gamma){
  sqrt(gamma*nu / N) + gamma / N
}

# UCB (nu is a vector)
ucb = function(nu, k_im, gamma){
    min(1, nu + 2*rad(nu, k_im, gamma))
}

# LCB
lcb = function(nu, k_im, gamma){
    max(0, nu - 2*rad(nu, k_im, gamma))
}

# Update intervals for selected arm
interval_function = function(arm_vector, k_im, gamma){
  ucb = sapply(arm_vector, ucb, k_im = k_im, gamma = gamma)
  lcb = sapply(arm_vector, lcb, k_im = k_im, gamma = gamma)
  list("up" = ucb, "low" = lcb)
}

# 3. Optimization functions
#---------------------------

# 3.1. Distance minimization
#---------------------------
# Create a function which measures L_inf distance to set
# Bc we are using L_inf distance we can simply compute the distance from point p
# to every vertex defining S and take the minimum
# Additionally, we use ptinpoly package for determining if a point is inside the set

inf_distance = function(point, vertices){

  # Bc hard bounds, we just care about being inside S
  
  # Compute L_\infty distance for points outside  
  # substract = function(b){point-b}
  # distances = apply(vertices, 1, substract) |> t() |>
  #            abs() |> apply(1, max)
  distance = 100
  inside = pip3d(Vertices = vertices, Faces = faces, Queries = t(point))
  
  if(inside == 0 | inside == 1){
    distance = 0 # distance is zero if inside the Set
  }
  distance
}

# For a given prob_vector and a given vector_v encoding the information of
# V_matrix, it computes the distance
obj_distance = function(vector_v, vertices, prob_vector){
  V = matrix(vector_v, nrow = m)
  inf_distance(point = t(prob_vector%*%V), vertices)
}

# Finds the vector_v (aka the V_matrix) with minimum distance provided a prob_vector
optim_distance = function(upper_bound, lower_bound, prob_vector, vertices){
  
  # Vectorize matrix of bounds
  lower_vector = as.vector(lower_bound)
  upper_vector = as.vector(upper_bound)
 
  # Optim
  min_distance = optim(par = rep(0, d*m), # vector_V
                 fn = obj_distance,
                 method = "L-BFGS-B",
                 vertices = vertices,
                 prob_vector = prob_vector,
                 lower = lower_vector,
                 upper = upper_vector)
  min_distance$value
}

# 3.2. Production function maximization
#---------------------------------------

# Production function of the firm (takes as input the average of r_t)
production_function = function(avg_reward){
  sqrt(avg_reward)
}

# Scaled g-function (and vector transformations)
production_function_nice = function(vector_v, prob_vector){
  V = matrix(vector_v, nrow = m)
  argument_vector = t(prob_vector%*%V)
  v_1 = argument_vector[1]
  v_2 = argument_vector[2]
  -(1/2*(1 + production_function(v_1) - v_2))
}

# Optimization function for the production level given a prob_vector
optim_value = function(upper_bound, lower_bound, prob_vector){
  
  # Vectorize matrix bounds
  lower_vector = as.vector(lower_bound)
  upper_vector = as.vector(upper_bound)
  
  # Optim
  prod_value = optim(par = rep(0, d*m), # vector_V
                     fn = production_function_nice,
                     method = "L-BFGS-B",
                     prob_vector = prob_vector,
                     lower = lower_vector,
                     upper = upper_vector)
  prod_value$value
}

# 4. OPT_LP()
#--------------------

dumb = function(x){ #Update if change g()
  if(x<1/2){
    (x*sqrt(2) - 2*x^2 + 1)*1/2
  }else{
    (sqrt(1/2) - x + 1)*1/2
  }
}

mu = c(1/2, dumb(1/4), dumb(1/2), dumb(3/4), dumb(1))

A = rbind(c(0, 1/8, 1/2, 3/4, 1),
          c(0, 1/2, 1, 1, 1),
          c(1, 1, 1, 1, 1))

# The rest of the LP problem is parameter dependent, so it is completed within BwW

# 5. BwW (Algorithm)
#---------------------

# Variable description:
#-----------------------
# K: Number of periods
# m: Number of arms
# d: Number of elements in \boldsymbol{v} (3 for the standard firm problem)
# delta: (1 - \delta) Probability guarantee of our results
# budget: Vector of non-standardized budget constraints for NON-reward resources
# prob_gap: Sensitivity of the probability simplex


BwW = function(K, delta, m, budget, prob_gap = 0.2, d = 3){
  
  # Arms in rows, resources in columns
  
  # Initialize
  arms = seq(from = 0, to = 1, length.out = m)
  
  # Theoretical parameters
  gamma = log(m*K*d / delta) #sth of this order
  eps = sqrt(gamma * m / budget[1]) + log(K)*gamma*m / budget[1]
  
  # Empirically useful parameters
  gamma = gamma/10
  eps = 0.05

  # Others
  k = rep(0, length(arms)) # Number of times arm i has been played up to i
  v_hat_matrix = matrix(0, nrow = m, ncol = d)
  H_i = list("low" = matrix(0, nrow = m, ncol = d),
             "up" = matrix(1, nrow = m, ncol = d))
  
  # Auxiliary initializations (store)
  selected_mixture = c()
  selected_x = c()
  selected_x_opt_st = c()
  selected_x_opt_dyn = c()
  
  g_value = rep(0, K)
  g_value_opt_st = rep(0, K)
  g_value_opt_dyn = rep(0, K)
  
  rem_wage = c(budget[1], rep(0, K-1))
  rem_space = c(budget[2], rep(0, K-1))
  rem_wage_opt_st = rem_wage
  rem_space_opt_st = rem_space
  rem_wage_opt_dyn = rem_wage
  rem_space_opt_dyn = rem_space
  
  nu_1 = rep(0, K)
  nu_2 = rep(0, K)
  nu_3 = rep(0, K)
  nu_1_dyn = rep(0, K)
  nu_2_dyn = rep(0, K)
  nu_3_dyn = rep(0, K)
  nu_1_st = rep(0, K)
  nu_2_st = rep(0, K)
  nu_3_st = rep(0, K)
  
  policy_index = 1
  dynamic_index = 1
  static_index = 1
  store_policy_index = K
  store_static_index = K
  store_dynamic_index = K
  total_index = K
  
  regret = c()
  regret_st = c()
  regret_dyn = c()
  
  # Generate iid data
  # U follows a Uniform distribution
  # V = 1/2*U
  u = runif(K, 0, 1)
  v = 1/2 * u
  
  # Solve OPT_LP
  b_vector = c(budget[1] / K, budget[2] / K, 1) # below 50 mixed arms solution
  
  opt_lp = lp(direction = "max",
              objective.in = mu,
              const.mat = A,
              const.dir = c("<=", "<=", "="),
              const.rhs = b_vector)
  
  opt_lp_value = opt_lp$objval
  opt_lp_probs = opt_lp$solution
  
  # Construct vertices (with hard bounds eps>0)
  vertices = create_set2(budget = budget, K = K, d = d, eps = eps)
  
  # Simplex partition
  all_P_clean = simplex_partition(m, prob_gap = prob_gap)
  
  for (i in 1:K){
    
    if(rem_space[i]>0 & rem_wage[i]>=0){ # Conditional on resource availability

      # Remove clearly suboptimal p (those who are not in set even with max interval)
      if(i ==1){
        
        # Some context:
        # (Technically) bounds can increase across time
        # So loop cannot run on a succession of subsets
        # However, if there are no V_matrix in [0, 1] then,
        # We can certainly delete those simplexes
        # We consider three different datasets: 
        # i = 0 (all_P_clean) - i = 1 (all_P_updated) - i > 1 (all_P_i)
        # with all_P_i = all_P_updated for i = 1
        
        S_distance = c(rep(0), nrow(all_P_clean)) # Vector for storing distances
        
        for (p in 1:nrow(all_P_clean)){
          
          prob_vector = all_P_clean[p,]
          S_distance[p] = optim_distance(lower_bound = H_i[["low"]], 
                                         upper_bound = H_i[["up"]],
                                         prob_vector = prob_vector,
                                         vertices = vertices)
        }
        
        all_P_updated = all_P_clean[S_distance <=1e-10, ]
        all_P_i = all_P_updated
  
      }else{
      
        S_distance = c(rep(0), nrow(all_P_updated)) # Vector for storing distances
        
        for (p in 1:nrow(all_P_updated)){
          prob_vector = all_P_updated[p,]
          S_distance[p] = optim_distance(lower_bound = H_i[["low"]],
                                         upper_bound = H_i[["up"]],
                                         prob_vector = prob_vector,
                                         vertices = vertices)
        }
      
      # Select as prob distribution candidates only those within the set
      all_P_i = all_P_updated[S_distance <=1e-10, ]
      }
      
      # Initialize store vector for production values
      production = c()
      
      # If data is empty
      if (length(all_P_i)==0 | length(all_P_i)==ncol(all_P_clean)){
        p_row = sample(1:nrow(all_P_clean), size = 1)
        p_i = all_P_clean[p_row,]
        print(paste("Be carefull, no prob vector satisfied the distance req",
                    "You are sampling at random",
                    "Iteration", i))
      }else{
        for (p in 1:nrow(all_P_i)){
          prob_vector = all_P_i[p,]
          production[p] = optim_value(upper_bound = H_i[["up"]], 
                                      lower_bound = H_i[["low"]],
                                      prob_vector = prob_vector)
        }
        # Select all vectors with maximum production
        opt_pool = which(production == min(production))

        # Sample uniformly from this pool
        opt_prod_arm = sample(opt_pool, size = 1)
        p_i = all_P_i[opt_prod_arm,] # Select probability distribution
      }  
  
      # Sample x according to p_i
      x_i = sample(arms, size = 1, prob = p_i)
      
      # Which arm?
      arm_number_i = which(arms == x_i)
      
      # vector nu is realized
      nu_1[i] = 1*(x_i>=v[i])*u[i]
      nu_2[i] = 1*(x_i>=v[i])*x_i
      nu_3[i] = 1*(x_i>=v[i])
      
      nu_i = c(nu_1[i], nu_2[i], nu_3[i])
      
      # Reward function
      g_value[i] = 1/2 * (1 + sqrt(1/i * sum(nu_1)) - 1/i * sum(nu_2))
  
      # Store
      selected_mixture[i] = paste(p_i, collapse = "")
      selected_x[i] = x_i
  
      # Update v_hat (updates empirical averages of resources)
      v_hat_matrix[arm_number_i,] = v_hat_matrix[arm_number_i,] * k[arm_number_i] /
                                      (k[arm_number_i] + 1) + nu_i / (k[arm_number_i] + 1)
      
      # Update bounds
      bounds = interval_function(arm_vector = v_hat_matrix[arm_number_i,], 
                                 k_im = k[arm_number_i] + 1,
                                 gamma = gamma)
      
      H_i[["low"]][arm_number_i,] = bounds[["low"]]
      H_i[["up"]][arm_number_i,] = bounds[["up"]]
      
      # Update k
      k[arm_number_i] = k[arm_number_i] + 1
      
      # Update resources depletion
      rem_wage[i+1] = rem_wage[i] - nu_2[i]
      rem_space[i+1] = rem_space[i] - nu_3[i]
      regret_end = opt_lp_value - g_value[i-1]
    }else{ # If resources are depleted then g_value[i] = 0
      
      if(policy_index == 1){
        store_policy_index = i # Time resources were depleted under BwW
        }
      policy_index = policy_index + 1
      g_value[i] = 0
    }
    
    # Regret wrt to OPL_LP
    regret[i] = opt_lp_value - g_value[i]
    
    
    # Best static arm (cost unaware)
    if (rem_wage_opt_st[i] >= 0 & rem_space_opt_st[i] > 0){
      x_opt_st = 1/4 #sqrt(2)/4
      selected_x_opt_st[i] = x_opt_st
      
      nu_1_st[i] = 1*(x_opt_st>=v[i])*u[i]
      nu_2_st[i] = 1*(x_opt_st>=v[i])*x_opt_st
      nu_3_st[i] = 1*(x_opt_st>=v[i])
      
      g_value_opt_st[i] = 1/2 * (1 + sqrt(1/i * sum(nu_1_st)) -
                                   1/i * sum(nu_2_st))
      
      # Update resources for optimal policy
      rem_wage_opt_st[i+1] = rem_wage_opt_st[i] - nu_2_st[i]
      rem_space_opt_st[i+1] = rem_space_opt_st[i] - nu_3_st[i]
      regret_end_st = opt_lp_value - g_value_opt_st[i-1]
      
    }else{
      if(static_index == 1){
        store_static_index = i
        }
      static_index = static_index + 1
      g_value_opt_st[i] = 0
    }
    
    regret_st[i] = opt_lp_value - g_value[i]
    
    # OPT_LP induced policy
    if (rem_wage_opt_dyn[i] >= 0 & rem_space_opt_dyn[i] > 0){
      
      x_opt_dyn = sample(arms, size = 1, prob = opt_probs)
      selected_x_opt_dyn[i] = x_opt_dyn
      
      nu_1_dyn[i] = 1*(x_opt_dyn>=v[i])*u[i]
      nu_2_dyn[i] = 1*(x_opt_dyn>=v[i])*x_opt_dyn
      nu_3_dyn[i] = 1*(x_opt_dyn>=v[i])
      
      g_value_opt_dyn[i] = 1/2 * (1 + sqrt(1/i * sum(nu_1_dyn)) - 
                                    1/i * sum(nu_2_dyn))
      
      # Update resources for optimal policy
      rem_wage_opt_dyn[i+1] = rem_wage_opt_dyn[i] - nu_2_dyn[i]
      rem_space_opt_dyn[i+1] = rem_space_opt_dyn[i] - nu_3_dyn[i]
      
      regret_end_dyn = opt_lp_value - g_value_opt_dyn[i-1]
      
    }else{
      if(dynamic_index == 1){
        store_dynamic_index = i
      }
      dynamic_index = dynamic_index + 1
      g_value_opt_dyn[i] = 0
    }
    
    regret_dyn[i] = opt_lp_value - g_value_opt_dyn[i]
    
    # Terminate code if resources are depleted in all 3 policies
    if((rem_wage_opt_dyn[i] <= 0 | rem_space_opt_dyn[i] == 0) &
       (rem_wage_opt_st[i] <= 0 | rem_space_opt_st[i] == 0) &
       (rem_wage[i] <= 0 | rem_space[i] == 0)){
      total_index = i
      break
    }
    
  }  
  
  # Return
  list("selected_mixture" = selected_mixture,
       "selected_x" = selected_x,
       "selected_x_opt_st" = selected_x_opt_st,
       "selected_x_opt_dyn" = selected_x_opt_dyn,
       "g_value" = g_value,
       "g_value_opt_st" = g_value_opt_st,
       "g_value_opt_dyn" = g_value_opt_dyn,
       "policy_index" = store_policy_index,
       "static_index" = store_static_index,
       "dynamic_index" = store_dynamic_index,
       "period_end" = total_index,
       "rem_wage" = rem_wage,
       "rem_space" = rem_space,
       "rem_wage_opt_st" = rem_wage_opt_st,
       "rem_space_opt_st" = rem_space_opt_st,
       "rem_wage_opt_dyn" = rem_wage_opt_dyn,
       "rem_space_opt_dyn" = rem_space_opt_dyn,
       "regret" = regret,
       "regret_st" = regret_st,
       "regret_dyn" = regret_dyn,
       "regret_end" = regret_end,
       "regret_end_st" = regret_end_st,
       "regret_end_dyn" = regret_end_dyn)
}

# intermediate_list = BwW(K = 1000, delta = 0.05, m = 5,
#                          budget = c(200, 150), d = 3, prob_gap = 0.2)

replicate_BwW = function(R, K, delta, m, budget, prob_gap = 0.2, d = 3){
  plan(multisession, workers = 6)
  future_replicate(R, BwW(K = K, delta = delta, m = m, 
                          budget = budget, d,
                          prob_gap = prob_gap),
                   future.seed = TRUE)
}

# 6. Data Preparation and Analysis
#-----------------------------------

results_list = replicate_BwW(R = 100, K = 1000, delta = 0.05, m = 5,
                             budget = c(200, 150), d = 3, prob_gap = 0.2)

# Regret and g-value
g_values = tibble(g_value = 
                  c(g_value_policy = 
                      do.call(rbind, results_list["g_value",]) |> colMeans(),
                  g_value_st = 
                    do.call(rbind, results_list["g_value_opt_st",]) |> colMeans(),
                  g_value_dyn = 
                    do.call(rbind, results_list["g_value_opt_dyn",]) |> colMeans()),
                  period = rep(1:K, 3),
                  policy = rep(c("BwW", "Static", "OPT_LP"), each = K))

g_values |> filter(g_value > 0) |> ggplot()+
            geom_smooth(aes(x = period, y = g_value, color = policy), se = FALSE) +
            theme_minimal() +
            ylab("g-value of action x_i induced by policy") +
            xlab("Period i")

ggsave(filename="../plots/regret.jpeg", bg = "white")

# Resource depletion rate

resources = tibble(rem_wage = 
                    c(rem_wage_policy = 
                        do.call(rbind, results_list["rem_wage",]) |> colMeans(),
                      rem_wage_st = 
                        do.call(rbind, results_list["rem_wage_opt_st",]) |> colMeans(),
                      rem_wage_dyn = 
                        do.call(rbind, results_list["rem_wage_opt_dyn",]) |> colMeans()),
                   rem_space = 
                     c(rem_space_policy = 
                         do.call(rbind, results_list["rem_space",]) |> colMeans(),
                       rem_space_st = 
                         do.call(rbind, results_list["rem_space_opt_st",]) |> colMeans(),
                       rem_space_dyn = 
                         do.call(rbind, results_list["rem_space_opt_dyn",]) |> colMeans()),
                  period = rep(1:K, 3),
                  policy = rep(c("BwW", "Static", "OPT_LP"), each = K))

resources |> filter(rem_wage > 0) |> ggplot()+
  geom_smooth(aes(x = period, y = rem_wage, color = policy), se = FALSE) +
  theme_minimal() +
  ylab("Remaining wage at the beginning of period i") +
  xlab("Period i")

ggsave(filename="../plots/rem_wage.jpeg", bg = "white")

resources |> filter(rem_space > 0) |> ggplot()+
  geom_smooth(aes(x = period, y = rem_space, color = policy), se = FALSE) +
  theme_minimal() +
  ylab("Remaining space at the beginning of period i") +
  xlab("Period i")

ggsave(filename="../plots/rem_space.jpeg", bg = "white")