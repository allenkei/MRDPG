source("utility.R")


# f1: Directed Dirichlet, T = 150 (3 change points)
# f2: Flip layers and change block numbers, T = 200 (5 change points) 
# f3: Block size change, first layer only, change to T = 200 (3 change points 1-2-3-1)
# f4: Sparsity fluctuates, T = 200 (5 change points 1-2-3-2-1-2)
# f5: Sparsity, weak difference, T = 200 (3 change points 1-2-3-1)


generate <- function(scenario, num_node = 50, num_seq = 1, save = FALSE) {
  
  num_time <- 200
  num_layer <- 4
  
  if(scenario == "f1"){
    
    cp_truth <- c(70, 140)
    d <- 5
    
    A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
    
    # begin simulate data
    for(seq_iter in 1:num_seq){
      params <- get_dirichlet_params(num_node, num_node, num_layer, d)
      X <- params[[1]]; Y <- params[[2]] # FIXED latent position
      W_1 <- params[[3]]; W_2 <- params[[4]]; # DIFFERENT weights
      
      A.all_seq[seq_iter,,,,] <- get_data_cp_dirichlet(num_time, cp_truth, num_node, num_node, num_layer, 
                                                       X, Y, W_1, W_2, directed = TRUE)
    }
    
  } else if(scenario == "f2"){
    
    # BLOCK NUMBER K CHANGED
    num_block_before_K <- 4 # different
    num_block_after_K <- 3 # different
    FL_K = FALSE # same
    
    # LAYER L FLIPPED
    num_block_before_L <- 4 # same
    num_block_after_L <- 4 # same 
    FL_L = TRUE # different
    
    # BLOCK
    sbm_params_K <- get_sbm_params(n=num_node, L=num_layer, n_c=c(num_block_before_K, num_block_after_K), FL_K)
    probability_K_1 = sbm_params_K[[1]]
    probability_K_2 = sbm_params_K[[2]]
    
    # LAYER
    sbm_params_L <- get_sbm_params(n=num_node, L=num_layer, n_c=c(num_block_before_L, num_block_after_L), FL_L)
    probability_L_1 = sbm_params_L[[1]]
    probability_L_2 = sbm_params_L[[2]]
    
    A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
    
    # begin simulate data
    for(seq_iter in 1:num_seq){
      
      A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
      
      for(t_iter in 1:20) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_1)
      for(t_iter in 21:60) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_2)
      for(t_iter in 61:80) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_L_2)
      for(t_iter in 81:160) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_L_1)
      for(t_iter in 161:180) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_2)
      for(t_iter in 181:200) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_1)
      
      A.all_seq[seq_iter,,,,] <- A.tensor
    }; rm(seq_iter, A.tensor)
    
  }else if(scenario == "f3"){

    block_size1 <- floor(c(3, 4, 3) / 10 * num_node) # fixed ratio
    block_size2 <- floor(c(4, 3, 3) / 10 * num_node) # fixed ratio
    block_size3 <- floor(c(5, 3, 2) / 10 * num_node) # fixed ratio
    
    sbm_params1 <- get_sbm_VS_FL_params(n=num_node, L=num_layer, block_size1, block_size2)
    sbm_params2 <- get_sbm_VS_FL_params(n=num_node, L=num_layer, block_size2, block_size3)
    probability_1 = sbm_params1[[1]]
    probability_2 = sbm_params1[[2]]
    probability_3 = sbm_params2[[2]] 
    
    A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
    
    # begin simulate data
    for(seq_iter in 1:num_seq){
      
      A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
      
      for(t_iter in 1:50) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
      for(t_iter in 51:100) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
      for(t_iter in 101:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_3)
      for(t_iter in 151:200) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
      
      A.all_seq[seq_iter,,,,] <- A.tensor
    }; rm(seq_iter, A.tensor)
    
  } else if(scenario == "f4"){

    epsilon <- 0.1
    
    sbm_params <- get_sbm_params_spa_inc(n=num_node, L=num_layer, epsilon = epsilon)
    probability_1 = sbm_params[[1]]
    probability_2 = sbm_params[[2]]
    probability_3 = sbm_params[[3]]
    
    A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
    
    # begin simulate data
    for(seq_iter in 1:num_seq){
      
      A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
      
      for(t_iter in 1:20) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
      for(t_iter in 21:60) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
      for(t_iter in 61:80) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_3)
      for(t_iter in 81:160) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
      for(t_iter in 161:180) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
      for(t_iter in 181:200) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
      
      A.all_seq[seq_iter,,,,] <- A.tensor
    }; rm(seq_iter, A.tensor)
    
  } else if(scenario == "f5") {

    epsilon <- 0.05
    
    sbm_params <- get_sbm_params_spa_inc(n=num_node, L=num_layer, epsilon = epsilon)
    probability_1 = sbm_params[[1]]
    probability_2 = sbm_params[[2]]
    probability_3 = sbm_params[[3]]
    
    A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
    
    # begin simulate data
    for(seq_iter in 1:num_seq){
      
      A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
      
      for(t_iter in 1:50) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
      for(t_iter in 51:100) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
      for(t_iter in 101:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_3)
      for(t_iter in 151:200) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
      
      A.all_seq[seq_iter,,,,] <- A.tensor
    }; rm(seq_iter, A.tensor)
  } 
  
  if (save == TRUE) {save(A.all_seq, file = paste0("data/seq",num_seq,"n",num_node,scenario,".RData"))}
  return(A.all_seq)
}




