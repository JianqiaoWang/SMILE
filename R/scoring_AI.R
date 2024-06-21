scoring_conjugate_AI <- function(X, y, kernel_fun, theta0 =c(0.5, 0.5), W = NULL, pre.condition= T,
                              max_iter = 100, eps = 1e-4,  verbose = FALSE,  nsim = 50,
                              scale = NULL, center = NULL, ind.col = cols_along(X),
                              ind.row = rows_along(X), MC.adaptive = T, ncores = 1){

  theta <- theta0
  converged <- FALSE
  iter <- 0
  n<- length(y)
  by = y/(1 + theta0[1]/theta0[2])

  while(!converged && iter < max_iter){

    theta2 = theta[1]
    sigma2 = theta[2]
    gamma = theta2/sigma2

    # by = (I + gamma S)^{-1}y
    by= cg(kernel_fun, b = y, mat = X, x0 = by, gamma = gamma, W = W, eps = 1e-5,
           verbose = verbose, scale = scale, center = center, ind.col = ind.col, ind.row = ind.row, ncores = ncores)

    my= cg(kernel_fun, b = by, mat = X, x0 = by, gamma = gamma, W = W, eps = 1e-5,
           verbose = verbose, scale = scale, center = center, ind.col = ind.col, ind.row = ind.row, ncores = ncores)


    # t1 = trace((I + gamma S)^{-1})
    nsim2 = min(5*(iter + 1), nsim)
    trace_res = trace_cg(nsim2,  X,  kernel_fun, W = W,  eps = 0.005/(iter+1),  verbose = verbose, gamma = gamma, n=n,
                         approx = "MC", scale = scale, center = center, ind.col = cols_along(X),
                         ind.row = rows_along(X), MC.adaptive = MC.adaptive, ncores = ncores)

    trace1= trace_res[1]
    trace2= trace_res[2]


    # calculate the score function and the expected gradient
    g <- gradient_var(y, by, theta2, sigma2, trace1, trace2)
    I <- AI(y, by, my, theta2, sigma2)

    # update the parameter
    theta_new <- theta + solve(I) %*% g
    converged <- max(abs(theta_new - theta)) < eps*max(abs(theta), 1)
    theta <- theta_new
    iter <- iter + 1
    if (verbose) {
      cat("Iteration", iter, ": theta =", theta, "\n")
      cat("simu trace: trance1", trace1, "trace2", trace2)
    }
  }


  if (iter >= max_iter) {
    warning("Maximum number of iterations reached without convergence")
  }

  return(var_cal(theta, I ))
}


scoring_conjugate_geno_AI <- function(X, y, kernel_fun, theta0 =c(0.5, 0.5), W = NULL, pre.condition= T,
                                      max_iter = 100, eps = 1e-4,  verbose = FALSE,  nsim = 50,
                                      scale = NULL, center = NULL, ind.col = cols_along(X),
                                      ind.row = rows_along(X), MC.adaptive = T, ncores = 1){

  # chage the trace calculation and adaptive accuracy
  stopifnot(ncol(W) == length(ind.col) ,  length(y) == length(ind.row) )

  theta <- theta0
  converged <- FALSE
  iter <- 0
  n<- length(y)
  by = y/(1 + theta0[1]/theta0[2])
  #v_rand_list = lapply(1:nsim, function(i){2* rbinom(n, 1, 0.5) - 1  } )
  v_rand_list = lapply(1:nsim, function(i){rnorm(n) } )
  bv_rand_list = v_rand_list

  while(!converged && iter < max_iter){

    theta2 = theta[1]
    sigma2 = theta[2]
    gamma = theta2/sigma2

    #eps_adaptive = max(10^(-(iter+1)), 1e-5)
    eps_adaptive2 = max(0.0001/(iter+1), 1e-6)


    # by = (I + gamma S)^{-1}y
    by= cg(kernel_fun, b = y, mat = X, x0 = by, gamma = gamma, W = W, eps = 1e-5,
           verbose = verbose, scale = scale, center = center, ind.col = ind.col, ind.row = ind.row, ncores = ncores)

    my= cg(kernel_fun, b = by, mat = X, x0 = by, gamma = gamma, W = W, eps = 1e-5,
           verbose = verbose, scale = scale, center = center, ind.col = ind.col, ind.row = ind.row, ncores = ncores)

    # t1 = trace((I + gamma S)^{-1})
    if(pre.condition){
      nsim2 = min(5*(iter + 1), nsim)
    }else{
      nsim2 = nsim
    }
    trace1 = 0;
    trace2 = 0

    for(m in 1:nsim2){

      v_rand = v_rand_list[[m]]

      bv0 = bv_rand_list[[m]]

      bv=cg(kernal_FBM, b = v_rand, mat = X, x0 = bv0, gamma = gamma, W = W,
                   eps = 1e-5, verbose = F, scale = scale, center = center, ind.col = ind.col, ind.row = ind.row,  ncores = ncores)

      bv_rand_list[[m]] = bv

      trace1 = trace1 +  sum(v_rand*bv)/nsim2
      trace2 = trace2 + sum(bv*bv)/nsim2
    }


    # calculate the score function and the expected gradient
    g <- gradient_var(y, by, theta2, sigma2, trace1, trace2)
    I <- AI(y, by, my, theta2, sigma2)


    # update the parameter
    theta_new <- theta + solve(I) %*% g
    converged <- max(abs(theta_new - theta)) < eps*max(abs(theta), 1)
    theta <- theta_new
    iter <- iter + 1
    if (verbose) {
      cat("Iteration", iter, ": theta =", theta, "\n")
      cat("simu trace: trance1", trace1, "trace2", trace2)
    }
  }


  if (iter >= max_iter) {
    warning("Maximum number of iterations reached without convergence")
  }

  return(var_cal(theta, I ))
}


# Evaluate the fisher informatrix of the likelihood function
AI = function(y, by, my, theta2, sigma2){
  n = length(y)
  fisher = matrix(NA, 2, 2)
  gamma = theta2/sigma2
  fisher[1,1] = 0.5*sum( (y - by)*(by-my) )/(gamma^2 * sigma2^3)
  fisher[1,2] = 0.5*sum(by*by -by*my)/(theta2*sigma2^2)
  fisher[2,1] = fisher[1,2]
  fisher[2,2] = 0.5*sum(by*my)/sigma2^3

  return(fisher)
}
