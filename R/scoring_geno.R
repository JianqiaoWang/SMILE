scoring_conjugate_geno <- function(X, y, kernel_fun, theta0 =c(0.5, 0.5), W = NULL, pre.condition= T,
                              max_iter = 100, eps = 1e-3, eps2 = 5e-4, verbose = T,  max_nsim = 50,
                              scale = NULL, center = NULL, ind.col = cols_along(X),
                              ind.row = rows_along(X),ncores = 1){

  # Input check: chage the trace calculation and adaptive accuracy
  stopifnot(ncol(W) == length(ind.col) ,  length(y) == length(ind.row) )


  #initiate the parameter for the scoring method
  theta <- theta0
  theta2 = theta[1]
  sigma2 = theta[2]
  gamma = theta2/sigma2
  converged <- FALSE
  iter <- 0
  n<- length(y)
  by = y/(1 + theta0[1]/theta0[2])

  #initiate mixed random-vectors for monte-carlo trace estimation
  #v_rand_list = lapply(1:10, function(i){
  #  if(i %% 2 == 1) {
  #    2 * rbinom(n, 1, 0.5) - 1  # Generate binomial values for odd i
  #  }else {
  #    rnorm(n)  # Generate normal values for even i
  #  }
  #}
  # )
  #v_rand_list = lapply(1:10, function(i){ 2 * rbinom(n, 1, 0.5) - 1})
  #print("normal")
  v_rand_list = lapply(1:max_nsim, function(i){  2 * rbinom(n, 1, 0.5) - 1 })
  bv_rand_list = v_rand_list


  # start the iteration procedure
  while(!converged && iter < max_iter){

    iter <- iter + 1
    theta2 = theta[1]
    sigma2 = theta[2]
    gamma = theta2/sigma2

    global_eps_ctrl = eps2/iter

    # step1 : calculate by
    by= cg(kernel_fun, b = y, mat = X, x0 = by, gamma = gamma, W = W, eps = 1e-5,
           verbose = verbose, scale = scale, center = center, ind.col = ind.col, ind.row = ind.row, ncores = ncores)


    # step2 : calculate the trace bottleneck and we need more discussion
    trace1 = 0; trace2 = 0; tr_nsim = min( floor(max_nsim/2) + 5*iter , max_nsim)

    for(m in 1:tr_nsim){

      v_rand = v_rand_list[[m]]

      bv0 = bv_rand_list[[m]]

      bv=cg_custom(kernal_FBM, b = v_rand, mat = X, x0 = bv0, gamma = gamma, W = W,
                   eps = global_eps_ctrl, verbose = T, scale = scale, center = center, ind.col = ind.col, ind.row = ind.row,  ncores = ncores)

      bv_rand_list[[m]] = bv #update

      trace1 = trace1 +  sum(v_rand*bv)/tr_nsim

      trace2 = trace2 + sum(bv*bv)/tr_nsim
    }

    # calculate the score function and the expected gradient
    g <- gradient_var(y, by, theta2, sigma2, trace1, trace2)
    I <- fisher_mat_var(y, by, theta2, sigma2, trace1, trace2)
    #var_theta = var_cal(theta, I )


    # update the parameter
    theta_new <- theta + solve(I) %*% g
    converged <- (max(abs(theta_new - theta)) < eps*max(abs(theta), 1))  #max( eps*max(abs(theta), 1), 0.01*(var_theta[4]))
    theta <- theta_new
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

scoring_conjugate_geno_ev <- function(X, y, kernel_fun, theta0 =c(0.5, 0.5), W = NULL, S_ev, pre.condition= T,
                                   max_iter = 100, eps = 1e-3, eps2 = 5e-4, verbose = T,  max_nsim = 50,
                                   scale = NULL, center = NULL, ind.col = cols_along(X),
                                   ind.row = rows_along(X),ncores = 1){

  # Input check: chage the trace calculation and adaptive accuracy
  stopifnot(ncol(W) == length(ind.col) ,  length(y) == length(ind.row) )


  #initiate the parameter for the scoring method
  theta <- theta0
  theta2 = theta[1]
  sigma2 = theta[2]
  gamma = theta2/sigma2
  converged <- FALSE
  iter <- 0
  n<- length(y)
  by = y/(1 + theta0[1]/theta0[2])

  # start the iteration procedure
  while(!converged && iter < max_iter){

    iter <- iter + 1
    theta2 = theta[1]
    sigma2 = theta[2]
    gamma = theta2/sigma2

    # step1 : calculate by
    by= cg(kernel_fun, b = y, mat = X, x0 = by, gamma = gamma, W = W, eps = 1e-5,
           verbose = verbose, scale = scale, center = center, ind.col = ind.col, ind.row = ind.row, ncores = ncores)

    # step2 : calculate the trace bottleneck and we need more discussion
    trace1= sum( 1/(1 + gamma* S_ev ) )
    trace2=sum( 1/(1 + gamma* S_ev )^2 )

    # calculate the score function and the expected gradient
    g <- gradient_var(y, by, theta2, sigma2, trace1, trace2)
    I <- fisher_mat_var(y, by, theta2, sigma2, trace1, trace2)
    #var_theta = var_cal(theta, I )


    # update the parameter
    theta_new <- theta + solve(I) %*% g
    converged <- (max(abs(theta_new - theta)) < eps*max(abs(theta), 1))  #max( eps*max(abs(theta), 1), 0.01*(var_theta[4]))
    theta <- theta_new
    if (T) {
      cat("Iteration", iter, ": theta =", theta, "\n")
      cat("simu trace: trance1", trace1, "trace2", trace2)
    }
  }


  if (iter >= max_iter) {
    warning("Maximum number of iterations reached without convergence")
  }

  return(var_cal(theta, I ))
}


kernal_FBM = function(x, mat, W = NULL, center = NULL, scale = NULL,  gamma = 1, ind.col = cols_along(mat), ind.row = rows_along(mat), ncores = 1 )
{

  if( class(mat) == "FBM.code256" |  class(mat) == "FBM"){
    prodVec_fun = bigstatsr::big_prodVec
    cprodVec_fun = bigstatsr::big_cprodVec

  }

  if(class(mat) == "bed"){
    prodVec_fun = bigsnpr::bed_prodVec
    cprodVec_fun = bigsnpr::bed_cprodVec
  }


  if(is.null(W)){

    gamma * as.numeric(prodVec_fun(mat,  cprodVec_fun(mat, x, ind.row = ind.row, ind.col = ind.col,  center = center, scale = scale, ncores = ncores),
                                   ind.row = ind.row, ind.col = ind.col, center = center, scale = scale, ncores = ncores)) + x

  }else{

    gamma *  as.numeric(prodVec_fun(mat, as.vector(W %*% cprodVec_fun(mat, x,ind.row = ind.row, ind.col = ind.col, center = center,scale = scale, ncores = ncores)),
                                    ind.row = ind.row, ind.col = ind.col, center = center, scale = scale, ncores = ncores )) +  x

  }

}

#
# Trace_estimation_geno = function(){
#
#   # there is no way to avoids first order trace
#   # the second order trace, AI for second order derivative
#
#   my = cg(kernel_fun, b = by, mat = X, x0 = by, gamma = gamma, W = W, eps = 1e-5,
#               verbose = verbose, scale = scale, center = center, ind.col = ind.col, ind.row = ind.row, ncores = ncores)
#
# }


