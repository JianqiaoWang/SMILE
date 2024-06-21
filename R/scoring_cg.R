
scoring_cg_S <- function(S, y, S_ev, theta0 =c(0.5, 0.5),
                         max_iter = 100, eps = 1e-6,
                         verbose = FALSE){

  theta <- theta0
  converged <- FALSE
  iter <- 0
  n<- length(y)
  theta2 = theta[1]
  sigma2 = theta[2]
  gamma = theta2/sigma2
  by = y/(1 + theta0[1]/theta0[2])

  while (!converged && iter < max_iter){

    theta2 = theta[1]
    sigma2 = theta[2]
    gamma = theta2/sigma2

    # by = (I + gamma S)^{-1}y
    #by= cg(Sx, b = y, S = S, gamma = gamma,  x0 = by, eps = 1e-6, verbose = verbose)
    by = solve( diag(n) + gamma * S, y)

    trace1= sum( 1/(1 + gamma* S_ev ) )
    trace2=sum( 1/(1 + gamma* S_ev )^2 )

    # calculate the score function and the expected gradient
    g <- gradient_var(y, by, theta2, sigma2, trace1, trace2)
    I <- fisher_mat_var(y, by, theta2, sigma2, trace1, trace2)

    # update the parameter
    theta_new <- theta + solve(I) %*% g
    converged <- max(abs(theta_new - theta)) < eps*max(abs(theta), 1)
    theta <- theta_new
    iter <- iter + 1
    if (T) {
      cat("Iteration", iter, ": theta =", theta, "\n")
     # cat("simu trace: trance1", trace1, "trace2", trace2)
    }
  }


  if (iter >= max_iter) {
    warning("Maximum number of iterations reached without convergence")
  }

  return(var_cal(theta, I ))
}


scoring_cg_S2 <- function(S, y, theta0 =c(0.5, 0.5), verbose = T,  max_nsim = 50, max_iter = 100, eps = 1e-6, eps2 = 5e-4){

  theta <- theta0
  converged <- FALSE
  iter <- 0
  n<- length(y)
  theta2 = theta[1]
  sigma2 = theta[2]
  gamma = theta2/sigma2
  by = y/(1 + theta0[1]/theta0[2])
  v_rand_list = lapply(1:max_nsim, function(i){  2 * rbinom(n, 1, 0.5) - 1 })
  bv_rand_list = v_rand_list

  while (!converged && iter < max_iter){

    theta2 = theta[1]
    sigma2 = theta[2]
    gamma = theta2/sigma2

    # by = (I + gamma S)^{-1}y
    by= cg(Sx, b = y, S = S, gamma = gamma, eps = 1e-6, x0 = by,  verbose = verbose)


    # step2 : calculate the trace bottleneck and we need more discussion
    trace1 = 0; trace2 = 0; tr_nsim = min( floor(max_nsim/2) + 5*iter , max_nsim)

    for(m in 1:tr_nsim){

      v_rand = v_rand_list[[m]]

      bv0 = bv_rand_list[[m]]

      bv=  cg(Sx, b = v_rand, S = S, gamma = gamma, eps = eps2, x0 = bv0, verbose = verbose)
      #cg_custom(kernal_FBM, b = v_rand, mat = X, x0 = bv0, gamma = gamma, W = W,
      #             eps = global_eps_ctrl, verbose = T, scale = scale, center = center, ind.col = ind.col, ind.row = ind.row,  ncores = ncores)

      bv_rand_list[[m]] = bv #update

      trace1 = trace1 +  sum(v_rand*bv)/tr_nsim

      trace2 = trace2 + sum(bv*bv)/tr_nsim
    }


    # calculate the score function and the expected gradient
    g <- gradient_var(y, by, theta2, sigma2, trace1, trace2)
    I <- fisher_mat_var(y, by, theta2, sigma2, trace1, trace2)

    # update the parameter
    theta_new <- theta + solve(I) %*% g
    converged <- max(abs(theta_new - theta)) < eps*max(abs(theta), 1)
    theta <- theta_new
    iter <- iter + 1
    if(T) {
      cat("Iteration", iter, ": theta =", theta, "\n")
      cat("simu trace: trance1", trace1, "trace2", trace2)
    }
  }


  if (iter >= max_iter) {
    warning("Maximum number of iterations reached without convergence")
  }

  return(var_cal(theta, I ))
}


scoring_cg_AI <- function(S, y, theta0 =c(0.5, 0.5), verbose = T,  max_nsim = 50, max_iter = 100, eps = 1e-6, eps2 = 5e-4){

  theta <- theta0
  converged <- FALSE
  iter <- 0
  n<- length(y)


  #initialization
  by = y/(1 + theta0[1]/theta0[2])
  my = by
  v_rand_list = lapply(1:max_nsim, function(i){  2 * rbinom(n, 1, 0.5) - 1 })
  bv_rand_list = v_rand_list

  while (!converged && iter < max_iter){

    theta2 = theta[1]
    sigma2 = theta[2]
    gamma = theta2/sigma2

    # by = (I + gamma S)^{-1}y
    by= cg(Sx, b = y, S = S, gamma = gamma, eps = 1e-6, x0 = by,  verbose = verbose)
    my= cg(Sx, b = by, S = S, gamma = gamma, eps = 1e-6, x0 = my, verbose = verbose)


    # step2 : calculate the trace bottleneck and we need more discussion
    trace1 = 0; trace2 = 0; tr_nsim = min( floor(max_nsim/2) + 5*iter , max_nsim)

    for(m in 1:tr_nsim){

      v_rand = v_rand_list[[m]]

      bv0 = bv_rand_list[[m]]

      bv=  cg(Sx, b = v_rand, S = S, gamma = gamma, eps = eps2, x0 = bv0, verbose = verbose)

      bv_rand_list[[m]] = bv #update

      trace1 = trace1 +  sum(v_rand*bv)/tr_nsim

      trace2 = trace2 + sum(bv*bv)/tr_nsim
    }


    # calculate the score function and the expected gradient
    g <- gradient_var(y, by, theta2, sigma2, trace1, trace2)
    I <- AI(y, by, my, theta2, sigma2)

    # update the parameter
    theta_new <- theta + solve(I) %*% g
    converged <- max(abs(theta_new - theta)) < eps*max(abs(theta), 1)
    theta <- theta_new
    iter <- iter + 1
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


Sx <- function(x, gamma, S) {
  gamma * S %*% x + x
}



cg = function(Ax, b, x0 = rep(0, length(b)), eps = 1e-6,
              verbose = TRUE, ...)
{
  m = length(b)
  x = x0
  r = b - Ax(x0, ...)
  p = r
  r2 = sum(r^2)
  for(i in 1:m)
  {
    Ap = Ax(p, ...)
    alpha = r2 / sum(p * Ap)
    x = x + alpha * p
    r = r - alpha * Ap
    r2_new = sum(r^2)
    err = sqrt(r2_new)

    if(verbose)
      cat(sprintf("Iteration %d, err = %.8f\n", i, err))

    if(err < eps)
      break;
    theta = r2_new / r2
    p = r + theta * p
    r2 = r2_new
  }
  x
}

cg_custom = function(Ax, b, x0 = rep(0, length(b)), eps = 1e-6,
                     verbose = TRUE,eval_fun, ...)
{ # customized relative error
  m = length(b)
  x = x0
  r = b - Ax(x0, ...)
  p = r
  r2 = sum(r^2)
  for(i in 1:m)
  {
    Ap = Ax(p, ...)
    alpha = r2 / sum(p * Ap)
    x0 = x
    x = x + alpha * p
    r = r - alpha * Ap


    diff_t1 = abs(sum(b*x) - sum(b*x0))/sum(b*x)
    diff_t2 = abs(sum(x*x) - sum(x0*x0))/sum(x*x)
    r2_new = sum(r^2)
    err = max( c(sqrt(r2_new), sqrt(diff_t1), sqrt(diff_t2) ))


    if(verbose)
      cat(sprintf("Iteration %d, err = %.8f\n", i, err))

    if(err < eps)
      break;
    theta = r2_new / r2
    p = r + theta * p
    r2 = r2_new
  }
  x
}

kernal_base = function(x, mat, W = NULL, center = NULL, scale = NULL,
                       gamma = 1, ind.col = cols_along(mat), ind.row = rows_along(mat), ncores = 1)
{
  # as.numeric(crossprod(mat, mat %*% x)) + lambda * x

  if(is.null(W)){

    gamma * mat %*% crossprod(mat, x) + x

  }else{

    gamma * mat %*% (W %*% crossprod(mat, x)) + x
  }

}



trace_cg = function( nsim2,  X,  kernel_fun, W = NULL, gamma, eps = 0.005/(iter+1), n,  verbose = FALSE,
                     approx = "MC", scale = NULL, center = NULL, ind.col = cols_along(X),
                     ind.row = rows_along(X), MC.adaptive, ncores = 1){

  trace1 = 0; trace2 = 0

  for(m in 1:nsim2){

    set.seed(m+gamma)
    v_rand = 2* rbinom(n, 1, 0.5) - 1
    bv0 = v_rand
    bv=cg_custom(kernel_fun, b = v_rand, mat = X, x0 = bv0, gamma = gamma, W = W,
                 eps = eps, verbose = F, scale = scale, center = center,
                 ind.col = ind.col, ind.row = ind.row,  ncores = ncores)

    trace1 = trace1 +  sum(v_rand*bv)/nsim2
    trace2 = trace2 + sum(bv*bv)/nsim2
  }

  return( c(trace1, trace2) )
}
