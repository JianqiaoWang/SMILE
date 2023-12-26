scoring_cg <- function(X, y, W = NULL, theta0 =c(0.5, 0.5),
                           max_iter = 100, eps = 1e-6,
                           verbose = FALSE){

  theta <- theta0
  converged <- FALSE
  iter <- 0
  n<- length(y)

  Trace_XWX =  if(is.null(W)){
    sum(X^2)
  }else{
    sum( W * (t(X) %*% X)  )
  }

  while (!converged && iter < max_iter){

    theta2 = theta[1]
    sigma2 = theta[2]
    gamma = theta2/sigma2

    # by = (I + gamma S)^{-1}y
    by= S_eigen$vectors %*% ( 1/(1 + gamma* S_eigen$values ) * ( t(S_eigen$vectors) %*% y) )
    # t1 = trace((I + gamma S)^{-1})
    trace1= sum( 1/(1 + gamma* S_eigen$values ) )
    trace2=sum( 1/(1 + gamma* S_eigen$values )^2 )

    # calculate the score function and the expected gradient
    g <- gradient_var(y, by, theta2, sigma2, trace1, trace2)
    I <- fisher_mat_var(y, by, theta2, sigma2, trace1, trace2)

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


