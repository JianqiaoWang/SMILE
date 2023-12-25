scoring_vanila <- function(X, y, W = NULL, theta0 =c(1, 0.5),
                                 max_iter = 100, eps = 1e-6,
                                 verbose = FALSE){

  theta <- theta0
  converged <- FALSE
  iter <- 0
  n<- length(y)

  S = ifelse(is.null(W), X %*% t(X) , X %*% W %*% t(X))
  S = (n/trace(S)) * S
  S_eigen = eigen(S)

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


var_cal = function(est, fisher_mat){
  heri.est = est[1]/(est[1] + est[2])
  var.mat = solve(fisher_mat)
  heri.gr = numDeriv::grad(heri, est)
  heri.var = t(heri.gr) %*% var.mat %*% heri.gr
  return(c(heri.est = heri.est, heri.se = sqrt(heri.var), sigma_g2.est = est[1], sigma_g2.se = sqrt(var.mat[1,1]),
           sigma_e2.est = est[2], sigma_e2.se = sqrt(var.mat[2,2]) ))
}

# Evaluate the gradient of the likelihood function
gradient_var = function(y, by, theta2, sigma2, trace1, trace2){

  n = length(y)
  gamma = theta2/sigma2
  gr_theta2 = -1/(2*gamma*sigma2)*(n -  trace1  ) + 1/(2*sigma2*sigma2*gamma)*( sum(y*by) - sum(by^2) )
  gr_sigma2 = -trace1/(2*sigma2) + 1/(2*sigma2^2) * sum(by*by)
  return(c(gr_theta2, gr_sigma2))
}

# Evaluate the fisher informatrix of the likelihood function
fisher_mat_var = function(y, by, theta2, sigma2, trace1, trace2){
  n = length(y)
  fisher = matrix(NA, 2, 2)
  gamma = theta2/sigma2
  fisher[1,1] = 0.5*(trace2 - 2*trace1 + n)/(gamma^2 * sigma2^2)
  fisher[1,2] = 0.5*(trace1 - trace2)/(theta2*sigma2)
  fisher[2,1] = fisher[1,2]
  fisher[2,2] = 0.5*trace2/sigma2^2

  return(fisher)
}
