arma2ma_grad <- function(a, b, h, grad = F)
{
  # do the polynomial division using matrices
  A <- matrix(a, ncol = 1);
  
  # construct a matrix containing "a" coefficients one row at a time
  for (j in 1:h)
  {
    A <- cbind(rbind(A, matrix(0, 1, j)),
               rbind(matrix(0, j, 1), matrix(a, ncol = 1)))
  }
  A <- A[1:(1 + h),]
  
  # if "bb" is a column containing coefficients in "b" and padded with zeros,
  # the division operation can then be cast as a solution to the linear system
  # of equations: A * theta = bb
  theta <- solve(A, matrix(c(b, numeric(1 + h - length(b))), ncol = 1 ))
  
  if (!grad)
  {
    # if a gradient was not requested, just return the solution theta
    return(theta)
  }
  else
  {
    # compute the gradient of theta with respect to the coefficients in "a"
    # and "b", where a0 = 1 is assumed
    p <- length(a) - 1
    q <- length(b)
    dtheta <- matrix(nrow = length(theta), ncol = p)

    L <- diag(1 + h)
    
    # the gradient with respect to "a_j" is also a solution to a system of
    # linear equations, but with "theta" adjusted
    if (p > 0)
    {
      for (i in 1:p)
      {
        dtheta[, i] <- -solve(A, matrix(c(numeric(i),
                                          theta[1:(1 + h - i)]), nrow = 1 + h))
      }
    }
    
    # finally, compute the gradient with respect to "b" as a solution to a
    # system of linear equations
    dtheta <- cbind(dtheta, solve(A, rbind(diag(q),
                                     matrix(0, nrow = 1 + h - q, ncol = q))))
    
    return(list(theta = theta, dtheta = dtheta))
  }
}

ardl_irfs_ci <- function(ardl_est, h = 40, cumirf = T, conf = 0.95)
{
  # obtain the desired percentile of the normal distribution
  z <- qnorm(1 - (1 - conf) / 2)
  
  # construct a lower-triangular matrix of ones to process cumulative IRFs
  L <- matrix(1, nrow = 1 + h, ncol = 1 + h)
  L[upper.tri(L)] <- 0
  
  # get the variance-covariance matrix of the estimated ARDL parameters
  Sig <- vcov(ardl_est)

  # setup indexing to access appropriate ARDL coefficient estimates  
  j <- 1 + ardl_est$order[1]
  a <- c(1, -ardl_est$coefficients[2:j])
  
  # get the number of exogenous variables
  k <- length(ardl_est$order) - 1
  
  # allocate space for the mid-interval IRF as the estimated IRF
  md <- matrix(nrow = 1 + h, ncol = k)
  colnames(md) <- rep("", k)

  # allocate space for the lower-bound of the CI
  lb <- matrix(nrow = 1 + h, ncol = k)
  colnames(lb) <- rep("", k)

  # allocate space for the upper-bound of the CI
  ub <- matrix(nrow = 1 + h, ncol = k)
  colnames(ub) <- rep("", k)
  
  for (i in 1:k)
  {
    j0 <- 1 + j
    j <- j0 + ardl_est$order[1 + i]
    b <- ardl_est$coefficients[j0:j]
    colnames(md)[i] <- names(ardl_est$coefficients[j0])
    colnames(lb)[i] <- names(ardl_est$coefficients[j0])
    colnames(ub)[i] <- names(ardl_est$coefficients[j0])
    
    # get the IRFs and gradients
    irfs <- arma2ma_grad(a, b, h, T)
    
    # compute the variance-covariance of IRFs using the delta method
    idx <- c(1 + c(1:ardl_est$order[1]), j0:j)
    irfs_var <- irfs$dtheta %*% Sig[idx, idx] %*% t(irfs$dtheta)
    
    if (cumirf)
    {
      # compute cumulative IRFs
      irfs_sd <- sqrt(diag(L %*% irfs_var %*% t(L)))
      md[, i] <- cumsum(irfs$theta)
    }
    else
    {
      # compute standard IRFs
      irfs_sd <- sqrt(diag(irfs_var))
      md[, i] <- irfs$theta
    }

    # compute the lower and upper bounds from the estimate and std errors    
    lb[, i] <- md[, i] - z * irfs_sd
    ub[, i] <- md[, i] + z * irfs_sd
  }
  
  return(list(md = md, lb = lb, ub = ub))
}