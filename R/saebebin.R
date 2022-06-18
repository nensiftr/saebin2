#' @title EB Binomial without Covariates

saebebin <- function(y,  n){
  result <- list(theta_i_hat_EB2 = NA, Parameter = list(alpha = NA, beta = NA))

  if (any(is.na(y)))
    stop("Argument y=", response, " contains NA values.")
  if (any(is.na(n)))
    stop("Argument n=", sampel, " contains NA values.")

  m <- length(y) #jumlah area
  nT <- sum(n) #jumlah seluruh sampel
  w <- n/nT #bobot

  #menghitung dugaan langsung proporsi dan variansnya
  theta_i <- y/n #dugaan langsung proporsi
  mse.theta_i <- theta_i * (1-theta_i) / n #varians dugaan langsung
  result$direct$est <- theta_i
  result$direct$mse <- mse.theta_i

  #menghitung proporsi dan ragam proporsi
  theta_ib <- w * theta_i
  theta_hat <- sum(theta_ib) #rataan terboboti
  s_theta2 <- w * (theta_i-theta_hat)^2
  sum_s_theta2 <- sum(s_theta2) #ragam terboboti

  #menduga parameter sebaran beta-binomial alpha dan beta dengan metode momen Kleinman
  n2nT <- (n^2) / nT
  sum_n2nT <- sum(n2nT)
  k11 <- (nT*sum_s_theta2) - theta_hat*(1-theta_hat) * (m-1)
  k12 <- theta_hat*(1-theta_hat) * (nT-sum_n2nT-(m-1))
  k1 <- k11/k12

  #nilai alpha beta
  alpha <- theta_hat*(1-k1)/k1
  beta <- alpha*(1/theta_hat-1)

  #pendugaan bayes dan ragam posterior bagi theta_i
  k21 <- (y+alpha)*(n-y+beta)
  k22 <- (n+alpha+beta+1)*(n+alpha+beta)^2
  theta_i_hat_EB1 <- (y+alpha)/(n+alpha+beta) #penduga bayes
  var_theta_i_hat_EB <- k21/k22 #ragam posterior

  #pendugaan empirical bayes bagi theta_i
  gamma <- n/(n+alpha+beta) #nilai gamma
  theta_i_hat_EB2 = gamma*theta_i+(1-gamma)*theta_hat #penduga empirical bayes
  result$Parameter$alpha <- alpha
  result$Parameter$beta <- beta
  result$theta_i_hat_EB2 <- theta_i_hat_EB2

  return(result)
}
