library(Rmosek)

lmkkmeans_train <- function(Km, parameters) {
  state <- list()
  state$time <- system.time({
    N <- dim(Km)[2]
    P <- dim(Km)[3]
    Theta <- matrix(1 / P, N, P)
    K_Theta <- matrix(0, nrow(Km), ncol(Km))
    for (m in 1:P) {
      K_Theta <- K_Theta + (Theta[,m,drop = FALSE] %*% t(Theta[,m,drop = FALSE])) * Km[,,m]  
    }

    objective <- rep(0, parameters$iteration_count)
    for (iter in 1:parameters$iteration_count) {
      print(sprintf("running iteration %d...", iter))
      H <- eigen(K_Theta, symmetric = TRUE)$vectors[, 1:parameters$cluster_count]
      HHT <- H %*% t(H)

      Q <- matrix(0, N * P, N * P)
      for (m in 1:P) {
        start_index <- (m - 1) * N + 1
        end_index <- m * N
        Q[start_index:end_index, start_index:end_index] <- diag(1, N, N) * Km[,,m] - HHT * Km[,,m]
      }
      problem <- list()
      problem$sense <- "min"
      problem$c <- rep(0, N * P)
      problem$A <- Matrix(rep(diag(1, N, N), P), nrow = N, ncol = N * P, sparse = TRUE)
      problem$bc <- rbind(blc = rep(1, N), buc = rep(1, N))
      problem$bx <- rbind(blx = rep(0, N * P), bux = rep(1, N * P))
      I <- matrix(1:(N * P), N * P, N * P, byrow = FALSE)
      J <- matrix(1:(N * P), N * P, N * P, byrow = TRUE)
      problem$qobj <- list(i = I[lower.tri(I, diag = TRUE)],
                           j = J[lower.tri(J, diag = TRUE)],
                           v = Q[lower.tri(Q, diag = TRUE)])
      opts <- list()
      opts$verbose <- 0
      result <- mosek(problem, opts)
      Theta <- matrix(result$sol$itr$xx, N, P, byrow = FALSE)
      K_Theta <- matrix(0, nrow(Km), ncol(Km))
      for (m in 1:P) {
        K_Theta <- K_Theta + (Theta[,m,drop = FALSE] %*% t(Theta[,m,drop = FALSE])) * Km[,,m]  
      }

      objective[iter] <- sum(diag(t(H) %*% K_Theta %*% H)) - sum(diag(K_Theta))
    }
    H_normalized <- H / matrix(sqrt(rowSums(H^2, 2)), nrow(H), parameters$cluster_count, byrow = FALSE)

    set.seed(NULL)
    state$clustering <- kmeans(H_normalized, centers = parameters$cluster_count, iter.max = 1000, nstart = 10)$cluster
    state$objective <- objective
    state$parameters <- parameters
    state$Theta <- Theta
  })
  return(state)
}
