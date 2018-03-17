kkmeans_train <- function(K, parameters) {
  state <- list()
  state$time <- system.time({
    H <- eigen(K, symmetric = TRUE)$vectors[, 1:parameters$cluster_count]
    objective <- sum(diag(t(H) %*% K %*% H)) - sum(diag(K))
    H_normalized <- H / matrix(sqrt(rowSums(H^2, 2)), nrow(H), parameters$cluster_count, byrow = FALSE)

    set.seed(NULL)
    state$clustering <- kmeans(H_normalized, centers = parameters$cluster_count, iter.max = 1000, nstart = 10)$cluster
    state$objective <- objective
    state$parameters <- parameters
  })    
  return(state)
}
