#initalize the parameters of the algorithm
parameters <- list()

#set the number of clusters
parameters$cluster_count <- 2

#set the number of iterations
parameters$iteration_count <- 10

#initialize the kernels
K <- ?? #should be an N x N x P matrix containing similarity values between samples

#perform training
state <- mkkmeans_train(K, parameters)

#display the clustering
print(state$clustering)

#display the kernel weights
print(state$theta)
