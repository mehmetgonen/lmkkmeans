%initalize the parameters of the algorithm
parameters = struct();

%set the number of clusters
parameters.cluster_count = 2;

%set the number of iterations
parameters.iteration_count = 10;

%initialize the kernels
K = ??; %should be an N x N X P matrix containing similarity values between samples

%perform training
state = lmkkmeans_train(K, parameters);

%display the clustering
display(state.clustering);

%display the kernel weights
display(state.Theta);
