%initalize the parameters of the algorithm
parameters = struct();

%set the number of clusters
parameters.cluster_count = 2;

%initialize the kernel
K = ??; %should be an N x N matrix containing similarity values between samples

%perform training
state = kkmeans_train(K, parameters);

%display the clustering
display(state.clustering);
