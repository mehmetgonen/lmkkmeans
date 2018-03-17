function state = kkmeans_train(K, parameters)
    tic;
    opt.disp = 0;
    [H, ~] = eigs(K, parameters.cluster_count, 'la', opt);
    objective = trace(H' * K * H) - trace(K);
    H_normalized = H ./ repmat(sqrt(sum(H.^2, 2)), 1, parameters.cluster_count);

    stream = RandStream.getGlobalStream;
    reset(stream);
    state.clustering = kmeans(H_normalized, parameters.cluster_count, 'MaxIter', 1000, 'Replicates', 10);
    state.objective = objective;
    state.parameters = parameters;
    state.time = toc;
end
