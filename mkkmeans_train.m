% Mehmet Gonen (mehmet.gonen@gmail.com)

function state = mkkmeans_train(Km, parameters)
    tic;
    P = size(Km, 3);
    theta = ones(P, 1) / P;
    K_theta = calculate_kernel_theta(Km, theta.^2);

    opt.disp = 0;
    objective = zeros(parameters.iteration_count, 1);
    for iter = 1:parameters.iteration_count
        fprintf(1, 'running iteration %d...\n', iter);
        [H, ~] = eigs(K_theta, parameters.cluster_count, 'la', opt);

        Q = zeros(P, P);
        for m = 1:P
            Q(m, m) = trace(Km(:, :, m)) - trace(H' * Km(:, :, m) * H);
        end
        res = mskqpopt(Q, zeros(P, 1), ones(1, P), 1, 1, zeros(P, 1), ones(P, 1), [], 'minimize echo(0)');
        theta = res.sol.itr.xx;
        K_theta = calculate_kernel_theta(Km, theta.^2);

        objective(iter) = trace(H' * K_theta * H) - trace(K_theta);
    end
    H_normalized = H ./ repmat(sqrt(sum(H.^2, 2)), 1, parameters.cluster_count);

    stream = RandStream.getGlobalStream;
    reset(stream);
    state.clustering = kmeans(H_normalized, parameters.cluster_count, 'MaxIter', 1000, 'Replicates', 10);
    state.objective = objective;
    state.parameters = parameters;
    state.theta = theta;
    state.time = toc;
end

function K_theta = calculate_kernel_theta(K, theta)
    K_theta = zeros(size(K(:, :, 1)));
    for m = 1:size(K, 3)
        K_theta = K_theta + theta(m) * K(:, :, m);
    end
end