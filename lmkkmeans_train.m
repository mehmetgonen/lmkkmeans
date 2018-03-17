function state = lmkkmeans_train(Km, parameters)
    tic;
    N = size(Km, 2);
    P = size(Km, 3);
    Theta = ones(N, P) / P;
    K_Theta = calculate_localized_kernel_theta(Km, Theta);

    opt.disp = 0;
    objective = zeros(parameters.iteration_count, 1);
    for iter = 1:parameters.iteration_count
        fprintf(1, 'running iteration %d...\n', iter);
        [H, ~] = eigs(K_Theta, parameters.cluster_count, 'la', opt);
        HHT = H * H';
        
        Q = zeros(N * P, N * P);
        for m = 1:P
            start_index = (m - 1) * N + 1;
            end_index = m * N;
            Q(start_index:end_index, start_index:end_index) = eye(N, N) .* Km(:, :, m) - HHT .* Km(:, :, m);
        end
        res = mskqpopt(Q, zeros(N * P, 1), repmat(eye(N, N), 1, P), ones(N, 1), ones(N, 1), zeros(N * P, 1), ones(N * P, 1), [], 'minimize echo(0)');
        Theta = reshape(res.sol.itr.xx, N, P);
        K_Theta = calculate_localized_kernel_theta(Km, Theta);

        objective(iter) = trace(H' * K_Theta * H) - trace(K_Theta);
    end
    H_normalized = H ./ repmat(sqrt(sum(H.^2, 2)), 1, parameters.cluster_count);

    stream = RandStream.getGlobalStream;
    reset(stream);
    state.clustering = kmeans(H_normalized, parameters.cluster_count, 'MaxIter', 1000, 'Replicates', 10);
    state.objective = objective;
    state.parameters = parameters;
    state.Theta = Theta;
    state.time = toc;
end

function K_Theta = calculate_localized_kernel_theta(K, Theta)
    K_Theta = zeros(size(K(:, :, 1)));
    for m = 1:size(K, 3)
        K_Theta = K_Theta + (Theta(:, m) * Theta(:, m)') .* K(:, :, m);
    end
end
