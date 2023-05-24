function loglik = loglik(theta, N0, N1, sum_qtd, n_dgaps, q_u, D)
    loglik = 0;
    if N1 > 0
        loglik = loglik + 2 * N1 * log(theta) - sum_qtd * theta;
    end
    if N0 > 0
        loglik = loglik + N0 * log(1 - theta * exp(-theta * q_u * D));
    end
end

% fini