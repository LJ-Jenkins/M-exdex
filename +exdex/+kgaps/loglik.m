function loglik = loglik(theta, N0, N1, sum_qs, n_kgaps)
    if theta < 0 || theta > 1
        loglik = -Inf;
        return
    end
    loglik = 0;
    if N1 > 0
        loglik = loglik + 2 * N1 * log(theta) - sum_qs * theta;
    end
    if N0 > 0
        loglik = loglik + N0 * log(1 - theta);
    end
end

% fini