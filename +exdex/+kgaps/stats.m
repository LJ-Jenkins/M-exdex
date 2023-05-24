function out = stats(data, u, q_u, k, inc_cens)
    %' Sufficient statistics for the \eqn{K}-gaps model
    %'
    %' Calculates sufficient statistics for the \eqn{K}-gaps model for the extremal
    %' index \eqn{\theta}. Called by \code{\link{kgaps}}.
    arguments
        data (:,:) double
        u (1,1) double {mustBeNumeric}
        q_u double
        k (1,1) double {mustBePositive,mustBeNumeric} = 1
        inc_cens (1,1) logical = true;
    end

    if isempty(q_u)
        q_u = mean(data > u, [], 'omitnan');
    end
    data = rmmissing(data);  
    % If all the data are smaller than the threshold then return null results
    if u >= max(data, [], 'omitnan')
        out = table(0, 0, 0, 0,'VariableNames',[ "N0", "N1", "sum_qs", "n_kgaps"]);
        return
    end
    % Sample size, positions, number and proportion of exceedances
    nx = numel(data);
    exc_u = find(data > u);
    N_u = numel(exc_u);
    % Inter-exceedances times and K-gaps
    T_u = diff(exc_u);
    S_k = max(T_u - k, 0);
    % N0, N1, sum of scaled K-gaps
    N1 = sum(S_k > 0);
    N0 = N_u - 1 - N1;
    sum_qs = sum(q_u * S_k);
    % Store the number of K-gaps, for use by nobs.kgaps()
    n_kgaps = N0 + N1;
    % Include right-censored inter-exceedance times?
    if true(inc_cens)
        % Right-censored inter-exceedance times and K-gaps
        T_u_cens = [exc_u(1) - 1, nx - exc_u(N_u)];
        S_k_cens = max(T_u_cens - k, 0);
        % N0, N1, sum of scaled K-gaps
        % S_k_cens = 0 adds no information, because P(S >= 0) = 1
        N1_cens = sum(S_k_cens > 0);
        n_kgaps = n_kgaps + N1_cens;
        % Remove the right-censored K-gaps that are equal to zero
        % (This is is not necessary here, but we do it for clarity)
        S_k_cens = S_k_cens(S_k_cens > 0);
        sum_s_cens = sum(q_u * S_k_cens);
        % Add contributions.
        % Note: we divide N1_cens by two because a right-censored non-zero K-gap
        % S_c contributes theta exp(-theta q_u S_c) to the K-gaps likelihood,
        % whereas a non-censored non-zero K-gap contributes
        % theta^2 exp(-theta q_u S_c).
        % See equation (4.3) of Attalides (2015)
        N1 = N1 + N1_cens / 2;
        sum_qs = sum_qs + sum_s_cens;
    end
    out = table(N0, N1, sum_qs, n_kgaps);
    clearvars -except out
end

% fini    