function out = stats(data, u, q_u, D, inc_cens)
    %' Sufficient statistics for the left-censored inter-exceedances time model
    %'
    %' Calculates sufficient statistics for the the left-censored inter-exceedances
    %' time \eqn{D}-gaps model for the extremal index \eqn{\theta}.
    arguments
        data (:,:) double
        u (1,1) double {mustBeNumeric}
        q_u double
        D (1,1) double {mustBePositive,mustBeNumeric} = 1
        inc_cens (1,1) logical = true;
    end

    if isempty(q_u)
        q_u = mean(data > u, [], 'omitnan');
    end
    data = rmmissing(data);  
    % If all the data are smaller than the threshold then return null results
    if u >= max(data, [], 'omitnan')
        out = table(0, 0, 0, 0,'VariableNames',[ "N0", "N1", "sum_qtd", "n_dgaps"]);
        return
    end
    % Sample size, positions, number and proportion of exceedances
    nx = numel(data);
    exc_u = find(data > u);
    N_u = numel(exc_u);
    % Inter-exceedances times and left-censoring indicator
    T_u = diff(exc_u);
    left_censored = T_u <= D;
    % N0, N1, sum of scaled inter-exceedance times that are greater than d,
    % that is, not left-censored
    N1 = sum(~left_censored);
    N0 = N_u - 1 - N1;
    T_gt_D = T_u(~left_censored);
    sum_qtd = sum(q_u * T_gt_D);
    % Store the number of D-gaps, for use by nobs.dgaps()
    n_dgaps = N0 + N1;
    % Include censored inter-exceedance times?
    if true(inc_cens)
        % censored inter-exceedance times and K-gaps
        T_u_cens = [exc_u(1) - 1, nx - exc_u(N_u)];
        % T_u_cens <= d adds no information, because we have no idea to which part
        % of the log-likelihood they would contribute
        left_censored_cens = T_u_cens <= D;
        % N0, N1, sum of scaled inter-exceedance times that are greater than D,
        % that is, not left-censored
        N1_cens = sum(~left_censored_cens);
        n_dgaps = n_dgaps + N1_cens;
        T_gt_D_cens = T_u_cens(~left_censored_cens);
        sum_qtd_cens = sum(q_u * T_gt_D_cens);
        % Add contributions.
        % Note: we divide N1_cens by two because a right-censored inter-exceedance
        % times that is not left-censored at d (i.e. is greater than d) contributes
        % theta exp(-theta q_u T_u) to the D-gaps likelihood, but an uncensored
        % observation contributes theta^2 exp(-theta q_u T_u).
        N1 = N1 + N1_cens / 2;
        sum_qtd = sum_qtd + sum_qtd_cens;
    end
    out = table(N0, N1, sum_qtd, n_dgaps);
end

% fini    