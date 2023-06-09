function res = kgaps(data, u, k, inc_cens, nv)
    %' Maximum likelihood estimation for the \eqn{K}-gaps model
    %'
    %' Calculates maximum likelihood estimates of the extremal index \eqn{\theta}
    %' based on the \eqn{K}-gaps model for threshold inter-exceedances times of
    %' Suveges and Davison (2010).

    arguments
        data (:,:) double
        u (1,1) double {mustBeNumeric}
        k (1,1) double {mustBePositive,mustBeNumeric} = 1
        inc_cens (1,1) logical = true
        nv.disp char {mustBeMember(nv.disp,{'y','n'})} = 'y'
    end

    if u >= max(data, [], 'omitnan')
        error("'u' must be less than 'the max of 'data'")
    end
    % If there are missing values then use split_by_NAs to extract sequences
    % of non-missing values
    if any(isnan(data))
        data = exdex.int.split_by_nans(data);
    end
    % Estimate the marginal exceedance probability q_u
    q_u = mean(data(~isnan(data)) > u);
    % Calculate sufficient statistics for each column in data and then sum
    stats_list =  arrayfun(@(i) exdex.kgaps.stats(data(:,i), u, q_u, k, inc_cens), ...
        1:size(data,2), 'un', 0);
    ss = sum(vertcat(stats_list{:}), 1);
    % If N0 = 0 then all exceedances occur singly (all K-gaps are positive)
    % and the likelihood is maximized at theta = 1.
    N0 = ss.N0;
    % If N1 = 0 then we are in the degenerate case where there is one cluster
    % (all K-gaps are zero) and the likelihood is maximized at theta = 0.
    N1 = ss.N1;
    if N1 == 0
        theta_mle = 0;
    elseif N0 == 0 
        theta_mle = 1;
    else
        sum_qs = ss.sum_qs;
        theta_mle = exdex.kgaps.quad_solve(N0, N1, sum_qs);
    end
    % Estimate standard error
    % For completeness add an estimate based on the expected information
    % If N1 = 0 then the estimate of theta is 0 and we return NA for se_exp
    % If N0 = 0 then the estimate of theta is 1 and the expected information is
    % not defined and we return NA for se_exp
    if N1 > 0 && N0 > 0
        exp_info = exdex.kgaps.exp_info(theta_mle, ss, inc_cens);
    else
        exp_info = nan;
    end
    se_exp = 1 / sqrt(exp_info);
    % Based on the observed information
    obs_info = 0;
    if N0 > 0
        obs_info = obs_info + N0 / (1 - theta_mle) ^ 2;
    end
    if N1 > 0
        obs_info = obs_info + 2 * N1 / theta_mle ^ 2;
    end
    theta_se = sqrt(1 / obs_info);
    % If K = 0 then the estimate of theta is 1 by default
    % We return a SE equal to 0 so that the estimation of SEs for return level
    % estimates in the lite package work in this case
    if k == 0
        theta_se = 0;
    end
    max_loglik = exdex.kgaps.loglik(theta_mle, N0, N1, ss.sum_qs, ss.n_kgaps);
    sum_qs = ss.sum_qs;
    n_kgaps = ss.n_kgaps;
    res = table(theta_mle, theta_se, se_exp, N0, N1, sum_qs, n_kgaps, k, u, inc_cens, max_loglik);
    if strcmpi(nv.disp, 'y')
        disp(res)
    end
end

% fini