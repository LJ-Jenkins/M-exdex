function res = dgaps(data, u, D, inc_cens, nv)
    %' Maximum likelihood estimation using left-censored inter-exceedances times
    %'
    %' Calculates maximum likelihood estimates of the extremal index \eqn{\theta}
    %' based on a model for threshold inter-exceedances times of
    %' Holesovsky and Fusek (2020).  We refer to this as the \eqn{D}-gaps model,
    %' because it uses a tuning parameter \eqn{D}, whereas the related \eqn{K}-gaps
    %' model of Suveges and Davison (2010) has a tuning parameter \eqn{K}.

    arguments
        data (:,:) double
        u (1,1) double {mustBeNumeric}
        D (1,1) double {mustBePositive,mustBeNumeric} = 1
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
    stats_list = arrayfun(@(i) exdex.dgaps.stats(data(:,i), u, q_u, D, inc_cens), ...
        1:size(data,2), 'un', 0);
    ss = sum(vertcat(stats_list{:}), 1);
    % Add q_u and D to the list of arguments to be passed to functions that
    % calculated quantities based on the log-likelihood
    ss.q_u = q_u;
    ss.D = D;
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
        f = @(theta_input) -exdex.dgaps.loglik(theta_input, N0, N1, ss.sum_qtd, ss.n_dgaps, ...
            q_u, D);
        theta_mle = fminbnd(f,0,1);
    end
    % Estimate standard error.  In some cases there will be problems with the
    % observed information.  Therefore, also calculate estimated standard errors
    % based on the expected information, using a modified (for inc_cens= TRUE)
    % version of equation (11) on page 202 of Holesovsky and Fusek (2020).
    % If N1 = 0 then the estimate of theta is 0 and we return NA for se_exp
    if N1 > 0
        exp_info = exdex.dgaps.exp_info(theta_mle, ss, inc_cens);
    else
        exp_info = nan;
    end
    % If N1 = 0 then the estimate of theta is 0. The contribution to obs_info
    % from the N0 > 0 case is not constrained to be positive unless D = 0 (when
    % the calculation is the same as K-gaps).  Therefore, we set the SE to NA
    % if N1 = 0 unless D = 0. Note: at least one of N0 and N1 must be positive.
    obs_info = 0;
    if N0 > 0
        if N1 > 0 || D == 0
            obs_info = obs_info - N0 * exdex.int.gdd_theta(theta_mle, q_u, D);
        else 
            obs_info = NA;
        end
    end
    if N1 > 0
        obs_info = obs_info + 2 * N1 / theta_mle ^ 2;
    end
    % The observed information is not guaranteed to be positive
    % If it is not positive then return NA for the estimated SE
    % Similarly for the expected information
    if ~isnan(obs_info) && obs_info <= 0
        theta_se = NA;
        se_exp = NA;
    else
        theta_se = sqrt(1 / obs_info);
        se_exp = 1 / sqrt(exp_info);
    end
    max_loglik = exdex.dgaps.loglik(theta_mle, N0, N1, ss.sum_qtd, ss.n_dgaps, ss.q_u, D);
    sum_qtd = ss.sum_qtd;
    n_dgaps = ss.n_dgaps;
    q_u = ss.q_u;
    res = table(theta_mle, theta_se, se_exp, N0, N1, sum_qtd, n_dgaps, q_u, D, u, ...
        inc_cens, max_loglik);
    if strcmpi(nv.disp, 'y')
        disp(res)
    end
end

% fini