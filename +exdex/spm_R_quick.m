function res = spm_R_quick(data, b, bias_adjust, constrain, varN, which_dj, nv)
    %' Semiparametric maxima estimator of the extremal index
    %'
    %' Estimates the extremal index \eqn{\theta} using a semiparametric block
    %' maxima estimator of Northrop (2015) (\code{"N2015"}) and a variant of this
    %' estimator studied by Berghaus and Bucher (2018) (\code{"BB2018"}), using
    %' both sliding (overlapping) block maxima and disjoint (non-overlapping) block
    %' maxima.  A simple modification (subtraction of \eqn{1 / b}, where \eqn{b} is
    %' the block size) of the Berghaus and Bucher (2018) estimator
    %' (\code{"BB2018b"}) is also calculated.  Estimates of uncertainty are
    %' calculated using the asymptotic theory developed by Berghaus and Bucher
    %' (2018).

    arguments
        data (:,1) double
        b (1,1) double {mustBeNumeric,mustBePositive}
        bias_adjust char {mustBeMember(bias_adjust,{'BB3', 'BB1', 'N', 'none'})} = 'N'
        constrain (1,1) logical = true
        varN (1,1) logical = true
        which_dj char {mustBeMember(which_dj,{'last','first'})} = 'first'
        nv.disp char {mustBeMember(nv.disp,{'y','n'})} = 'y'
    end

    if any(isnan(data))
        error("Data contains nan values")
    elseif any(isinf(data))
        error("Data contains Inf values")
    end

    % Find the number of (disjoint) blocks
    k_n = floor(length(data) / b);
    if k_n < 1
        error("b is too large: it is larger than length(data)")
    end
    %
    % Estimate sigma2_dj based on Section 4 of Berghaus and Bucher (2018)
    % We require the disjoint maxima to do this.  If sliding = TRUE then
    % pass these to spm_sigmahat_dj using the dj_maxima argument
    % Only do this is b_ok = TRUE.
    % Otherwise, just calculate point estimates of theta
    % At this point these estimates have not been bias-adjusted, unless
    % bias_adjust = "N".
    %
    % Find all sets of maxima of disjoint blocks of length b
    all_max = exdex.spm.all_max(data, b, which_dj);
    res = exdex.spm.ests_sigmahat_dj(all_max, b, which_dj, bias_adjust);
    % Sliding maxima
    Fhaty = exdex.int.ecdf2(all_max.xs, all_max.ys);
    % Avoid over-writing the `disjoint' sample size k_n: it is needed later
    k_n_sl = length(all_max.ys);
    % Now, m is the number of xs that contribute to the sliding maxima
    m = length(all_max.xs);
    const = -log(m - b + k_n_sl);
    if strcmpi(bias_adjust, 'N')
        Fhaty = (m * Fhaty - b) / (m - b);
    end
    res.theta_sl = [-1 / mean(b * exdex.int.log0const(Fhaty, const)), 1 / (b * mean(1 - Fhaty))];
    % N2015 and BB2018
    %
    % Add the values of the Y-data and the Z-data to the output
    res.data_sl = [-b * log(Fhaty), b * (1 - Fhaty)];
    %
    % Estimate the sampling variances of the estimators
    %
    res.sigma2sl = res.sigma2dj_for_sl - (3 - 4 * log(2)) ./ res.theta_sl .^ 2;
    % res.sigma2sl could contain non-positive values
    % If it does then replace them with NA
    res.sigma2sl(res.sigma2sl <= 0) = nan;
    if true(varN)
        index = 1:2;
    else
        index = [2 2];
    end
    res.se_dj = res.theta_dj .^ 2 .* sqrt(res.sigma2dj(index) / k_n);
    res.se_sl = res.theta_sl .^ 2 .* sqrt(res.sigma2sl(index) / k_n);
    %
    % Perform BB2018 bias-adjustment if required
    %
    if strcmpi(bias_adjust, 'BB3')
        res.bias_dj = res.theta_dj / k_n + res.theta_dj .^ 3 .* res.sigma2dj / k_n;
        res.theta_dj = res.theta_dj - res.bias_dj;
        BB3adj_sl = res.theta_sl / k_n + res.theta_sl .^ 3 .* res.sigma2sl / k_n;
        BB1adj_sl = res.theta_sl / k_n;
        if isnan(res.se_sl)
            res.bias_sl = BB1adj_sl;
        else
            res.bias_sl = BB3adj_sl;
        end
        res.theta_sl = res.theta_sl - res.bias_sl;
        if isnan(res.se_sl(1))
            warning("'bias_adjust' has been changed to ''BB1'' for estimator N2015")
        end
        if isnan(res.se_sl(2))
            warning("'bias_adjust' has been changed to ''BB1'' for estimator BB2018")
        end
    elseif strcmpi(bias_adjust, 'BB1')
        res.bias_dj = res.theta_dj / k_n;
        res.theta_dj = res.theta_dj - res.bias_dj;
        res.bias_sl = res.theta_sl / k_n;
        res.theta_sl = res.theta_sl - res.bias_sl;
    else
        res.bias_dj = [0, 0];
        res.bias_sl = [0, 0];
    end
    %
    % Save the unconstrained estimates, so that they can be returned
    res.uncon_theta_dj = res.theta_dj;
    res.uncon_theta_sl = res.theta_sl;
    %
    % Constrain to (0, 1] if required
    if true(constrain)
        res.theta_dj = min(res.theta_dj, 1);
        res.theta_sl = min(res.theta_sl, 1);
    end
    %
    res.bias_adjust = bias_adjust;
    res.b = b;
    n = ["N2015, sliding", "BB2018, sliding", "N2015, disjoint", "BB2018, disjoint"]';
    res.summary = table(n, [res.theta_sl, res.theta_dj]', [res.se_sl, res.se_dj]', [res.bias_sl, ...
        res.bias_dj]', 'VariableNames', ["Method", "Estimate", "Std. Error", "Bias Adjustment"]);
    if strcmpi(nv.disp, 'y')
        disp(res.summary)
    end
end

% fini