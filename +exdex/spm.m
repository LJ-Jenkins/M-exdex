function res = spm(data, b, bias_adjust, constrain, varN, which_dj)
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
        bias_adjust
        constrain (1,1) logical = true
        varN (1,1) logical = true
        which_dj char {mustBeMember({'last','first'})}
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
    all_max = exdex.spm.all_max_rcpp(data, b, 'all');
    % Here, m is the number of xs that contribute to each set of disjoint block
    % maxima
    m = size(all_max.xd, 1);
    res = cpp_sigma2hat_dj(all_max = all_max, b = b, kn = k_n, m = m,
                            bias_adjust = bias_adjust, which_dj = which_dj);
    est_names = c("N2015", "BB2018");
    % In res theta_dj ... are 2x1 matrices.  Convert them to named vectors.
    res.theta_dj = as.vector(res.theta_dj);
    names(res.theta_dj) = est_names;
    res.sigma2dj = as.vector(res.sigma2dj);
    names(res.sigma2dj) = est_names;
    res.sigma2dj_for_sl = as.vector(res.sigma2dj_for_sl);
    names(res.sigma2dj_for_sl) = est_names;
    colnames(res.data_dj) = est_names;
    % Sliding maxima
    Fhaty = ecdf2(all_max.xs, all_max.ys);
    % Avoid over-writing the `disjoint' sample size k_n: it is needed later
    k_n_sl = length(all_max.ys);
    % Now, m is the number of xs that contribute to the sliding maxima
    m = length(all_max.xs);
    const = -log(m - b + k_n_sl);
    if bias_adjust == "N"
        Fhaty = (m * Fhaty - b) / (m - b);
    end
    res.theta_sl = c(-1 / mean(b * log0const(Fhaty, const)),
                    1 / (b * mean(1 - Fhaty)));
    names(res.theta_sl) = est_names;
    %
    % Add the values of the Y-data and the Z-data to the output
    res.data_sl = cbind(N2015 = -b * log(Fhaty), BB2018 = b * (1 - Fhaty))
    %
    % Estimate the sampling variances of the estimators
    %
    res.sigma2sl = res.sigma2dj_for_sl - (3 - 4 * log(2)) / res.theta_sl ^ 2
    % res.sigma2sl could contain non-positive values
    % If it does then replace them with NA
    res.sigma2sl[res.sigma2sl <= 0] = NA
    indexN = ifelse(varN, 2, 1)
    if (varN) {
    index = 1:2
    } else {
    index = c(2, 2)
    }
    res.se_dj = res.theta_dj ^ 2 * sqrt(res.sigma2dj[index] / k_n)
    res.se_sl = res.theta_sl ^ 2 * sqrt(res.sigma2sl[index] / k_n)
    %
    % Store the raw (not bias-adjusted) estimates
    %
    res.raw_theta_dj = res.theta_dj
    res.raw_theta_sl = res.theta_sl
    %
    % Perform BB2018 bias-adjustment if required
    %
    if (bias_adjust == "BB3") {
    res.bias_dj = res.theta_dj / k_n + res.theta_dj ^ 3 * res.sigma2dj / k_n
    res.theta_dj = res.theta_dj - res.bias_dj
    BB3adj_sl = res.theta_sl / k_n + res.theta_sl ^ 3 * res.sigma2sl / k_n
    BB1adj_sl = res.theta_sl / k_n
    res.bias_sl = ifelse(is.na(res.se_sl), BB1adj_sl, BB3adj_sl)
    res.theta_sl = res.theta_sl - res.bias_sl
    } else if (bias_adjust == "BB1") {
    res.bias_dj = res.theta_dj / k_n
    res.theta_dj = res.theta_dj - res.bias_dj
    res.bias_sl = res.theta_sl / k_n
    res.theta_sl = res.theta_sl - res.bias_sl
    } else {
    res.bias_dj = res.bias_sl = c(N2015 = 0, BB2018 = 0)
    }
    %
    % Add estimates, bias and standard errors for an estimator that results from
    % applying a further bias-adjustment to the BB2018 estimator.  The adjustment
    % is to subtract 1/b.  This stems from the fact that the expected value of
    % the BB2018 variable Z is 1/(theta + 1/b), not 1/theta.
    %
    % Estimates
    res.theta_dj = c(res.theta_dj, res.theta_dj["BB2018"] - 1 / b)
    names(res.theta_dj)[3] = "BB2018b"
    res.theta_sl = c(res.theta_sl, res.theta_sl["BB2018"] - 1 / b)
    names(res.theta_sl)[3] = "BB2018b"
    % Bias adjustment
    if (bias_adjust == "BB3" || bias_adjust == "BB1") {
    res.bias_dj = c(res.bias_dj, res.bias_dj["BB2018"] + 1 / b)
    res.bias_sl = c(res.bias_sl, res.bias_sl["BB2018"] + 1 / b)
    } else {
    res.bias_dj = c(res.bias_dj, 1 / b)
    res.bias_sl = c(res.bias_sl, 1 / b)
    }
    names(res.bias_dj)[3] = "BB2018b"
    names(res.bias_sl)[3] = "BB2018b"
    % Standard errors
    res.se_dj = c(res.se_dj, res.se_dj["BB2018"])
    names(res.se_dj)[3] = "BB2018b"
    res.se_sl = c(res.se_sl, res.se_sl["BB2018"])
    names(res.se_sl)[3] = "BB2018b"
    % Save the unconstrained estimates, so that they can be returned
    res.uncon_theta_dj = res.theta_dj
    res.uncon_theta_sl = res.theta_sl
    % Constrain to (0, 1] if required
    if (constrain) {
    res.theta_dj = pmin(res.theta_dj, 1L)
    res.theta_sl = pmin(res.theta_sl, 1L)
    }
    % Constrain estimates to be non-negative. BB3 bias-adjustment could result in
    % negative estimates as could subtracting 1/b from BB2018 to form BB2018b.
    res.theta_dj = pmax(res.theta_dj, 0L)
    res.theta_sl = pmax(res.theta_sl, 0L)
    %
    res.bias_adjust = bias_adjust
    res.b = b
    res.call = Call
    %
    class(res) = c("spm", "exdex")
    return(res)
}

end






















