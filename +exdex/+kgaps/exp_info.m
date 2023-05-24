function val = exp_info(theta, ss, inc_cens)
    % How many observations are not right-censored?
    % Note: if inc_cens = TRUE then ss$N1 and ss$n_dgaps are inflated by
    % right-censored observations.  ss$n_dgaps is inflated by 1 for each
    % right-censored observation and ss$N1 by 1/2.
    % Subtract the excess from N0 + N1 to get the total number of observations
    % that are not right-censored
    not_right_censored = ss.N0 + ss.N1 - (ss.n_kgaps - ss.N0 - ss.N1);
    term1 = 1 / (1 - theta);
    term2 = 2 / theta;
    % Add expected contributions from right-censored inter-exceedance times that
    % are not left-censored
    if true(inc_cens)
        term3 = 2 / theta;
    else
        term3 = 0;
    end
    val = not_right_censored * (term1 + term2) + term3;
end

% fini