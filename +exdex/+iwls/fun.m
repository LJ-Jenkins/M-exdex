function out = iwls_fun(n_wls, N, S_1_sort, exp_qs, ws, nx) 
    %
    % This function implements the algorithm on page 46 of Suveges (2007).
    % [In step (1) there is a typo in the paper: in x_i the N_C+1 should be N.]
    %
    % Returns: A list with components
    %    theta : A numeric scalar.  The new estimate of theta.
    %    n_wls : A numeric scalar.  The new value of n_wls.
    %
    % Extract the values corresponding to the largest n_wls 1-gaps
    % Extract the largest n_wls scaled 1-gaps (ordered largest to smallest)
    chi_i = S_1_sort(1:n_wls);
    % Standard exponential quantiles, based on N 1-gaps (largest to smallest)
    x_i = exp_qs(1:n_wls);
    % Extract the weights for the values in chi_i
    ws = ws(1:n_wls);
    % Weighted least squares for (chi_i, x_i) (for 'fitlm' this order is reversed)
    temp = fitlm(x_i, chi_i, 'Weights', ws);
    ab = temp.Coefficients.Estimate;
    % Estimate theta
    theta = min(exp(ab(1) / ab(2)), 1);
    % Update n_wls
    n_wls = floor(theta * (N - 1));
    out = table(theta, n_wls);
end

% fini