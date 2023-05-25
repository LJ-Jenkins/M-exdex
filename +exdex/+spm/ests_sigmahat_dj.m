function res = ests_sigmahat_dj(all_max, b, which_dj, bias_adjust)
    % Which of the raw values in x are <= each of the values in y?
    % For each of the block maxima in y calculate the numbers of the raw
    % values in each block that are <= the block maximum
    % k_n is the number of blocks
    k_n = size(all_max.yd, 1);
    % m is the number of raw observations
    m = size(all_max.xd, 1);
    % Value to replace log(0), in the unlikely event that this happens
    const = -log(m - b + k_n);
    block = repelem(1:k_n, b)';
    % Use all sets of disjoint maxima to estimate sigmahat2_dj for sliding maxima
    which_vals = 1:size(all_max.yd, 2);
    temp = arrayfun(@(in) UsumN_fn(in), which_vals, 'un', 0);
    temp = horzcat(temp{:});
    Nhat = horzcat(temp{1, :});
    % BB2018
    Zhat = b * (1 - Nhat);
    That = mean(Zhat, 1)';
    Usum = horzcat(temp{3, :});
    Usum =  (k_n * That - (k_n - 1) * (Usum)')';
    Bhat = ((Zhat + Usum)' - 2 * That)';
    % N2015
    ZhatN = -b * log(Nhat);
    ThatN = mean(ZhatN, 1)';
    UsumN = horzcat(temp{2, :});
    UsumN =  (k_n * ThatN - (k_n - 1) * UsumN')';
    % Bhat was mean-centred, but BhatN isn't (quite)
    BhatN = ((ZhatN + UsumN)' - 2 * ThatN)';
    BhatN = ((BhatN)' - mean(BhatN, 1)')';
    % Estimate sigma2_dj.
    % First calculate an estimate for each set of disjoint block maxima
    sigmahat2_dj = sum(Bhat .^ 2, 1) ./ length(Bhat);
    sigmahat2_djN = sum(BhatN .^ 2, 1) ./ length(BhatN);
    % Then, for sliding maxima, take the mean value over all sets of maxima
    sigmahat2_dj_for_sl = sum(sigmahat2_dj) / length(sigmahat2_dj);
    sigmahat2_dj_for_slN = sum(sigmahat2_djN) / length(sigmahat2_djN);
    sigma2dj_for_sl = [sigmahat2_dj_for_slN, sigmahat2_dj_for_sl];
    % For disjoint maxima pick either the first or last value, based on which_dj
    if strcmpi(which_dj, 'first')
        j = 1;
    elseif strcmpi(which_dj, 'last')
        j = length(sigmahat2_dj);
    end
    sigma2dj = [sigmahat2_djN(j), sigmahat2_dj(j)];
    %
    % Point estimates: disjoint maxima. Component i of ThatN (N) and That (BB)
    % contains (the reciprocal of) point estimates of theta based on set of
    % disjoint maxima i.  Use the mean of these estimates as an overall estimate.
    % Perform the Northrop (2015) `bias-adjustment' of Fhaty (Nhat here) if
    % requested and recalculate That and ThatN
    %
    if strcmpi(bias_adjust, 'N') 
        Nhat = (m * Nhat - b) / (m - b);
        That = mean(b * (1 - Nhat), 1)';
        ThatN = mean(-b * exdex.int.log0const(Nhat, const), 1)';
    end
    theta_dj = 1 ./ [ThatN(j), That(j)];
    res.sigma2dj = sigma2dj;
    res.sigma2dj_for_sl = sigma2dj_for_sl;
    res.theta_dj = theta_dj;
    res.data_dj = [-b * log(Nhat(:, j)), b * (1 - Nhat(:, j))]; % N2015 and BB2018

    function out = UsumN_fn(i)
        % Function to calculate Fhaty and UsumN for each set of disjoint block maxima
        y = all_max.yd(:,i);
        x = all_max.xd(:,i);
        % This returns a list with k_n elements.  The ith element of the list
        % (a of length vector k_n) contains the numbers of values in the ith block
        % that are <= each of the block maxima in y
        sum_fun_single = @(x_input) {sum_fun(x_input, y)};
        nums_list = splitapply(sum_fun_single, x, block);
        % Make this into a matrix
        % Column j contains the numbers of values in the ith block that are <= each
        % of the block maxima in y
        % Row i contains the numbers of values in each block that are <=
        % block maximum i in y
        nums_mat = horzcat(nums_list{:});
        % Therefore, the row sums contain the total number of values that are <=
        % each block maximum in y
        % The corresponding proportion is Fhaty in spm(): ecdf of x evaluated at y,
        % in the disjoint (sliding = FALSE) case
        Fhaty = sum(nums_mat, 2) / m;
        % For each block, we want an equivalent vector obtained when we delete that
        % block
        % Column j contains the numbers of values outside of block j that are <= each
        % of the block maxima in y
        % Row i contains the number of values that are outside block 1, ..., k_n
        % and <= block maximum i in y
        % The proportions are Fhat_{-j}(M_{ni}), i, j = 1, ..., k_n
        logm = ~eye(k_n);
        FhatjMni = cell2mat(arrayfun(@(i) sum(nums_mat(:, logm(i, :)), 2), 1:k_n, ...
            'un', 0)) / (m - b);
        % Column j enables us to calculate Yhatni(j) and/or Zhatni(j)
        UsumN = (-b * mean(exdex.int.log0const(FhatjMni, const), 1))';
        Usum = (b * (1 - mean(FhatjMni, 1)))';
        out = {Fhaty ; UsumN ; Usum};
    end

    function out = sum_fun(x, y)
        out = cell2mat(arrayfun(@(y) sum(x <= y), y, 'un', 0));
    end

end

% fini