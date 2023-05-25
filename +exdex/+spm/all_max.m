function out = all_max(x, b, which_dj)
    %' Sliding and disjoint block maxima
    %'
    %' Calculates the (sliding) maxima of all blocks of \code{b} contiguous values
    %' and all sets of the maxima of disjoint blocks of \code{b} contiguous values
    %' in the vector \code{x}.  This provides the first step of computations in
    %' \code{\link{spm}}.

    arguments
        x (:,1) double
        b (1,1) double {mustBePositive,mustBeNumeric}
        which_dj char {mustBeMember(which_dj,{'all','first','last'})}
    end

    % First calculate the sliding block maxima.  All the disjoint maxima that
    % we need are contained in s_max, and we need the sliding maxima anyway
    ys = movmax(x, [0 b - 1]); % get movmax to behave the same as RcppRoll:roll_max and zoo:rollapply
    ys(end + 2 - b:end) = [];
    % The number of maxima of blocks of length b
    n = length(x);
    n_max = floor(n / b);
    % Set the possible first indices
    if strcmpi(which_dj,'all')
        first_value = 1:(n - n_max * b + 1);
    elseif strcmpi(which_dj,'first')
        first_value = 1;
    elseif strcmpi(which_dj,'last')
        first_value = n - n_max * b + 1;
    end
    temp = cell2mat(arrayfun(@(z) get_maxima(z, b, n_max, ys, x), ...
        first_value, 'un', 0));
    i = false(height(temp),1);
    i(1:n_max) = true;
    yd = temp(i, :);
    xd = temp(~i, :);
    out = struct('ys', ys, 'xs', x, 'yd', yd, 'xd', xd);
end

function m = get_maxima(first, b, n_max, ys, x)
    % A function to return block maxima and contributing values starting from
    % the first value first
    % s_ind = seq.int(from = first, by = b, length.out = n_max);
    s_ind = first:b:(b * n_max);
    m = [ys(s_ind) ; x(first:(first + n_max * b - 1))];
end

% fini