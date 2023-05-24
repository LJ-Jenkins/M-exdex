function theta_mle = quad_solve(N0, N1, sum_qs)
    aa = sum_qs;
    bb = -(N0 + 2 * N1 + sum_qs);
    cc = 2 * N1;
    qq = -(bb - sqrt(bb ^ 2 - 4 * aa * cc)) / 2;
    theta_mle = cc / qq;
end

% fini