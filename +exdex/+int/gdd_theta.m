function val = gdd_theta(theta, q_u, D)
  d = q_u * D;
  etd = exp(theta * d);
  val = -(theta * d ^ 2 * etd - 2 * d * etd + 1) / (etd - theta) ^ 2;
end

% fini