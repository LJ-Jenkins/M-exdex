function out = log0const(x, const)
    out = log(x + ~x) + const * ~x;
end

% fini