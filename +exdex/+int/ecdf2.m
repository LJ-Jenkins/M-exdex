function out = ecdf2(x, y)
    out = cell2mat(arrayfun(@(y) sum(x <= y) / length(x), y, 'un', 0));
end

% fini