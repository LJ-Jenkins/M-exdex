function newx = split_by_nans(x)
    %' Divides data into parts that contain no missing values
    %'
    %' Splits the values in a numeric matrix column-wise into sequences of
    %' non-missing values.
    arguments
        x(:,:) double
    end

    x(end + 1,:) = nan;
    x = reshape(x, [numel(x) ,1]);
    i = ~isnan(x);
    r = regionprops(i ,x ,'PixelValues');
    c = struct2cell(r);
    n = cellfun(@numel, c); 
    m = max(n);
    padcells = @(k) [c{k} ; nan(m - n(k), 1)];
    c = arrayfun(padcells, 1:numel(c), 'un', 0);
    newx = horzcat(c{:});
end

% fini