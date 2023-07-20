function [res] = finite_diffs(in)
    res = zeros(1, length(in)); 
    for i = 2 : length(in)
        res(1, i) = in(1, i) - in(1, i-1); 
    end
end