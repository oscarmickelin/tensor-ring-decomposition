function res = permuteTR(tr, k)
%%% shifts the indices of a TR-representation
%%% k-1 places to the left

[d, ~] = size(tr);
res = cell(d,1);

order = [k:d, 1:(k-1)];
for k = 1:d
    res{k} = tr{order(k)};
end
