function res = ipermuteTR(tr, start_ind)
%%% shifts the indices of a TR-representation
%%% start_ind-1 places to the right

[d, ~] = size(tr);
res = cell(d,1);
perm = [start_ind:d, 1:(start_ind-1)];
order(perm) = 1:length(perm);
for k = 1:d
    res{k} = tr{order(k)};
end
end