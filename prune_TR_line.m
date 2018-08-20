function res = prune_TR_line(tr)
%prunes the output of ex_find_best_graph back to TR representation

[d, ~] = size(tr);

res = cell(d,1);
res{1} = tr{1};

curr_ind = 2;
for k = 3:2:d-2
    s = size(tr{k});
    res{curr_ind} = reshape(tr{k}, [s(1), s(3), s(4)]);
    curr_ind = curr_ind + 1;
end

s = size(tr{d});
res{curr_ind} = reshape(tr{d}, [s(1), s(3), 1]);
curr_ind = curr_ind + 1;

for k = 2:2:d-1
    res{curr_ind} = tr{k};
    curr_ind = curr_ind + 1;
end

end