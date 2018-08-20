function [res, i_ranks] = choose_starting_index(T)
%%% chooses the starting index k minimizing
%%% rank(IM_k)

d = ndims(T);
i_ranks = zeros(1,d);
for k = 1:d
    i_ranks(k) = interaction_rank(T, k);
end
[~, res] = min(i_ranks);


