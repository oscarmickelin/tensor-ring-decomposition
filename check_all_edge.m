function res = check_all_edge(tr, prec)
%%% Inserts exactly one edge between the first index
%%% of the TR-representation and a different index.
%%% The index resulting in the lowest storage cost
%%% is chosen.

[d, ~] = size(tr);

best_store = Inf;
for k = 3:d
    tmp = insert_gen_edge(tr, k, prec);    
    curr_store = storage_size(tmp);
    if curr_store < best_store
        best_store = curr_store;
        res = tmp;
    end
end