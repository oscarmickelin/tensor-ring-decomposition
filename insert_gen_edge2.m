function res = insert_gen_edge2(tr, k, l, prec)
%%% Inserts an edge between indices k and l for a TR-representation tr.
%%% For the edge rank, the lowest storage result is chosen.

res = ipermuteTR(insert_gen_edge( permuteTR(tr, k), l-k+1, prec), k);
end