function res = CP2TRcheckall(V, n, maxn, prec)
%%% Transforms a canonical decomposition into TR-format with first index
%%% resulting in lowest storage cost.
%%% V is a cell array of dimension d \times 1 containing matrices of dimension nk x r.

[d, ~] = size(V);
best_cost = Inf;
best_k = -1;
best_tr = Inf;
done = 0;
for k = 1:n:maxn
    disp(k)
    W = permuteTR(V, k);
    V1 = W{1};
    [n1, r] = size(V1);

    divs = divisors(r);

    for l = 1:min(length(divs)-1, 500) %no need to check the last one (perm.)
        div = divs(l);
        tr = roundingTR(CP2TR(W, div), prec);
        curr_cost = storage_size(tr);
        if curr_cost < best_cost
            best_cost = curr_cost;
            best_k = k;
            best_tr = tr;
            if all_ones(tr) == 1    %no need to check the rest
                disp('all ones')
                done = 1;
                break
            end            
        end
    end
    if done == 1
        break
    end
end
res = ipermuteTR(best_tr, best_k);
end

function yes = all_ones(tr)
%returns one if r_k = 1 for all k; 0 otherwise.
[d, ~] = size(tr);
yes = 1;
for k = 1:d
    [r0, nk, r1] = size(tr{k});
    if r0 ~= 1
        yes = 0;
        break
    end
end
end