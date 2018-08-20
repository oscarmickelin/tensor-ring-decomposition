function res1 = normTR(tr)
%%% computes the norm of tr

[d, ~] = size(tr);
best_r = Inf;
best_k = -1;

for k = 1:d
    [r0, ~, ~] = size(tr{k});
    if r0 < best_r
        best_r = r0;
        best_k = k;
        if r0 == 1
            break
        end
    end
end

res1 = norm_no_perm(permuteTR(tr, best_k));

end
