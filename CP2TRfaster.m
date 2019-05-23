function tr = CP2TRfaster(V, r0, prec)
%%% Transforms a canonical decomposition into TR-format with first index
%%% r0. V is a cell array of dimension d \times 1 containing matrices of dimension nk x r. r0 must
%%% be a divisor of r.
%%% Performs rounding on subsystems to speed up runtime.

[d, ~] = size(V);
if r0 == 1
    tr = roundingTR(CP2TR(V, r0), prec);
    return
end

tmptr = permuteTR(CP2TR(V, r0), 2);
[p0, n1, p1] = size(tmptr{1});

tr = cell(d,1);

for k = 1:r0
    tr_add = cell(d,1);
    tmp = tmptr{1};
    tr_add{1} = tmp(:,:,p1/r0*(k-1)+1:p1/r0*k);
    for l = 2:d-1
        tmp = tmptr{l};
        [a0, nk, a1] = size(tmp);
        tr_add{l} = tmp(1:a0/r0,:, 1:a1/r0);
    end
    tmp = tmptr{d};
    [a0, nk, a1] = size(tmp);
    tr_add{d} = tmp(a0/r0*(k-1)+1:a0/r0*k, :, :);
    if k == 1
        tr = tr_add;
    else
        tr = roundingTR(add(tr, tr_add), prec/r0);
    end
end

tr = ipermuteTR(tr, 2);
end