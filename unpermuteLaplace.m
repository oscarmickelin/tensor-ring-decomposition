function res = unpermuteLaplace(tr, n)
%%% Converts Laplacian in QTT-format to TT-format

[dn, ~] = size(tr);
d = dn/n;
res = cell(d, 1);

%merge cores
for k = 1:(d)
    tmp = tr{(k-1)*n + 1};
    for l = 1:n-1
        [r0, nk, r1] = size(tmp);
        tmpnext = tr{(k-1)*n+1+l};
        [p0, mk, p1] = size(tmpnext);
        tmp = reshape(tmp, [r0*nk, r1]) * reshape(tmpnext, [p0, mk*p1]);
        tmp = reshape(tmp, [r0, numel(tmp)/(r0*p1), p1]);
    end
    res{k} = tmp;
end

for k = 1:(d)
    tmp = res{k};
    [r0, prodn, r1] = size(tmp);
    tmp = reshape(tmp, [r0, 2*ones(1,2*n), r1]);
    %unshuffle i,j
    dtmp = ndims(tmp);
    tmp = permute(tmp, [1, 2:2:dtmp-1, 3:2:dtmp-1, dtmp]);
    %reshape into i_1i_2 ..., j_1j_2...
    tmp = reshape(tmp, [r0, 2^n*2^n, r1]);
    res{k} = tmp;
end