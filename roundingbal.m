function res = roundingbal(tr, prec)
%%% Truncates a TR-representation with precision prec.
%%% Balances the norms of the cores.

[d, ~] = size(tr);
Gold = tr{1};
norm_error = prec/sqrt(d-1);
nrm=zeros(d,1);

[r0, ~, ~] = size(Gold);
[~, nlast, ~] = size(tr{d});

for k = 1:d-1
    [rold, n, rnew] = size(Gold);
    Gold = reshape(Gold, [rold*n, rnew]);
    [Q, R] = qr(Gold, 0);
    nrm(k+1) = norm(R, 'fro');
    if (nrm(k+1)~=0)
        R=R./nrm(k+1);
    end
    tr{k} = reshape(Q, [rold, n, numel(Q)/(n*rold)]);
    Gnew = tr{k+1};
    [rnew, n, rnext] = size(Gnew);
    Gnew = reshape(Gnew, [rnew, numel(Gnew)/rnew]);
    tmp = R*Gnew;
    Gnew = reshape(tmp, [numel(tmp)/(n*rnext), n, rnext]);
    Gold = Gnew;
end

for k = d:-1:2
    Gnew = tr{k-1};
    [rnew, n, rold] = size(Gold);
    Gold = reshape(Gold, [rnew, n*rold]);
    [U,S,V]=svd(Gold, 'econ');

    rtrunc = truncation_index(diag(S),norm_error*norm(diag(S)));
    U=U(:,1:rtrunc);
    S=S(1:rtrunc, 1:rtrunc);
    V = V(:, 1:rtrunc);
    [rnext, n, rnew] = size(Gnew);
    Gnew = reshape(Gnew, [rnext*n, rnew])*(U*S);
    Gnew = reshape(Gnew, [rnext, n, rtrunc]);
    tr{k} = reshape(V', [rtrunc, n, rold]);
    tr{k-1} = Gnew;
    Gold = Gnew;
end
tr{1} = Gnew;

nrm(1) = norm(reshape(Gnew, [numel(Gnew), 1]), 'fro');
if (nrm(1) ~= 0)
    tr{1} = tr{1}./nrm(1);
end

nrm0=sum(log(abs(nrm))); 
nrm0=nrm0/d; nrm0=exp(nrm0);
for k=1:d
    tr{k}=nrm0*tr{k};
end
res = tr;
    

