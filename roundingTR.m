function res = roundingTR(tr, prec)
%%% Truncates a TR-representation with relative precision prec.

[d, ~] = size(tr);
norm_not_extended = normTR(tr);

Gold = tr{1};
[r0, ~, ~] = size(Gold);
[~, nlast, ~] = size(tr{d});
nrm=zeros(d,1);
if r0 ~= 1
    norm_error = prec*norm_not_extended/sqrt(r0*(d));
else
    norm_error = prec*norm_not_extended/sqrt(r0*(d-1));
end


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


[rold, n, rnew] = size(Gold);
Gold = reshape(Gold, [rold*n, rnew]);
[Q, R] = qr(Gold, 0);

nrm(1) = norm(R, 'fro');
if (nrm(1) ~= 0)
    R = R/nrm(1);
end
norm_error = norm_error/prod(nrm);

tr{d} = reshape(Q, [rold, n, numel(Q)/(n*rold)]);
Gnew = tr{d};

[U,S,V]=svd(R, 'econ');

rtrunc = truncation_index(diag(S),norm_error*norm(diag(S)));

if rtrunc < rnew    %only do this if the rank is reduced
    U=U(:,1:rtrunc);
    S=S(1:rtrunc, 1:rtrunc);
    V = V(:, 1:rtrunc);

    Gnew = reshape(Gnew, [rold*n, numel(Gnew)/(rold*n)])*(U*S);
    Gnew = reshape(Gnew, [rold, n, rtrunc]);
    tr{d} = Gnew;
    Gtmp = tr{1};
    [r0, n1, r1] = size(Gtmp);
    Gtmp = V'*reshape( Gtmp, [r0, n1*r1]);
    tr{1} = reshape(Gtmp, [rtrunc, n1, r1]);

else
    tr{d} = reshape(Q*R, [rold, n, rnew]);
    Gnew = tr{d};
end


Gold = Gnew;
for k = d:-1:2
    Gnew = tr{k-1};
    [rnew, n, rold] = size(Gold);
    Gold = reshape(Gold, [rnew, n*rold]);
    [U,S,V]=svd(Gold, 'econ');

    rtrunc = truncation_index(diag(S),norm_error*norm(diag(S)));
    U=U(:,1:rtrunc);
    S=S(1:rtrunc, 1:rtrunc);

    V = V(:, 1:rtrunc);
    [rnext, nnew, rnew] = size(Gnew);
    Gnew = reshape(Gnew, [rnext*nnew, rnew])*(U*S);
    Gnew = reshape(Gnew, [rnext, nnew, rtrunc]);
    tr{k} = reshape(V', [rtrunc, n, rold]);
    Gold = Gnew;
end
tr{1} = Gnew;
nrm0=sum(log(abs(nrm))); 
nrm0=nrm0/d; nrm0=exp(nrm0);
for k=1:d
    tr{k}=nrm0*tr{k};
end
res = tr;

