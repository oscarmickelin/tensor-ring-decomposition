function tr = CP2TR(V, r0)
%%% Transforms a canonical decomposition into TR-format with first index
%%% r0. V is a cell array of dimension d \times 1 containing matrices of dimension nk x r. r0 must
%%% be a divisor of r.

[d, ~] = size(V);

tr = cell(d,1);

V1 = V{1};
[n1, r] = size(V1);
tr{1} = permute(reshape(reshape(V1', [n1*r, 1])', [1, r, n1]), [1,3,2]);

for k = 2:(d-1)
    Vk = V{k};
    [nk, r] = size(Vk);
    g = zeros(r, nk, r);
    tmp = reshape(Vk', [nk*r, 1]);
    for l = 1:nk
        g(:,l,:) = diag(tmp( (l-1)*r + 1 : l*r));
    end
    tr{k} = g;
end

Vd = V{d};
[nd, r] = size(Vd);

tr{d} = reshape(reshape(Vd', [nd*r, 1]), [r, nd, 1]);
tr = TT2TR(tr, r0);