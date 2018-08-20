function tr = TT2TR(tt, r0)
%%% converts a TT-representation into a TR-representation with left end
%%% index r0 by reshaping, stacking and rounding.
%%% r0 must be a divisor of r1 of tt.

[d, ~] = size(tt);
tr = cell(d,1);

%%% stack the first two cores
tmp = tt{1};
[r0old, n1, r1] = size(tmp);
G = zeros(r0*r0old, n1, r1/r0);
for k = 1:n1
    a = reshape(tmp(:,k,:), [r0old, r1]);
    for l = 1:r0
        G(((l-1)*r0old + 1):l*r0old, k, :) = a(:, ((l-1)*r1/r0 + 1):l*r1/r0);
    end
end

tr{1} = G;

tmp = tt{2};
[r1, n2, r2] = size(tmp);
G = zeros(r1/r0, n2, r2*r0);
for k = 1:n2
    a = reshape(tmp(:,k,:), [r1, r2]);
    for l = 1:r0
        G(:, k, ((l-1)*r2 + 1):l*r2 ) = a(((l-1)*r1/r0 + 1):l*r1/r0, :);
    end
end

tr{2} = G;

%%% put remaining cores on diagonal
for k = 3:d
    tmp = tt{k};
    [rb, nk, rn] = size(tmp);
    G = zeros(rb*r0*r0old, nk, rn*r0*r0old);
    for l = 1:nk
        G(:,l, :) = kron(eye(r0*r0old), reshape(tmp(:,l,:), [rb, rn]));
    end
    tr{k} = G;
end
        

end