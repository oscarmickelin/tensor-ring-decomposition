function [res] = fullTR(tr)
%%% Converts a TR-representation into full format.
%%% Algorithm: multiply together the reshaped cores, then make this into
%%% an \alpha_0 , . , \alpha_0 tensor and sum over alpha_0 in a loop. Then
%%% reshape this into an n1 x ... x nd tensor.

a = tr{1};
[d, ~] = size(tr);
[r0, n, r1] = size(a);
ns = zeros(1,d);
ns(1) = n;

for k = 2:d
    [rold, n, rnew] = size(a);
    a = reshape(a, [rold*n, rnew]);
    b = tr{k};
    [rnew, n, rnext] = size(b);
    ns(k) = n;
    b = reshape(b, [rnew, numel(b)/rnew]);
    tmp = a*b;
    a = reshape(tmp, [rold, numel(tmp)/(rold*rnext), rnext]);
end

a = reshape(a, [r0, numel(a)/(r0^2), r0]);

res = zeros(size(numel(a)/r0^2));
for k = 1:r0
    res = res + a(k, :, k);
end

res = reshape(res, ns);

end