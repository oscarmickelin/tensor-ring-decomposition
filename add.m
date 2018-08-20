function tr = add(tr1, tr2)
%%% adds two tr decompositions of order-d-tensors

[d, ~] = size(tr1);
tr = cell(d,1);

[r01, n1, r11] = size(tr1{1});
[r02, n1, r12] = size(tr2{1});

g = zeros( max(r01, r02), n1, r11 + r12);
g(1:r01, :, 1:r11) = tr1{1};
g(1:r02, :, (r11+1):(r11+r12)) = tr2{1};
tr{1} = g;

for k = 2:d-1
    [r01, n1, r11] = size(tr1{k});
    [r02, n1, r12] = size(tr2{k});

    g = zeros( r01 + r02, n1, r11 + r12);
    g(1:r01, :, 1:r11) = tr1{k};
    g((r01+1):(r01+r02), :, (r11+1):(r11+r12)) = tr2{k};
    tr{k} = g;
end

[r01, n1, r11] = size(tr1{d});
[r02, n1, r12] = size(tr2{d});

g = zeros( r01 + r02, n1, max(r11,r12));
g(1:r01, :, 1:r11) = tr1{d};
g((r01+1):(r01+r02), :, 1:r12) = tr2{d};
tr{d} = g;