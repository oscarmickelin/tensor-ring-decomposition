function res = mult_mv2(B, x)

[d, ~] = size(x);
res = cell(d,1);
for k = 1:(d)
    [z1, n, z2] = size(x{k});
    v = x{k};
    [s1, n1, s2] = size(B{k});
    Y = zeros(s1*z1, n, s2*z2);
    Breshaped = reshape(B{k}, [s1, n, n, s2]);
    Y = reshape( permute(Breshaped, [4, 1, 2 ,3]), [s2*s1*n, n])...
        *reshape(permute(v, [2, 3, 1]), [n, z2*z1]);
    Y = reshape(Y, [s2, s1, n, z2, z1]);
    Y = permute(Y, [5, 2, 3, 4, 1]);
    Y = reshape(Y, [z1*s1, n, z2*s2]);
    res{k} = Y;
end
% pause
% 
% [z1, n, z2] = size(x{d});
% v = x{d};
% [s1, n1, s2] = size(B{d});
% Y = zeros(s1*z1,n, s2*z2);
% 
% Breshaped = reshape(B{d}, [s1, n, n, s2]);
% for i = 1:n
%     tmp = 0;
%     for j = 1:n    
% squeeze(v(:, j, :))
% 	v(:, j, :)
%         tmp = tmp + kron(squeeze(Breshaped(:, i, j, :)), squeeze(v(:, j, :)));
%     end
%     size(tmp)
%     size(Y(:, i, :))
%     Y(:, i, :) = tmp;
% end
% res{d} = Y;
