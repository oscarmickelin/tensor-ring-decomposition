function res = norm_no_perm(tr)
%%% Computes the norm of a TR-representation,
%%% without using cyclic shifts.

[d, ~] = size(tr);
%balance to prevent inaccuracies
rs = zeros(d,1);
for k = 1:d
    tmp = tr{k};
    [r0, nk, r1] = size(tmp);
    rs(k) = norm(reshape(tmp, [r0*nk, r1]), 'fro');
    tr{k} = tr{k}/rs(k);
end

tmp = tr{1};
[r0, nk, r1] = size(tmp);
mat = zeros(r0^2, r1^2);
for j = 1:nk  
    mat = mat + kron(reshape(tmp(:,j,:), [r0,r1]),reshape(tmp(:,j,:), [r0,r1]));
end

for k = 2:d
    tmp = tr{k};
    [r0, nk, r1] = size(tmp);
    [ro, col] = size(mat);
    tmpres = zeros(ro, r1^2);
    for j = 1:nk
        for l = 1:ro
            A = reshape(mat(l,:)', [r0, r0]);
            B = reshape(tmp(:,j,:), [r0, r1]);
            ntmp = B'*A*B;
            tmpres(l,:) = tmpres(l,:) + ntmp(:)';
        end
    end
    mat = tmpres;
end
res = sqrt(abs(trace(mat)));    %sometimes arithmetic errors
                                %give a small inaccuracy. In case it is
                                %negative, we put abs here.
res = res*prod(rs);                                
end
