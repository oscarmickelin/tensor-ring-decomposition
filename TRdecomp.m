function cores = TRdecomp(A, prec, r0)
%%% Converts a tensor A in full format into TR-format,
%%% with first TR-rank r0 and using relative accuracy prec.
%%% r0 must be a divisor of the first unfolding rank of A.

    s = size(A);
    d = length(s);
    norm_error = prec/sqrt(d-1);
    cores = cell(d,1);
    C = A;
    rold = 1;
    
    %initial step
    n = s(1);
    C = reshape(C, [rold*n, numel(C)/(rold*n)]);
    [U,S,V] = svd(C,'econ');
    rnew = truncation_index(diag(S),norm_error*norm(diag(S)));
    %%% choose r0 dividing rnew
    factors = divisors(rnew);
    if nargin == 2
        if length(factors) > 1
            r0 = randsample(factors, 1);
        else
            r0 = 1;
        end
    end

    U=U(:,1:rnew);
    S=S(1:rnew, 1:rnew);
    V = V(:, 1:rnew);

    Unew = zeros([rold*r0, n, rnew/r0]);
    for index = 1:r0
        Unew(index, :, :) = U(:, ( (index-1)*rnew/r0 + 1):(index*rnew/r0));
    end
    cores{1} = Unew;
    C = (S*V');

    %reshape C into a 3-tensor
    Cnew = zeros([rnew/r0, numel(C)/(rnew), r0]);
    for index = 1:r0
        Cnew(:, :, index) = C(( (index-1)*(rnew/r0) + 1):(index*(rnew/r0)), :);
    end
    C = Cnew;
    %%%reshape to merge i_d and alpha_0
    %make it full
    C = reshape(C, [rnew/r0 s(2:end) r0]);
    %merge last two indices
    C = reshape(C, [rnew/r0 s(2:end-1) s(end)*r0]);
    %and take first unfolding of this
    C = reshape(C, [rnew/r0, numel(C)/(rnew/r0)]);
    rold = rnew/r0;    
    %%%
    %intermediate steps
    for k =2:d-1
        n = s(k);
        C = reshape(C, [rold*n, numel(C)/(rold*n)]);        
        [U,S,V] = svd(C,'econ');
        rnew = truncation_index(diag(S),norm_error*norm(diag(S)));
        U=U(:,1:rnew);
        S=S(1:rnew, 1:rnew);
        V = V(:, 1:rnew);
        cores{k} = reshape(U, [rold, n, rnew]);
        C = (S*V');
        rold = rnew;
    end
    
    C = reshape(C, [rold, s(end), r0]);
    cores{d} = C;
   
end