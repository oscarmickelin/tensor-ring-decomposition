function final_cores = TRprealt2(A, prec)
%%% Converts a tensor A in full format into TR-format,
%%% using relative accuracy prec.
%%% Cyclic shift and divisor of first unfolding rank of A
%%% are chosen heuristically.

    d = ndims(A);
    norm_error = prec/sqrt(d-1);
    cores = cell(d,1);
    final_cores = cell(d,1);

    %%% pick initial permutation
    [start_ind, i_ranks] = choose_starting_index(A);
    A = permute(A, [start_ind:d, 1:(start_ind-1)]);
    s = size(A);

    C = A;
    rold = 1;
    %initial step
    n = s(1);
    C = reshape(C, [rold*n, numel(C)/(rold*n)]);
    [U,S,V] = svd(C,'econ');
    rnew = truncation_index(diag(S),norm_error*norm(S));

    factors = divisors(rnew);
    n_ind = start_ind -1;
    if n_ind == 0
        n_ind = d;
    end

    [~, m_ind] = min( abs(rnew./factors - i_ranks(n_ind)) + abs(factors-i_ranks(start_ind)));
    r0 = factors(m_ind);

    U=U(:,1:rnew);
    S=S(1:rnew, 1:rnew);
    V = V(:, 1:rnew);
    %%% reshaping U by stacking columns into 3-tensor
    Unew = zeros([rold*r0, n, rnew/r0]);
    for index = 1:r0
        Unew(index, :, :) = U(:, ( (index-1)*rnew/r0 + 1):(index*rnew/r0));
    end
    cores{1} = Unew;
    C = (S*V');

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
        rnew = truncation_index(diag(S),norm_error*norm(S));
        U=U(:,1:rnew);
        S=S(1:rnew, 1:rnew);
        V = V(:, 1:rnew);
        cores{k} = reshape(U, [rold, n, rnew]);
        C = (S*V');
        rold = rnew;
    end
    
    C = reshape(C, [rold, s(end), r0]);
    cores{d} = C;
    
    %unpermute the cores
    perm = [start_ind:d, 1:(start_ind-1)];
    order(perm) = 1:length(perm);
    for k = 1:d
        final_cores{k} = cores{order(k)};
    end
        
end