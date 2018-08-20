function res = interaction_rank(T, k)
%%% computes the k:th interaction rank of T
    d = ndims(T);
    s = size(T);
    n_ind = k+1;
    if k == d
        n_ind = 1;
    end
    T_r = reshape( permute( T, [k:d, 1:(k-1)]), ...
        [ s(k)*s( n_ind ), numel(T)/(s(k)*s( n_ind )) ] );
    res = rank( T_r );
end
    