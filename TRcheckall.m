function res = TRcheckall(A, prec, m)
%%% Converts a tensor A in full format into TR-format,
%%% using relative accuracy prec.
%%% Loops over all cyclic shifts (with step size m)
%%% and divisors of the first unfolding rank of A
%%% and chooses the one with smallest storage cost.

s = size(A);
d = ndims(A);
norm_error = prec/sqrt(d-1);
curr_small = Inf;
res = cell(d,1);
curr_ind = -1;

for start_ind = 1:m:d
        B = permute(A, [start_ind:d, 1:(start_ind-1)]);
        s = size(B);
        unf = reshape(B, [s(1), prod(s(2:end))]);
        [U,S,V] = svd(unf,'econ');
        rnew = truncation_index(diag(S),norm_error*norm(diag(S)));
        factors = divisors(rnew);
        l = length(factors);
        for i = 1:l
            r0 = factors(i);
            trtmp = TRdecomp(B, prec, r0);
            
            a = storage_size(trtmp);
            if a < curr_small
                tr = trtmp;
                curr_ind = start_ind;
                curr_small = a;
            end
        end
end
%unpermute the cores
res = ipermuteTR(tr, curr_ind);