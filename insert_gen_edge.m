function res = insert_gen_edge(tr, l, prec)
%%% Inserts an edge between indices 1 and l for a TR-representation tr.
%%% For the edge rank, the lowest storage result is chosen.

[d, ~] = size(tr);
tr_orig = tr;
normie = normTR(tr);
trsub = tr(1:l);

[r0, n, r1] = size(trsub{1});
trsub{1} = reshape( trsub{1}, [1, r0*n, r1]);
[rdlast, nl, rd] = size(trsub{l});
trsub{l} = reshape( trsub{l}, [rdlast, nl*rd, 1]);
divs = divisors(r1);

if l < d
    tr_rest = tr(l+1:d);
    [r0r, nr, r1r] = size(tr_rest{1});
    tr_rest{1} = reshape( tr_rest{1}, [1, r0r*nr, r1r]);
    [rdlastr, nr, rdr] = size(tr_rest{end});
    tr_rest{end} = reshape( tr_rest{end}, [rdlastr, nr*rdr, 1]);
    norm_v = normTR(tr_rest);
    ntot = (normTR(trsub)*norm_v)/normie;
else
    ntot = 1;
end

if l > 2
    tr_other_rest = tr(2:(l-1));
    [r0r, nr, r1r] = size(tr_other_rest{1});
    tr_other_rest{1} = reshape( tr_other_rest{1}, [1, r0r*nr, r1r]);
    [rdlastr, nr, rdr] = size(tr_other_rest{end});
    tr_other_rest{end} = reshape( tr_other_rest{end}, [rdlastr, nr*rdr, 1]);
    norm_other_v = normTR(tr_other_rest);
else
    norm_other_v = 1;
end
best_store = Inf;

for k = 1:length(divs)
    tr_tmp = tr;

    % round first cycle
    zsub_ins = roundingTR( TT2TR(trsub, divs(k)), prec/ntot);

    [a0, nprod, a1] = size( zsub_ins{1});
    zsub_ins{1} = reshape( zsub_ins{1}, [a0, r0, n, a1]);
    [adlast, nprod, ad] = size( zsub_ins{l});
    zsub_ins{l} = reshape( zsub_ins{l}, [adlast, nl, rd, ad]);

    tr_tmp(1:l) = zsub_ins(1:l);
    
    %round second cycle
    tr_2nd = cell(d-l+2, 1);
    tr_2nd(1:d-l+1) = tr_tmp(l:d);
    tr_2nd{d-l+2} = tr_tmp{1};

    [adlast, n2nd, rd, ad] = size( tr_2nd{1});
    tr_2nd{1} = permute( reshape(tr_2nd{1}, [adlast*n2nd, rd, ad]), [3,1,2]);
    
    [a0, r0, nl2nd, a1] = size(tr_2nd{d-l+2});
    tr_2nd{d-l+2} = permute( reshape(tr_2nd{d-l+2}, [a0, r0, nl2nd*a1]), [2,3,1]);

    tr_2nd = roundingTR(tr_2nd, prec/(normTR(tr_2nd)*norm_other_v*normie));
    
    [b0, nprod, b1] = size( tr_2nd{1});
    
    tr_2nd{1} = reshape( tr_2nd{1}, [b0, adlast, n2nd, b1]);
    [bdlast, nprod, bd] = size( tr_2nd{d-l+2});
    
    tr_2nd{d-l+2} =reshape( tr_2nd{d-l+2}, [bdlast, nl2nd, a1, bd]);
    
    tr_tmp(l:d) = tr_2nd(1:(d-l+1));
    tr_tmp{1} = tr_2nd{d-l+2};
    
    curr_store = storage_size(tr_tmp);
    if curr_store < best_store
        best_store = curr_store;
        res = tr_tmp;
    end
end