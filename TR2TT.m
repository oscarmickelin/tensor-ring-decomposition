function tt = TR2TT(tr)
%%% converts a TR-representation into a TT-representation by summation and
%%% truncation

[d, ~] = size(tr);
[r0, ~, ~] = size(tr{1});
tt = cell(d,1);

tmp = tr{1};
tt{1} = tmp(1,:,:);
for k = 2:(d-1)
    tt{k} = tr{k};
end
tmp = tr{d};
tt{d} = tmp(:,:,1);

for l = 2:r0
    tt2 = cell(d,1);
    tmp = tr{1};
    tt2{1} = tmp(l,:,:);
    for k = 2:(d-1)
        tt2{k} = tr{k};
    end
    tmp = tr{d};
    tt2{d} = tmp(:,:,l);
    
    tt = add(tt, tt2);
end