function tt = newTT2oldTTm(t)
%%% converts the TT-representation in the TT-toolbox
%%% by Oseledets et al into the format used for
%%% the TR-format

d = length(t.ps)-1;
tt = cell(d,1);
ranks = t.r;
ns = t.n;
ms = t.m;

for k = 1:d
    cr=t.core; ps=t.ps; corek=cr(ps(k):ps(k+1)-1);
    tt{k} = reshape(corek, [ranks(k), ns(k)*ms(k), ranks(k+1)]);
end
    
