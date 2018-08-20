function [normL, normR] = normSweep(tr)
%can be done faster by reusing the G^{<k}
%as in the norm calculation

[d, ~] = size(tr);
normL = zeros(d,1);
normR = zeros(d,1);
for k = 1:d
    A = tr(1:k);
    s = size(A{k});
    if length(s) == 2
        s = [s 1];
    end
    A{k} = reshape(A{k}, [s(1), s(2)*s(3), 1]);
    normL(k) = normTR(A);
end

for k = d:-1:1
    A = tr(k:d);
    s = size(A{1});
    if length(s) == 2
        s = [s 1];
    end
    A{1} = reshape(A{1}, [1, s(1)*s(2), s(3)]);
    A
    normR(k) = normTR(A);
end
end
    