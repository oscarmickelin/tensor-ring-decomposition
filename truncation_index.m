function rnew = truncation_index(vector, bound)

rnew = length(vector);
norm_vector = cumsum(vector(end:-1:1).^2);
l = length(norm_vector);
for i = l:-1:1
    if norm_vector(i) < bound^2
        rnew = l - i;
        break
    end
end
end