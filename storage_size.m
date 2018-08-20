function total_size = storage_size(tr)
%%% Computes the number of parameters in the TR-representation
    d = length(tr);
    total_size = 0;
    for i = 1:d
        s = size(tr{i});
        total_size = total_size + prod(s);
    end
end