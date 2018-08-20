function res = mult(tr, number)
%%% multiplies a TR-representation by a number

res = tr;
res{1} = number*tr{1};
end