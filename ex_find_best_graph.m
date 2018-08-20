%%% Compares storage costs of TR-representation
%%% to representation with edges inserted every b_sz steps

% clc; clf; clear all; close all
ds = 5:2:30;
store_TR = zeros(length(ds),1);
store_TT = zeros(length(ds),1);
store_edge = zeros(length(ds),1);
b_sz = 3;  %batch size
trs = cell(length(ds));

for d_ind = 1:length(ds)
    d = ds(d_ind);
    disp(d)
    a = 0; b = 1; n = 5; dx = (b-a)/(2^n-1);
    e = ones(1,n);
    x = a:+dx:b; 
    V = cell(d,1);
    prec = 1e-9;
    maxn = d;
    %x1*x3 + x3*x5 + ... 
    r = floor(d/2);
    for k = 1:d
        V{k} = ones(2^n, r);
    end
    
    prod_ind = 1:2:d-2;
    not_prod_ind = 2:2:d-1;
    curr_ind = 1;
    for l = 1:d
        if ismember(l, prod_ind)
            tmp = V{l};
            tmp(:,curr_ind) = x;
            V{l} = tmp;            
            tmp = V{l+2};
            tmp(:,curr_ind) = x;
            V{l+2} = tmp;  
            curr_ind = curr_ind + 1;
        end
    end
    tr = roundingTR(CP2TR(V,1), prec/(d+1));

    [normL, normR] = normSweep(tr);
    n_tr = normTR(tr);
    for b_ind = 1:(b_sz-1):d-1
        A = tr(b_ind:d);
        tmp = A{1};
        s = size(tmp);
        A{1} = reshape(A{1}, [1, prod(s(1:end-2))*s(end-1), s(end)]);     
        tmp2 = A{d-b_ind+1};
        s2 = size(tmp2);
        try
            A{d-b_ind+1} = reshape(A{d-b_ind+1}, [s2(1), s2(2)*s2(3)]);  
        catch
            s2 = [s2 1];
            A{d-b_ind+1} = reshape(A{d-b_ind+1}, [s2(1), s2(2)*s2(3)]);
        end
        prec_ind = prec/(d+1)*n_tr/(normR(b_ind)*normL(b_ind));
        Anew = insert_gen_edge2(A, 1, b_sz, prec_ind);     

        tmp = Anew{1};
        safter = size(tmp);
        Anew{1} = reshape(tmp, [s(1:end-1), safter(end)]);
        tmp2 = Anew{d-b_ind+1};
        s2after = size(tmp2);
        Anew{d-b_ind+1} = reshape(tmp2, [s2after(1:end-1), s2(end-1), s2(end)]);
        tr(b_ind:d) = Anew;
    end
    
    trs{d_ind} = tr;

    tr2 = CP2TRcheckallfaster(V,1,maxn,prec);
    tr3 = roundingTR(CP2TR(V,1), prec);
    store_edge(d_ind) = storage_size(tr);
    store_TR(d_ind) = storage_size(tr2);
    store_TT(d_ind) = storage_size(tr3);
end


figure(1)
hold on

plot(ds, store_TR./store_TT, 'kv-');
plot(ds, store_edge./store_TT, 'g^-');

set(gca,...
'FontUnits','points',...
'FontSize',18,...
'TickLabelInterpreter','latex',...
'FontName','Times')

xlabel('$d$',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',18,...
    'FontName','Times')

yticks(0:0.2:1)
ylabel('Storage quotient',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',18,...
    'FontName','Times')

legend({'$g:$ TR-format', '$g:$ Alg.~$4.4$'}, 'location', 'BestOutside',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',18,...
        'FontName','Times')
    grid()
    pbaspect([2 1 1])
    
    
%     print -dpdf insert_edge_structure_storage.pdf

    