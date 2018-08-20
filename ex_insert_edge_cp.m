%%% Compares storage costs in TR-format to TR with edge inserted,
%%% when converting from canonical format.

clc; clf; clear all; close all
ds = 5:1:30;
store_TR = zeros(length(ds),1);
store_TT = zeros(length(ds),1);
store_edge = zeros(length(ds),1);
    
    
for d_ind = 1:length(ds)
    d = ds(d_ind);
    disp(d)
    a = 0; b = 1; n = 5; dx = (b-a)/(2^n-1);
    x = a:+dx:b; 
    V = cell(d,1);
    prec = 1e-9;
    maxn = d;

    %\sum_{k=1}^{r} \sin(x_1)\sin(x_{d/2})
    r = 20;
    for l = 1:d
        tmp = V{l};
        for k = 1:r
             tmp(:,k) = ones(size(x));
        end
        V{l} = tmp;
    end


    for k = 1:r
        for l = [1 floor(d/2)]
            tmp = V{l};
            tmp(:,k) = sin(k*x);
            V{l} = tmp;
        end
    end

    tr_help = CP2TRcheckallfaster(V,1,maxn,prec/2);
    tr2 = roundingbal(CP2TR(V, 1), prec);
    tr3 = check_all_edge(tr_help, prec/2);
    tr = CP2TRcheckallfaster(V,1,maxn,prec);

    store_TR(d_ind) = storage_size(tr);
    store_TT(d_ind) = storage_size(tr2);
    store_edge(d_ind) = storage_size(tr3);
end


figure(1)
hold on

plot(ds, store_TR./store_TT, 'bo-');
plot(ds, store_edge./store_TT, 'rs-');

set(gca,...
'FontUnits','points',...
'FontSize',24,...
'TickLabelInterpreter','latex',...
'FontName','Times')

xlim([5, 30]);
xticks([10, 20, 30]);
xlabel('$d$',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',24,...
    'FontName','Times')

yticks(0:0.2:1)
ylabel('Storage quotient',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',24,...
    'FontName','Times')
    grid()

legend({'$f:$ TR-format', '$f:$ Alg.~$4.4$'},...
        'location', 'BestOutside',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',24,...
        'FontName','Times')
    pbaspect([2 1 1])

% print -dpdf insert_edge_fig_journal.pdf