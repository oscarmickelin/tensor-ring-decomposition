%%% Compares storage costs in TT- and TR-formats,
%%% when converting from canonical format.
%%% Performs rounding on subsystems to speed up runtime.

clc; clf; clear all; close all

ds = 3:26
storages_TR = zeros(length(ds),1);
storages_TT = zeros(length(ds),1);
norms_TR = zeros(length(ds),1);
norms_TT = zeros(length(ds),1);

for dind = 1:length(ds)
    d = ds(dind);
    disp(d)
    a = 0; b = 1; n = 5; dx = (b-a)/(2^n-1);
    x = a:+dx:b; 
    prec = 1e-8;
    V = cell(d,1);

    %x1*xd + x2 + ... + x_(d-1)
%     r = d-1;
%     for k = 1:d
%         V{k} = ones(2^n, r);
%     end
%     for l = [1 d]
%         tmp = V{l};
%         tmp(:,1) = x;
%         V{l} = tmp;
%     end
% 
%     for l = 2:(d-1)
%         tmp = V{l};
%         tmp(:,l) = x;
%         V{l} = tmp;
%     end

r=20;
for l = 1:d
    tmp = V{l};
    for k = 1:r
         tmp(:,k) = 1-exp(-x.^2/k);
    end
    V{l} = tmp;
end


for k = 1:r
    for l = [1 d]
    tmp = V{l};
    tmp(:,k) = sin(k*x);
    V{l} = tmp;
    end
end


% % \sum_k sin(kx1)*sin(kxd)
% r=20;
% for l = 1:d
%     tmp = V{l};
%     for k = 1:r
%          tmp(:,k) = ones(size(x));
%     end
%     V{l} = tmp;
% end
% 
% 
% for k = 1:r
%     for l = [1 d]
%     tmp = V{l};
%     tmp(:,k) = sin(k*x);
%     V{l} = tmp;
%     end
% end

% r=20;
% for l = 1:d
%     tmp = V{l};
%     for k = 1:r
%          tmp(:,k) = x.^k;
%     end
%     V{l} = tmp;
% end
% 
% tmp = V{1};
% for k = 1:r  
%      tmp(:,k) = k*(r-k)*tmp(:,k);
% end
% V{1} = tmp;


    x = CP2TRcheckallfaster(V,1,d,prec);
    x2 = roundingbal(CP2TR(V,1), prec);
    storages_TR(dind) = storage_size(x);
    storages_TT(dind) = storage_size(x2);
    norms_TR(dind) = normTR(x);
    norms_TT(dind) = normTR(x2);
end

figure(1)
subplot(2,1,1)
semilogy(ds, storages_TR, 'ro-')
hold on
semilogy(ds, storages_TT, 'bs-')
grid()
set(gca,...
'FontUnits','points',...
'FontSize',18,...
'TickLabelInterpreter','latex',...
'FontName','Times')

legend({'TR', 'TT'}, 'location', 'northwest',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',18,...
        'FontName','Times')
xlabel('$d$',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',18,...
    'FontName','Times')

ylabel('Storage cost',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',18,...
    'FontName','Times')
subplot(2,1,2)
hold on
plot(ds, storages_TR./storages_TT, 'bo-')
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

ylabel('Storage cost quotient',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',18,...
    'FontName','Times')
grid()