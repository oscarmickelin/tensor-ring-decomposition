%%% Simulates a d-dimensional wave equation with Dirichlet boundary
%%% conditions and initial conditions specified in canonical format.
%%% NOTE: requires Oseledet et al's TT-toolbox in the current path (https://github.com/oseledets/TT-Toolbox).

clc; clf; clear all; close all

d = 5;
a = 0; b = 1; n = 5; dx = (b-a)/(2^n-1);
x = a:+dx:b; 
V = cell(d,1);
prec = 1e-12;
prec_inv = 1e-12;
maxn = d;

%     %x1*xd + x2 + ... + x_(d-1)
% r = d-1;
% for k = 1:d
%     V{k} = ones(2^n, r);
% end
% for l = [1 d]
%     tmp = V{l};
%     tmp(:,1) = x;
%     V{l} = tmp;
% end
% 
% for l = 2:(d-1)
%     tmp = V{l};
%     tmp(:,l) = x;
%     V{l} = tmp;
% end
% 
% for k = 1:d
%     tmp = V{k};
%     tmp(end, :) = 0;
%     tmp(1, :) = 0;
%     V{k} = tmp;
% end

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
    for l = [1 d]
    tmp = V{l};
    tmp(:,k) = sin(k*x);
    V{l} = tmp;
    end
end
for k = 1:d
    tmp = V{k};
    tmp(end, :) = 0;
    tmp(1, :) = 0;
    V{k} = tmp;
end


disp('finding best starting TR')
x2 = CP2TRcheckallfaster(V,1,maxn,prec);    
x2old = x2;

xTT = roundingbal(CP2TR(V, 1), prec);
xTTold = xTT;
dt = 1e-3; 
maxT = 500*dt;
maxTiter = ceil(maxT/dt)


N = n*ones(d, 1)';
D2 = tt_qlaplace_dd(N);
B2 = D2*(-1/dx^2);


storages2 = zeros(maxTiter+1, 1);
storages2(1) = storage_size(x2);
storagesTT = zeros(maxTiter+1, 1);
storagesTT(1) = storage_size(xTT);

time = zeros(maxTiter, 1);
time2 = zeros(maxTiter, 1);
timeTT = zeros(maxTiter, 1);

norm_TR = zeros(maxTiter+1, 1);
norm_TT = zeros(maxTiter+1, 1);

disp('computing inverse matrix')
m = tt_eye(2, n*d) - (dt)^2*B2;
invm=tt_minres_selfprec2(m,  prec_inv);
invm = newTT2oldTTm(invm);
invm = unpermuteLaplace(invm, n);

disp('Starting BE now')
for k = 1:maxTiter
        if mod(k, 20) == 0
            disp(k)
        end
    
    tic;
    x2tmp = roundingTR(add( mult(x2, 2), mult(x2old, -1)), prec );
    x2next = roundingTR( mult_mv2(invm, x2tmp), prec);
    x2old = x2;
    x2 = x2next;
    storages2(k+1) = storage_size(x2);
    t = toc;
    time2(k) = t;
    
    tic;
    xTTtmp = roundingbal(add( mult(xTT, 2), mult(xTTold, -1)), prec );
    xTTnext = roundingbal( mult_mv2(invm, xTTtmp), prec);
    xTTold = xTT;
    xTT = xTTnext;
    storagesTT(k+1) = storage_size(xTT);
    t = toc;
    timeTT(k) = t;
end

figure(1)
semilogy(1:maxTiter+1, storages2, '-r')
hold on
semilogy(1:maxTiter+1, storagesTT, '--b')

legend('TR', 'TT')
ylabel('Storage cost')
xlabel('Iteration number')

figure(2)
hold on
plot(1:maxTiter, cumsum(time2), '-r')
plot(1:maxTiter, cumsum(timeTT), '--b')

legend('TR', 'TT')
ylabel('Runtime (s)')
xlabel('Iteration number')

storages2(500)/storagesTT(500)
sum(time2(1:500))/sum(timeTT(1:500))
