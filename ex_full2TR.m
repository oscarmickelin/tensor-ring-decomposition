%%% Compares storage costs in TT- and TR-formats, using heuristic
%%% and non-heuristic algorithms when converting from full format.

clc; clf; clear all; close all
d = 5; a = 0; b = 1; n = 20; N = n^d; h = (b-a)/(n-1);
e = ones(1,n);
x1 = a:+h:b;
x = zeros(d,N);
x(1,:) = kron(kron(kron(kron(e,e),e),e),x1);
x(2,:) = kron(kron(kron(kron(e,e),e),x1),e);
x(3,:) = kron(kron(kron(kron(e,e),x1),e),e);
x(4,:) = kron(kron(kron(kron(e,x1),e),e),e);
x(5,:) = kron(kron(kron(kron(x1,e),e),e),e);

% A = exp(cos(x(1,:).*x(5,:) + x(2,:)+x(3,:)+ x(4,:)));
A = 1./sqrt(1 + x(1,:).^2+ x(2,:).^2+ x(3,:).^2 + x(4,:).^2 + x(5,:).^2);
% A = exp(cos(x(1,:).*x(5,:) + x(2,:).*x(1,:) + x(4,:)+x(3,:)));
% A = exp(((x(1, :).*x(2,:) + x(5,:).*x(4,:).*x(3,:).*x(1,:).^10)));

% d = 4; a = 1e-10; b = 1; n = 20; N = n^d; h = (b-a)/(n-1);
% e = ones(1,n);
% x1 = a:+h:b;
% x = zeros(d,N);
% x(1,:) = kron(kron(kron(e,e),e),x1);
% x(2,:) = kron(kron(kron(e,e),x1),e);
% x(3,:) = kron(kron(kron(e,x1),e),e);
% x(4,:) = kron(kron(kron(x1,e),e),e);
% 
% A = x(1,:).*(sqrt(1 + (x(2,:) + x(3,:).^2).*x(4,:)./(x(1,:).^2))-1) ...
%     + (x(1,:) + 3*x(4,:)).*exp(1+sin(x(3,:)));

A = reshape(A, n*ones(1,d));

prec = 1e-12;
max_iter = 1;

tic
for k = 1:max_iter
zTT = TRdecomp(A, prec, 1);
end
t1=toc/max_iter;

[~,~,r1] = size(zTT{1});
divs = divisors(r1);
[~, r0] = min( abs(divs - r1./divs));

tic
for k = 1:max_iter
zBal = TRdecomp(A, prec, divs(r0));
end
tbal = toc/max_iter;

tic
for k = 1:max_iter
zTR = TRcheckall(A, prec, 1);
end
tall = toc/max_iter;

tic
for k = 1:max_iter
zPrealt2 = TRprealt2(A, prec);
end
tprealt2 = toc/max_iter;

disp('Storage quotients')
storage_size(zBal)/storage_size(zTT)
storage_size(zTR)/storage_size(zTT)
storage_size(zPrealt2)/storage_size(zTT)

disp('Runtime quotients')
tbal/t1
tall/t1
tprealt2/t1