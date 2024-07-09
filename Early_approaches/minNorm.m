clear all;
close all;

%% Min Norm
% Set up
m = 5;
n = 20;
rng(1);
%
A = randn(m,n);
x0 = zeros(n,1);
%
x0(5) = 2;
x0(10) = -1;
x0(18) = 1;
%
x0 = cumsum(x0);

y = A*x0;

% Figure out x using least norm
x_ln = (A'*A)'*((A'*A)*(A'*A)')\(A'*y);

figure;
stem(x0);
figure;
stem(x_ln)