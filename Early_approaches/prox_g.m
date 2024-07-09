clear all;
close all;
%% Proximal gradient
% Set up
m = 10;
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

% Figure out x using proximal gradient
f = @(x) norm(A*x -y)^2 ;
grad = @(x) A'*(A*x-y);

sigma = @(s,mu) max(abs(s)-mu,0).*sign(s);

% tunning
L = norm(A'*A);
tau = 1/L;
lambda = 0.1;

%
x_k = zeros(n,1);
cost = f(x_k);      %initialize cost array
tol = 1e-6;   % Convergence tolerance

% main loop
for k = 1:1000
    z_k = x_k-tau*grad(x_k);
    % second min problem
    mu = lambda*tau;
    x_k_new = sigma(z_k,mu);

    if norm(x_k_new - x_k)^2 < tol
        fprintf('Converged after %d iterations\n', k);
        break;
    end
    x_k = x_k_new;
    cost = [cost,f(x_k)];
end

figure;
stem(x0);
figure;
stem(x_k)

figure;
plot(cost, 'linewidth',2);
ylabel('cost');
xlabel('iteration');
set(gca,'fontsize',18);
set(gca,'linewidth',2);
