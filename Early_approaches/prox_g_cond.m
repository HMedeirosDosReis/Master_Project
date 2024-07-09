clear all;
close all;
%% Proximal gradient Large condition number
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
% Rescale singular values of A
A = A';
[U,S,V] = svd(A);
S(10,10) = S(10,10)*0.003;
A = U*S*V';
A = A';
%
y = A*x0;
lambda = 0.1;
% Figure out x using proximal gradient
f = @(x) 0.5*norm(A*x -y)^2 +lambda*norm(x,1);
grad = @(x) A'*(A*x-y);

sigma = @(s,mu) max(abs(s)-mu,0).*sign(s);

[x_k,cost] = proximal_grad(A,x0,n,lambda,f,grad,sigma,1e-3);

%fprintf("%.10f\n",cond(A))

figure;
stem(x0);
figure;
stem(x_k)

figure;
semilogy(cost, 'linewidth',2);
ylabel('cost');
xlabel('iteration');
set(gca,'fontsize',18);
set(gca,'linewidth',2);

function [x,cost] = proximal_grad(A,x0,n,lambda,f,grad,sigma, tolerance)
    
    % tunning
    L = norm(A'*A);
    tau = 1/L;
    % initialization
    x_k = zeros(n,1);
    cost = f(x_k);      %initialize cost array
    tol = tolerance;   % Convergence tolerance

    % main loop
    for k = 1:100000
        z_k = x_k-tau*grad(x_k);
        % second min problem
        mu = lambda*tau;
        x_k_new = sigma(z_k,mu);
    
        if norm(x_k_new - x0)^2 < tol
            fprintf('Converged after %d iterations\n', k);
            break;
        end
        x_k = x_k_new;
        cost = [cost,f(x_k)];
    end
    x = x_k;
end
