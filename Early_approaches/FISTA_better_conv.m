clear all;
close all;
%% FISTA Large condition number
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
x0 = cumsum(x0);

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

[x_k,cost] = FISTA(A,x0,n,lambda,f,grad,sigma, 1e-3);

fprintf("%.10f\n",cond(A))

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


function [x,cost] = FISTA(A,x0,n,lambda,f,grad,sigma, tolerance)
    
    % tunning
    L = norm(A'*A);
    tau = 1/L;
    
    %
    x_k = zeros(n,1);
    cost = f(x_k);      %initialize cost array
    tol = tolerance;   % Convergence tolerance
    % x_0, t, y_k
    y_k = x_k;
    t_k = 1;
    % main loop
    for k = 1:10000
        z_k = y_k-tau*grad(y_k);
        % second min problem
        mu = lambda*tau;
        x_k_new = sigma(z_k,mu); %4.1
        t_k_new = (1+sqrt(1+4*t_k^2))/2; %4.2
        y_k = x_k_new+(t_k-1)/t_k_new*(x_k_new-x_k); %4.3
        if norm(x_k_new - x0)^2 < tol %primal dual gap? %gradient small %use proximal grad t get optimal andcompare
            fprintf('Converged after %d iterations\n', k);
            break;
        end
        x_k = x_k_new;
        t_k = t_k_new;
        cost = [cost,f(x_k)];
        norm(grad(x_k))
    end
    x=x_k;
end