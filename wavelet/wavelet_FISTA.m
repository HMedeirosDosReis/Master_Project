clear all;
close all;
clc;
%% Haar wavelet FISTA Large condition
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
% A = A';
% [U,S,V] = svd(A);
% S(10,10) = S(10,10)*0.003;
% A = U*S*V';
% A = A';
% Haar wavelet transform --------WHERE SHOULD THIS GO
% need A_hat = A*inv(W) --- A is 10x20 -- how to get W?
% need c = W*x
%A = inverse_wavelet_transform(A);
%x0 = wavelet_transform(x0);

%
y = A*x0;
lambda = 0.1;

% Figure out x using proximal gradient
f = @(x) 0.5*norm(A*x -y)^2 +lambda*norm(x,1);
grad = @(x) A'*(A*x-y);

sigma = @(s,mu) max(abs(s)-mu,0).*sign(s);

% tunning
L = norm(A'*A);
tau = 1/L;

%
x_k = zeros(n,1);
cost = f(x_k);      %initialize cost array
tol = 1e-3;   % Convergence tolerance
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
%     if norm(x_k_new - x0)^2 < tol %primal dual gap? %gradient small %use proximal grad t get optimal andcompare
%         fprintf('Converged after %d iterations\n', k);
%         break;
%     end
    x_k = x_k_new;
    t_k = t_k_new;
    cost = [cost,f(x_k)];
    norm(grad(x_k))
end

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

function [X, Y_l, Y_h, Z_l, Z_h, Z] = wavelet_transform(X)
    
    [m,n] = size(X); 
    
    Y_l = zeros(m,n);
    Y_h = zeros(m,n);
    
        for i = 1:n
            x = X(:,i);
            
            y_l = conv(x,[1/sqrt(2),1/sqrt(2),0],"same"); 
            y_h = conv(x,[1/sqrt(2),-1/sqrt(2),0],"same");
            
            Y_l(:,i) = y_l;
            Y_h(:,i) = y_h;
        end
    
    Z_l = Y_l(1:2:end,:); 
    Z_h = Y_h(1:2:end,:);
    
    Z = [Z_l; Z_h];

end

%inverse 

function F = inverse_wavelet_transform(Z)

    [n, numCol] = size(Z);
    
    F = zeros(n, numCol); 
    
        for i = 1:numCol
            
            z = Z(:,i); 
            z_l = z(1:n/2);
            z_h = z(n/2+1:end);
            
            w_l = zeros(n,1); 
            w_h = zeros(n,1);
            
            w_l(1:2:end) = z_l;
            w_h(1:2:end) = z_h;
        
            f = conv(w_l,[0,1/sqrt(2),1/sqrt(2)],"same") + ...
                 conv(w_h,[0,-1/sqrt(2),1/sqrt(2)],"same");
                 
            F(:,i) = f;
        end

end
