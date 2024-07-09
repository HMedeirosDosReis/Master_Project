tic
%%
clear all;
close all;
%%
% Load the measurement matrix and the sinogram from
% file Data164.mat
load Data328 A m

%% Proximal Gradient Descent 
lambda = 0.1;
% Figure out x using proximal gradient
f = @(x) 0.5*norm(A*x(:)-m(:))^2+lambda*norm(x(:),1);
grad = @(x) A'*(A*x(:)-m(:));
sigma = @(s,mu) max(abs(s)-mu,0).*sign(s);

% x0 is ground truth, how do I get it in this dataset?
x0 = [];
% n = sqrt(size(A,2)); in the example they used this, why the sqrt?
n = size(A);
n = n(2);
tolerance = 1e-3;
L = 625.6242; % normest(A'*A) %very expensive 
[x,cost] = proximal_grad(A,x0,n,lambda,f,grad,sigma, tolerance,L); 

%%

h=0.1;
% Take a look at the sinograms and the reconstructions
figure;


% 
imagesc(reshape(x,sqrt(size(A,2)),sqrt(size(A,2))),[0,h]);
colormap gray;
axis square;
axis off;
toc
%% My function
% This is where you should implement the algo
function [x,cost] = proximal_grad(A,x0,n,lambda,f,grad,sigma, tolerance,L)
    % tunning
  
    tau = 1/L;
    % initialization
    x_k = zeros(n,1);
    cost = f(x_k);      %initialize cost array
    tol = tolerance;   % Convergence tolerance

    % main loop
    for k = 1:100
        z_k = x_k-tau*grad(x_k);
        % second min problem
        mu = lambda*tau;
        x_k_new = sigma(z_k,mu);
    
%         if norm(x_k_new - x0)^2 < tol
%             fprintf('Converged after %d iterations\n', k);
%             break;
%         end
        x_k = x_k_new;
        cost = [cost,f(x_k)];
    end
    x = x_k;
end