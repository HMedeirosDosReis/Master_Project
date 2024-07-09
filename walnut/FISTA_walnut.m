tic
%%
clear all;
close all;
rng(1)

%%
% Load the measurement matrix and the sinogram from
% file Data164.mat
load Data328 A m

% estimates from powermethod using wavelet, no memory to compute A'A
% 164 = 1.2513e+03
% 82 = 2.5027e+03
% 328 = 625.6242
%% FISTA

lambda = 0.1;
% Figure out x FISTA
% prev  0.5*norm(A*x -y)^2 +lambda*norm(x,1);
f = @(x) 0.5*norm(A*x(:)-m(:))^2+lambda*norm(x(:),1);
grad = @(x) A'*(A*x(:)-m(:));
sigma = @(s,mu) max(abs(s)-mu,0).*sign(s);

n = size(A);
n = n(2);
tolerance = 1e-3;

L = 625.6242; % normest(A'*A) %very expensive 
[x, cost] = FISTA(A,n,lambda,f,grad,sigma, tolerance,L); 

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
function [x,cost] = FISTA(A,n,lambda,f,grad,sigma, tolerance,L)
    
    % tunning

    tau = 1/L;
    
    %
    x_k = zeros(n,1);
    cost = f(x_k);      %initialize cost array
    tol = tolerance;   % Convergence tolerance
    % x_0, t, y_k
    y_k = x_k;
    t_k = 1;
    % main loop
    for k = 1:1000
        z_k = y_k-tau*grad(y_k);
        % second min problem
        mu = lambda*tau;
        x_k_new = sigma(z_k,mu); %4.1
        t_k_new = (1+sqrt(1+4*t_k^2))/2; %4.2
        y_k = x_k_new+(t_k-1)/t_k_new*(x_k_new-x_k); %4.3
%         if norm(x_k_new - x0)^2 < tol %primal dual gap? %gradient small %use proximal grad t get optimal andcompare
%             fprintf('Converged after %d iterations\n', k);
%             break;
%         end
        x_k = x_k_new;
        t_k = t_k_new;
        cost = [cost,f(x_k)];
    end
    x=x_k;
end