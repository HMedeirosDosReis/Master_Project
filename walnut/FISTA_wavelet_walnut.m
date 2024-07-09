tic
%%
clear all;
close all;
rng(1)
%%
% Load the measurement matrix and the sinogram from
% file Data164.mat
load Data328 A m

%% FISTA 

% Figure out x using wavelet and FISTA
% prev  0.5*norm(A*x -y)^2 +lambda*norm(x,1);

%[x, y_l, y_h, z_l, z_h, Z] = w_t(X);
%[x, y_l, y_h, z_l, z_h, Z] = w_t(Z');
%Z = Z';

%F = inv_t(Z);
%F = inv_t(F');

%W_inv = the inverse wavelet transform of A?
% c = W*x -> wavelet tranform of x?

% Need to call it twice inside the funcion? because of the transpose?

B = @(c) A*inv_t(c);

Bt = @(y) w_t(A'*y);
lambda = 0.1;
%f = @(x) 0.5*norm(A*W_inv*c -m(:))^2 +lambda*norm(c,1);
f = @(c) 0.5*norm(B(c) -m(:))^2 +lambda*norm(c,1);
grad = @(c) Bt(B(c)-m(:)); % is this gradient correct?
sigma = @(s,mu) max(abs(s)-mu,0).*sign(s);


% n = sqrt(size(A,2)); in the example they used this, why the sqrt?
n = size(A);
n = n(2);
tolerance = 1e-3;
%L = normest(A'*A);
%L = pm(B,Bt,n);
L = 625.6242; % normest(A'*A) %very expensive 
[c,cost] = FISTA(n,lambda,f,grad,sigma, tolerance,L); 

x = inv_t(c);

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
function [x,cost] = FISTA(n,lambda,f,grad,sigma, tolerance,L)
    
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

function nc = pm(B,Bt,n)
    c= rand(n,1);
    for k= 1:100
        nc = norm(c);
        c = c/nc;
        c = Bt(B(c));
    end
end

%% Wavelet

function Z = w_t(X)
    
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

function F = inv_t(Z)

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