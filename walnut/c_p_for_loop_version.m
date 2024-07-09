%%
clear all;
close all;
%%
% Load the measurement matrix and the sinogram from
% file Data164.mat
load Data82 A m
%%

% Chose parameters sigma, tau >0
% such that sigma*tau*norm(K)^2 =1
sigma = 1; 

n = size(A,2);
KtK = @(x) conv2(x, [0,1,0;1,-2,1;0,1,0],'same');
tau = 1/(pm(KtK,n)^2);
lambda = 0.5;
%initialize x and y
y_k = zeros(n,1);
x_k = zeros(n,1);
x_bar_k = x_k;

% function for pcg -> doubt here, it this function correct?
%A'*A is too big to calculate in the large dataset
f = @(x) (eye(n)+tau/lambda*(A'*A)); 
% argmin x norm(x-x0)^2+tau/lambda(norm(Ax-b)^2) 
prox_tauG = @(x,A) pcg(f, x+tau/lambda*A'*m(:));

for k = 1:10
    prev_x_k = x_k;
    prev_y_k = y_k;
    y_k = prox_sigF(prev_y_k+sigma*K(x_bar_k)); 
    x_k = prox_tauG(prev_x_k-tau*Kt(y_k),A);
    x_bar_k = 2*x_k-prev_x_k;
end
% output x_k

%%
% argmin of y 1/2*norm(y-y0)^2+F*(y)
function clipped = prox_sigF(y0)
    for i=1:size(y0,2) % 2 or 1, make sure later
        if y0(i) > 1
            y0(i) = 1;
        end
        if y0(i) <= 1
            y0(i) = y0(i);
        end
        if y0(i) < -1
            y0(i) = -1;
        end
    end
    clipped = y0;
end

function Kx = K(X)

    [m,n] = size(X); 
    
    Y_l = zeros(m,n);
    Y_h = zeros(m,n);
    
        for i = 1:n
            x = X(:,i);
            
            y_l = conv(x,[1,-1,0],"same"); 
            y_h = conv(x,[1;-1;0],"same");
            
            Y_l(:,i) = y_l;
            Y_h(:,i) = y_h;
        end
    
    Z_l = Y_l(1:2:end,:); 
    Z_h = Y_h(1:2:end,:);
    
    Kx = [Z_l; Z_h];

end

function Kty = Kt(Z)
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
        
            f = conv(w_l,[0,-1,1],"same") + ...
                 conv(w_h,[0,-1,1],"same");
                 
            F(:,i) = f;
        end
    Kty = F;
end


function nc = pm(KtK,n)
    c= rand(n,1);
    for k= 1:100
        nc = norm(c);
        c = c/nc;
        c = KtK(c);
    end
end