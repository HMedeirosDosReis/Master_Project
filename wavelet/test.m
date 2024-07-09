close all;
clear all;
clc;
rng(1);

n = 8;
x = rand(n,1);

[x, y_l, y_h, z_l, z_h, z] = wavelet_transform(x);
f = inverse_wavelet_transform(z);

Fx = [y_l;y_h];

% Haar wavelet transform
A = rand(n,n);
[x, y_l, y_h, z_l, z_h,WA] = wavelet_transform(A);
A_hat = inverse_wavelet_transform(WA);

hold on
plot(A)
plot(A_hat)
hold off

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

