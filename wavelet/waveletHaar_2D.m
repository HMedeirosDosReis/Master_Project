close all;
clear all;
clc;
rng(1);

% Need to modify the data. 
n = 8;
X = imread("luffy.jpg");
X = im2double(X);
X = rgb2gray(X);
X = imresize(X,4*[64 64]);

Z = (wavelet_transform(wavelet_transform(X)'))';

F = (inverse_wavelet_transform(inverse_wavelet_transform(Z)'))';

figure;
imagesc(Z);
axis image
colorbar;

figure;
imagesc(F);
axis image
colorbar;

function Z = wavelet_transform(X)
    
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
    disp("he")

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

