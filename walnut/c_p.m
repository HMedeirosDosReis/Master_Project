%PLOT COSTS
%%
clear all;
close all;
%%
% Load the measurement matrix and the sinogram from
% file Data164.mat
% load Data82 A m
load Data328 A m
load output_data.mat
%%

% Chose parameters sigma, tau >0
% such that sigma*tau*norm(K)^2 =1
sigma = 1; 

global N
global n
N = size(A,2); %vectorized image length
n = round(sqrt(N)); %image is n-by-n array
KtK = @(x) Kt(K(x));
% pm returns 7.9603
tau = 1/7.9603;%1/pm(KtK);
lambda = 0.01;
%initialize x and y
y_k = zeros(2*N,1);
x_k = zeros(N,1);
x_bar_k = x_k;


f = @(x) x+(tau/lambda)*A'*(A*x); %M*x
% argmin x norm(x-x0)^2+tau/lambda(norm(Ax-b)^2) 
prox_tauG = @(x,x0) pcg(f, x+tau/lambda*A'*m(:),1e-6,100,[],[],x0);
cost = @(x) 0.5*norm(A*x-m(:))^2+lambda*norm(K(x),1);
my_cost = cost(x_k);

tic
for k = 1:1000
    prev_x_k = x_k;
    prev_y_k = y_k;
    y_k = prox_sigF(prev_y_k+sigma*K(x_bar_k)); 
    x_k = prox_tauG(prev_x_k-tau*Kt(y_k),prev_x_k);
    x_bar_k = 2*x_k-prev_x_k;
    my_cost = [my_cost,cost(x_k)];
    if (abs(cost(x_k) -ground_truth_admm) < 0.1*ground_truth_admm)
        break
    end
    %disp(cost(x_k))
end
% output x_k
elapsed_time = toc;
save_cost = my_cost(end);
save_image = x_k;
save('cp_within10.mat', 'elapsed_time', 'save_cost', 'save_image');
%%
figure;
imshow(reshape(x_k,[n n]),[])
figure;
plot(my_cost)
%%
% argmin of y 1/2*norm(y-y0)^2+F*(y)
function clipped = prox_sigF(y0)
     y0(y0>1) = 1; 
     y0(y0<-1) = -1;
     clipped = y0;
end

function Kx = K(x)
    global n
    x = reshape(x,[n n]);
    Y_l = conv2(x, [1,-1,0],'same');
    Y_h = conv2(x, [1;-1;0],'same');

    Kx = [Y_l(:); Y_h(:)];
end

function Kty = Kt(y)
    global N
    global n
 
    y_l = reshape(y(1:N),[n n]);
    y_h = reshape(y(N+1:end),[n n]);
    Kty = conv2(y_l,[0,-1,1],'same')+conv2(y_h, [0;-1;1],'same');
    Kty = Kty(:);

end


function nc = pm(KtK)
    global N
    c= rand(N,1);
    for k= 1:100
        nc = norm(c);
        c = c/nc;
        c = KtK(c);
    end
    disp(nc)
end