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
% load('G:\My Drive\MarquetteDr\Summer2023\walnut\results_cp_vs_admm/admm_within10.mat')


%%

% Chose parameters sigma, tau >0
% such that sigma*tau*norm(K)^2 =1

global N
global n
global A_global
A_global = A;
N = size(A,2); %vectorized image length
n = round(sqrt(N)); %image is n-by-n array
KtK = @(x) Kt(K(x));
pm(KtK)
%%
% pm returns 1.2513e+03, 625.6285 for 328
Lsqu = 625.6285;
sigma = 1;
tau =1/Lsqu;%1/pm(KtK);%
 
lambda = 0.01;
theta = 1;
epis = 0.01;
%initialize x and y
y_k = zeros(2*N,1);
x_k = zeros(N,1);
p_k = zeros(size(m(:),1),1);
x_bar_k = x_k;


% f = @(x) x+(tau/lambda)*A'*(A*x); %M*x
% argmin x norm(x-x0)^2+tau/lambda(norm(Ax-b)^2) 
% prox_tauG = @(x,x0) pcg(f, x+tau/lambda*A'*m(:),1e-6,100,[],[],x0);
cost = @(x) 0.5*norm(A*x-m(:))^2+lambda*norm(K(x),1);
my_cost = cost(x_k);
tic
for k = 1:1000
    prev_x_k = x_k;
    p_k = (p_k+sigma*(A*x_bar_k-m(:)))/(1+sigma);
    y_k = prox_sigF(y_k+sigma*K(x_bar_k),lambda); 
    x_k = prev_x_k-tau*A'*p_k-tau*Kt(y_k);
    x_bar_k = x_k+theta*(x_k-prev_x_k);
    %disp(k)
     my_cost = [my_cost,cost(x_k)];
     if (abs(cost(x_k) -ground_truth_admm) < 0.1*ground_truth_admm)
         break
     end
    %disp(cost(x_k))
end
elapsed_time = toc;
save_cost = my_cost(end);
save_image = x_k;
%save('admm_within10.mat', 'elapsed_time', 'save_cost', 'save_image');
% output x_k
%%
figure;
imshow(reshape(x_k,[n n]),[])
figure;
plot(my_cost)

%%
% argmin of y 1/2*norm(y-y0)^2+F*(y)
function clipped = prox_sigF(y0,lambda)
     y0(y0>lambda) = lambda; 
     y0(y0<-lambda) = -lambda;
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
    global A_global
    c= rand(N,1);
    for k= 1:100
        nc = norm(c);
        c = c/nc;
        c =  A_global'*(A_global*c) + KtK(c);
    end
    disp(nc)
end