%%
clear all;
close all;

%%
% Load the measurement matrix and the sinogram from
% file Data164.mat
load Data82 A m

%%
% This is where you should implement the algo
% y = A*x0;
% 
% % Figure out x using least norm
x_ln = (A'*A)'*((A'*A)*(A'*A)')\(A'*m(:));

%x     = pcg(fun,b);


%%


% Take a look at the sinograms and the reconstructions
figure;

% First subplot
subplot(1,2,1);
imagesc(m);
colormap gray;
axis square;
axis off;
title('Sinogram, 120 projections');

% Second subplot
subplot(1,2,2);
imagesc(reshape(x, N, N));
colormap gray;
axis square;
axis off;
title({'Tikhonov reconstruction,'; '120 projections'});
