close all;
clc;
clear all;

% Wavelet Transform

% RGB image
img = imread('luffy.jpg');
figure;
imshow(img)

% Retrieve approximation and detail coefficients for each channel R, G, B

% ------- Not sure what this do
[xaR, xhR, xvR, xdR] = dwt2(img(:,:,1), 'haar');
[xaG, xhG, xvG, xdG] = dwt2(img(:,:,2), 'haar');
[xaB, xhB, xvB, xdB] = dwt2(img(:,:,3), 'haar');

xa(:,:,1) = xaR; xa(:,:,2) = xaG; xa(:,:,3) = xaB;
xh(:,:,1) = xhR; xa(:,:,2) = xhG; xa(:,:,3) = xhB;
xv(:,:,1) = xvR; xa(:,:,2) = xvG; xa(:,:,3) = xvB;
xd(:,:,1) = xdR; xd(:,:,2) = xdG; xd(:,:,3) = xdB;

xA = xa/255;

figure, imshow(xA*0.3);
figure, imshow(log10(xh)*0.3);
figure, imshow(log10(xv)*0.3);
figure, imshow(log10(xd)*0.3);

% Apply DWT on each of the partial components above
[xaaR, xhhR, xvvR, xddR] = dwt2(xa(:,:,1), 'haar');
[xaaG, xhhG, xvvG, xddG] = dwt2(xa(:,:,2), 'haar');
[xaaB, xhhB, xvvB, xddB] = dwt2(xa(:,:,3), 'haar');

xaa(:,:,1) = xaaR; xaa(:,:,2) = xaaG; xaa(:,:,3) = xaaB;
xhh(:,:,1) = xhhR; xaa(:,:,2) = xhhG; xhh(:,:,3) = xhhB;
xvv(:,:,1) = xvvR; xvv(:,:,2) = xvvG; xvv(:,:,3) = xvvB;
xdd(:,:,1) = xddR; xdd(:,:,2) = xddG; xdd(:,:,3) = xddB;

xAA = xaa/255;

figure, imshow(xAA*0.3);
figure, imshow(log10(xhh)*0.5);
figure, imshow(log10(xvv)*0.5);
figure, imshow(log10(xdd)*0.5);

