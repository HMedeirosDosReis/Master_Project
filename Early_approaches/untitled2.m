n = 20;
rng(1);
%
x0 = zeros(n,1);
%
x0(5) = 2;
x0(10) = -1;
x0(18) = 1;
x0new = cumsum(x0);

figure;
stem(x0new)