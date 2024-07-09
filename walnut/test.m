load Data328 A m
n = size(A);
n = n(2);
x = zeros(n,1);
A'*(A*x(:)-m(:))