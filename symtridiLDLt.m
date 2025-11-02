% LDL^T decomposition for a symmetric tridiagonal matrix

a = 2*(71+2)^2 * ones(71,1);
[n,~] = size(a);
c = -1*(n+2)^2 * ones(70,1);

u = zeros(n,1);
for i=1:n
    u(i) = i / (n+1);
end

b = zeros(n,1);
% Building b
for i=1:n
    b(i) = 1 / (1-u(i));
end

l = []; % we know L is unitriangular
d = [];

d(1) = a(1); % the first diag element is just a_11
for i=2:n
    lx = c(i-1) / d(i-1);
    l = [l; lx]; % append the next element of l vector
    dx = a(i) - (l(i-1)^2 * d(i-1));
    d = [d; dx]; % append the next element of d vector
end

% Forw substitution (Lz = b) where z=Dy
z = [];
z(1) = b(1);
for i=2:n
    zx = b(i) - l(i-1) * z(i-1);
    z = [z; zx];
end
D = diag(d);
y = D \ z;

% Back substitution (L^Tx = y)
x = zeros(n,1);
x(n) = y(n);
for i=2:n
    k = n-i+1; % reverse the mapping
    x(k) = y(k) - x(k+1) * l(k);
end

% exact solution of f
x_true = zeros(n,1);
for i=1:n
    x_true(i) = (u(i)-1) * log(1-u(i));
end

rel_errors = zeros(n,1);
for i=1:n
    rel_errors(i) = abs(x(i) - x_true(i)) / x_true(i);
end
