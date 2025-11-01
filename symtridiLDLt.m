% LDL^T decomposition for a symmetric tridiagonal matrix

a = [2;12;45];
c = [4;12];
[n,~] = size(a);

l = []; % we know L is unitriangular
d = [];

d(1) = a(1); % the first diag element is just a_11
for i=2:n
    lx = c(i-1) / d(i-1);
    l = [l; lx];
    dx = a(i) - (l(i-1)^2 * d(i-1));
    d = [d; dx];
end

% Forw substitution