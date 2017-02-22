function [U,S,V] = randsvd2(A,p,q)

[m,n] = size(A);
nr = min(max(2*p,p+5),n);
r = randn(n,nr);
Y = A*r;
[QY,~,~] = qr(Y,0);

if (q > 0)
    Ap = A';
end

% perform subspace iteration

for j = 1:q
    Y = Ap*QY;
    [QY,~,~] = qr(Y,0);
    Y = A*QY;
    [QY,~,~] = qr(Y,0);
end

AY = QY'*A;
[U,S,V] = svd(AY,0);
U = QY*U;
U = U(:,1:p);
S = S(1:p,1:p);
V = V(:,1:p);

end

